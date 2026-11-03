"""
fuentes.py - Lectura de ejes de carretera y puntos kilométricos.

Los datos pueden venir de varios archivos a la vez. Cada origen es un dict:
    {"tipo": "mdb" | "shp", "rol": "ejes" | "pk", "ruta": ..., "capa": ..., "epsg_origen": ...}
y cargar_fuentes() los junta en un único objeto `Datos`:
  - tipo "mdb": almacén .mdb de GeoMedia. Lector principal: mdbtools (comandos
                mdb-export / mdb-tables) con decodificación propia de la geometría
                binaria de GeoMedia; no necesita GDAL. Si mdbtools no está
                instalado, se usa el controlador "Geomedia" de GDAL/OGR (ODBC).
                Las coordenadas se toman tal cual (se asume EPSG:25830).
  - tipo "shp": shapefile exportado desde GeoMedia (pyshp); ignora "capa".

Todas las coordenadas salen en EPSG:25830 (ETRS89 / UTM 30N), en metros.
"""
import csv
import io
import math
import os
import re
import shutil
import struct
import subprocess
from collections import Counter
from pathlib import Path

import shapefile

csv.field_size_limit(10 ** 9)   # geometrías largas exportadas en hexadecimal

CAMPO_EJE_VIA = "NOMBRE"   # capa de ejes: nombre de la vía
CAMPO_PK_VIA = "ID_VIA"    # capa de PK: nombre de la vía
CAMPO_PK_NUM = "VIA_PK"    # capa de PK: número de PK
CAMPO_MUN_NOMBRE = "MUNICIPIO"   # capa de municipios (fondo del mapa)
CAMPO_MUN_COD = "COD_MUN"
CAMPO_MUN_PROV = "PROVINCIA"
CAMPO_PROV_NOMBRE = "PROVINCIA"   # capa de provincias (fondo del mapa)
CAMPO_PROV_COD = "COD_PROV"
EPSG_DESTINO = 25830


def norm(txt):
    return str(txt).strip().upper()


def a_float(txt):
    return float(str(txt).strip().replace(",", "."))


def _punto_etiqueta(m):
    """
    Dónde poner el nombre de un municipio: un punto garantizado dentro de su parte principal
    (el centro del conjunto puede caer fuera: p. ej. Almería incluye la isla de Alborán).
    Devuelve (punto, caja de la parte principal, área de la parte principal).
    """
    from shapely.geometry import Polygon
    mejor = None
    for a in m["anillos"]:
        try:
            p = Polygon(a)
            area = p.area
        except Exception:
            continue
        if mejor is None or area > mejor[0]:
            mejor = (area, a, p)
    if mejor is None:
        return m["centro"], m["bbox"], 0.0
    area, anillo, poligono = mejor
    xs, ys = [q[0] for q in anillo], [q[1] for q in anillo]
    caja = (min(xs), min(ys), max(xs), max(ys))
    try:
        r = poligono.representative_point()
        punto = (r.x, r.y)
    except Exception:
        punto = ((caja[0] + caja[2]) / 2, (caja[1] + caja[3]) / 2)
    return punto, caja, area


class Datos:
    """Tramos de eje y PK agrupados por vía (clave = nombre en mayúsculas)."""

    def __init__(self):
        self.tramos = {}    # clave -> [[(x, y), ...], ...]
        self.pks = {}       # clave -> {pk: (x, y)}
        self.nombres = {}   # clave -> nombre tal como viene en los datos
        self.avisos = []
        self.duplicados = []   # (vía, pk) repetidos; se resumen en un solo aviso
        self.pks_alt = {}      # clave -> {pk: [otras posiciones]} de los PK repetidos
        self.municipios = []   # fondo del mapa: {cod, nombre, provincia, anillos, bbox, centro}
        self.provincias = []   # fondo del mapa: {cod, nombre, anillos}
        self._redes = {}       # clave -> (aristas, índice espacial), se calculan al usarse

    def add_tramo(self, via, coords):
        k = norm(via)
        if k and len(coords) >= 2:
            self.nombres.setdefault(k, str(via).strip())
            self.tramos.setdefault(k, []).append(coords)

    def add_pk(self, via, pk, punto):
        k = norm(via)
        if not k:
            return
        self.nombres.setdefault(k, str(via).strip())
        d = self.pks.setdefault(k, {})
        if pk in d:
            self.duplicados.append((str(via).strip(), pk))
            self.pks_alt.setdefault(k, {}).setdefault(pk, []).append(punto)
        else:
            d[pk] = punto

    def add_municipio(self, cod, nombre, provincia, anillos, centro=None):
        anillos = [a for a in anillos if len(a) >= 3]
        if not anillos:
            return
        xs = [p[0] for a in anillos for p in a]
        ys = [p[1] for a in anillos for p in a]
        bbox = (min(xs), min(ys), max(xs), max(ys))
        if centro is None:
            centro = ((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2)
        self.municipios.append({"cod": str(cod or "").strip(), "nombre": str(nombre or "").strip(),
                                "provincia": str(provincia or "").strip(), "anillos": anillos,
                                "bbox": bbox, "centro": centro})

    def add_provincia(self, cod, nombre, anillos):
        anillos = [a for a in anillos if len(a) >= 3]
        if anillos:
            self.provincias.append({"cod": str(cod or "").strip(), "nombre": str(nombre or "").strip(),
                                    "anillos": anillos})

    def _red(self, clave, capa):
        """
        Red única de fronteras de una capa: se unen los contornos de todos sus polígonos, de modo
        que cada frontera compartida queda una sola vez, y se funden las aristas entre nudos.
        Se calcula la primera vez que se necesita y se reutiliza.
        """
        if clave not in self._redes:
            from shapely.geometry import LineString
            from shapely.ops import linemerge, unary_union
            from shapely.strtree import STRtree
            anillos = [a for e in capa for a in e["anillos"]]
            if not anillos:
                self._redes[clave] = ([], None)
            else:
                fus = linemerge(unary_union([LineString(a) for a in anillos]))
                bordes = list(getattr(fus, "geoms", [fus]))
                self._redes[clave] = (bordes, STRtree(bordes))
        return self._redes[clave]

    def _poligonos_tierra(self):
        """Construye y cachea los polígonos de provincias (o municipios) para pintar el relleno crema."""
        if "poligonos_tierra" not in self._redes:
            from shapely.geometry import Polygon
            from shapely.strtree import STRtree
            origen = self.provincias if self.provincias else self.municipios
            polis = []
            for item in origen:
                for a in item["anillos"]:
                    if len(a) >= 3:
                        try:
                            p = Polygon(a)
                            if not p.is_valid:
                                p = p.buffer(0)
                            if not p.is_empty:
                                for sub in getattr(p, "geoms", [p]):
                                    if sub.geom_type == "Polygon" and sub.area > 1.0:
                                        polis.append(sub)
                        except Exception:
                            continue
            self._redes["poligonos_tierra"] = (polis, STRtree(polis) if polis else None)
        return self._redes["poligonos_tierra"]

    def _rellenos(self, bbox, tolerancia):
        from shapely.geometry import box
        from shapely import clip_by_rect
        polis, arbol = self._poligonos_tierra()
        if arbol is None:
            return []
        nd = 0 if tolerancia >= 5 else 1
        res = []
        for i in arbol.query(box(*bbox)):
            simple = polis[i].simplify(tolerancia, preserve_topology=True)
            if simple.is_empty:
                continue
            try:
                recorte = clip_by_rect(simple, *bbox)
            except Exception:
                recorte = simple.intersection(box(*bbox))
            if recorte.is_empty:
                continue
            for parte in getattr(recorte, "geoms", [recorte]):
                if parte.geom_type == "Polygon" and len(parte.exterior.coords) >= 3:
                    res.append([(round(x, nd), round(y, nd)) for x, y in parte.exterior.coords])
        return res

    def _lineas(self, clave, capa, bbox, tolerancia):
        from shapely.geometry import box
        bordes, arbol = self._red(clave, capa)
        if arbol is None:
            return []
        from shapely import clip_by_rect
        nd = 0 if tolerancia >= 5 else 1
        res = []
        for i in arbol.query(box(*bbox)):
            simple = bordes[i].simplify(tolerancia, preserve_topology=False)
            recorte = clip_by_rect(simple, *bbox)
            if recorte.is_empty:
                continue
            for parte in getattr(recorte, "geoms", [recorte]):
                if parte.geom_type == "LineString" and len(parte.coords) >= 2:
                    res.append([(round(x, nd), round(y, nd)) for x, y in parte.coords])
        return res

    def fondo_en(self, bbox, tolerancia):
        """
        Fondo de mapa para el rectángulo `bbox`: (lineas_municipios, etiquetas, lineas_provincias, rellenos_tierra).
        """
        lineas = self._lineas("municipios", self.municipios, bbox, tolerancia)
        provincias = self._lineas("provincias", self.provincias, bbox, tolerancia)
        rellenos = self._rellenos(bbox, tolerancia)
        x0, y0, x1, y1 = bbox
        etiquetas = []
        for m in self.municipios:
            a0, b0, a1, b1 = m["bbox"]
            if a1 < x0 or a0 > x1 or b1 < y0 or b0 > y1:
                continue
            if "etq" not in m:
                m["etq"] = _punto_etiqueta(m)
            centro, caja, area = m["etq"]
            etiquetas.append({"nombre": m["nombre"], "provincia": m["provincia"],
                              "centro": centro, "bbox": caja, "area": area})
        return lineas, etiquetas, provincias, rellenos

    def lista_tramos(self, via):
        return self.tramos.get(norm(via), [])

    def pks_multi(self, via):
        """[(pk, [posiciones])] ordenado por pk; un PK repetido en los datos tiene varias posiciones."""
        k = norm(via)
        alt = self.pks_alt.get(k, {})
        return [(pk, [p] + alt.get(pk, [])) for pk, p in sorted(self.pks.get(k, {}).items())]

    def lista_pks(self, via):
        return sorted(self.pks.get(norm(via), {}).items())

    def rangos_pk_via(self, via):
        """
        Devuelve la lista de intervalos [(pk_ini, pk_fin), ...] continuos de PK de una vía.
        Si entre dos PKs consecutivos hay un salto > 2 km Y además la distancia en línea recta
        entre sus coordenadas supera los 2.500 m (vía discontinua como la N-323A), abre un nuevo intervalo.
        """
        k = norm(via)
        pks_ord = self.pks_multi(k)
        if not pks_ord:
            return []
        if len(pks_ord) == 1:
            v = pks_ord[0][0]
            return [(v, v)]

        rangos = []
        ini_actual = pks_ord[0][0]
        prev_pk, prev_pts = pks_ord[0]

        for curr_pk, curr_pts in pks_ord[1:]:
            salto_pk = curr_pk - prev_pk
            dist_min_m = min(
                math.hypot(p1[0] - p2[0], p1[1] - p2[1])
                for p1 in prev_pts for p2 in curr_pts
            )
            if salto_pk > 2.0 and dist_min_m > 2500.0:
                rangos.append((ini_actual, prev_pk))
                ini_actual = curr_pk
            prev_pk, prev_pts = curr_pk, curr_pts

        rangos.append((ini_actual, prev_pk))
        return rangos

    def resumen_vias(self):
        res = []
        for k in sorted(self.tramos):
            pks = self.pks.get(k, {})
            rangos = self.rangos_pk_via(k)
            res.append({
                "via": self.nombres[k],
                "tramos": len(self.tramos[k]),
                "n_pk": len(pks),
                "pk_min": min(pks) if pks else None,
                "pk_max": max(pks) if pks else None,
                "rangos_pk": [[r[0], r[1]] for r in rangos]
            })
        return res


# ------------------------------------------------------------- Shapefile
def leer_ejes_shp(datos, ruta, enc="latin-1"):
    r = shapefile.Reader(str(ruta), encoding=enc)
    _exigir_campos([f[0] for f in r.fields[1:]], [CAMPO_EJE_VIA], "ejes")
    n = 0
    for sr in r.iterShapeRecords():
        via = sr.record[CAMPO_EJE_VIA]
        pts = sr.shape.points
        partes = list(sr.shape.parts) + [len(pts)]
        for a, b in zip(partes[:-1], partes[1:]):
            datos.add_tramo(via, [(p[0], p[1]) for p in pts[a:b]])
            n += 1
    return n


def leer_pk_shp(datos, ruta, enc="latin-1"):
    r = shapefile.Reader(str(ruta), encoding=enc)
    _exigir_campos([f[0] for f in r.fields[1:]], [CAMPO_PK_VIA, CAMPO_PK_NUM], "PK")
    n = 0
    for sr in r.iterShapeRecords():
        if not sr.shape.points:
            continue
        try:
            pk = a_float(sr.record[CAMPO_PK_NUM])
        except (TypeError, ValueError):
            continue
        p = sr.shape.points[0]
        datos.add_pk(sr.record[CAMPO_PK_VIA], pk, (p[0], p[1]))
        n += 1
    return n


def leer_municipios_shp(datos, ruta, enc="latin-1"):
    r = shapefile.Reader(str(ruta), encoding=enc)
    campos = [f[0] for f in r.fields[1:]]
    _exigir_campos(campos, [CAMPO_MUN_NOMBRE], "municipios")
    n = 0
    for sr in r.iterShapeRecords():
        pts = sr.shape.points
        partes = list(sr.shape.parts) + [len(pts)]
        anillos = [[(p[0], p[1]) for p in pts[a:b]] for a, b in zip(partes[:-1], partes[1:])]
        rec = sr.record
        datos.add_municipio(rec[CAMPO_MUN_COD] if CAMPO_MUN_COD in campos else "", rec[CAMPO_MUN_NOMBRE],
                            rec[CAMPO_MUN_PROV] if CAMPO_MUN_PROV in campos else "", anillos)
        n += 1
    return n


def leer_provincias_shp(datos, ruta, enc="latin-1"):
    r = shapefile.Reader(str(ruta), encoding=enc)
    campos = [f[0] for f in r.fields[1:]]
    _exigir_campos(campos, [CAMPO_PROV_NOMBRE], "provincias")
    n = 0
    for sr in r.iterShapeRecords():
        pts = sr.shape.points
        partes = list(sr.shape.parts) + [len(pts)]
        anillos = [[(p[0], p[1]) for p in pts[a:b]] for a, b in zip(partes[:-1], partes[1:])]
        datos.add_provincia(sr.record[CAMPO_PROV_COD] if CAMPO_PROV_COD in campos else "",
                            sr.record[CAMPO_PROV_NOMBRE], anillos)
        n += 1
    return n


def _exigir_campos(existentes, necesarios, capa):
    faltan = [c for c in necesarios if c not in existentes]
    if faltan:
        raise RuntimeError(f"La capa de {capa} no tiene el campo {', '.join(faltan)}. "
                           f"Campos disponibles: {', '.join(existentes)}")


CAPA_DEFECTO = {"ejes": "CARRETERAS_BCN", "pk": "PK_VIA", "municipios": "MUNICIPIOS", "provincias": "PROVINCIAS"}


# ------------------------------------------- MDB con mdbtools (sin GDAL)
def mdbtools_disponible():
    return bool(shutil.which("mdb-export") and shutil.which("mdb-tables"))


def _ejecuta(cmd):
    try:
        return subprocess.run(cmd, capture_output=True)
    except FileNotFoundError:
        raise RuntimeError(f"No se encuentra el comando {cmd[0]} (instale mdbtools).")


def _error_cmd(p):
    return p.stderr.decode("latin-1", errors="ignore").strip()


def _tablas_mdb(ruta):
    if not os.path.exists(str(ruta)):
        raise RuntimeError(f"No existe el archivo '{ruta}'.")
    p = _ejecuta(["mdb-tables", "-1", str(ruta)])
    if p.returncode != 0:
        raise RuntimeError(f"mdb-tables no pudo abrir '{ruta}': {_error_cmd(p)}")
    return [t.strip() for t in p.stdout.decode("latin-1", errors="ignore").splitlines()
            if t.strip() and not t.startswith("MSys")]


def _cabecera_mdb(ruta, tabla):
    """Nombres de columna de una tabla (lee solo la primera línea del export)."""
    p = subprocess.Popen(["mdb-export", str(ruta), tabla], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL)
    try:
        linea = p.stdout.readline()
    finally:
        p.kill()
        p.stdout.close()
        p.wait()
    return next(csv.reader([linea.decode("latin-1", errors="ignore")]), [])


def _filas_mdb(ruta, tabla):
    """Filas de una tabla como dicts. Geometría binaria en hexadecimal si mdbtools lo permite."""
    p = _ejecuta(["mdb-export", "-b", "hex", str(ruta), tabla])
    hex_ok = p.returncode == 0
    if not hex_ok:
        p = _ejecuta(["mdb-export", str(ruta), tabla])
        if p.returncode != 0:
            raise RuntimeError(f"mdb-export falló en la tabla '{tabla}': {_error_cmd(p)}")
    texto = None
    if hex_ok:
        try:
            texto = p.stdout.decode("utf-8")
        except UnicodeDecodeError:
            pass
    if texto is None:
        texto = p.stdout.decode("latin-1", errors="ignore")
    rd = csv.DictReader(io.StringIO(texto))
    return list(rd.fieldnames or []), list(rd)


def _col(campos, nombre):
    for c in campos:
        if c.lower() == nombre.lower():
            return c
    return None


def _col_geom(campos):
    return _col(campos, "geometry") or next((c for c in campos if "geom" in c.lower()), None)


_HEX = re.compile(r"^[0-9A-Fa-f]+$")


def _raw_a_bytes(s):
    """Salida cruda de mdb-export: escapes \\xNN literales o caracteres latin-1."""
    out, i, n = [], 0, len(s)
    while i < n:
        if s[i] == "\\" and i + 3 < n and s[i + 1] == "x":
            try:
                out.append(int(s[i + 2:i + 4], 16))
                i += 4
                continue
            except ValueError:
                pass
        out.append(ord(s[i]) & 0xFF)
        i += 1
    return bytes(out)


def _a_bytes(txt):
    if not txt:
        return b""
    t = txt.strip()
    if t and len(t) % 2 == 0 and _HEX.match(t):
        return bytes.fromhex(t)
    return _raw_a_bytes(txt)


def _interpolar_arco_3pt(p1, p2, p3, n_segs=16):
    """Discretiza un arco circular de GeoMedia (0xC4) definido por 3 puntos en una polilínea."""
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    d = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-9:
        return [p1, p2, p3]
    ux = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / d
    uy = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    if r > 1e7 or not _finito(ux, uy):
        return [p1, p2, p3]
    a1 = math.atan2(y1 - uy, x1 - ux)
    a2 = math.atan2(y2 - uy, x2 - ux)
    a3 = math.atan2(y3 - uy, x3 - ux)

    def _norm_ang(a):
        return a % (2.0 * math.pi)

    da2 = _norm_ang(a2 - a1)
    da3 = _norm_ang(a3 - a1)
    if da2 > da3:
        da3 -= 2.0 * math.pi
    return [(ux + r * math.cos(a1 + da3 * i / n_segs),
             uy + r * math.sin(a1 + da3 * i / n_segs)) for i in range(n_segs + 1)]


def _extraer_subgeometrias_gdo(b):
    """
    Escanea un bloque binario GDO buscando todas las cabeceras internas `XX FF D2 0F`
    o bloques secuenciales de líneas/anillos.
    """
    res = []
    pos = 0
    n_bytes = len(b)
    while pos <= n_bytes - 20:
        idx = b.find(b"\xFF\xD2\x0F", pos + 1)
        if idx < 1:
            break
        inicio = idx - 1
        t_sub = b[inicio]
        datos_sub = b[inicio + 16:]
        if t_sub == 0xC1 and len(datos_sub) >= 48:
            x1, y1, _, x2, y2, _ = struct.unpack("<dddddd", datos_sub[:48])
            if _finito(x1, y1) and _finito(x2, y2):
                res.append([(x1, y1), (x2, y2)])
            pos = inicio + 16 + 48
            continue
        elif t_sub == 0xC4 and len(datos_sub) >= 72:
            x1, y1, _, x2, y2, _, x3, y3, _ = struct.unpack("<ddddddddd", datos_sub[:72])
            if _finito(x1, y1) and _finito(x2, y2) and _finito(x3, y3):
                res.append(_interpolar_arco_3pt((x1, y1), (x2, y2), (x3, y3)))
            pos = inicio + 16 + 72
            continue
        elif t_sub in (0xC2, 0xC3) and len(datos_sub) >= 4:
            n_pts = struct.unpack("<i", datos_sub[:4])[0]
            if 1 <= n_pts <= 2_000_000 and len(datos_sub) >= 4 + n_pts * 24:
                pts = [struct.unpack_from("<dd", datos_sub, 4 + i * 24) for i in range(n_pts)]
                pts_v = [p for p in pts if _finito(*p)]
                if len(pts_v) >= 2:
                    res.append(pts_v)
                pos = inicio + 20 + n_pts * 24
                continue
        pos = idx + 3
    return res


def decodifica_geometria(b):
    """
    Decodifica la geometría binaria GDO de GeoMedia.
    Soporta:
      - Puntos: 0xC0, 0xC8, 0xC9
      - Segmento de 2 puntos: 0xC1 (LineGeometry)
      - Polilínea: 0xC2 (PolylineGeometry)
      - Polígono / Anillo: 0xC3 (PolygonGeometry)
      - Arco de 3 puntos: 0xC4 (ArcGeometry)
      - Compuestas / Colecciones / Multilíneas: 0xC5, 0xC6, 0xC7, 0xCA, 0xCB
    Devuelve ("punto", (x, y)), ("linea", [(x, y), ...]), ("anillos", [[(x, y), ...], ...]),
    ("otro", tipo) o None si es ilegible.
    """
    try:
        if not b or len(b) < 20:
            return None
        if b[0] == 0xC2 and len(b) > 1 and b[1] != 0xFF:
            b = b[1:]
        if len(b) < 20:
            return None
        if not (b[1] == 0xFF and b[2] == 0xD2 and b[3] == 0x0F):
            pos = b.find(b"\xFF\xD2\x0F")
            if pos < 1:
                return None
            b = b[pos - 1:]
        if len(b) < 20:
            return None

        tipo, datos = b[0], b[16:]

        if tipo in (0xC0, 0xC8, 0xC9):
            if len(datos) < 24:
                return None
            x, y, _z = struct.unpack("<ddd", datos[:24])
            return ("punto", (x, y)) if _finito(x, y) else None

        if tipo == 0xC1:
            if len(datos) < 48:
                return None
            x1, y1, _z1, x2, y2, _z2 = struct.unpack("<dddddd", datos[:48])
            if _finito(x1, y1) and _finito(x2, y2):
                return ("linea", [(x1, y1), (x2, y2)])
            return None

        if tipo == 0xC4:
            if len(datos) < 72:
                return None
            x1, y1, _z1, x2, y2, _z2, x3, y3, _z3 = struct.unpack("<ddddddddd", datos[:72])
            if _finito(x1, y1) and _finito(x2, y2) and _finito(x3, y3):
                return ("linea", _interpolar_arco_3pt((x1, y1), (x2, y2), (x3, y3)))
            return None

        if tipo == 0xC2:
            n = struct.unpack("<i", datos[:4])[0]
            datos = datos[4:]
            if n < 1 or n > 500000 or len(datos) < n * 24:
                return None
            pts = []
            for i in range(n):
                x, y, _z = struct.unpack("<ddd", datos[i * 24:i * 24 + 24])
                if _finito(x, y):
                    pts.append((x, y))
            if len(pts) >= 2:
                return ("linea", pts)
            if len(pts) == 1:
                return ("punto", pts[0])
            return None

        if tipo in (0xC3, 0xC5, 0xC6, 0xC7, 0xCA, 0xCB):
            anillos = _anillos(b)
            if not anillos:
                anillos = _extraer_subgeometrias_gdo(b)
            return ("anillos", anillos) if anillos else None

        sub_partes = _extraer_subgeometrias_gdo(b)
        if sub_partes:
            return ("anillos", sub_partes)

        return ("otro", tipo)
    except (struct.error, IndexError):
        return None


def _anillos(b):
    """
    Extrae las listas de coordenadas [[(x, y), ...], ...] de un polígono o multilínea de GeoMedia.
    """
    if len(b) < 20:
        return []
    t = b[0]
    if t == 0xC1 and len(b) >= 64:
        x1, y1, _, x2, y2, _ = struct.unpack_from("<dddddd", b, 16)
        return [[(x1, y1), (x2, y2)]] if (_finito(x1, y1) and _finito(x2, y2)) else []
    if t == 0xC4 and len(b) >= 88:
        x1, y1, _, x2, y2, _, x3, y3, _ = struct.unpack_from("<ddddddddd", b, 16)
        if _finito(x1, y1) and _finito(x2, y2) and _finito(x3, y3):
            return [_interpolar_arco_3pt((x1, y1), (x2, y2), (x3, y3))]
        return []
    if t in (0xC2, 0xC3):
        n = struct.unpack("<i", b[16:20])[0]
        if n < 1 or n > 2_000_000 or len(b) < 20 + n * 24:
            return []
        pts = [struct.unpack_from("<dd", b, 20 + i * 24) for i in range(n)]
        return [[p for p in pts if _finito(*p)]]
    if t == 0xC5:
        tam = struct.unpack("<i", b[16:20])[0]
        res = _anillos(b[20:20 + tam])
        return res if res else _extraer_subgeometrias_gdo(b[20:])
    if t in (0xC6, 0xC7, 0xCA, 0xCB):
        k = struct.unpack("<i", b[16:20])[0]
        off, res = 20, []
        if 1 <= k <= 100000:
            try:
                for _ in range(k):
                    if off + 4 > len(b):
                        break
                    tam = struct.unpack("<i", b[off:off + 4])[0]
                    if tam <= 0 or off + 4 + tam > len(b):
                        break
                    res += _anillos(b[off + 4:off + 4 + tam])
                    off += 4 + tam
            except struct.error:
                pass
        return res if res else _extraer_subgeometrias_gdo(b[16:])
    return []


def _finito(x, y):
    return x == x and y == y and abs(x) != float("inf") and abs(y) != float("inf")


def _resolver_tabla(ruta, capa, rol, datos):
    tablas = _tablas_mdb(ruta)
    for t in tablas:
        if _clave_capa(t) == _clave_capa(capa or ""):
            return t
    if capa and capa != CAPA_DEFECTO[rol]:
        raise RuntimeError(f"No encuentro la tabla '{capa}'. Tablas disponibles: {', '.join(tablas)}")
    necesarios = {"ejes": [CAMPO_EJE_VIA], "pk": [CAMPO_PK_VIA, CAMPO_PK_NUM],
                  "municipios": [CAMPO_MUN_NOMBRE], "provincias": [CAMPO_PROV_NOMBRE]}[rol]
    cands = []
    for t in tablas:
        campos = _cabecera_mdb(ruta, t)
        if _col_geom(campos) and all(_col(campos, n) for n in necesarios):
            cands.append(t)
    if len(cands) == 1:
        datos.avisos.append(f"{Path(str(ruta)).name}: no hay tabla «{capa}»; se usa «{cands[0]}», "
                            f"la única con geometría y los campos {', '.join(necesarios)}.")
        return cands[0]
    if not cands:
        raise RuntimeError(f"Ninguna tabla tiene geometría y los campos {', '.join(necesarios)}. "
                           f"Tablas disponibles: {', '.join(tablas)}")
    raise RuntimeError(f"Varias tablas sirven ({', '.join(cands)}); indique cuál en «Capa».")


def _escala_mdb(ruta, tabla, col_geom, datos):
    """
    Factor (sx, sy) que pasa las coordenadas almacenadas a metros. GeoMedia guarda en
    GCoordSystem la matriz almacenamiento -> sistema real (p. ej. 0,01 si se guarda en cm),
    y GeometryProperties/FieldLookup enlazan cada tabla con su sistema.
    """
    try:
        _, lookup = _filas_mdb(ruta, "FieldLookup")
        _, props = _filas_mdb(ruta, "GeometryProperties")
        _, sistemas = _filas_mdb(ruta, "GCoordSystem")
    except RuntimeError:
        return 1.0, 1.0
    por_id = {(r.get("IndexID") or "").strip(): r for r in lookup
              if _clave_capa(r.get("FeatureName") or "") == _clave_capa(tabla)}
    cand = [p for p in props if (p.get("IndexID") or "").strip() in por_id]
    if not cand:
        return 1.0, 1.0
    pref = [p for p in cand if (por_id[(p["IndexID"] or "").strip()].get("FieldName") or "").lower()
            == col_geom.lower()] or cand
    guid = (pref[0].get("GCoordSystemGUID") or "").strip().upper()
    cs = next((c for c in sistemas if (c.get("CSGUID") or "").strip().upper() == guid), None)
    if cs is None:
        return 1.0, 1.0
    m = []
    for i in range(1, 17):
        try:
            m.append(float(cs.get(f"Stor2CompMatrix{i}") or 0))
        except ValueError:
            m.append(0.0)
    if any(abs(m[i]) > 1e-12 for i in (1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14)):
        raise RuntimeError(f"El sistema de coordenadas «{cs.get('Name') or guid}» de '{tabla}' tiene una "
                           "matriz con giros o desplazamientos, que esta lectura no soporta.")
    sx, sy = m[0], m[5]
    if sx <= 0 or sy <= 0:
        return 1.0, 1.0
    if abs(sx - 1) > 1e-12 or abs(sy - 1) > 1e-12:
        datos.avisos.append(f"{Path(str(ruta)).name} / {tabla}: coordenadas guardadas con escala {sx:g} "
                            f"(«{cs.get('Name') or 'sistema'}»); se convierten a metros.")
    return sx, sy


def _leer_mdbtools(datos, ruta, capa, rol):
    tabla = _resolver_tabla(ruta, capa, rol, datos)
    campos, filas = _filas_mdb(ruta, tabla)
    cg = _col_geom(campos)
    if rol == "ejes":
        cv, cp = _col(campos, CAMPO_EJE_VIA), None
        necesarios = [CAMPO_EJE_VIA]
    else:
        cv, cp = _col(campos, CAMPO_PK_VIA), _col(campos, CAMPO_PK_NUM)
        necesarios = [CAMPO_PK_VIA, CAMPO_PK_NUM]
    if not cg or not cv or (rol == "pk" and not cp):
        raise RuntimeError(f"La tabla '{tabla}' necesita geometría y los campos {', '.join(necesarios)}. "
                           f"Campos disponibles: {', '.join(campos)}")
    sx, sy = _escala_mdb(ruta, tabla, cg, datos)
    n, omitidas, vistas = 0, Counter(), 0
    for row in filas:
        g = decodifica_geometria(_a_bytes(row.get(cg) or ""))
        if g is None:
            omitidas["geometría ilegible"] += 1
            continue
        via = row.get(cv)
        if not via:
            continue
        if rol == "ejes":
            if g[0] == "linea":
                partes_linea = [g[1]]
            elif g[0] == "anillos":
                partes_linea = [sub for sub in g[1] if len(sub) >= 2]
            else:
                detalle = f"tipo 0x{g[1]:02X}" if g[0] == "otro" and isinstance(g[1], int) else g[0]
                omitidas[f"no es una línea ({detalle})"] += 1
                continue
            if not partes_linea:
                omitidas["línea sin vértices suficientes"] += 1
                continue
            for sub in partes_linea:
                coords = [(x * sx, y * sy) for x, y in sub]
                datos.add_tramo(via, coords)
        else:
            if g[0] != "punto":
                omitidas["no es un punto"] += 1
                continue
            try:
                pk = a_float(row.get(cp))
            except (TypeError, ValueError):
                omitidas["PK no numérico"] += 1
                continue
            coords = [(g[1][0] * sx, g[1][1] * sy)]
            datos.add_pk(via, pk, coords[0])
        n += 1
        if vistas < 1:
            vistas += 1
            x, y = coords[0]
            if abs(x) < 1000 and abs(y) < 1000:
                raise RuntimeError(f"Las coordenadas de '{tabla}' parecen grados ({x:g}, {y:g}), no UTM en "
                                   "metros. Esta lectura no reproyecta: exporte en EPSG:25830.")
    for motivo, k in omitidas.items():
        datos.avisos.append(f"{Path(str(ruta)).name} / {tabla}: {k} entidades omitidas ({motivo}).")
    return n


def leer_municipios_mdbtools(datos, ruta, capa=None):
    capa = capa or CAPA_DEFECTO["municipios"]
    tabla = _resolver_tabla(ruta, capa, "municipios", datos)
    campos, filas = _filas_mdb(ruta, tabla)
    cg, cn = _col_geom(campos), _col(campos, CAMPO_MUN_NOMBRE)
    if not cg or not cn:
        raise RuntimeError(f"La tabla '{tabla}' necesita geometría y el campo {CAMPO_MUN_NOMBRE}. "
                           f"Campos disponibles: {', '.join(campos)}")
    cc, cp = _col(campos, CAMPO_MUN_COD), _col(campos, CAMPO_MUN_PROV)
    cen = _col(campos, "centroide")
    sx, sy = _escala_mdb(ruta, tabla, cg, datos)
    n = omitidas = 0
    for row in filas:
        g = decodifica_geometria(_a_bytes(row.get(cg) or ""))
        if g is None or g[0] != "anillos":
            omitidas += 1
            continue
        anillos = [[(x * sx, y * sy) for x, y in a] for a in g[1]]
        centro = None
        if cen:
            c = decodifica_geometria(_a_bytes(row.get(cen) or ""))
            if c and c[0] == "punto":
                centro = (c[1][0] * sx, c[1][1] * sy)
        datos.add_municipio(row.get(cc) if cc else "", row.get(cn), row.get(cp) if cp else "", anillos, centro)
        n += 1
    if omitidas:
        datos.avisos.append(f"{Path(str(ruta)).name} / {tabla}: {omitidas} municipios omitidos (geometría no reconocida).")
    return n


def leer_provincias_mdbtools(datos, ruta, capa=None):
    capa = capa or CAPA_DEFECTO["provincias"]
    tabla = _resolver_tabla(ruta, capa, "provincias", datos)
    campos, filas = _filas_mdb(ruta, tabla)
    cg, cn = _col_geom(campos), _col(campos, CAMPO_PROV_NOMBRE)
    if not cg or not cn:
        raise RuntimeError(f"La tabla '{tabla}' necesita geometría y el campo {CAMPO_PROV_NOMBRE}. "
                           f"Campos disponibles: {', '.join(campos)}")
    cc = _col(campos, CAMPO_PROV_COD)
    sx, sy = _escala_mdb(ruta, tabla, cg, datos)
    n = omitidas = 0
    for row in filas:
        g = decodifica_geometria(_a_bytes(row.get(cg) or ""))
        if g is None or g[0] != "anillos":
            omitidas += 1
            continue
        datos.add_provincia(row.get(cc) if cc else "", row.get(cn),
                            [[(x * sx, y * sy) for x, y in a] for a in g[1]])
        n += 1
    if omitidas:
        datos.avisos.append(f"{Path(str(ruta)).name} / {tabla}: {omitidas} provincias omitidas (geometría no reconocida).")
    return n


def leer_ejes_mdbtools(datos, ruta, capa=None):
    return _leer_mdbtools(datos, ruta, capa or CAPA_DEFECTO["ejes"], "ejes")


def leer_pk_mdbtools(datos, ruta, capa=None):
    return _leer_mdbtools(datos, ruta, capa or CAPA_DEFECTO["pk"], "pk")


def _diagnostico_mdbtools(ruta):
    res = []
    for t in _tablas_mdb(ruta):
        campos = _cabecera_mdb(ruta, t)
        n = None
        if shutil.which("mdb-count"):
            p = _ejecuta(["mdb-count", str(ruta), t])
            try:
                n = int(p.stdout.decode("latin-1", errors="ignore").strip())
            except ValueError:
                pass
        g = _col_geom(campos)
        res.append({"capa": t, "geometria": f"columna {g}" if g else "sin geometría",
                    "entidades": n, "campos": campos, "srs": None})
    return res


def diagnostico_mdb(ruta):
    """Tablas del .mdb con sus campos; usa mdbtools si está instalado y, si no, GDAL."""
    return _diagnostico_mdbtools(ruta) if mdbtools_disponible() else _diagnostico_gdal(ruta)


# ------------------------------------------------------------------- MDB con GDAL
def _ogr():
    try:
        from osgeo import ogr, osr
    except ImportError:
        raise RuntimeError(
            "Falta GDAL para Python (módulo 'osgeo'). Hace falta una instalación de GDAL "
            "con soporte ODBC (por ejemplo la de OSGeo4W/QGIS). Mientras tanto, la app "
            "puede trabajar con shapefiles exportados desde GeoMedia.")
    ogr.UseExceptions()
    osr.UseExceptions()
    return ogr, osr


def abrir_mdb(ruta):
    ogr, _ = _ogr()
    drv = ogr.GetDriverByName("Geomedia")
    if drv is None:
        raise RuntimeError(
            "Esta instalación de GDAL no incluye el controlador 'Geomedia' (necesita "
            "soporte ODBC). Alternativas: instalar mdbtools (comando mdb-export), que la "
            "app usa automáticamente y no necesita GDAL, o exportar a shapefile.")
    try:
        ds = drv.Open(str(ruta), 0)
    except Exception as e:
        ds = None
        causa = str(e)
    else:
        causa = ""
    if ds is None:
        raise RuntimeError(
            f"No se pudo abrir '{ruta}' como almacén de GeoMedia. {causa} Compruebe que la "
            "ruta es correcta y que el controlador ODBC de Access (32 o 64 bits) coincide "
            "con la arquitectura de Python.")
    return ds


def _diagnostico_gdal(ruta):
    ogr, _ = _ogr()
    ds = abrir_mdb(ruta)
    capas = []
    for i in range(ds.GetLayerCount()):
        lyr = ds.GetLayer(i)
        defn = lyr.GetLayerDefn()
        srs = lyr.GetSpatialRef()
        try:
            n = lyr.GetFeatureCount()
        except Exception:
            n = None
        capas.append({
            "capa": lyr.GetName(),
            "geometria": ogr.GeometryTypeToName(lyr.GetGeomType()),
            "entidades": n,
            "campos": [defn.GetFieldDefn(j).GetName() for j in range(defn.GetFieldCount())],
            "srs": srs.GetName() if srs is not None else None,
        })
    return capas


def _clave_capa(s):
    return s.strip().upper().replace(" ", "_")


def _buscar_capa(ds, nombre):
    disponibles = []
    for i in range(ds.GetLayerCount()):
        lyr = ds.GetLayer(i)
        disponibles.append(lyr.GetName())
        if _clave_capa(lyr.GetName()) == _clave_capa(nombre):
            return lyr
    raise RuntimeError(f"No encuentro la capa '{nombre}'. Capas disponibles: {', '.join(disponibles)}")


def _campos_capa(lyr):
    defn = lyr.GetLayerDefn()
    return [defn.GetFieldDefn(j).GetName() for j in range(defn.GetFieldCount())]


def _transformador(lyr, epsg_origen, datos, nombre):
    _, osr = _ogr()
    destino = osr.SpatialReference()
    destino.ImportFromEPSG(EPSG_DESTINO)
    destino.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if epsg_origen:
        origen = osr.SpatialReference()
        origen.ImportFromEPSG(int(epsg_origen))
    else:
        origen = lyr.GetSpatialRef()
    if origen is None:
        datos.avisos.append(f"Capa {nombre}: no declara sistema de coordenadas; se asume EPSG:{EPSG_DESTINO}. "
                            "Si no lo es, indique el EPSG de origen.")
        return None
    origen = origen.Clone()
    origen.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if origen.IsSame(destino):
        return None
    datos.avisos.append(f"Capa {nombre}: sistema de origen «{origen.GetName()}», "
                        f"reproyectada a EPSG:{EPSG_DESTINO}. Compruebe que es el correcto.")
    return osr.CoordinateTransformation(origen, destino)


def _lineas(geom):
    ogr, _ = _ogr()
    t = ogr.GT_Flatten(geom.GetGeometryType())
    if t == ogr.wkbLineString:
        pts = geom.GetPoints() or []
        yield [(p[0], p[1]) for p in pts]
    elif t in (ogr.wkbMultiLineString, ogr.wkbGeometryCollection):
        for i in range(geom.GetGeometryCount()):
            yield from _lineas(geom.GetGeometryRef(i))


def _puntos(geom):
    ogr, _ = _ogr()
    t = ogr.GT_Flatten(geom.GetGeometryType())
    if t == ogr.wkbPoint:
        yield (geom.GetX(), geom.GetY())
    elif t in (ogr.wkbMultiPoint, ogr.wkbGeometryCollection):
        for i in range(geom.GetGeometryCount()):
            yield from _puntos(geom.GetGeometryRef(i))


def _geometria(feat, ct):
    g = feat.GetGeometryRef()
    if g is None:
        return None
    if ct is not None:
        g = g.Clone()
        g.Transform(ct)
    return g


def leer_ejes_ogr(datos, ds, capa, epsg_origen=None):
    lyr = _buscar_capa(ds, capa)
    _exigir_campos(_campos_capa(lyr), [CAMPO_EJE_VIA], f"ejes ({lyr.GetName()})")
    ct = _transformador(lyr, epsg_origen, datos, lyr.GetName())
    n = 0
    for feat in lyr:
        g = _geometria(feat, ct)
        via = feat.GetField(CAMPO_EJE_VIA)
        if g is None or via is None:
            continue
        for linea in _lineas(g):
            datos.add_tramo(via, linea)
            n += 1
    return n


def leer_pk_ogr(datos, ds, capa, epsg_origen=None):
    lyr = _buscar_capa(ds, capa)
    _exigir_campos(_campos_capa(lyr), [CAMPO_PK_VIA, CAMPO_PK_NUM], f"PK ({lyr.GetName()})")
    ct = _transformador(lyr, epsg_origen, datos, lyr.GetName())
    n = 0
    for feat in lyr:
        g = _geometria(feat, ct)
        via = feat.GetField(CAMPO_PK_VIA)
        try:
            pk = a_float(feat.GetField(CAMPO_PK_NUM))
        except (TypeError, ValueError):
            continue
        if g is None or via is None:
            continue
        for p in _puntos(g):
            datos.add_pk(via, pk, p)
            n += 1
            break
    return n


def leer_municipios_ogr(datos, ds, capa, epsg_origen=None):
    from shapely import wkt
    lyr = _buscar_capa(ds, capa)
    campos = _campos_capa(lyr)
    _exigir_campos(campos, [CAMPO_MUN_NOMBRE], f"municipios ({lyr.GetName()})")
    ct = _transformador(lyr, epsg_origen, datos, lyr.GetName())
    n = 0
    for feat in lyr:
        g = _geometria(feat, ct)
        if g is None:
            continue
        geo = wkt.loads(g.ExportToWkt())
        polis = [geo] if geo.geom_type == "Polygon" else list(getattr(geo, "geoms", []))
        anillos = []
        for p in polis:
            if p.geom_type == "Polygon":
                anillos.append([(x, y) for x, y, *_ in p.exterior.coords])
                anillos += [[(x, y) for x, y, *_ in i.coords] for i in p.interiors]
        datos.add_municipio(feat.GetField(CAMPO_MUN_COD) if CAMPO_MUN_COD in campos else "",
                            feat.GetField(CAMPO_MUN_NOMBRE),
                            feat.GetField(CAMPO_MUN_PROV) if CAMPO_MUN_PROV in campos else "", anillos)
        n += 1
    return n


def leer_provincias_ogr(datos, ds, capa, epsg_origen=None):
    from shapely import wkt
    lyr = _buscar_capa(ds, capa)
    campos = _campos_capa(lyr)
    _exigir_campos(campos, [CAMPO_PROV_NOMBRE], f"provincias ({lyr.GetName()})")
    ct = _transformador(lyr, epsg_origen, datos, lyr.GetName())
    n = 0
    for feat in lyr:
        g = _geometria(feat, ct)
        if g is None:
            continue
        geo = wkt.loads(g.ExportToWkt())
        polis = [geo] if geo.geom_type == "Polygon" else list(getattr(geo, "geoms", []))
        anillos = []
        for p in polis:
            if p.geom_type == "Polygon":
                anillos.append([(x, y) for x, y, *_ in p.exterior.coords])
                anillos += [[(x, y) for x, y, *_ in i.coords] for i in p.interiors]
        datos.add_provincia(feat.GetField(CAMPO_PROV_COD) if CAMPO_PROV_COD in campos else "",
                            feat.GetField(CAMPO_PROV_NOMBRE), anillos)
        n += 1
    return n


def cargar_fuentes(fuentes, enc_shp="latin-1"):
    """Lee todos los orígenes indicados y los junta en un único `Datos`."""
    d = Datos()
    abiertos = {}
    cuenta = {"ejes": 0, "pk": 0, "municipios": 0, "provincias": 0}
    fondo = ("municipios", "provincias")
    for i, f in enumerate(fuentes, start=1):
        rol = f.get("rol")
        tipo = f.get("tipo", "mdb")
        ruta = str(f.get("ruta", "")).strip()
        if rol not in cuenta:
            raise RuntimeError(f"Origen {i}: el rol debe ser 'ejes' o 'pk'.")
        if not ruta and rol in fondo:
            continue
        if not ruta:
            raise RuntimeError(f"Origen {i} ({rol}): falta la ruta del archivo.")
        capa = (f.get("capa") or "").strip() or CAPA_DEFECTO[rol]
        epsg = f.get("epsg_origen") or None
        try:
            if tipo == "mdb" and mdbtools_disponible():
                if epsg and int(epsg) != EPSG_DESTINO:
                    raise RuntimeError("El lector mdbtools no reproyecta: deje vacío «EPSG origen» "
                                       f"(se asume EPSG:{EPSG_DESTINO}) o use un origen ya en ese sistema.")
                leer = {"ejes": leer_ejes_mdbtools, "pk": leer_pk_mdbtools,
                        "municipios": leer_municipios_mdbtools, "provincias": leer_provincias_mdbtools}[rol]
                n = leer(d, ruta, capa)
            elif tipo == "mdb":
                if ruta not in abiertos:
                    abiertos[ruta] = abrir_mdb(ruta)
                leer = {"ejes": leer_ejes_ogr, "pk": leer_pk_ogr, "municipios": leer_municipios_ogr,
                        "provincias": leer_provincias_ogr}[rol]
                n = leer(d, abiertos[ruta], capa, epsg)
            elif tipo == "shp":
                leer = {"ejes": leer_ejes_shp, "pk": leer_pk_shp, "municipios": leer_municipios_shp,
                        "provincias": leer_provincias_shp}[rol]
                n = leer(d, ruta, f.get("enc") or enc_shp)
            else:
                raise RuntimeError(f"tipo de origen desconocido: {tipo}")
        except (RuntimeError, OSError, shapefile.ShapefileException) as e:
            if rol in fondo:
                d.avisos.append(f"{rol.capitalize()} ({Path(ruta).name}): no se pudo cargar; el mapa irá sin esa capa de fondo. {e}")
                continue
            raise RuntimeError(f"Origen {i} ({rol}, {Path(ruta).name}): {e}")
        if n == 0:
            d.avisos.append(f"Origen {i} ({rol}, {Path(ruta).name}): no se leyó ninguna entidad.")
        cuenta[rol] += n
    if d.duplicados:
        vias = sorted({v for v, _ in d.duplicados})
        ej = ", ".join(f"{v} PK {pk:g}" for v, pk in d.duplicados[:3])
        d.avisos.append(f"{len(d.duplicados)} PK repetidos en {len(vias)} vías (p. ej. {ej}); "
                        "en cada uno se elige la posición más cercana al eje de la vía.")
    if not cuenta["ejes"] or not cuenta["pk"]:
        raise RuntimeError("Hace falta al menos un origen de ejes y uno de PK con datos.")
    return d


def cargar_shp(ruta_ejes, ruta_pk, enc="latin-1"):
    return cargar_fuentes([{"tipo": "shp", "rol": "ejes", "ruta": ruta_ejes},
                           {"tipo": "shp", "rol": "pk", "ruta": ruta_pk}], enc)


def cargar_mdb(ruta, capa_ejes="CARRETERAS_BCN", capa_pk="PK_VIA", epsg_origen=None):
    return cargar_fuentes([{"tipo": "mdb", "rol": "ejes", "ruta": ruta, "capa": capa_ejes, "epsg_origen": epsg_origen},
                           {"tipo": "mdb", "rol": "pk", "ruta": ruta, "capa": capa_pk, "epsg_origen": epsg_origen}])


def cargar_ogr(ds, capa_ejes, capa_pk, epsg_origen=None):
    d = Datos()
    leer_ejes_ogr(d, ds, capa_ejes, epsg_origen)
    leer_pk_ogr(d, ds, capa_pk, epsg_origen)
    return d
