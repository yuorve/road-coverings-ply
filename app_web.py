#!/usr/bin/env python3
"""
app_web.py - Interfaz web local para generar archivos PLY de zonas de influencia.

Uso:   python app_web.py            (abre http://127.0.0.1:5000 en el navegador)
Datos: se guardan en proyecto_ply.json (carpeta actual, o variable PLY_PROYECTO).

Requisitos: pip install flask pyshp "shapely>=2"
            (GDAL con soporte ODBC solo si se lee directamente el .mdb)
"""
import csv as _csv
import io as _io
import json
import os
import threading
import webbrowser
from pathlib import Path

from flask import Flask, Response, jsonify, request
from werkzeug.exceptions import HTTPException

import fuentes as fu
import ply_zonas as pz

AQUI = Path(__file__).resolve().parent
PROYECTO = Path(os.environ.get("PLY_PROYECTO", "proyecto_ply.json"))

app = Flask(__name__)
estado = {"datos": None, "gen": None, "gen_clave": None}

DEFECTO = {
    "fuentes": [
        {"tipo": "mdb", "rol": "ejes", "ruta": "", "capa": "CARRETERAS_BCN", "epsg_origen": ""},
        {"tipo": "mdb", "rol": "pk", "ruta": "", "capa": "PK_info_PCA", "epsg_origen": ""},
        {"tipo": "mdb", "rol": "municipios", "ruta": "", "capa": "da02_term_munic_mod_AF_SHP", "epsg_origen": ""},
        {"tipo": "mdb", "rol": "provincias", "ruta": "", "capa": "da03_provincia_AF_SHP", "epsg_origen": ""},
    ],
    "opciones": {"ancho": 20, "quad_segs": 18, "esz_inicio": 1,
                 "clase": "Cobertura MVIA.F6_Coberturas", "saltos": "crlf", "enc_ply": "latin-1",
                 "mostrar_municipios": True},
    "mantenedores": [],
    "asignaciones": [],
    "ultimo_esz": 0,
}


class ErrorUsuario(Exception):
    pass


@app.errorhandler(Exception)
def _error(e):
    if isinstance(e, HTTPException):
        return e
    if isinstance(e, (ErrorUsuario, ValueError, RuntimeError, FileNotFoundError, KeyError)):
        return jsonify(error=str(e).strip("'\"") if isinstance(e, KeyError) else str(e)), 400
    raise e


# ------------------------------------------------------------- proyecto
def limpia_ruta(s):
    return str(s or "").strip().strip('"').strip("'")


def _parsear_esz_positivo(val):
    """Devuelve el entero positivo si `val` es un número > 0, o None en caso contrario."""
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        n = int(float(s.replace(",", ".")))
        return n if n > 0 else None
    except ValueError:
        return None


def _prioridad_esz(esz_val):
    """
    Prioridad para desempatar dos coberturas con el mismo tramo exacto:
    1º La que tenga número de ESZ positivo asignado.
    2º Si ambas tienen ESZ, la de menor valor numérico.
    """
    esz_pos = _parsear_esz_positivo(esz_val)
    return (0, esz_pos) if esz_pos is not None else (1, float("inf"))


def _formatear_rangos_pk(rangos):
    return " · ".join(f"{a:g} a {b:g}" if a != b else f"{a:g}" for a, b in rangos)


def _ajustar_y_validar_pks_en_tramos(pk0, pk1, rangos_pk, lim_inf_global, lim_sup_global):
    """
    Valida que el intervalo [pk0, pk1] caiga íntegramente dentro de UN ÚNICO tramo continuo
    de la vía (rangos_pk).
    - Si cae en un hueco entre tramos (ej. 40 a 66 en N-323A), lanza ValueError.
    - Si abarca dos o más tramos desconectados (ej. 25 a 110 en N-323A), lanza ValueError
      indicando qué tramos está intentando cruzar.
    - Si pertenece a un único tramo y solo sobresale por decimales (<= 1,5 km), lo recorta
      a los límites exactos de ese tramo.
    """
    p_ini, p_fin = min(pk0, pk1), max(pk0, pk1)

    if not rangos_pk:
        p_ini = max(lim_inf_global, p_ini)
        p_fin = min(lim_sup_global, p_fin)
        if p_ini >= p_fin:
            raise ValueError(f"PK fuera del rango de la vía ({lim_inf_global:g} a {lim_sup_global:g}).")
        return p_ini, p_fin

    n_r = len(rangos_pk)
    intervalos = []
    for idx, (r_a, r_b) in enumerate(rangos_pk):
        a_ef = min(float(r_a), lim_inf_global) if idx == 0 else float(r_a)
        b_ef = max(float(r_b), lim_sup_global) if idx == n_r - 1 else float(r_b)
        intervalos.append((a_ef, b_ef, float(r_a), float(r_b)))

    if n_r == 1:
        a_ef, b_ef, _, _ = intervalos[0]
        corte_a = max(p_ini, a_ef)
        corte_b = min(p_fin, b_ef)
        if corte_a >= corte_b:
            raise ValueError(f"El rango PK {p_ini:g} a {p_fin:g} está fuera del tramo de la vía ({a_ef:g} a {b_ef:g}).")
        return corte_a, corte_b

    TOL_DECIMAL_BORDE = 1.5

    tramos_ini = []
    tramos_fin = []
    tramos_tocados = []

    for idx, (a_ef, b_ef, r_a, r_b) in enumerate(intervalos):
        ini_tol = (a_ef - TOL_DECIMAL_BORDE) if idx == 0 else (r_a - TOL_DECIMAL_BORDE)
        fin_tol = (b_ef + TOL_DECIMAL_BORDE) if idx == n_r - 1 else (r_b + TOL_DECIMAL_BORDE)

        if ini_tol <= p_ini <= fin_tol:
            tramos_ini.append(idx)
        if ini_tol <= p_fin <= fin_tol:
            tramos_fin.append(idx)

        sol_a = max(p_ini, ini_tol)
        sol_b = min(p_fin, fin_tol)
        if sol_b > sol_a + 1e-4:
            tramos_tocados.append(idx)

    txt_r = _formatear_rangos_pk(rangos_pk)

    comunes = set(tramos_ini).intersection(tramos_fin)
    if len(comunes) == 1:
        idx_unico = next(iter(comunes))
        a_ef, b_ef, r_a, r_b = intervalos[idx_unico]
        lim_a = a_ef if idx_unico == 0 else r_a
        lim_b = b_ef if idx_unico == n_r - 1 else r_b
        corte_a = max(p_ini, lim_a)
        corte_b = min(p_fin, lim_b)
        if corte_a >= corte_b:
            raise ValueError(
                f"El rango PK {p_ini:g} a {p_fin:g} queda fuera del tramo ({r_a:g} a {r_b:g}). Tramos válidos: {txt_r}."
            )
        return corte_a, corte_b

    if len(tramos_tocados) >= 2 or (tramos_ini and tramos_fin and tramos_ini[0] != tramos_fin[-1]):
        idx_1 = tramos_ini[0] if tramos_ini else tramos_tocados[0]
        idx_2 = tramos_fin[-1] if tramos_fin else tramos_tocados[-1]
        r1 = rangos_pk[idx_1]
        r2 = rangos_pk[idx_2]
        raise ValueError(
            f"PK desde ({p_ini:g}) y PK hasta ({p_fin:g}) pertenecen a tramos discontinuos distintos "
            f"({r1[0]:g} a {r1[1]:g} y {r2[0]:g} a {r2[1]:g}). Cree una asignación por cada tramo ({txt_r})."
        )

    raise ValueError(
        f"El rango PK {p_ini:g} a {p_fin:g} se sale de los tramos continuos de la vía ({txt_r})."
    )


def normaliza_opciones(o):
    base = dict(DEFECTO["opciones"])
    base.update({k: v for k, v in (o or {}).items() if k in base})
    try:
        base["ancho"] = float(str(base["ancho"]).replace(",", "."))
        base["quad_segs"] = int(base["quad_segs"])
        base["esz_inicio"] = int(base["esz_inicio"])
    except (TypeError, ValueError):
        raise ErrorUsuario("Opciones no válidas: ancho, segmentos y ESZ inicial deben ser números.")
    if base["ancho"] <= 0:
        raise ErrorUsuario("El ancho del buffer debe ser mayor que 0.")
    if base["saltos"] not in ("crlf", "lf"):
        base["saltos"] = "crlf"
    base["mostrar_municipios"] = base["mostrar_municipios"] not in (False, "false", "False", 0, "0", "", None)
    return base


def migra(g):
    f = g.pop("fuente", None)
    if f and "fuentes" not in g:
        if f.get("tipo") == "shp":
            g["fuentes"] = [
                {"tipo": "shp", "rol": "ejes", "ruta": f.get("ruta_ejes", ""), "capa": "", "epsg_origen": ""},
                {"tipo": "shp", "rol": "pk", "ruta": f.get("ruta_pk", ""), "capa": "", "epsg_origen": ""}]
        else:
            g["fuentes"] = [
                {"tipo": "mdb", "rol": "ejes", "ruta": f.get("ruta", ""),
                 "capa": f.get("capa_ejes", "CARRETERAS_BCN"), "epsg_origen": f.get("epsg_origen", "")},
                {"tipo": "mdb", "rol": "pk", "ruta": f.get("ruta", ""),
                 "capa": f.get("capa_pk", "PK_info_PCA"), "epsg_origen": f.get("epsg_origen", "")}]
    return g


def sanea_fuentes(lista):
    res = []
    for f in lista or []:
        res.append({"tipo": "shp" if f.get("tipo") == "shp" else "mdb",
                    "rol": f.get("rol") if f.get("rol") in ("ejes", "pk", "municipios", "provincias") else "ejes",
                    "ruta": limpia_ruta(f.get("ruta")),
                    "capa": str(f.get("capa", "")).strip(),
                    "epsg_origen": str(f.get("epsg_origen", "")).strip()})
    for rol, capa in (
        ("ejes", "CARRETERAS_BCN"),
        ("pk", "PK_info_PCA"),
        ("municipios", "da02_term_munic_mod_AF_SHP"),
        ("provincias", "da03_provincia_AF_SHP"),
    ):
        if not any(f["rol"] == rol for f in res):
            res.append({"tipo": "mdb", "rol": rol, "ruta": "", "capa": capa, "epsg_origen": ""})
    return res


def cargar_proyecto():
    p = json.loads(json.dumps(DEFECTO))
    if PROYECTO.exists():
        guardado = migra(json.loads(PROYECTO.read_text(encoding="utf-8")))
        for k in p:
            if k in guardado:
                if isinstance(p[k], dict):
                    p[k].update(guardado[k])
                else:
                    p[k] = guardado[k]
    p["fuentes"] = sanea_fuentes(p["fuentes"])
    return p


def _clave_orden_asignacion(a):
    val = (a.get("value") or "").strip().upper()
    es_vacia = 1 if not val and not a.get("via") else 0
    via = (a.get("via") or "").strip().upper()
    try:
        pk0 = float(str(a.get("pk_desde") or "0").replace(",", "."))
    except ValueError:
        pk0 = 0.0
    esz_fijo = _parsear_esz_positivo(a.get("esz")) or 10**9
    return (es_vacia, val, via, pk0, esz_fijo)


def _mismo_sector_continuo(rangos_via, pk_fin_anterior, pk_ini_siguiente, max_hueco_km=30.0):
    """
    Comprueba si entre pk_fin_anterior y pk_ini_siguiente la carretera es continua
    o si por el contrario pertenecen a dos bloques discontinuos separados (ej. N-323A).
    """
    if pk_ini_siguiente <= pk_fin_anterior + 1.5:
        return True
    if (pk_ini_siguiente - pk_fin_anterior) > max_hueco_km:
        return False
    if not rangos_via or len(rangos_via) <= 1:
        return True
    for r_ini, r_fin in rangos_via:
        if (r_ini - 2.0) <= pk_fin_anterior <= (r_fin + 2.0) and (r_ini - 2.0) <= pk_ini_siguiente <= (r_fin + 2.0):
            return True
    return False


def _consolidar_coberturas(lista_asig):
    """
    0) Si dos o más coberturas en la misma Vía tienen exactamente el mismo tramo (solape al 100%),
       conserva únicamente la que tenga número de ESZ asignado (y si ambas lo tienen, la de menor ESZ).
    1) Elimina únicamente micro-fragmentos parásitos de cruces/enlaces (<= 1,2 km) embebidos
       dentro de un tramo mayor.
    2) Si un tramo válido de otro mantenedor está dentro de un tramo más largo, divide
       automáticamente el tramo largo en dos para evitar cualquier solape.
    3) Une tramos del mismo Mantenedor SOLO cuando son consecutivos en la vía.
    """
    por_via = {}
    sin_pk_por_via = {}

    for idx_orden, a in enumerate(lista_asig):
        via_k = fu.norm(a.get("via", ""))
        val = (a.get("value") or "").strip()
        esz_pos = _parsear_esz_positivo(a.get("esz"))
        raw_d = str(a.get("pk_desde", "") or "").strip()
        raw_h = str(a.get("pk_hasta", "") or "").strip()

        if raw_d == "" and raw_h == "":
            cand_sin_pk = {
                "esz": esz_pos,
                "via": a.get("via", ""),
                "pk_desde": "",
                "pk_hasta": "",
                "value": val,
                "_prio": (*_prioridad_esz(esz_pos), idx_orden)
            }
            actual = sin_pk_por_via.get(via_k)
            if actual is None or cand_sin_pk["_prio"] < actual["_prio"]:
                sin_pk_por_via[via_k] = cand_sin_pk
            continue

        try:
            p0 = fu.a_float(raw_d)
            p1 = fu.a_float(raw_h)
            if p0 > p1:
                p0, p1 = p1, p0
            por_via.setdefault(via_k, []).append({
                "esz": esz_pos,
                "via": a.get("via", ""),
                "p0": p0,
                "p1": p1,
                "len": p1 - p0,
                "value": val,
                "val_k": val.upper(),
                "_orden": idx_orden
            })
        except ValueError:
            pass

    sin_pk_num = [
        {k: v for k, v in item.items() if k != "_prio"}
        for item in sin_pk_por_via.values()
    ]

    consolidadas = []
    MAX_LONGITUD_PARASITO_KM = 1.2
    d = estado["datos"]

    for via_k, tramos_brutos in por_via.items():
        rangos_via = d.rangos_pk_via(via_k) if d is not None else []

        mejores_por_rango = {}
        for t in tramos_brutos:
            clave_rango = (round(t["p0"], 3), round(t["p1"], 3))
            prio_t = (*_prioridad_esz(t["esz"]), t["_orden"])
            existente = mejores_por_rango.get(clave_rango)
            if existente is None or prio_t < existente[0]:
                mejores_por_rango[clave_rango] = (prio_t, t)

        tramos = [par[1] for par in mejores_por_rango.values()]
        tramos.sort(key=lambda x: (x["p0"], x["p1"]))

        long_por_mant = {}
        for t in tramos:
            long_por_mant[t["val_k"]] = long_por_mant.get(t["val_k"], 0.0) + t["len"]

        limpios = []
        for i, cand in enumerate(tramos):
            es_parasito = False
            if cand["len"] <= MAX_LONGITUD_PARASITO_KM and long_por_mant.get(cand["val_k"], 0.0) <= MAX_LONGITUD_PARASITO_KM:
                for j, otro in enumerate(tramos):
                    if i == j or cand["val_k"] == otro["val_k"]:
                        continue
                    if abs(cand["p0"] - otro["p0"]) <= 1e-3 and abs(cand["p1"] - otro["p1"]) <= 1e-3:
                        continue
                    if long_por_mant.get(otro["val_k"], 0.0) >= max(5.0, cand["len"] * 4.0):
                        solape_ini = max(cand["p0"], otro["p0"] - 0.2)
                        solape_fin = min(cand["p1"], otro["p1"] + 0.2)
                        if (solape_fin - solape_ini) >= cand["len"] * 0.6:
                            es_parasito = True
                            break
            if not es_parasito:
                limpios.append(dict(cand))

        sin_envolventes = []
        for i, t in enumerate(limpios):
            fragmentos_t = [dict(t)]
            for j, otro in enumerate(limpios):
                if i == j or t["val_k"] == otro["val_k"]:
                    continue
                if otro["len"] < t["len"] - 1e-3:
                    nuevos_frag = []
                    for frag in fragmentos_t:
                        if frag["p0"] + 1e-3 < otro["p0"] and otro["p1"] < frag["p1"] - 1e-3:
                            izq = dict(frag, p0=frag["p0"], p1=otro["p0"], len=otro["p0"] - frag["p0"])
                            der = dict(frag, esz=None, p0=otro["p1"], p1=frag["p1"], len=frag["p1"] - otro["p1"])
                            nuevos_frag.extend([izq, der])
                        elif otro["p0"] <= frag["p0"] + 1e-3 and frag["p0"] + 1e-3 < otro["p1"] < frag["p1"] - 1e-3:
                            der = dict(frag, p0=otro["p1"], p1=frag["p1"], len=frag["p1"] - otro["p1"])
                            nuevos_frag.append(der)
                        elif frag["p0"] + 1e-3 < otro["p0"] < frag["p1"] - 1e-3 and otro["p1"] >= frag["p1"] - 1e-3:
                            izq = dict(frag, p0=frag["p0"], p1=otro["p0"], len=otro["p0"] - frag["p0"])
                            nuevos_frag.append(izq)
                        else:
                            nuevos_frag.append(frag)
                    fragmentos_t = nuevos_frag
            sin_envolventes.extend(fragmentos_t)

        unicos_post = {}
        for idx_p, t in enumerate(sin_envolventes):
            clave_r = (round(t["p0"], 3), round(t["p1"], 3))
            prio_t = (*_prioridad_esz(t["esz"]), t.get("_orden", idx_p))
            ex = unicos_post.get(clave_r)
            if ex is None or prio_t < ex[0]:
                unicos_post[clave_r] = (prio_t, t)

        lista_post = [par[1] for par in unicos_post.values()]
        lista_post.sort(key=lambda x: (x["p0"], x["p1"]))

        fusionados = []
        for t in lista_post:
            if not fusionados:
                fusionados.append(t)
                continue
            prev = fusionados[-1]
            if prev["val_k"] == t["val_k"] and _mismo_sector_continuo(rangos_via, prev["p1"], t["p0"], max_hueco_km=30.0):
                prev["p1"] = max(prev["p1"], t["p1"])
                prev["len"] = prev["p1"] - prev["p0"]
                if prev["esz"] is None:
                    prev["esz"] = t["esz"]
                elif t["esz"] is not None:
                    prev["esz"] = min(prev["esz"], t["esz"])
            else:
                fusionados.append(t)

        for a_t, b_t in zip(fusionados[:-1], fusionados[1:]):
            if a_t["val_k"] != b_t["val_k"] and b_t["p0"] < a_t["p1"]:
                solape_km = a_t["p1"] - b_t["p0"]
                if solape_km <= 2.5 and a_t["p0"] < b_t["p0"] and b_t["p1"] > a_t["p1"]:
                    pk_medio = round((a_t["p1"] + b_t["p0"]) / 2.0, 3)
                    a_t["p1"] = pk_medio
                    b_t["p0"] = pk_medio

        for f in fusionados:
            if f["p1"] - f["p0"] > 1e-4:
                consolidadas.append({
                    "esz": f["esz"],
                    "via": f["via"],
                    "pk_desde": f"{f['p0']:g}".replace(".", ","),
                    "pk_hasta": f"{f['p1']:g}".replace(".", ","),
                    "value": f["value"]
                })

    return consolidadas + sin_pk_num


def sanea_proyecto(b):
    """
    Valida, consolida tramos y ordena por Mantenedor.
    - Si una vía no tiene PKs (recorrido completo), deja obligatoriamente en blanco pk_desde y pk_hasta.
    - Si las asignaciones ya traen un ESZ positivo (> 0), se mantiene ese ESZ.
    """
    p = json.loads(json.dumps(DEFECTO))
    b = migra(dict(b))
    p["fuentes"] = sanea_fuentes(b["fuentes"]) if "fuentes" in b else p["fuentes"]
    p["opciones"] = normaliza_opciones(b.get("opciones"))
    vistos = set()
    p["mantenedores"] = [m for m in (" ".join(str(x).split()) for x in b.get("mantenedores", []))
                         if m and not (m in vistos or vistos.add(m))]
    asig = []
    d = estado["datos"]
    cache_limites = {}

    for a in b.get("asignaciones", []):
        esz_pos = _parsear_esz_positivo(a.get("esz"))
        via_str = str(a.get("via", "")).strip()
        pk_d = str(a.get("pk_desde", "")).strip()
        pk_h = str(a.get("pk_hasta", "")).strip()
        val = " ".join(str(a.get("value", "")).split())

        if val and val not in vistos:
            vistos.add(val)
            p["mantenedores"].append(val)

        if d is not None and via_str:
            k_via = fu.norm(via_str)
            if k_via not in cache_limites:
                pks_via = d.pks_multi(via_str)
                tramos_via = d.lista_tramos(via_str)
                if len(pks_via) >= 2 and tramos_via:
                    min_pk = float(pks_via[0][0])
                    max_pk = float(pks_via[-1][0])
                    min_ef, max_ef = pz.rango_efectivo_via(tramos_via, pks_via)
                    lim_inf = min_pk if min_pk <= min_ef + 0.05 else min_ef
                    lim_sup = max_pk if max_ef <= max_pk + 0.05 else max_ef
                    rangos_v = d.rangos_pk_via(via_str)
                    cache_limites[k_via] = ("CON_PK", min_pk, max_pk, lim_inf, lim_sup, rangos_v)
                elif tramos_via:
                    cache_limites[k_via] = ("SIN_PK",)
                else:
                    cache_limites[k_via] = None

            info_lim = cache_limites[k_via]
            if info_lim is not None:
                if info_lim[0] == "SIN_PK":
                    pk_d = ""
                    pk_h = ""
                elif info_lim[0] == "CON_PK":
                    _, min_pk, max_pk, lim_inf, lim_sup, rangos_v = info_lim
                    if len(rangos_v) == 1 or (pk_d != "" and pk_h != ""):
                        try:
                            v0 = min_pk if pk_d == "" else fu.a_float(pk_d)
                            v1 = max_pk if pk_h == "" else fu.a_float(pk_h)
                            adj0, adj1 = _ajustar_y_validar_pks_en_tramos(v0, v1, rangos_v, lim_inf, lim_sup)
                            pk_d = f"{adj0:g}".replace(".", ",")
                            pk_h = f"{adj1:g}".replace(".", ",")
                        except ValueError:
                            pass

        asig.append({"esz": esz_pos,
                     "via": via_str,
                     "pk_desde": pk_d,
                     "pk_hasta": pk_h,
                     "value": val})

    asig = _consolidar_coberturas(asig)

    p["mantenedores"].sort(key=lambda x: x.upper())
    asig.sort(key=_clave_orden_asignacion)

    esz_ini = p["opciones"]["esz_inicio"]
    esz_usados = {a["esz"] for a in asig if a.get("esz") is not None and a["esz"] > 0}
    siguiente_esz = max([esz_ini] + list(esz_usados)) if esz_usados else esz_ini

    for a in asig:
        if a.get("esz") is not None and a["esz"] > 0:
            a["esz"] = int(a["esz"])
        else:
            while siguiente_esz in esz_usados:
                siguiente_esz += 1
            a["esz"] = siguiente_esz
            esz_usados.add(siguiente_esz)
            siguiente_esz += 1

    ultimo = max(esz_usados) if esz_usados else (esz_ini - 1)
    p["asignaciones"], p["ultimo_esz"] = asig, ultimo
    return p


def guardar_proyecto(p):
    PROYECTO.write_text(json.dumps(p, ensure_ascii=False, indent=2), encoding="utf-8")


# --------------------------------------------------------------- datos
def conectar(fuentes):
    d = fu.cargar_fuentes(sanea_fuentes(fuentes))
    estado.update(datos=d, gen=None, gen_clave=None)
    return d


def generador(o):
    if estado["datos"] is None:
        raise ErrorUsuario("Primero conecte con los datos.")
    clave = (o["ancho"], o["quad_segs"])
    if estado["gen"] is None or estado["gen_clave"] != clave:
        estado["gen"] = pz.Generador(estado["datos"], o["ancho"], o["quad_segs"])
        estado["gen_clave"] = clave
    return estado["gen"]


def parsea_fila(a):
    via = str(a.get("via", "")).strip()
    if not via:
        raise ValueError("Falta la vía.")
    d = estado["datos"]
    if d is not None:
        via = d.nombres.get(fu.norm(via), via)
        pks_via = d.pks_multi(via)
        tramos_via = d.lista_tramos(via)
        rangos_via = d.rangos_pk_via(via)
    else:
        pks_via, tramos_via, rangos_via = [], [], []

    raw_desde = str(a.get("pk_desde", "") or "").strip()
    raw_hasta = str(a.get("pk_hasta", "") or "").strip()

    if len(pks_via) < 2:
        if not tramos_via:
            raise ValueError(f"La vía {via} no existe en la capa de carreteras.")
        a["pk_desde"] = ""
        a["pk_hasta"] = ""
        long_km = pz.longitud_total_km(tramos_via)
        return via, 0.0, long_km

    if len(rangos_via) > 1 and (not raw_desde or not raw_hasta):
        txt_r = _formatear_rangos_pk(rangos_via)
        raise ValueError(
            f"La vía {via} es discontinua ({txt_r}); indique PK desde y PK hasta dentro de uno de sus tramos."
        )

    min_pk = float(pks_via[0][0])
    max_pk = float(pks_via[-1][0])
    min_ef, max_ef = pz.rango_efectivo_via(tramos_via, pks_via)

    try:
        pk0 = min_pk if not raw_desde else fu.a_float(raw_desde)
        pk1 = max_pk if not raw_hasta else fu.a_float(raw_hasta)
    except ValueError:
        raise ValueError("PK desde / PK hasta deben ser números (admite coma decimal).")

    if pk0 >= pk1:
        raise ValueError(f"PK desde ({pk0:g}) debe ser menor que PK hasta ({pk1:g}).")

    lim_inf = min_pk if min_pk <= min_ef + 0.05 else min_ef
    lim_sup = max_pk if max_ef <= max_pk + 0.05 else max_ef

    pk0_adj, pk1_adj = _ajustar_y_validar_pks_en_tramos(pk0, pk1, rangos_via, lim_inf, lim_sup)
    a["pk_desde"] = f"{pk0_adj:g}".replace(".", ",")
    a["pk_hasta"] = f"{pk1_adj:g}".replace(".", ",")

    return via, pk0_adj, pk1_adj


def _mapear_columnas_csv(cabeceras):
    """Detecta las columnas VIA, PK inicio, PK fin, Mantenedor y opcionalmente ESZ."""
    mapa = {}
    limpias = {c.strip().lower(): c for c in cabeceras if c}
    for cand in ("esz", "id_esz", "cod_esz", "codigo_esz", "id"):
        if cand in limpias:
            mapa["esz"] = limpias[cand]
            break
    for cand in ("via", "vía", "carretera", "nombre", "id_via", "road"):
        if cand in limpias:
            mapa["via"] = limpias[cand]
            break
    for cand in ("pk inicio", "pk ini", "pk_inicio", "pk_ini", "pk_desde", "pk desde", "inicio", "desde"):
        if cand in limpias:
            mapa["pk_inicio"] = limpias[cand]
            break
    for cand in ("pk fin", "pk final", "pk_fin", "pk_hasta", "pk hasta", "fin", "hasta"):
        if cand in limpias:
            mapa["pk_fin"] = limpias[cand]
            break
    for cand in ("mantenedor", "value", "valor", "organismo", "titular", "conservacion"):
        if cand in limpias:
            mapa["value"] = limpias[cand]
            break
    return mapa


# ---------------------------------------------------------------- rutas
@app.get("/")
def index():
    return Response((AQUI / "index.html").read_text(encoding="utf-8"), mimetype="text/html")


@app.get("/api/proyecto")
def get_proyecto():
    d = estado["datos"]
    return jsonify(proyecto=cargar_proyecto(),
                   vias=d.resumen_vias() if d else None,
                   n_municipios=len(d.municipios) if d else 0,
                   n_provincias=len(d.provincias) if d else 0,
                   avisos=d.avisos if d else [])


@app.put("/api/proyecto")
def put_proyecto():
    p = sanea_proyecto(request.get_json(force=True))
    guardar_proyecto(p)
    return jsonify(proyecto=p)


@app.post("/api/capas")
def capas():
    ruta = limpia_ruta((request.get_json(force=True) or {}).get("ruta"))
    if not ruta:
        raise ErrorUsuario("Indique la ruta del archivo .mdb.")
    return jsonify(capas=fu.diagnostico_mdb(ruta))


@app.post("/api/conectar")
def post_conectar():
    d = conectar((request.get_json(force=True) or {}).get("fuentes", []))
    return jsonify(vias=d.resumen_vias(), n_municipios=len(d.municipios), n_provincias=len(d.provincias),
                   avisos=d.avisos)


@app.post("/api/importar_csv")
def importar_csv():
    d = estado["datos"]
    if d is None:
        raise ErrorUsuario("Primero conecte y lea los datos de origen (Paso 1) antes de importar el CSV.")

    if "archivo" not in request.files:
        raise ErrorUsuario("No se ha enviado ningún archivo CSV.")

    raw_bytes = request.files["archivo"].read()
    if not raw_bytes:
        raise ErrorUsuario("El archivo CSV está vacío.")

    texto = None
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            texto = raw_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    if texto is None:
        texto = raw_bytes.decode("latin-1", errors="ignore")

    primera_linea = next((l for l in texto.splitlines() if l.strip()), "")
    sep = ";"
    for s in (";", ",", "\t", "|"):
        if primera_linea.count(s) >= 3:
            sep = s
            break

    lector = _csv.DictReader(_io.StringIO(texto), delimiter=sep)
    if not lector.fieldnames:
        raise ErrorUsuario("El CSV no tiene cabeceras válidas.")

    col_map = _mapear_columnas_csv(lector.fieldnames)
    faltan = [k for k in ("via", "pk_inicio", "pk_fin", "value") if k not in col_map]
    if faltan:
        raise ErrorUsuario(
            f"Faltan columnas en el CSV ({', '.join(faltan)}). "
            f"Se esperaban: VIA, PK inicio, PK fin, Mantenedor (y opcionalmente ESZ). Detectadas: {list(lector.fieldnames)}"
        )

    tiene_col_esz = "esz" in col_map
    validas = []
    descartes = []
    mantenedores_nuevos = []
    cache_limites = {}
    es_vacio = lambda s: s == "" or s.lower() in ("nan", "none", "null", "-")

    for idx, fila in enumerate(lector, start=2):
        via_raw = str(fila.get(col_map["via"]) or "").strip()
        raw_ini = str(fila.get(col_map["pk_inicio"]) or "").strip()
        raw_fin = str(fila.get(col_map["pk_fin"]) or "").strip()
        mantenedor = " ".join(str(fila.get(col_map["value"]) or "").split())
        esz_csv = _parsear_esz_positivo(fila.get(col_map["esz"])) if tiene_col_esz else None

        if not via_raw and not raw_ini and not raw_fin and not mantenedor:
            continue

        clave_via = fu.norm(via_raw)
        if clave_via not in d.tramos:
            descartes.append({"fila": idx, "via": via_raw or "(vacía)", "pk_ini": raw_ini, "pk_fin": raw_fin,
                              "motivo": "La vía no existe en la capa CARRETERAS"})
            continue

        if not mantenedor:
            descartes.append({"fila": idx, "via": via_raw, "pk_ini": raw_ini, "pk_fin": raw_fin,
                              "motivo": "Falta el Mantenedor"})
            continue

        nombre_canonico = d.nombres.get(clave_via, via_raw)
        pks_via = d.pks_multi(via_raw)
        tramos_via = d.lista_tramos(nombre_canonico)

        if len(pks_via) < 2:
            validas.append({
                "esz": esz_csv,
                "via": nombre_canonico,
                "pk_desde": "",
                "pk_hasta": "",
                "value": mantenedor
            })
            if mantenedor not in mantenedores_nuevos:
                mantenedores_nuevos.append(mantenedor)
            continue

        if clave_via not in cache_limites:
            min_pk = float(pks_via[0][0])
            max_pk = float(pks_via[-1][0])
            min_ef, max_ef = pz.rango_efectivo_via(tramos_via, pks_via)
            lim_inf = min_pk if min_pk <= min_ef + 0.05 else min_ef
            lim_sup = max_pk if max_ef <= max_pk + 0.05 else max_ef
            rangos_v = d.rangos_pk_via(nombre_canonico)
            cache_limites[clave_via] = (min_pk, max_pk, lim_inf, lim_sup, rangos_v)

        min_pk, max_pk, lim_inf, lim_sup, rangos_v = cache_limites[clave_via]

        try:
            pk_ini_orig = min_pk if es_vacio(raw_ini) else fu.a_float(raw_ini)
            pk_fin_orig = max_pk if es_vacio(raw_fin) else fu.a_float(raw_fin)
        except ValueError:
            descartes.append({"fila": idx, "via": nombre_canonico, "pk_ini": raw_ini, "pk_fin": raw_fin,
                              "motivo": "PK no numérico"})
            continue

        if pk_ini_orig >= pk_fin_orig:
            descartes.append({"fila": idx, "via": nombre_canonico, "pk_ini": f"{pk_ini_orig:g}", "pk_fin": f"{pk_fin_orig:g}",
                              "motivo": f"PK inicio ({pk_ini_orig:g}) >= PK fin ({pk_fin_orig:g})"})
            continue

        try:
            pk_ini, pk_fin = _ajustar_y_validar_pks_en_tramos(
                pk_ini_orig, pk_fin_orig, rangos_v, lim_inf, lim_sup
            )
        except ValueError as err_rango:
            descartes.append({
                "fila": idx, "via": nombre_canonico, "pk_ini": f"{pk_ini_orig:g}", "pk_fin": f"{pk_fin_orig:g}",
                "motivo": str(err_rango)
            })
            continue

        validas.append({
            "esz": esz_csv,
            "via": nombre_canonico,
            "pk_desde": f"{pk_ini:g}".replace(".", ","),
            "pk_hasta": f"{pk_fin:g}".replace(".", ","),
            "value": mantenedor
        })
        if mantenedor not in mantenedores_nuevos:
            mantenedores_nuevos.append(mantenedor)

    validas = _consolidar_coberturas(validas)
    validas.sort(key=_clave_orden_asignacion)
    return jsonify(validas=validas, descartes=descartes, mantenedores=mantenedores_nuevos)


def fondo_ventana(cx, cy, escala, o, holgura=1.15):
    vacio = {"lineas": [], "provincias": [], "etiquetas": [], "rellenos": []}
    d = estado["datos"]
    if not o["mostrar_municipios"] or d is None or not (d.municipios or d.provincias):
        return vacio
    vw, vh = 900 / escala * holgura, 380 / escala * holgura
    bbox = (cx - vw / 2, cy - vh / 2, cx + vw / 2, cy + vh / 2)
    tol = max(1.0, 0.6 / escala)
    lineas, etiquetas, provincias, rellenos = d.fondo_en(bbox, tol)
    m_por_px = 1 / escala
    etiquetas = [e for e in etiquetas
                 if max(e["bbox"][2] - e["bbox"][0], e["bbox"][3] - e["bbox"][1]) >= 25 * m_por_px]
    for e in etiquetas:
        e["centro"] = [round(e["centro"][0]), round(e["centro"][1])]
        e["bbox"] = [round(v) for v in e["bbox"]]
        e["area"] = round(e["area"])
    return {"lineas": [[list(p) for p in l] for l in lineas],
            "provincias": [[list(p) for p in l] for l in provincias],
            "rellenos": [[list(p) for p in r] for r in rellenos],
            "etiquetas": etiquetas}


def municipios_fondo(res, o):
    pts = [(x, y) for r in res if r.get("ok") for pol in r["polis"] for x, y in pol]
    if not pts:
        return {"lineas": [], "provincias": [], "etiquetas": [], "rellenos": []}
    x0, x1, y0, y1 = min(p[0] for p in pts), max(p[0] for p in pts), min(p[1] for p in pts), max(p[1] for p in pts)
    mg = 0.3 * max(x1 - x0, y1 - y0, 500.0)
    escala = min(860 / (x1 - x0 + 2 * mg), 340 / (y1 - y0 + 2 * mg))
    return fondo_ventana((x0 + x1) / 2, (y0 + y1) / 2, escala, o)


@app.post("/api/fondo")
def post_fondo():
    b = request.get_json(force=True) or {}
    v = b.get("vista") or {}
    try:
        cx, cy, escala = float(v["cx"]), float(v["cy"]), float(v["escala"])
    except (KeyError, TypeError, ValueError):
        raise ErrorUsuario("Vista no válida.")
    if not escala > 1e-7:
        raise ErrorUsuario("Vista no válida: escala demasiado pequeña.")
    if estado["datos"] is None:
        raise ErrorUsuario("Primero conecte con los datos.")
    return jsonify(fondo=fondo_ventana(cx, cy, escala, normaliza_opciones(b.get("opciones"))))


@app.post("/api/previsualizar")
def previsualizar():
    b = request.get_json(force=True)
    p_saneado = sanea_proyecto(b)
    guardar_proyecto(p_saneado)
    o = p_saneado["opciones"]
    gen = generador(o)
    pz.recoger_avisos()
    filas = p_saneado["asignaciones"]
    res, rangos, items_geom = [], [], []
    for a in filas:
        try:
            via, pk0, pk1 = parsea_fila(a)
            res_aud = gen.poligonos_con_auditoria(via, pk0, pk1)
            polis = res_aud[0]
            aviso_trazado = res_aud[1] if len(res_aud) > 1 else None
            pts_desvio = res_aud[2] if len(res_aud) > 2 else []
            res.append({"ok": True, "n_poligonos": len(polis),
                        "n_vertices": 0,
                        "polis": [],
                        "huecos": False,
                        "fuera_trazado": bool(aviso_trazado),
                        "pts_desvio": pts_desvio,
                        "aviso": aviso_trazado})
            rangos.append((via, pk0, pk1))
            items_geom.append({"via": via, "value": a.get("value", ""), "polis": polis})
        except ValueError as e:
            res.append({"ok": False, "error": str(e), "polis": [], "fuera_trazado": False, "pts_desvio": [], "aviso": None})
            rangos.append(None)
            items_geom.append({"via": "", "value": "", "polis": []})

    polis_sin_solape = pz.descontar_solapes_geometricos(items_geom)
    for k, polis_k in enumerate(polis_sin_solape):
        if res[k].get("ok"):
            res[k]["n_poligonos"] = len(polis_k)
            res[k]["n_vertices"] = sum(len(pz.anillo_cm(p)) for p in polis_k)
            res[k]["polis"] = [
                [[round(x, 1), round(y, 1)] for x, y in p.simplify(2.5).exterior.coords]
                for p in polis_k
            ]

    for i, j, msg in pz.solapes(rangos):
        pz.aviso(msg)
        for k in (i, j):
            if res[k].get("ok"):
                res[k]["aviso"] = msg if not res[k].get("aviso") else f"{res[k]['aviso']} | {msg}"

    guardar_proyecto(p_saneado)
    return jsonify(proyecto=p_saneado, resultados=res, fondo=municipios_fondo(res, o), avisos=pz.recoger_avisos())


def construir(p):
    o = p["opciones"]
    gen = generador(o)
    asig = p["asignaciones"]
    if not asig:
        raise ErrorUsuario("No hay asignaciones que exportar.")
    rangos = []
    for i, a in enumerate(asig, start=1):
        try:
            rangos.append(parsea_fila(a))
        except ValueError as e:
            raise ErrorUsuario(f"Asignación {i} (ESZ {a.get('esz')}): {e}")
        if not a["value"]:
            raise ErrorUsuario(f"Asignación {i} (ESZ {a.get('esz')}): falta el mantenedor.")

    pz.recoger_avisos()
    for _, _, msg in pz.solapes(rangos):
        pz.aviso(msg)

    filas = sorted(((a["esz"], r[0], r[1], r[2], a["value"]) for a, r in zip(asig, rangos)),
                   key=lambda f: f[0])
    return pz.generar(filas, gen)


@app.get("/api/exportar/<tipo>")
def exportar(tipo):
    if tipo not in ("ply", "csv"):
        raise ErrorUsuario("Tipo de exportación desconocido.")
    p = cargar_proyecto()
    bloques, tabla = construir(p)
    if tipo == "csv":
        datos, mime, nombre = pz.tabla_csv_bytes(tabla), "text/csv", "zonas_asignaciones.csv"
    else:
        o = p["opciones"]
        try:
            datos = pz.escribir_ply(None, bloques, o["clase"], o["enc_ply"],
                                    "\r\n" if o["saltos"] == "crlf" else "\n")
        except UnicodeEncodeError:
            raise ErrorUsuario(f"Hay caracteres que no existen en la codificación {o['enc_ply']}. "
                               "Revise los nombres de los mantenedores.")
        mime, nombre = "text/plain", "zonas.ply"
    return Response(datos, mimetype=mime,
                    headers={"Content-Disposition": f'attachment; filename="{nombre}"'})


def main():
    url = "http://127.0.0.1:5000"
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    print(f"Interfaz en {url}   (Ctrl+C para salir)")
    app.run(host="127.0.0.1", port=5000, debug=False)


if __name__ == "__main__":
    main()
