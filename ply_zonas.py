#!/usr/bin/env python3
"""
ply_zonas.py - Núcleo y línea de comandos para generar archivos PLY (formato
GeoMedia) con zonas de influencia de carreteras entre dos puntos kilométricos.
"""
import argparse
import bisect
import csv
import heapq
import math
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, MultiLineString, Point, Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import linemerge, split, substring, unary_union
from shapely.strtree import STRtree

import fuentes
from fuentes import a_float, norm

_avisos = []


def aviso(msg):
    _avisos.append(msg)
    print(f"AVISO: {msg}", file=sys.stderr)


def recoger_avisos():
    """Devuelve los avisos acumulados (sin repetir) y los borra."""
    res = list(dict.fromkeys(_avisos))
    _avisos.clear()
    return res


# ---------------------------------------------------- motor geométrico
def _centroide_pk(pts):
    """Si un PK tiene varios puntos (ej. doble calzada), devuelve su punto medio."""
    if len(pts) == 1:
        return Point(pts[0])
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return Point(sum(xs) / len(xs), sum(ys) / len(ys))


def _fusionar_tramos_via(tramos_coords):
    """Convierte las listas de coordenadas de CARRETERAS_BCN en LineStrings fusionados."""
    lineas = [LineString(t) for t in tramos_coords if len(t) >= 2]
    if len(lineas) > 1:
        try:
            m = linemerge(lineas)
            lineas = [m] if isinstance(m, LineString) else list(m.geoms)
        except Exception:
            pass
    return lineas


def longitud_total_km(tramos_coords):
    """Calcula la longitud total en km de los tramos de una vía cuando no tiene PKs."""
    lineas = _fusionar_tramos_via(tramos_coords)
    if not lineas:
        return 0.0
    total_m = sum(l.length for l in lineas)
    return round(max(total_m / 1000.0, 0.001), 3)


def _buscar_puente_en_red_vial(pt_salida, pt_entrada, red_global_tree, red_global_lineas):
    """
    Cuando una vía tiene una discontinuidad en CARRETERAS_BCN porque comparte variante
    con otra denominación (ej. N-340 con A-7), busca el camino asfaltado real en la red global
    que conecta pt_salida con pt_entrada.
    """
    dist_recta = math.dist(pt_salida, pt_entrada)
    if red_global_tree is None or dist_recta < 10.0 or dist_recta > 25000.0:
        return None

    margen = max(2500.0, dist_recta * 0.6)
    minx = min(pt_salida[0], pt_entrada[0]) - margen
    miny = min(pt_salida[1], pt_entrada[1]) - margen
    maxx = max(pt_salida[0], pt_entrada[0]) + margen
    maxy = max(pt_salida[1], pt_entrada[1]) + margen

    idxs = red_global_tree.query(box(minx, miny, maxx, maxy))
    if len(idxs) == 0:
        return None

    cand = [red_global_lineas[i] for i in idxs]
    if len(cand) > 1:
        try:
            m = linemerge(cand)
            cand = [m] if isinstance(m, LineString) else list(m.geoms)
        except Exception:
            pass

    p_sal_geom = Point(pt_salida)
    p_ent_geom = Point(pt_entrada)

    mejor_l = None
    mejor_len = math.inf
    for l in cand:
        if l.distance(p_sal_geom) <= 45.0 and l.distance(p_ent_geom) <= 45.0:
            d1, d2 = l.project(p_sal_geom), l.project(p_ent_geom)
            tramo_l = abs(d2 - d1)
            if dist_recta * 0.75 <= tramo_l <= dist_recta * 3.0 and tramo_l < mejor_len:
                rec = substring(l, min(d1, d2), max(d1, d2))
                coords_rec = list(rec.coords)
                if d1 > d2:
                    coords_rec.reverse()
                mejor_l = coords_rec
                mejor_len = tramo_l

    if mejor_l is not None:
        return mejor_l

    SNAP_RED = 25.0
    m = len(cand)
    ext = [(l.coords[0], l.coords[-1]) for l in cand]

    cola = []
    dist = {}
    prev = {}

    for i, l in enumerate(cand):
        d_sal = l.distance(p_sal_geom)
        if d_sal <= 45.0:
            proy_s = l.project(p_sal_geom)
            coste_fw = (l.length - proy_s) + d_sal * 5.0
            estado_fw = (i, True)
            if coste_fw < dist.get(estado_fw, math.inf):
                dist[estado_fw] = coste_fw
                prev[estado_fw] = ("INICIO", proy_s)
                heapq.heappush(cola, (coste_fw, estado_fw))

            coste_bw = proy_s + d_sal * 5.0
            estado_bw = (i, False)
            if coste_bw < dist.get(estado_bw, math.inf):
                dist[estado_bw] = coste_bw
                prev[estado_bw] = ("INICIO", proy_s)
                heapq.heappush(cola, (coste_bw, estado_bw))

    if not cola:
        return None

    mejor_destino = None
    mejor_coste_dest = dist_recta * 3.2

    while cola:
        d, u = heapq.heappop(cola)
        if d > mejor_coste_dest or d > dist.get(u, math.inf):
            continue

        i_u, fw_u = u
        l_u = cand[i_u]
        d_ent = l_u.distance(p_ent_geom)
        if d_ent <= 45.0:
            proy_e = l_u.project(p_ent_geom)
            if prev[u][0] == "INICIO":
                proy_s = prev[u][1]
                valido = (proy_e >= proy_s) if fw_u else (proy_e <= proy_s)
                if valido:
                    coste_real = abs(proy_e - proy_s) + d_ent * 5.0
                    if coste_real < mejor_coste_dest:
                        mejor_coste_dest = coste_real
                        mejor_destino = (u, proy_e)
            else:
                resto = (l_u.length - proy_e) if fw_u else proy_e
                coste_real = d - resto + d_ent * 5.0
                if coste_real < mejor_coste_dest:
                    mejor_coste_dest = coste_real
                    mejor_destino = (u, proy_e)

        pt_fin_u = ext[i_u][1] if fw_u else ext[i_u][0]
        for j in range(m):
            if j == i_u:
                continue
            for fw_v in (True, False):
                pt_ini_v = ext[j][0] if fw_v else ext[j][1]
                salto = math.dist(pt_fin_u, pt_ini_v)
                if salto <= SNAP_RED:
                    nd = d + cand[j].length + salto * 4.0
                    v = (j, fw_v)
                    if nd < dist.get(v, math.inf) and nd < mejor_coste_dest:
                        dist[v] = nd
                        prev[v] = ("PASO", u)
                        heapq.heappush(cola, (nd, v))

    if mejor_destino is None:
        return None

    u_fin, proy_fin = mejor_destino
    pasos = []
    curr = u_fin
    while True:
        tipo_p, dato_p = prev[curr]
        if tipo_p == "INICIO":
            pasos.append((curr[0], curr[1], dato_p, proy_fin if curr == u_fin else None))
            break
        else:
            pasos.append((curr[0], curr[1], None, proy_fin if curr == u_fin else None))
            curr = dato_p
    pasos.reverse()

    coords_puente = []
    for idx_c, fw_c, p_ini_c, p_fin_c in pasos:
        l_c = cand[idx_c]
        d_a = p_ini_c if p_ini_c is not None else (0.0 if fw_c else l_c.length)
        d_b = p_fin_c if p_fin_c is not None else (l_c.length if fw_c else 0.0)
        if abs(d_b - d_a) < 0.05:
            continue
        sub = substring(l_c, min(d_a, d_b), max(d_a, d_b))
        c_sub = list(sub.coords)
        if d_a > d_b:
            c_sub.reverse()
        for pt in c_sub:
            if not coords_puente or math.dist(coords_puente[-1], pt) > 0.05:
                coords_puente.append(pt)

    return coords_puente if len(coords_puente) >= 2 else None


def _vector_unitario(p_a, p_b):
    dx, dy = p_b[0] - p_a[0], p_b[1] - p_a[1]
    n = math.hypot(dx, dy)
    return (dx / n, dy / n) if n > 1e-6 else (0.0, 0.0)


def _construir_cadena_principal(lineas, pks_ordenados, red_global_tree=None, red_global_lineas=None, pk_rango_objetivo=None):
    """
    Construye el eje longitudinal de referencia de la vía siguiendo los PKs.
    Si la vía es físicamente discontinua (como la N-323A dividida en 4 sectores), devuelve
    un MultiLineString con cada bloque continuo por separado sin trazar rectas en el aire.
    """
    if not lineas:
        return None

    pks_sector = pks_ordenados
    if pk_rango_objetivo is not None and pks_ordenados and len(pks_ordenados) >= 2:
        p0_obj, p1_obj = min(pk_rango_objetivo), max(pk_rango_objetivo)
        cercanos = [p for p in pks_ordenados if (p0_obj - 60.0) <= p[0] <= (p1_obj + 60.0)]
        if len(cercanos) >= 2:
            pks_sector = cercanos

    if len(lineas) == 1:
        l = lineas[0]
        if pks_sector and len(pks_sector) >= 2:
            p0 = _centroide_pk(pks_sector[0][1])
            p1 = _centroide_pk(pks_sector[-1][1])
            if l.project(p0) > l.project(p1):
                l = LineString(list(l.coords)[::-1])
        return l

    n = len(lineas)
    extremos = [(l.coords[0], l.coords[-1]) for l in lineas]

    def _camino_dijkstra(idx_ini, dir_ini, idx_fin, dir_fin):
        origen = (idx_ini, dir_ini)
        destino = (idx_fin, dir_fin)
        if origen == destino:
            return [origen], 0.0

        dist = {origen: 0.0}
        prev = {}
        cola = [(0.0, origen)]

        while cola:
            d, u = heapq.heappop(cola)
            if u == destino:
                break
            if d > dist.get(u, math.inf):
                continue
            i_u, fw_u = u
            pt_salida = extremos[i_u][1] if fw_u else extremos[i_u][0]

            for j in range(n):
                if j == i_u:
                    continue
                for fw_v in (True, False):
                    pt_entrada = extremos[j][0] if fw_v else extremos[j][1]
                    salto = math.dist(pt_salida, pt_entrada)
                    if salto > 15000.0:
                        continue
                    coste_salto = salto if salto <= 5.0 else (salto * 25.0 + 50.0)
                    nd = d + lineas[j].length + coste_salto
                    v = (j, fw_v)
                    if nd < dist.get(v, math.inf):
                        dist[v] = nd
                        prev[v] = u
                        heapq.heappush(cola, (nd, v))

        if destino not in dist:
            return None, math.inf
        camino = []
        curr = destino
        while curr != origen:
            camino.append(curr)
            curr = prev[curr]
        camino.append(origen)
        camino.reverse()
        return camino, dist[destino]

    anclas = []
    if pks_sector:
        ultimo_idx = None
        for pk_val, pts in pks_sector:
            candidatos_pk = []
            for p_xy in pts:
                pt = Point(p_xy)
                for i in range(n):
                    d_i = lineas[i].distance(pt)
                    if d_i <= 500.0:
                        penalizacion = 0.0 if (ultimo_idx is not None and i == ultimo_idx) else 25.0
                        candidatos_pk.append((d_i + penalizacion, d_i, pt, i))
            if candidatos_pk:
                candidatos_pk.sort(key=lambda x: x[0])
                _, _, pt_elegido, idx_mejor = candidatos_pk[0]
                ultimo_idx = idx_mejor
                proy = lineas[idx_mejor].project(pt_elegido)
                anclas.append((pk_val, pt_elegido, idx_mejor, proy))

    if len(anclas) < 2:
        usados = {0}
        bloques = [list(lineas[0].coords)]
        while len(usados) < n:
            pt_actual = bloques[-1][-1]
            mejor_j, mejor_inv, mejor_d = None, False, math.inf
            for j in range(n):
                if j in usados:
                    continue
                d_ini = math.dist(pt_actual, extremos[j][0])
                d_fin = math.dist(pt_actual, extremos[j][1])
                if d_ini < mejor_d:
                    mejor_d, mejor_j, mejor_inv = d_ini, j, False
                if d_fin < mejor_d:
                    mejor_d, mejor_j, mejor_inv = d_fin, j, True
            if mejor_j is None or mejor_d > 5000.0:
                break
            usados.add(mejor_j)
            c_j = list(lineas[mejor_j].coords)
            if mejor_inv:
                c_j.reverse()
            if mejor_d <= 250.0:
                bloques[-1].extend(c_j)
            else:
                bloques.append(c_j)
        geoms = [LineString(b) for b in bloques if len(b) >= 2]
        return geoms[0] if len(geoms) == 1 else MultiLineString(geoms)

    tramos_ancla = []
    for k, (pk_val, pt, idx_t, proy) in enumerate(anclas):
        if not tramos_ancla or tramos_ancla[-1][0] != idx_t:
            mismo_tramo = [a for a in anclas[k + 1:] if a[2] == idx_t]
            if mismo_tramo:
                fw = mismo_tramo[0][3] >= proy
            elif k + 1 < len(anclas):
                pt_sig = anclas[k + 1][1]
                d_fin = math.dist(extremos[idx_t][1], (pt_sig.x, pt_sig.y))
                d_ini = math.dist(extremos[idx_t][0], (pt_sig.x, pt_sig.y))
                fw = d_fin <= d_ini
            elif tramos_ancla:
                prev_idx, prev_fw = tramos_ancla[-1]
                pt_prev_sal = extremos[prev_idx][1] if prev_fw else extremos[prev_idx][0]
                fw = math.dist(pt_prev_sal, extremos[idx_t][0]) <= math.dist(pt_prev_sal, extremos[idx_t][1])
            else:
                fw = True
            tramos_ancla.append((idx_t, fw))

    secuencia_tramos = [tramos_ancla[0]]
    for (idx_a, fw_a), (idx_b, fw_b) in zip(tramos_ancla[:-1], tramos_ancla[1:]):
        if idx_a == idx_b:
            continue
        ya_recorrido = unary_union([lineas[i] for i, _ in secuencia_tramos])
        if lineas[idx_b].distance(ya_recorrido) < 5.0 and lineas[idx_b].intersection(ya_recorrido.buffer(800.0)).length > lineas[idx_b].length * 0.7:
            continue

        camino, _ = _camino_dijkstra(secuencia_tramos[-1][0], secuencia_tramos[-1][1], idx_b, fw_b)
        if camino and len(camino) > 1:
            for paso in camino[1:]:
                if paso[0] != secuencia_tramos[-1][0]:
                    secuencia_tramos.append(paso)
        else:
            if idx_b != secuencia_tramos[-1][0]:
                secuencia_tramos.append((idx_b, fw_b))

    MAX_SALTO_EXT = 4500.0
    usados_idx = {idx for idx, _ in secuencia_tramos}

    while True:
        idx_ult, fw_ult = secuencia_tramos[-1]
        c_ult = list(lineas[idx_ult].coords)
        if not fw_ult:
            c_ult.reverse()
        pt_sal = c_ult[-1]
        p_ant = c_ult[max(0, len(c_ult) - 5)]
        vec_avance = _vector_unitario(p_ant, pt_sal)

        mejor_j, mejor_fw, mejor_d = None, True, MAX_SALTO_EXT
        for j in range(n):
            if j in usados_idx:
                continue
            for fw_cand, p_ent in (
                (True, extremos[j][0]),
                (False, extremos[j][1]),
            ):
                d_salto = math.dist(pt_sal, p_ent)
                if d_salto < mejor_d:
                    if d_salto > 25.0:
                        vec_salto = _vector_unitario(pt_sal, p_ent)
                        if (vec_avance[0] * vec_salto[0] + vec_avance[1] * vec_salto[1]) < -0.2:
                            continue
                    mejor_d, mejor_j, mejor_fw = d_salto, j, fw_cand
        if mejor_j is not None:
            secuencia_tramos.append((mejor_j, mejor_fw))
            usados_idx.add(mejor_j)
        else:
            break

    while True:
        idx_pri, fw_pri = secuencia_tramos[0]
        c_pri = list(lineas[idx_pri].coords)
        if not fw_pri:
            c_pri.reverse()
        pt_ent = c_pri[0]
        p_post = c_pri[min(4, len(c_pri) - 1)]
        vec_retro = _vector_unitario(p_post, pt_ent)

        mejor_j, mejor_fw, mejor_d = None, True, MAX_SALTO_EXT
        for j in range(n):
            if j in usados_idx:
                continue
            for fw_cand, p_sal_cand in (
                (True, extremos[j][1]),
                (False, extremos[j][0]),
            ):
                d_salto = math.dist(pt_ent, p_sal_cand)
                if d_salto < mejor_d:
                    if d_salto > 25.0:
                        vec_salto = _vector_unitario(pt_ent, p_sal_cand)
                        if (vec_retro[0] * vec_salto[0] + vec_retro[1] * vec_salto[1]) < -0.2:
                            continue
                    mejor_d, mejor_j, mejor_fw = d_salto, j, fw_cand
        if mejor_j is not None:
            secuencia_tramos.insert(0, (mejor_j, mejor_fw))
            usados_idx.add(mejor_j)
        else:
            break

    bloques_coords = []
    actual_coords = []

    for idx_t, fw in secuencia_tramos:
        c_seg = list(lineas[idx_t].coords)
        if not fw:
            c_seg.reverse()

        if actual_coords:
            salto = math.dist(actual_coords[-1], c_seg[0])
            if salto > 15.0:
                puente = None
                if red_global_tree is not None and salto <= 12000.0:
                    puente = _buscar_puente_en_red_vial(
                        actual_coords[-1], c_seg[0], red_global_tree, red_global_lineas
                    )
                if puente:
                    for pt_p in puente:
                        if math.dist(actual_coords[-1], pt_p) > 0.05:
                            actual_coords.append(pt_p)
                elif salto > 250.0:
                    if len(actual_coords) >= 2:
                        bloques_coords.append(LineString(actual_coords))
                    actual_coords = []

        for pt in c_seg:
            if not actual_coords or math.dist(actual_coords[-1], pt) > 0.05:
                actual_coords.append(pt)

    if len(actual_coords) >= 2:
        bloques_coords.append(LineString(actual_coords))

    if not bloques_coords:
        return lineas[0]
    return bloques_coords[0] if len(bloques_coords) == 1 else MultiLineString(bloques_coords)


def _calibrar_pks_sobre_eje(eje, pks_ordenados):
    """
    Proyecta los PKs existentes sobre un LineString (o sobre el componente continuo correspondiente)
    garantizando que la distancia a lo largo del eje sea estrictamente creciente.
    """
    if isinstance(eje, MultiLineString):
        coords_cat = [pt for g in eje.geoms for pt in g.coords]
        eje = LineString(coords_cat) if len(coords_cat) >= 2 else eje.geoms[0]

    proyecciones = []
    for pk_val, pts in pks_ordenados:
        p = min((Point(q) for q in pts), key=lambda pt: eje.distance(pt))
        if eje.distance(p) <= 600.0:
            proyecciones.append((float(pk_val), eje.project(p)))

    if len(proyecciones) < 2:
        return [(0.0, 0.0), (round(eje.length / 1000.0, 3), eje.length)]

    valores_s = [s for _, s in proyecciones]
    cola, idx_cola, padre = [], [], [-1] * len(valores_s)
    for i, v in enumerate(valores_s):
        j = bisect.bisect_left(cola, v)
        if j == len(cola):
            cola.append(v)
            idx_cola.append(i)
        else:
            cola[j] = v
            idx_cola[j] = i
        padre[i] = idx_cola[j - 1] if j else -1

    keep, curr = [], idx_cola[-1] if idx_cola else -1
    while curr >= 0:
        keep.append(curr)
        curr = padre[curr]
    keep.reverse()

    tabla = [proyecciones[i] for i in keep]
    return tabla if len(tabla) >= 2 else proyecciones


def rango_efectivo_via(tramos_coords, pks_ordenados):
    """
    Devuelve (pk_min_efectivo, pk_max_efectivo) incluyendo el trazado real de CARRETERAS_BCN
    que continúa antes del primer PK o después del último PK para poder extrapolarlo.
    """
    lineas = _fusionar_tramos_via(tramos_coords)
    if not lineas:
        return 0.0, 0.0
    if not pks_ordenados or len(pks_ordenados) < 2:
        return 0.0, longitud_total_km(tramos_coords)

    eje = _construir_cadena_principal(lineas, pks_ordenados)
    partes = [eje] if isinstance(eje, LineString) else list(eje.geoms)

    p0 = _centroide_pk(pks_ordenados[0][1])
    p_ult = _centroide_pk(pks_ordenados[-1][1])
    comp_ini = min(partes, key=lambda c: c.distance(p0))
    comp_fin = min(partes, key=lambda c: c.distance(p_ult))

    s_min = comp_ini.project(p0)
    s_max = comp_fin.project(p_ult)

    pk_min = float(pks_ordenados[0][0])
    pk_max = float(pks_ordenados[-1][0])

    len_total_vias = sum(l.length for l in lineas)
    len_eje = sum(c.length for c in partes)
    extra_no_cosido = max(0.0, len_total_vias - len_eje)

    metros_antes = max(0.0, s_min) + extra_no_cosido
    metros_despues = max(0.0, comp_fin.length - s_max) + extra_no_cosido

    pk_min_ef = max(0.0, pk_min - (metros_antes / 1000.0))
    pk_max_ef = pk_max + (metros_despues / 1000.0)
    return round(pk_min_ef, 3), round(pk_max_ef, 3)


def _distancia_en_eje_para_pk(eje, tabla_calibrada, target_pk):
    """
    Convierte un PK en su distancia exacta en metros a lo largo del trazado real `eje`.
    """
    pks = [t[0] for t in tabla_calibrada]
    dists = [t[1] for t in tabla_calibrada]

    if target_pk <= pks[0]:
        delta_m = (pks[0] - target_pk) * 1000.0
        return max(0.0, dists[0] - delta_m)

    if target_pk >= pks[-1]:
        delta_m = (target_pk - pks[-1]) * 1000.0
        return min(eje.length, dists[-1] + delta_m)

    idx = bisect.bisect_left(pks, target_pk)
    if idx < len(pks) and abs(pks[idx] - target_pk) < 1e-6:
        return dists[idx]

    pk_low, pk_high = pks[idx - 1], pks[idx]
    d_low, d_high = dists[idx - 1], dists[idx]
    frac = (target_pk - pk_low) / (pk_high - pk_low)
    return d_low + frac * (d_high - d_low)


def _extraer_tramos_ambos_sentidos(lineas, eje_principal, tabla_calibrada, pks_ordenados, p_ini, p_fin):
    """
    Extrae el subtramo (o subtramos si la vía es discontinua, ej. N-323A) entre p_ini y p_fin
    sin trazar rectas en el aire entre bloques desconectados, y añade las calzadas paralelas.
    """
    componentes = [eje_principal] if isinstance(eje_principal, LineString) else list(eje_principal.geoms)
    resultado = []

    for comp in componentes:
        if len(componentes) == 1 and tabla_calibrada is not None:
            tabla_c = tabla_calibrada
            pk_min_c, s_min_c = tabla_c[0]
            pk_max_c, s_max_c = tabla_c[-1]
            min_ef_c = pk_min_c - (s_min_c / 1000.0)
            max_ef_c = pk_max_c + ((comp.length - s_max_c) / 1000.0)
            if p_fin <= min_ef_c + 1e-3 or p_ini >= max_ef_c - 1e-3:
                continue
            d1 = _distancia_en_eje_para_pk(comp, tabla_c, max(p_ini, min_ef_c))
            d2 = _distancia_en_eje_para_pk(comp, tabla_c, min(p_fin, max_ef_c))
        else:
            pks_comp = [
                (pk_val, pts) for pk_val, pts in pks_ordenados
                if min(comp.distance(Point(q)) for q in pts) <= 500.0
            ]
            if len(pks_comp) >= 2:
                tabla_c = _calibrar_pks_sobre_eje(comp, pks_comp)
                pk_min_c, s_min_c = tabla_c[0]
                pk_max_c, s_max_c = tabla_c[-1]
                min_ef_c = pk_min_c - (s_min_c / 1000.0)
                max_ef_c = pk_max_c + ((comp.length - s_max_c) / 1000.0)
                if p_fin <= min_ef_c + 1e-3 or p_ini >= max_ef_c - 1e-3:
                    continue
                d1 = _distancia_en_eje_para_pk(comp, tabla_c, max(p_ini, min_ef_c))
                d2 = _distancia_en_eje_para_pk(comp, tabla_c, min(p_fin, max_ef_c))
            elif len(pks_comp) == 1:
                pk_unico, pts_u = pks_comp[0]
                s_u = comp.project(_centroide_pk(pts_u))
                min_ef_c = pk_unico - (s_u / 1000.0)
                max_ef_c = pk_unico + ((comp.length - s_u) / 1000.0)
                if p_fin <= min_ef_c + 1e-3 or p_ini >= max_ef_c - 1e-3:
                    continue
                d1 = max(0.0, s_u + (max(p_ini, min_ef_c) - pk_unico) * 1000.0)
                d2 = min(comp.length, s_u + (min(p_fin, max_ef_c) - pk_unico) * 1000.0)
            else:
                continue

        s_min, s_max = min(d1, d2), max(d1, d2)
        if s_max - s_min < 0.1:
            continue

        tramo_ppal = substring(comp, s_min, s_max)
        resultado.append(tramo_ppal)

        pt_ini_ppal = comp.interpolate(s_min)
        pt_fin_ppal = comp.interpolate(s_max)
        corredor_ppal = tramo_ppal.buffer(2.5)
        zona_busqueda = tramo_ppal.buffer(1500.0)

        for l in lineas:
            if not l.intersects(zona_busqueda):
                continue
            dif = l.difference(corredor_ppal)
            if dif.is_empty:
                continue
            partes_dif = [dif] if isinstance(dif, LineString) else list(getattr(dif, "geoms", []))
            for sub in partes_dif:
                if sub.length < 15.0 or sub.distance(comp) > 1500.0:
                    continue
                pr_0 = comp.project(Point(sub.coords[0]))
                pr_1 = comp.project(Point(sub.coords[-1]))
                a_sub, b_sub = min(pr_0, pr_1), max(pr_0, pr_1)

                if b_sub <= s_min + 5.0 or a_sub >= s_max - 5.0:
                    continue
                if pr_0 > pr_1:
                    sub = LineString(list(sub.coords)[::-1])

                corte_ini = sub.project(pt_ini_ppal) if a_sub < s_min - 5.0 else 0.0
                corte_fin = sub.project(pt_fin_ppal) if b_sub > s_max + 5.0 else sub.length

                if corte_fin - corte_ini > 5.0:
                    rec_sec = substring(sub, max(0.0, corte_ini), min(sub.length, corte_fin))
                    if rec_sec.length > 5.0:
                        resultado.append(rec_sec)

    return resultado


def _cortar_poligono_con_huecos(poligono):
    """
    Si dos sentidos separados de una carretera se bifurcan y más adelante se vuelven a unir,
    la unión de sus buffers genera un polígono con hueco interior (`interiors`).
    Cortamos el polígono a través del hueco para que cada sentido quede como un polígono tubular
    independiente sin rellenar jamás la montaña/mediana central.
    """
    if not poligono.interiors:
        return [poligono]

    pendientes = [poligono]
    finales = []

    while pendientes:
        actual = pendientes.pop(0)
        if not actual.interiors:
            finales.append(actual)
            continue

        hueco = actual.interiors[0]
        pt_centro = Polygon(hueco).representative_point()
        minx, miny, maxx, maxy = actual.bounds
        hx0, hy0, hx1, hy1 = Polygon(hueco).bounds
        if (hy1 - hy0) >= (hx1 - hx0):
            linea_corte = LineString([(minx - 10.0, pt_centro.y), (maxx + 10.0, pt_centro.y)])
        else:
            linea_corte = LineString([(pt_centro.x, miny - 10.0), (pt_centro.x, maxy + 10.0)])

        try:
            partes = split(actual, linea_corte)
            nuevos = [g for g in getattr(partes, "geoms", [partes]) if g.geom_type == "Polygon" and g.area > 0.5]
            if len(nuevos) > 1:
                pendientes.extend(nuevos)
            else:
                finales.append(actual)
        except Exception:
            finales.append(actual)

    return finales


def _prioridad_jerarquia_via(via, area_poligonos):
    """
    Determina qué vía tiene prioridad en las bifurcaciones/enlaces (ej. A-4 tiene prioridad
    sobre N-IVa) para que la vía secundaria recorte su cola sin pisar la vía principal.
    """
    v = norm(via)
    es_ramal = any(p in v for p in ("ENLACE", "RAMAL", "VARIANTE", "ACCESO")) or (
        len(v) >= 4 and v[-1].isalpha()
    )
    if not es_ramal and (v.startswith("AP-") or (v.startswith("A-") and len(v) <= 5)):
        jerarquia = 0
    elif not es_ramal and v.startswith("N-"):
        jerarquia = 1
    else:
        jerarquia = 2
    return (jerarquia, -area_poligonos)


def descontar_solapes_geometricos(lista_items):
    """
    Recorta geométricamente las colas de bifurcación donde una carretera secundaria/ramal
    se despega de una vía principal, usando un índice espacial STRtree para máxima velocidad.
    """
    if len(lista_items) <= 1:
        return [it["polis"] for it in lista_items]

    geoms_item = []
    for idx, it in enumerate(lista_items):
        polis = it.get("polis") or []
        u = unary_union(polis) if len(polis) > 1 else (polis[0] if polis else Polygon())
        geoms_item.append({
            "idx": idx,
            "via": norm(it.get("via", "")),
            "val": (it.get("value") or "").strip().upper(),
            "geom": u,
            "area": u.area if not u.is_empty else 0.0
        })

    geoms_validas = [it["geom"] for it in geoms_item if not it["geom"].is_empty]
    if not geoms_validas:
        return [[] for _ in lista_items]

    arbol = STRtree(geoms_validas)
    items_validos = [it for it in geoms_item if not it["geom"].is_empty]
    prioridades = [
        (*_prioridad_jerarquia_via(it["via"], it["area"]), pos)
        for pos, it in enumerate(items_validos)
    ]

    resultado_por_idx = [[] for _ in lista_items]

    for pos, item in enumerate(items_validos):
        g = item["geom"]
        mi_prio = prioridades[pos]

        idxs_vecinos = arbol.query(g)
        superiores = [
            items_validos[j]["geom"]
            for j in idxs_vecinos
            if j != pos and prioridades[j] < mi_prio
        ]

        if superiores:
            try:
                union_sup = unary_union(superiores) if len(superiores) > 1 else superiores[0]
                if g.intersects(union_sup):
                    recortada = g.difference(union_sup)
                    partes = [
                        p for p in getattr(recortada, "geoms", [recortada])
                        if p.geom_type == "Polygon" and p.area > 5.0
                    ]
                    g = unary_union(partes) if partes else Polygon()
            except Exception:
                pass

        if not g.is_empty:
            brutos = [g] if g.geom_type == "Polygon" else [
                p for p in getattr(g, "geoms", []) if p.geom_type == "Polygon" and p.area > 5.0
            ]
            resultado_por_idx[item["idx"]] = [
                sub for p in brutos for sub in _cortar_poligono_con_huecos(p)
            ]

    return resultado_por_idx


# ------------------------------------------------------------- geometría
class Generador:
    """Construye las zonas de influencia de una vía siguiendo siempre el trazado de CARRETERAS_BCN."""

    def __init__(self, datos, ancho=20.0, quad_segs=18):
        self.datos, self.ancho, self.quad_segs = datos, ancho, quad_segs
        self._cache_via = {}
        self._red_tree = None
        self._red_lineas = None

    def _obtener_red_global(self):
        """Construye un índice espacial STRtree con todos los tramos de CARRETERAS_BCN."""
        if self._red_tree is None:
            todas = []
            for lista_t in self.datos.tramos.values():
                for coords in lista_t:
                    if len(coords) >= 2:
                        todas.append(LineString(coords))
            self._red_lineas = todas
            self._red_tree = STRtree(todas) if todas else None
        return self._red_tree, self._red_lineas

    def _obtener_eje_via(self, via, pk_desde=None, pk_hasta=None):
        k = norm(via)
        tramos = self.datos.lista_tramos(via)
        pks = self.datos.pks_multi(via)
        if not tramos:
            raise ValueError(f"{via}: no hay tramos de esa vía en la capa de ejes.")

        sector = None
        if pk_desde is not None and pk_hasta is not None and pks and (pks[-1][0] - pks[0][0] > 250):
            sector = (int(min(pk_desde, pk_hasta) // 100), int(max(pk_desde, pk_hasta) // 100))

        clave_cache = (k, sector)
        if clave_cache not in self._cache_via:
            lineas = _fusionar_tramos_via(tramos)
            rtree, rlineas = self._obtener_red_global()
            rango_obj = (pk_desde, pk_hasta) if sector is not None else None
            eje = _construir_cadena_principal(lineas, pks, rtree, rlineas, rango_obj)
            tabla = _calibrar_pks_sobre_eje(eje, pks) if (len(pks) >= 2 and isinstance(eje, LineString)) else None
            self._cache_via[clave_cache] = (lineas, eje, tabla, pks)
        return self._cache_via[clave_cache]

    def _verificar_fidelidad_trazado(self, lista_tramos_eje, via, pk_desde, pk_hasta, tol_m=None):
        """
        Verifica que los tramos generados sigan en todo momento algún tramo real
        de CARRETERAS_BCN. Ignora micro-saltos puntuales en rotondas o enlaces (< 75 m)
        y devuelve (mensaje_aviso, [[x, y, desv_m], ...]) con las coordenadas exactas
        donde se producen las desviaciones reales para señalarlas en el mapa.
        """
        rtree, rlineas = self._obtener_red_global()
        if rtree is None or not lista_tramos_eje:
            return None, []

        if tol_m is None:
            tol_m = max(25.0, self.ancho * 1.25)

        PASO_MUESTREO = 40.0
        metros_fuera = 0.0
        max_desv = 0.0
        max_racha_m = 0.0
        pts_desvio = []

        for geom in lista_tramos_eje:
            partes = [geom] if isinstance(geom, LineString) else list(getattr(geom, "geoms", []))
            for lin in partes:
                L = lin.length
                if L < 1.0:
                    continue
                n_pasos = max(2, int(math.ceil(L / PASO_MUESTREO)))
                paso_real = L / n_pasos
                racha_actual_m = 0.0
                peor_pt_racha = None
                peor_d_racha = 0.0

                for i in range(n_pasos + 1):
                    pt = lin.interpolate(i * paso_real)
                    ventana = max(tol_m * 2.0, 60.0)
                    idxs = rtree.query(box(pt.x - ventana, pt.y - ventana, pt.x + ventana, pt.y + ventana))
                    if len(idxs) > 0:
                        d_min = min(rlineas[idx].distance(pt) for idx in idxs)
                    else:
                        idxs_amplio = rtree.query(box(pt.x - 2000.0, pt.y - 2000.0, pt.x + 2000.0, pt.y + 2000.0))
                        if len(idxs_amplio) > 0:
                            d_min = min(rlineas[idx].distance(pt) for idx in idxs_amplio)
                        else:
                            d_min = 2000.0

                    if d_min > tol_m:
                        racha_actual_m += paso_real
                        metros_fuera += paso_real
                        if d_min > max_desv:
                            max_desv = d_min
                        if racha_actual_m > max_racha_m:
                            max_racha_m = racha_actual_m
                        if d_min > peor_d_racha:
                            peor_d_racha = d_min
                            peor_pt_racha = (round(pt.x, 1), round(pt.y, 1), round(d_min))
                    else:
                        if peor_pt_racha is not None and (racha_actual_m >= 75.0 or peor_d_racha >= 35.0):
                            pts_desvio.append(list(peor_pt_racha))
                        racha_actual_m = 0.0
                        peor_pt_racha = None
                        peor_d_racha = 0.0

                if peor_pt_racha is not None and (racha_actual_m >= 75.0 or peor_d_racha >= 35.0):
                    pts_desvio.append(list(peor_pt_racha))

        if max_desv > tol_m and (max_racha_m >= 75.0 or max_desv >= 35.0) and pts_desvio:
            msg = (
                f"{via} (PK {pk_desde:g}-{pk_hasta:g}): se aparta del trazado de CARRETERAS_BCN "
                f"en ~{round(metros_fuera)} m (desviación máx. {round(max_desv)} m)."
            )
            aviso(msg)
            return msg, pts_desvio
        return None, []

    def poligonos_con_auditoria(self, via, pk_desde, pk_hasta):
        """
        Devuelve (lista_poligonos_sin_huecos, aviso_trazado, pts_desvio) respetando vías
        discontinuas y creando la cobertura individual para cada sentido cuando las calzadas se separan.
        """
        lineas, eje, tabla, pks = self._obtener_eje_via(via, pk_desde, pk_hasta)

        if not pks or len(pks) < 2:
            aviso(f"{via}: sin registros en la tabla de PK; se toma el recorrido completo existente.")
            buffers = [
                l.buffer(self.ancho, cap_style="flat", join_style="round", quad_segs=self.quad_segs)
                for l in lineas
            ]
            z = unary_union(buffers)
            brutos = [z] if z.geom_type == "Polygon" else [p for p in getattr(z, "geoms", []) if p.geom_type == "Polygon"]
            polis = [sub for p in brutos for sub in _cortar_poligono_con_huecos(p)]
            return polis, None, []

        p_ini, p_fin = min(pk_desde, pk_hasta), max(pk_desde, pk_hasta)
        if p_fin - p_ini < 1e-6:
            raise ValueError(f"{via}: el rango PK {pk_desde:g}-{pk_hasta:g} tiene longitud nula.")

        tramos_sentidos = _extraer_tramos_ambos_sentidos(lineas, eje, tabla, pks, p_ini, p_fin)
        if not tramos_sentidos:
            raise ValueError(
                f"{via}: el tramo PK {p_ini:g}-{p_fin:g} cae en un salto discontinuo o fuera del trazado de CARRETERAS_BCN."
            )

        aviso_trazado, pts_desvio = self._verificar_fidelidad_trazado(tramos_sentidos, via, p_ini, p_fin)

        buffers_individuales = [
            t.buffer(self.ancho, cap_style="flat", join_style="round", quad_segs=self.quad_segs)
            for t in tramos_sentidos if not t.is_empty and t.length >= 0.1
        ]
        if not buffers_individuales:
            raise ValueError(f"{via}: el buffer generado está vacío.")

        z = unary_union(buffers_individuales)
        brutos = [z] if z.geom_type == "Polygon" else [p for p in getattr(z, "geoms", []) if p.geom_type == "Polygon"]
        polis = [sub for p in brutos for sub in _cortar_poligono_con_huecos(p)]
        return polis, aviso_trazado, pts_desvio

    def poligonos(self, via, pk_desde, pk_hasta):
        """Lista de polígonos Shapely (remates planos) para el rango de PK o recorrido completo."""
        res_aud = self.poligonos_con_auditoria(via, pk_desde, pk_hasta)
        return res_aud[0]


def anillo_cm(poligono):
    """Anillo exterior en sentido horario, en cm enteros, sin repetir el cierre."""
    pts = []
    for x, y in orient(poligono, sign=-1.0).exterior.coords[:-1]:
        p = (round(x * 100), round(y * 100))
        if not pts or p != pts[-1]:
            pts.append(p)
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts.pop()
    return pts


def solapes(rangos):
    """
    rangos: lista de (via, pk_desde, pk_hasta) o None (fila no válida).
    Devuelve [(i, j, mensaje)] con índices base 0 de los pares que se solapan.
    """
    por_via = {}
    for i, r in enumerate(rangos):
        if r is not None:
            por_via.setdefault(norm(r[0]), []).append((min(r[1], r[2]), max(r[1], r[2]), i))
    res = []
    for via, lista in por_via.items():
        lista.sort()
        for a, b in zip(lista[:-1], lista[1:]):
            if b[0] < a[1] - 1e-9:
                res.append((a[2], b[2],
                            f"{via}: las asignaciones {a[2] + 1} (PK {a[0]:g}-{a[1]:g}) y "
                            f"{b[2] + 1} (PK {b[0]:g}-{b[1]:g}) se solapan."))
    return res


# ----------------------------------------------------------------- salida
def escribir_ply(ruta_o_none, bloques, clase, enc="latin-1", nl="\r\n"):
    """
    bloques: lista de (esz, value, [(x, y), ...]). Separadores: tabulador.
    Con ruta escribe el archivo; con None devuelve los bytes.
    """
    L = [f"REM:\tPolygons for Feature Class:\t{clase}", "REM:"]
    for n, (esz, value, pts) in enumerate(bloques, start=1):
        L += ["", f"REM:\tPolygon#\t{n}", "", f"ESZ:\t{esz}", f"VALUE:\t{value}", f"NUM:\t{len(pts)}"]
        L += [f"XY:\t{x}\t{y}" for x, y in pts]
        L.append("END")
    datos = (nl.join(L) + nl).encode(enc)
    if ruta_o_none is None:
        return datos
    with open(ruta_o_none, "wb") as f:
        f.write(datos)


def generar(asignaciones, gen):
    """
    asignaciones: lista de (esz, via, pk_desde, pk_hasta, value).
    Devuelve (bloques, tabla): bloques para escribir_ply y filas para el CSV,
    descontando solapes geométricos en bifurcaciones/enlaces.
    """
    items = []
    for esz, via, pk0, pk1, value in asignaciones:
        polis = gen.poligonos(via, pk0, pk1)
        items.append({"esz": esz, "via": via, "pk0": pk0, "pk1": pk1, "value": value, "polis": polis})

    polis_recortados = descontar_solapes_geometricos(items)

    bloques, tabla = [], []
    for it, polis in zip(items, polis_recortados):
        for p in polis:
            bloques.append((it["esz"], it["value"], anillo_cm(p)))
        tabla.append((it["esz"], it["via"], it["pk0"], it["pk1"], it["value"], len(polis)))
    return bloques, tabla


def tabla_csv_bytes(tabla):
    import io
    buf = io.StringIO()
    w = csv.writer(buf, delimiter=";", lineterminator="\r\n")
    w.writerow(["ESZ", "VIA", "PK_DESDE", "PK_HASTA", "VALUE", "N_POLIGONOS"])
    w.writerows(tabla)
    return buf.getvalue().encode("utf-8-sig")


# -------------------------------------------------------------------- CLI
def leer_asignaciones(ruta):
    with open(ruta, newline="", encoding="utf-8-sig") as f:
        muestra = f.read(2048)
        f.seek(0)
        delim = ";" if muestra.splitlines()[0].count(";") else ","
        filas = []
        for i, fila in enumerate(csv.DictReader(f, delimiter=delim), start=2):
            fila = {k.strip().lower(): (v or "").strip() for k, v in fila.items() if k}
            try:
                filas.append((fila["via"], a_float(fila["pk_desde"]),
                              a_float(fila["pk_hasta"]), fila["value"]))
            except (KeyError, ValueError) as e:
                raise ValueError(f"Asignaciones, línea {i}: dato no válido ({e}).")
    return filas


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mdb")
    ap.add_argument("--mdb-ejes", action="append", default=[])
    ap.add_argument("--mdb-pk", action="append", default=[])
    ap.add_argument("--capa-ejes", default="CARRETERAS_BCN")
    ap.add_argument("--capa-pk", default="PK_VIA")
    ap.add_argument("--epsg-origen")
    ap.add_argument("--ejes", action="append", default=[])
    ap.add_argument("--pk", action="append", default=[])
    ap.add_argument("--asignaciones", required=True)
    ap.add_argument("--salida", required=True)
    ap.add_argument("--ancho", type=float, default=20.0, help="metros a cada lado (20)")
    ap.add_argument("--quad-segs", type=int, default=18)
    ap.add_argument("--esz-inicio", type=int, default=1)
    ap.add_argument("--clase", default="Cobertura MVIA.F6_Coberturas")
    ap.add_argument("--enc-shp", default="latin-1")
    ap.add_argument("--enc-ply", default="latin-1")
    ap.add_argument("--saltos", choices=["crlf", "lf"], default="crlf")
    a = ap.parse_args()

    mdb = lambda rol, ruta, capa: {"tipo": "mdb", "rol": rol, "ruta": ruta, "capa": capa,
                                   "epsg_origen": a.epsg_origen}
    shp = lambda rol, ruta: {"tipo": "shp", "rol": rol, "ruta": ruta}
    orig = ([mdb("ejes", r, a.capa_ejes) for r in a.mdb_ejes] + [mdb("pk", r, a.capa_pk) for r in a.mdb_pk]
            + [shp("ejes", r) for r in a.ejes] + [shp("pk", r) for r in a.pk])
    if a.mdb:
        orig += [mdb("ejes", a.mdb, a.capa_ejes), mdb("pk", a.mdb, a.capa_pk)]
    if not orig:
        sys.exit("ERROR: indique los orígenes de datos (--mdb-ejes/--mdb-pk, --ejes/--pk o --mdb).")
    datos = fuentes.cargar_fuentes(orig, a.enc_shp)
    for m in datos.avisos:
        aviso(m)

    asig = leer_asignaciones(a.asignaciones)
    sol = solapes([(v, p0, p1) for v, p0, p1, _ in asig])
    if sol:
        raise ValueError(sol[0][2])

    gen = Generador(datos, a.ancho, a.quad_segs)
    filas = [(esz, v, p0, p1, val) for esz, (v, p0, p1, val) in enumerate(asig, start=a.esz_inicio)]
    bloques, tabla = generar(filas, gen)
    escribir_ply(a.salida, bloques, a.clase, a.enc_ply, "\r\n" if a.saltos == "crlf" else "\n")
    ruta_csv = Path(a.salida).with_name(Path(a.salida).stem + "_asignaciones.csv")
    ruta_csv.write_bytes(tabla_csv_bytes(tabla))
    print(f"{len(bloques)} polígonos escritos en {a.salida}\nTabla de asignaciones: {ruta_csv}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, KeyError, RuntimeError) as e:
        sys.exit(f"ERROR: {e}")
