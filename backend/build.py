#!/usr/bin/env python3
"""
Nieve España · generador de datos
=================================
Este programa se ejecuta solo en GitHub Actions (cada 12 horas) y prepara
los archivos que luego lee la web. Así la web no pide nada a ningún servicio
meteorológico y nunca choca con límites de uso.

Fuentes:
  1. ECMWF Open Data (oficial, CC BY 4.0): malla de toda Europa para las capas
     del mapa y una previsión "directa" en cada estación usando el perfil
     vertical de temperatura del modelo.
  2. Open-Meteo (gratis para uso NO comercial): ECMWF 9 km, ICON-EU y
     Météo-France AROME ya corregidos a cada cota, y el conjunto ICON-EPS
     para las probabilidades.
  3. AEMET OpenData (oficial): observaciones reales para verificar cada día
     cuánto acierta cada fuente y ajustar sus pesos.

Si una fuente falla, las demás siguen funcionando y la web avisa.
"""
import datetime as dt
import json
import math
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import requests

# ---------------------------------------------------------------------------
# 1) CONFIGURACIÓN
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
DATA = SITE / "data"
GRID_DIR = DATA / "grid"
VERIF = ROOT / "verif"
TZ = ZoneInfo("Europe/Madrid")
UTC = dt.timezone.utc

OM_API = "https://api.open-meteo.com/v1/forecast"
OM_ENS = "https://ensemble-api.open-meteo.com/v1/ensemble"
OM_ELEV = "https://api.open-meteo.com/v1/elevation"
AEMET_OBS = "https://opendata.aemet.es/opendata/api/observacion/convencional/todas"

# Fuentes de la mezcla y su peso inicial (la verificación con AEMET los ajusta)
SOURCES = [
    {"id": "ecmwf_ifs",            "name": "ECMWF IFS 9 km",           "w": 0.40, "color": "#38bdf8"},
    {"id": "icon_eu",              "name": "ICON-EU 7 km",             "w": 0.30, "color": "#a78bfa"},
    {"id": "meteofrance_seamless", "name": "Météo-France AROME",       "w": 0.20, "color": "#f472b6"},
    {"id": "ecmwf_od",             "name": "ECMWF oficial (perfil)",   "w": 0.10, "color": "#fbbf24"},
]
OM_MODELS = [s for s in SOURCES if s["id"] != "ecmwf_od"]
OM_VARS = ("precipitation,temperature_2m,relative_humidity_2m,wind_speed_10m,"
           "wind_gusts_10m,wind_direction_10m,cloud_cover,freezing_level_height,snow_depth")
ELEV = ["top", "ref", "base"]          # cima, media, base

# Malla de Europa (0,25°): de Canarias/norte de África a Escandinavia
GRID = {"S": 34.0, "N": 71.5, "W": -25.0, "E": 45.0, "step": 0.25}
GRID["nLat"] = round((GRID["N"] - GRID["S"]) / GRID["step"]) + 1   # 151
GRID["nLon"] = round((GRID["E"] - GRID["W"]) / GRID["step"]) + 1   # 281
SFC_STEPS = list(range(0, 145, 3))        # cada 3 h hasta 6 días
PL_STEPS = list(range(0, 145, 6))         # niveles de presión cada 6 h
PL_LEVELS = [1000, 925, 850, 700, 600, 500]
NODATA = -32768

STATIONS = json.loads((ROOT / "backend" / "stations.json").read_text(encoding="utf-8"))
for s in STATIONS:
    s["ref"] = round((s["base"] + s["top"]) / 2)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def isnum(v):
    return v is not None and not (isinstance(v, float) and math.isnan(v))


def local_key(d):
    return d.astimezone(TZ).strftime("%Y-%m-%dT%H:00")


def key_to_utc(k):
    return dt.datetime.strptime(k, "%Y-%m-%dT%H:%M").replace(tzinfo=TZ).astimezone(UTC)


def get_json(url, params=None, label="", headers=None, tries=4, raw=False):
    """Petición con reintentos. Si salta un límite de uso, espera 65 s."""
    err = "?"
    for k in range(tries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=120)
            if raw:
                r.raise_for_status()
                return r
            j = r.json()
        except Exception as e:
            err, wait = f"{label}: {e}", 20 * (k + 1)
        else:
            if r.ok and not (isinstance(j, dict) and j.get("error")):
                return j
            reason = j.get("reason") if isinstance(j, dict) else r.status_code
            err = f"{label}: {reason}"
            if r.status_code == 400:
                raise RuntimeError(err)
            wait = 65 if "limit" in str(reason).lower() else 20 * (k + 1)
        log("  aviso:", err, f"→ reintento en {wait} s")
        time.sleep(wait)
    raise RuntimeError(err)


# ---------------------------------------------------------------------------
# 2) FÍSICA DE LA NIEVE (la misma que explica la web)
# ---------------------------------------------------------------------------
SLR_T = [-15, -10, -5, -2, 0, 1]
SLR_R = [1.8, 1.6, 1.3, 1.0, 0.8, 0.7]     # cm de nieve por mm de agua


def wet_bulb(T, RH):
    """Temperatura de bulbo húmedo (Stull 2011). Funciona con números o arrays."""
    RH = np.clip(RH, 5, 99)
    return (T * np.arctan(0.151977 * np.sqrt(RH + 8.313659)) + np.arctan(T + RH)
            - np.arctan(RH - 1.676331) + 0.00391838 * RH ** 1.5 * np.arctan(0.023101 * RH) - 4.686035)


def precip_to_snow(p, T, RH, wind):
    """mm de agua → (cm de nieve, mm de lluvia). Admite arrays de numpy."""
    p, T = np.asarray(p, float), np.asarray(T, float)
    RH = np.full_like(T, np.nan) if RH is None else np.asarray(RH, float)
    wind = np.zeros_like(T) if wind is None else np.nan_to_num(np.asarray(wind, float))
    tw = np.where(np.isnan(RH), T - 1, wet_bulb(T, np.nan_to_num(RH, nan=85)))
    fs = np.clip((1.5 - tw) / 1.5, 0, 1)
    r = np.interp(T, SLR_T, SLR_R)
    r = np.where(wind > 30, r * (1 - np.minimum(0.25, (wind - 30) / 200)), r)
    p0 = np.where(p > 0, p, 0)
    snow, rain = p0 * fs * r, p0 * (1 - fs)
    bad = np.isnan(p) | np.isnan(T)
    return np.where(bad, np.nan, snow), np.where(bad, np.nan, rain)


def arr(values):
    """Lista con None → array con NaN."""
    return np.array([np.nan if v is None else v for v in values], dtype=float)


# ---------------------------------------------------------------------------
# 3) OPEN-METEO: 3 modelos × 28 estaciones × 3 cotas
# ---------------------------------------------------------------------------
def fetch_openmeteo():
    lat, lon, ele = [], [], []
    for s in STATIONS:
        for k in ELEV:
            lat.append(s["lat"]); lon.append(s["lon"]); ele.append(s[k])
    out, errors = {}, {}
    for m in OM_MODELS:
        try:
            log("Open-Meteo:", m["name"])
            r = get_json(OM_API, {"latitude": ",".join(map(str, lat)), "longitude": ",".join(map(str, lon)),
                                  "elevation": ",".join(map(str, ele)), "hourly": OM_VARS, "models": m["id"],
                                  "timezone": "Europe/Madrid", "past_days": 2, "forecast_days": 7}, m["name"])
            r = r if isinstance(r, list) else [r]
            out[m["id"]] = [x["hourly"] for x in r]
        except Exception as e:
            errors[m["id"]] = str(e)
            log("  ERROR", e)
        time.sleep(10)          # sin prisas: muy por debajo de 600 llamadas/minuto
    return out, errors


def derive_om(h):
    """Respuesta horaria de Open-Meteo → nuestras variables."""
    p = arr(h["precipitation"]); T = arr(h["temperature_2m"])
    rh = arr(h.get("relative_humidity_2m", [None] * len(p)))
    w = arr(h.get("wind_speed_10m", [None] * len(p)))
    snow, rain = precip_to_snow(p, T, rh, w)
    return {"p": p, "snow": snow, "rain": rain, "temp": T, "wind": w,
            "gust": arr(h.get("wind_gusts_10m", [None] * len(p))),
            "dir": arr(h.get("wind_direction_10m", [None] * len(p))),
            "cloud": arr(h.get("cloud_cover", [None] * len(p))),
            "fl": arr(h.get("freezing_level_height", [None] * len(p))),
            # espesor de nieve en el suelo según el modelo (m → cm). No sabe de nieve artificial ni pisado.
            "depth": arr(h.get("snow_depth", [None] * len(p))) * 100}


# ---------------------------------------------------------------------------
# 4) ECMWF OPEN DATA: descarga GRIB2 y lectura
# ---------------------------------------------------------------------------
def download_ecmwf(tmp):
    from ecmwf.opendata import Client
    last = None
    for source in ("ecmwf", "aws", "azure"):
        try:
            log("ECMWF Open Data: fuente", source)
            c = Client(source=source, model="ifs", resol="0p25")
            run = c.latest(type="fc", stream="oper", step=144, param="2t")
            common = {"date": int(run.strftime("%Y%m%d")), "time": run.hour, "stream": "oper", "type": "fc"}
            log("  ejecución", run.strftime("%Y-%m-%d %HZ"), "· superficie…")
            c.retrieve(**common, step=SFC_STEPS, param=["2t", "2d", "10u", "10v", "tp"], target=str(tmp / "sfc.grib2"))
            log("  niveles de presión…")
            c.retrieve(**common, step=PL_STEPS, param=["t", "gh", "r"], levelist=PL_LEVELS, target=str(tmp / "pl.grib2"))
            return run.replace(tzinfo=UTC)
        except Exception as e:
            last = e
            log("  falló:", e)
    raise RuntimeError(f"ECMWF Open Data no disponible: {last}")


class Cropper:
    """Recorta un campo global a la malla de Europa (sur→norte, oeste→este)."""

    def __init__(self, Ni, Nj, lat1, lon1, dlat, dlon, jpos):
        lats = lat1 + np.arange(Nj) * dlat * (1 if jpos else -1)
        lons = (lon1 + np.arange(Ni) * dlon + 180) % 360 - 180
        eps = 1e-6
        j = np.where((lats >= GRID["S"] - eps) & (lats <= GRID["N"] + eps))[0]
        i = np.where((lons >= GRID["W"] - eps) & (lons <= GRID["E"] + eps))[0]
        self.j, self.i, self.Ni, self.Nj = j[np.argsort(lats[j])], i[np.argsort(lons[i])], Ni, Nj
        if len(self.j) != GRID["nLat"] or len(self.i) != GRID["nLon"]:
            raise RuntimeError(f"La malla de ECMWF no encaja: {len(self.j)}×{len(self.i)}")

    def __call__(self, values):
        return values.reshape(self.Nj, self.Ni)[np.ix_(self.j, self.i)].astype(np.float32)


def read_grib(path):
    """Devuelve {(nombre, nivel, paso): campo recortado}."""
    import eccodes
    fields, crop = {}, None
    with open(path, "rb") as f:
        while True:
            h = eccodes.codes_grib_new_from_file(f)
            if h is None:
                break
            try:
                if crop is None:
                    g = lambda k: eccodes.codes_get_double(h, k)
                    crop = Cropper(eccodes.codes_get_long(h, "Ni"), eccodes.codes_get_long(h, "Nj"),
                                   g("latitudeOfFirstGridPointInDegrees"), g("longitudeOfFirstGridPointInDegrees"),
                                   g("jDirectionIncrementInDegrees"), g("iDirectionIncrementInDegrees"),
                                   eccodes.codes_get_long(h, "jScansPositively"))
                name = eccodes.codes_get_string(h, "shortName")
                level = eccodes.codes_get_long(h, "level")
                step = eccodes.codes_get_long(h, "endStep")
                fields[(name, level, step)] = crop(eccodes.codes_get_values(h))
            finally:
                eccodes.codes_release(h)
    return fields


def sfc(fields, name, step):
    for (n, _lev, st), v in fields.items():
        if n == name and st == step:
            return v
    raise KeyError(f"{name} paso {step}")


def freezing_level(T_levels, gh_levels):
    """Altura de 0 °C a partir del perfil (niveles ordenados de abajo arriba, T en °C)."""
    fl = np.full(T_levels[0].shape, np.nan, dtype=np.float32)
    for l in range(len(T_levels) - 1):
        t0, t1, z0, z1 = T_levels[l], T_levels[l + 1], gh_levels[l], gh_levels[l + 1]
        cross = np.isnan(fl) & (t0 > 0) & (t1 <= 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = z0 + (0 - t0) * (z1 - z0) / (t1 - t0)
        fl = np.where(cross, z, fl)
    fl = np.where(np.isnan(fl) & (T_levels[0] <= 0), 0, fl)          # helada desde el suelo
    fl = np.where(np.isnan(fl), gh_levels[-1], fl)                    # todo por encima de 0 °C
    return fl


def rh_from_dew(T, Td):
    return 100 * np.exp(17.625 * Td / (243.04 + Td)) / np.exp(17.625 * T / (243.04 + T))


def bilinear(field, lat, lon):
    gy, gx = (lat - GRID["S"]) / GRID["step"], (lon - GRID["W"]) / GRID["step"]
    r0, c0 = int(math.floor(gy)), int(math.floor(gx))
    r1, c1 = min(r0 + 1, GRID["nLat"] - 1), min(c0 + 1, GRID["nLon"] - 1)
    fy, fx = gy - r0, gx - c0
    return float((field[r0, c0] * (1 - fx) + field[r0, c1] * fx) * (1 - fy)
                 + (field[r1, c0] * (1 - fx) + field[r1, c1] * fx) * fy)


def build_ecmwf(sfc_f, pl_f, run):
    """Genera los archivos de la malla y devuelve lo necesario para las estaciones."""
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    for old in GRID_DIR.glob("s*.bin"):
        old.unlink()
    # perfiles en los pasos de 6 h
    prof = {}
    for st in PL_STEPS:
        T = [pl_f[("t", lev, st)] - 273.15 for lev in PL_LEVELS]
        Z = [pl_f[("gh", lev, st)] for lev in PL_LEVELS]
        R = [pl_f[("r", lev, st)] for lev in PL_LEVELS]
        prof[st] = (T, Z, R, freezing_level(T, Z))
    times, cum = [], np.zeros((GRID["nLat"], GRID["nLon"]), np.float32)
    prev_tp = None
    for st in SFC_STEPS:
        T = sfc(sfc_f, "2t", st) - 273.15
        rh = rh_from_dew(T, sfc(sfc_f, "2d", st) - 273.15)
        u, v = sfc(sfc_f, "10u", st), sfc(sfc_f, "10v", st)
        tp = sfc(sfc_f, "tp", st) * 1000                          # m → mm, acumulado
        p = np.zeros_like(tp) if prev_tp is None else np.maximum(0, tp - prev_tp)
        prev_tp = tp
        snow, _ = precip_to_snow(p, T, rh, np.hypot(u, v) * 3.6)
        cum = cum + np.nan_to_num(snow)
        if st in prof:
            fl = prof[st][3]
        else:                                                      # interpola entre pasos de 6 h
            a, b = st - 3, st + 3
            fl = 0.5 * (prof[a][3] + prof[b][3]) if a in prof and b in prof else prof[max(k for k in prof if k <= st)][3]
        pack = np.stack([cum * 10, T * 10, u * 10, v * 10, fl])
        pack = np.where(np.isfinite(pack), np.clip(np.round(pack), -32767, 32767), NODATA).astype("<i2")
        pack.tofile(GRID_DIR / f"s{st:03d}.bin")
        times.append(local_key(run + dt.timedelta(hours=st)))
    meta = {"source": "ECMWF IFS 0,25° (Open Data, CC BY 4.0)", "run": run.strftime("%Y-%m-%dT%H:%MZ"),
            "S": GRID["S"], "W": GRID["W"], "N": GRID["N"], "E": GRID["E"], "step": GRID["step"],
            "nLat": GRID["nLat"], "nLon": GRID["nLon"], "steps": SFC_STEPS, "times": times,
            "vars": ["snowcum_cm_x10", "temp_c_x10", "u_ms_x10", "v_ms_x10", "freezing_level_m"],
            "nodata": NODATA}
    (GRID_DIR / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    log(f"  malla escrita: {len(times)} pasos de {GRID['nLat']}×{GRID['nLon']}")
    return prof


def ecmwf_station_series(sfc_f, prof, run, timeline):
    """Previsión en cada estación y cota usando el perfil vertical (como hace Snow-Forecast)."""
    hours = np.array([(key_to_utc(k) - run).total_seconds() / 3600 for k in timeline])
    valid = (hours >= 0) & (hours <= SFC_STEPS[-1])
    out = []
    for s in STATIONS:
        tp = np.array([bilinear(sfc(sfc_f, "tp", st), s["lat"], s["lon"]) * 1000 for st in SFC_STEPS])
        uu = np.array([bilinear(sfc(sfc_f, "10u", st), s["lat"], s["lon"]) for st in SFC_STEPS])
        vv = np.array([bilinear(sfc(sfc_f, "10v", st), s["lat"], s["lon"]) for st in SFC_STEPS])
        cumh = np.interp(hours, SFC_STEPS, tp)
        p = np.maximum(0, np.diff(cumh, prepend=cumh[0]))
        wind = np.interp(hours, SFC_STEPS, np.hypot(uu, vv) * 3.6)
        wdir = np.interp(hours, SFC_STEPS, (np.degrees(np.arctan2(-uu, -vv)) + 360) % 360)
        fl_pt = np.array([bilinear(prof[st][3], s["lat"], s["lon"]) for st in PL_STEPS])
        per_el = []
        for k in ELEV:
            z = s[k]
            Tz, Rz = [], []
            for st in PL_STEPS:
                T, Z, R, _ = prof[st]
                zs = [bilinear(Z[l], s["lat"], s["lon"]) for l in range(len(PL_LEVELS))]
                Tz.append(np.interp(z, zs, [bilinear(T[l], s["lat"], s["lon"]) for l in range(len(PL_LEVELS))]))
                Rz.append(np.interp(z, zs, [bilinear(R[l], s["lat"], s["lon"]) for l in range(len(PL_LEVELS))]))
            temp = np.interp(hours, PL_STEPS, Tz)
            rh = np.interp(hours, PL_STEPS, Rz)
            snow, rain = precip_to_snow(p, temp, rh, wind)
            nan = np.where(valid, 1.0, np.nan)
            per_el.append({"p": p * nan, "snow": snow * nan, "rain": rain * nan, "temp": temp * nan,
                           "wind": wind * nan, "gust": np.full(len(hours), np.nan), "dir": wdir * nan,
                           "cloud": np.full(len(hours), np.nan),
                           "fl": np.interp(hours, PL_STEPS, fl_pt) * nan,
                           "depth": np.full(len(hours), np.nan)})
        out.append(per_el)
    return out


# ---------------------------------------------------------------------------
# 5) MEZCLA, CONJUNTO, COMPROBACIÓN DE COORDENADAS
# ---------------------------------------------------------------------------
def blend(per, weights, n):
    out = {}
    for k in ["p", "snow", "rain", "temp", "wind", "gust", "cloud", "fl", "depth"]:
        sv, sw = np.zeros(n), np.zeros(n)
        for sid, d in per.items():
            v = d[k]
            ok = ~np.isnan(v)
            sv[ok] += weights[sid] * v[ok]
            sw[ok] += weights[sid]
        with np.errstate(invalid="ignore", divide="ignore"):
            out[k] = np.where(sw > 0, sv / np.where(sw > 0, sw, 1), np.nan)
    d = np.full(n, np.nan)
    for sid in per:                               # dirección: la de la primera fuente con dato
        v = per[sid]["dir"]
        d = np.where(np.isnan(d), v, d)
    out["dir"] = d
    return out


def fetch_ensemble(now_key):
    log("Open-Meteo: conjunto ICON-EPS (probabilidades)")
    r = get_json(OM_ENS, {"latitude": ",".join(str(s["lat"]) for s in STATIONS),
                          "longitude": ",".join(str(s["lon"]) for s in STATIONS),
                          "elevation": ",".join(str(s["ref"]) for s in STATIONS),
                          "hourly": "precipitation,temperature_2m", "models": "dwd_icon_seamless_eps",
                          "timezone": "Europe/Madrid", "forecast_days": 7}, "conjunto")
    r = r if isinstance(r, list) else [r]
    res = []
    for loc in r:
        h = loc["hourly"]
        now = h["time"].index(now_key) if now_key in h["time"] else 0
        pkeys = [k for k in h if k == "precipitation" or k.startswith("precipitation_member")]

        def total(pk, hours):
            tk = pk.replace("precipitation", "temperature_2m")
            p = arr(h[pk][now + 1: now + 1 + hours]); T = arr(h.get(tk, [None] * len(h[pk]))[now + 1: now + 1 + hours])
            s, _ = precip_to_snow(p, T, None, None)
            return float(np.nansum(s))

        def stats(hours):
            v = np.sort([total(k, hours) for k in pkeys])
            q = lambda x: float(v[min(len(v) - 1, int(x * len(v)))])
            pr = lambda x: round(100 * float(np.mean(v >= x)))
            return {"n": len(v), "p10": round(q(0.1), 1), "p50": round(q(0.5), 1), "p90": round(q(0.9), 1),
                    "p5": pr(5), "p20": pr(20), "p50cm": pr(50)}
        res.append({"from": now_key, "d3": stats(72), "d7": stats(168)})
    return res


def check_coords():
    r = get_json(OM_ELEV, {"latitude": ",".join(str(s["lat"]) for s in STATIONS),
                           "longitude": ",".join(str(s["lon"]) for s in STATIONS)}, "altitudes")
    r = r[0] if isinstance(r, list) else r
    for s, e in zip(STATIONS, r["elevation"]):
        s["dem"] = round(e)
        s["coordOk"] = s["base"] - 250 <= s["dem"] <= s["top"] + 250


# ---------------------------------------------------------------------------
# 6) VERIFICACIÓN CON AEMET
# Ventana de 24 h de 06 a 06 UTC. Se compara la última previsión emitida antes
# de que empiece la ventana con la lluvia/nieve (en agua) medida por AEMET en la
# estación oficial más cercana a cada estación de esquí.
# ---------------------------------------------------------------------------
def window_after(t):
    ws = t.replace(hour=6, minute=0, second=0, microsecond=0)
    return ws if ws > t else ws + dt.timedelta(days=1)


def archive_forecast(timeline, per_source_mid, blend_mid, now_utc):
    ws = window_after(now_utc)
    we = ws + dt.timedelta(hours=24)
    sel = np.array([ws < key_to_utc(k) <= we for k in timeline])
    rec = {"issued": now_utc.isoformat(), "window": [ws.isoformat(), we.isoformat()], "sources": {}}
    for sid in [s["id"] for s in SOURCES] + ["blend"]:
        vals = {}
        for i, st in enumerate(STATIONS):
            d = blend_mid[i] if sid == "blend" else per_source_mid[i].get(sid)
            if d is None:
                continue
            v = d["p"][sel]
            if len(v) and np.isfinite(v).sum() >= 20:
                vals[st["name"]] = round(float(np.nansum(v)), 1)
        if vals:
            rec["sources"][sid] = vals
    (VERIF / "fc").mkdir(parents=True, exist_ok=True)
    (VERIF / "fc" / f"{ws.date()}.json").write_text(json.dumps(rec, ensure_ascii=False), encoding="utf-8")
    log("  previsión archivada para la ventana", ws.date())


def haversine(a, b, c, d):
    R = 6371
    p1, p2 = math.radians(a), math.radians(c)
    dphi, dl = p2 - p1, math.radians(d - b)
    x = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(x))


def verify_with_aemet(now_utc):
    key = os.environ.get("AEMET_API_KEY", "").strip()
    if not key:
        return "sin clave de AEMET (verificación desactivada)"
    we = now_utc.replace(hour=6, minute=0, second=0, microsecond=0)
    if now_utc < we:
        we -= dt.timedelta(days=1)
    ws = we - dt.timedelta(hours=24)
    date = str(ws.date())
    scores_path = VERIF / "scores.jsonl"
    done = set()
    if scores_path.exists():
        done = {json.loads(l)["date"] for l in scores_path.read_text(encoding="utf-8").splitlines() if l.strip()}
    if date in done:
        return f"ventana {date} ya verificada"
    fc_path = VERIF / "fc" / f"{date}.json"
    if not fc_path.exists():
        return f"no hay previsión archivada para {date} (normal los primeros días)"
    fc = json.loads(fc_path.read_text(encoding="utf-8"))
    log("AEMET: observaciones de las últimas 24 h")
    meta = get_json(AEMET_OBS, headers={"api_key": key}, label="AEMET")
    datos = get_json(meta["datos"], label="AEMET datos", raw=True)
    obs = json.loads(datos.content.decode("latin-1"))
    agg = {}
    for o in obs:
        if o.get("prec") is None or "fint" not in o:
            continue
        t = dt.datetime.strptime(o["fint"][:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=UTC)
        if not (ws < t <= we):
            continue
        a = agg.setdefault(o["idema"], {"lat": o["lat"], "lon": o["lon"], "alt": o.get("alt", 0),
                                        "name": o.get("ubi", o["idema"]), "sum": 0.0, "n": 0})
        a["sum"] += float(o["prec"]); a["n"] += 1
    lines = []
    for st in STATIONS:
        best, score = None, 1e9
        for idema, a in agg.items():
            if a["n"] < 20:
                continue
            d = haversine(st["lat"], st["lon"], a["lat"], a["lon"])
            if d > 30:
                continue
            sc = d + abs(a["alt"] - st["base"]) / 100       # cerca y a una altitud parecida
            if sc < score:
                best, score = (idema, a, d), sc
        if not best:
            continue
        idema, a, d = best
        fcs = {sid: v[st["name"]] for sid, v in fc["sources"].items() if st["name"] in v}
        if fcs:
            lines.append({"date": date, "station": st["name"], "obs": round(a["sum"], 1), "fc": fcs,
                          "aemet": {"id": idema, "name": a["name"], "alt": a["alt"], "km": round(d, 1)}})
    with scores_path.open("a", encoding="utf-8") as f:
        for l in lines:
            f.write(json.dumps(l, ensure_ascii=False) + "\n")
        if not lines:   # marca el día como hecho aunque no haya pares válidos
            f.write(json.dumps({"date": date, "station": None}) + "\n")
    return f"ventana {date}: {len(lines)} estaciones comparadas"


def summarize_verification():
    """Error medio de cada fuente en los últimos 45 días y pesos ajustados."""
    path = VERIF / "scores.jsonl"
    prior = {s["id"]: s["w"] for s in SOURCES}
    if not path.exists():
        return {"weights": prior, "n": 0, "sources": {}, "pairs": {}}
    limit = str((dt.date.today() - dt.timedelta(days=45)))
    rows = [json.loads(l) for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
    rows = [r for r in rows if r.get("station") and r["date"] >= limit]
    stats = {}
    for sid in [s["id"] for s in SOURCES] + ["blend"]:
        err = [(r["fc"][sid] - r["obs"]) for r in rows if sid in r["fc"] and max(r["obs"], r["fc"][sid]) >= 0.5]
        if err:
            stats[sid] = {"n": len(err), "mae": round(float(np.mean(np.abs(err))), 2), "bias": round(float(np.mean(err)), 2)}
    # pesos: inverso del error, mezclado con el peso inicial según cuántos casos haya
    ids = [s["id"] for s in SOURCES]
    inv = {i: 1 / (stats[i]["mae"] + 0.5) for i in ids if i in stats and stats[i]["n"] >= 5}
    w = dict(prior)
    if inv:
        tot = sum(inv.values())
        for i in ids:
            n = stats.get(i, {}).get("n", 0) if i in inv else 0
            learned = inv[i] / tot * sum(prior[j] for j in inv) if i in inv else prior[i]
            w[i] = (prior[i] * 15 + learned * n) / (15 + n)
        s = sum(w.values())
        w = {i: round(v / s, 3) for i, v in w.items()}
    pairs = {}
    for r in rows:
        pairs[r["station"]] = r["aemet"]
    return {"weights": w, "n": len({r["date"] for r in rows}), "sources": stats, "pairs": pairs,
            "updated": dt.datetime.now(UTC).isoformat()}


# ---------------------------------------------------------------------------
# 7) PROGRAMA PRINCIPAL
# ---------------------------------------------------------------------------
def rnd(a, nd):
    return [None if not np.isfinite(x) else (round(float(x), nd) if nd > 0 else int(round(float(x)))) for x in a]


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    VERIF.mkdir(parents=True, exist_ok=True)
    now_utc = dt.datetime.now(UTC)
    now_key = local_key(now_utc)
    status = {"generated": now_utc.isoformat(), "sources": {}}

    # 7.1 Verificación del día anterior (antes de construir, para usar los pesos nuevos)
    try:
        status["sources"]["aemet"] = {"ok": True, "msg": verify_with_aemet(now_utc)}
    except Exception as e:
        traceback.print_exc()
        status["sources"]["aemet"] = {"ok": False, "msg": str(e)}
    summary = summarize_verification()
    weights = summary["weights"]
    (VERIF / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    (DATA / "verif.json").write_text(json.dumps(summary, ensure_ascii=False), encoding="utf-8")
    log("Pesos en uso:", weights)

    # 7.2 Open-Meteo
    om, om_err = {}, {}
    try:
        om, om_err = fetch_openmeteo()
    except Exception as e:
        om_err = {"todos": str(e)}
    status["sources"]["openmeteo"] = {"ok": bool(om), "errors": om_err}

    # 7.3 ECMWF Open Data
    ec_series, run = None, None
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            run = download_ecmwf(tmp)
            log("Leyendo GRIB…")
            sfc_f, pl_f = read_grib(tmp / "sfc.grib2"), read_grib(tmp / "pl.grib2")
            prof = build_ecmwf(sfc_f, pl_f, run)
            if om:
                timeline = next(iter(om.values()))[0]["time"]
            else:
                timeline = [local_key(run + dt.timedelta(hours=h)) for h in range(0, SFC_STEPS[-1] + 1)]
            ec_series = ecmwf_station_series(sfc_f, prof, run, timeline)
        status["sources"]["ecmwf"] = {"ok": True, "run": run.strftime("%Y-%m-%dT%H:%MZ")}
    except Exception as e:
        traceback.print_exc()
        status["sources"]["ecmwf"] = {"ok": False, "msg": str(e)}

    # 7.4 Estaciones
    if om or ec_series:
        timeline = next(iter(om.values()))[0]["time"] if om else \
            [local_key(run + dt.timedelta(hours=h)) for h in range(0, SFC_STEPS[-1] + 1)]
        n = len(timeline)
        used = [s for s in SOURCES if (s["id"] in om) or (s["id"] == "ecmwf_od" and ec_series)]
        try:
            check_coords()
        except Exception as e:
            log("  sin comprobación de coordenadas:", e)
        try:
            ens = fetch_ensemble(now_key)
            status["sources"]["ensemble"] = {"ok": True}
        except Exception as e:
            ens = None
            status["sources"]["ensemble"] = {"ok": False, "msg": str(e)}
        out_st, per_mid, blend_mid = [], [], []
        for i, s in enumerate(STATIONS):
            els = []
            for k in range(3):
                per = {}
                for src in used:
                    if src["id"] == "ecmwf_od":
                        per["ecmwf_od"] = ec_series[i][k]
                    else:
                        per[src["id"]] = derive_om(om[src["id"]][i * 3 + k])
                b = blend(per, weights, n)
                els.append({"b": {"snow": rnd(b["snow"], 2), "rain": rnd(b["rain"], 2), "temp": rnd(b["temp"], 1),
                                  "wind": rnd(b["wind"], 0), "gust": rnd(b["gust"], 0), "dir": rnd(b["dir"], 0),
                                  "cloud": rnd(b["cloud"], 0), "fl": rnd(b["fl"], 0), "depth": rnd(b["depth"], 0)},
                            "per": {sid: {"snow": rnd(d["snow"], 2)} for sid, d in per.items()}})
                if k == 1:
                    per_mid.append(per); blend_mid.append(b)
            out_st.append({**{k: s[k] for k in ("name", "zone", "lat", "lon", "base", "ref", "top", "infonieve", "aviso") if k in s},
                           "dem": s.get("dem"), "coordOk": s.get("coordOk"), "el": els,
                           "ens": ens[i] if ens else None})
        stations = {"generated": now_utc.isoformat(), "time": timeline,
                    "sources": [{**s, "w": weights.get(s["id"], s["w"])} for s in used], "stations": out_st}
        (DATA / "stations.json").write_text(json.dumps(stations, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
        log("stations.json escrito con", len(used), "fuentes")
        try:
            archive_forecast(timeline, per_mid, blend_mid, now_utc)
        except Exception as e:
            log("  no se pudo archivar la previsión:", e)
    else:
        log("ATENCIÓN: ninguna fuente de estaciones disponible; se mantiene el archivo anterior")

    (DATA / "meta.json").write_text(json.dumps(status, ensure_ascii=False, indent=1), encoding="utf-8")
    ok = status["sources"]["openmeteo"]["ok"] or status["sources"]["ecmwf"]["ok"]
    log("Terminado." if ok else "Terminado CON ERRORES.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
