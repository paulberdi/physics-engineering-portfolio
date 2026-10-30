# -*- coding: utf-8 -*-
"""
Plantilla de ajuste — Experiencia 5 (caída en fluido viscoso, Ley de Stokes)

Qué hace:
1) Para CADA corrida ajusta s(t) al modelo lineal con rozamiento viscoso (Stokes):
       s(t)= s0 + v_t (t - t0) - v_t τ (1 - exp(-(t - t0)/τ))
   y obtiene: s0, t0, v_t, τ, con sus errores (±).
2) Con v_t y τ estima la viscosidad η de la glicerina por dos vías (consistencia):
       η_(v_t) = 2 r^2 g (ρ_s - ρ_f) / (9 v_t)
       η_(τ)   = 2 ρ_s r^2 / (9 τ)
3) Resume por corrida y promedia por ESFERA (mismo 'sphere_id').
4) Exporta CSV y muestra resultados.

Cómo usar:
- Declara una sola vez RHO_FLUID (densidad del fluido, kg/m^3).
- Declara esferas en SPHERES con: id, radio_mm y rho_s (kg/m^3).
- Pega las corridas en DATASETS con: 'sphere_id', 't' (s), 's' (m).
  Opcional: 'name', 'fit_window'=[tmin, tmax] para recortar, 'invert_sign'=True
- Requiere: numpy, pandas. (Si tienes SciPy, mejor ajuste con curve_fit).
"""

import numpy as np
import pandas as pd

g = 9.80665  # m/s^2

# ============ EDITA AQUÍ: propiedades del fluido ============
RHO_FLUID = 1260.0   # kg/m^3  (ejemplo glicerina ~ 1.26 g/cm^3). Cambia por tu medición con densitómetro.

# ============ EDITA AQUÍ: esferas utilizadas ============
# ρ_esfera típica (acero) ~ 7850 kg/m^3. Ajusta tus radios y densidades reales.
SPHERES = [
    {"sphere_id": "esf_peq", "radio_mm": 3.0,  "rho_s": 7850.0},
    {"sphere_id": "esf_med", "radio_mm": 4.0,  "rho_s": 7850.0},
    {"sphere_id": "esf_grd", "radio_mm": 4.75, "rho_s": 7850.0},
]
# Acceso rápido por id
SPH_MAP = {s["sphere_id"]: s for s in SPHERES}

# ============ PEGA AQUÍ TUS CORRIDAS (t en s, s en m) ============
DATASETS = [
    # EJEMPLO de formato:
    # {
    #   "name": "corrida_1",
    #   "sphere_id": "esf_peq",
    #   "invert_sign": True,             # si tu 's' decrece en el tiempo y quieres reportar v_t > 0
    #   "fit_window": [0.3, 2.5],        # opcional: recorta datos antes de ajustar
    #   "t": [0.00, 0.03, 0.06, ...],
    #   "s": [0.210, 0.205, 0.199, ...]
    # },
]

# ============ Utilidades matemáticas ============

def _as_np_clean(x):
    a = np.asarray(x, dtype=float)
    return a[np.isfinite(a)]

def _window_mask(t, w):
    if w is None: 
        return np.ones_like(t, dtype=bool)
    tmin, tmax = w
    return (t >= tmin) & (t <= tmax)

def _linregress_xy(x, y):
    # y ~ a x + b, por mínimos cuadrados (numpy)
    X = np.column_stack([x, np.ones_like(x)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = beta
    # varianzas (estimador clásico)
    yhat = X @ beta
    res = y - yhat
    dof = max(len(y) - 2, 1)
    s2 = (res @ res) / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    se_a, se_b = np.sqrt(np.diag(cov))
    return a, b, se_a, se_b

# Modelo Stokes (lineal en v): s(t)= s0 + v_t Δt - v_t τ (1 - e^{-Δt/τ}), con Δt = (t - t0) >= 0
def s_model(t, s0, vt, tau, t0):
    dt = np.maximum(t - t0, 0.0)  # evita exponentes con dt < 0 si el ajuste pone t0 levemente mayor al primer t
    return s0 + vt*dt - vt*tau*(1.0 - np.exp(-dt/tau))

def guess_initial_params(t, s):
    # s0 ~ primer valor
    s0 = s[0]
    # vt ~ pendiente de la cola: usa 20–30% final
    n = len(t)
    k0 = int(n*0.7)
    a, b, *_ = _linregress_xy(t[k0:], s[k0:])
    vt = max(abs(a), 1e-4)  # magnitud
    # tau: orden ~ 0.1–1 s. Estímalo por curvatura usando un ajuste cuadrático temprano
    k1 = int(n*0.3)
    k1 = max(k1, 5)
    p2, p1, _p0 = np.polyfit(t[:k1], s[:k1], 2)
    # En el régimen temprano s ≈ s0 + (g_eff*tau)t - 0.5(g_eff)*t^2 + ... — difícil directo;
    # usa heurística segura:
    tau = max((t[-1]-t[0]) / 5.0, 0.05)
    # t0: inicio efectivo (pequeño)
    t0 = t[0]
    return s0, vt, tau, t0

def fit_stokes(t, s):
    """
    Ajusta s(t) al modelo Stokes. Si SciPy está disponible, usa curve_fit para estimar
    parámetros y covarianza. Si no, hace un ajuste mixto:
      1) Estima vt con regresión lineal de la cola.
      2) Optimiza tau por búsqueda 1D minimizando SSE con vt fijo (búsqueda gruesa + refinamiento).
      3) Estima s0 y t0 finamente por mínimo cuadrados linealizado local.
    Devuelve: params, errors, R2
    """
    t = _as_np_clean(t)
    s = _as_np_clean(s)
    if len(t) != len(s) or len(t) < 8:
        raise ValueError("Cada corrida requiere ≥8 puntos limpios.")
    
    # Intento con SciPy
    try:
        from scipy.optimize import curve_fit
        p0 = guess_initial_params(t, s)
        popt, pcov = curve_fit(s_model, t, s, p0=p0, maxfev=10000)
        s0, vt, tau, t0 = popt
        # Errores 1σ
        perr = np.sqrt(np.diag(pcov))
        yhat = s_model(t, *popt)
        R2 = 1 - np.sum((s - yhat)**2) / np.sum((s - np.mean(s))**2)
        return (s0, vt, tau, t0), (perr[0], perr[1], perr[2], perr[3]), R2
    except Exception:
        pass

    # ——— Fallback SIN SciPy ———
    # 1) vt de la cola
    n = len(t)
    k0 = int(n*0.7)
    a, b, se_a, se_b = _linregress_xy(t[k0:], s[k0:])
    vt0 = abs(a)
    # 2) tau por búsqueda simple
    tspan = t[-1] - t[0]
    taus = np.geomspace(max(tspan/50, 0.01), max(tspan*2, 0.2), 60)
    s0_0 = s[0]
    t0_0 = t[0]
    best = None
    for tau in taus:
        yhat = s_model(t, s0_0, vt0, tau, t0_0)
        sse = float(np.sum((s - yhat)**2))
        if (best is None) or (sse < best[0]):
            best = (sse, tau)
    tau0 = best[1]
    # 3) Refina s0 y t0 (pequeño ajuste lineal local): linealiza por derivadas numéricas en s0 y t0
    #   s ≈ s_model + (∂s/∂s0)Δs0 + (∂s/∂t0)Δt0, con vt0 y tau0 fijos
    eps = 1e-6
    y0 = s_model(t, s0_0, vt0, tau0, t0_0)
    # Derivadas numéricas
    ds_ds0 = (s_model(t, s0_0+eps, vt0, tau0, t0_0) - y0)/eps
    ds_dt0 = (s_model(t, s0_0, vt0, tau0, t0_0+eps) - y0)/eps
    A = np.column_stack([ds_ds0, ds_dt0])
    d = s - y0
    dtheta, *_ = np.linalg.lstsq(A, d, rcond=None)
    s0 = s0_0 + dtheta[0]
    t0 = t0_0 + dtheta[1]
    # Recalcula con estos s0,t0
    yhat = s_model(t, s0, vt0, tau0, t0)
    R2 = 1 - np.sum((s - yhat)**2) / np.sum((s - np.mean(s))**2)
    # Incertidumbres: aproxima con matriz normal
    # Parám: [s0, vt, tau, t0] ~ estima sólo diag con diferencias finitas
    def jac(t, p):
        s0, vt, tau, t0 = p
        eps = 1e-6
        base = s_model(t, *p)
        J = []
        for i in range(4):
            dp = list(p)
            dp[i] += eps
            J.append((s_model(t, *dp) - base)/eps)
        return np.column_stack(J)
    p = np.array([s0, vt0, tau0, t0])
    J = jac(t, p)
    res = s - yhat
    dof = max(len(t) - 4, 1)
    sigma2 = (res @ res) / dof
    try:
        cov = sigma2 * np.linalg.pinv(J.T @ J)
        perr = np.sqrt(np.diag(cov))
    except Exception:
        perr = np.array([np.nan, np.nan, np.nan, np.nan])
    return (s0, vt0, tau0, t0), tuple(perr.tolist()), R2

def viscosity_from_vt_tau(vt, tau, r_m, rho_s, rho_f):
    # Dos estimaciones físicas de η; devuelve ambas y su promedio simple si son compatibles
    eta_vt  = (2 * (r_m**2) * g * (rho_s - rho_f)) / (9 * max(vt, 1e-12))
    eta_tau = (2 * rho_s * (r_m**2)) / (9 * max(tau, 1e-12))
    return eta_vt, eta_tau

# ============ Pipeline principal ============

rows = []
for d in DATASETS:
    name = d.get("name", "")
    sid  = d["sphere_id"]
    sph  = SPH_MAP.get(sid)
    if sph is None:
        raise ValueError(f"sphere_id '{sid}' no está declarado en SPHERES.")
    r_m   = float(sph["radio_mm"]) * 1e-3
    rho_s = float(sph["rho_s"])
    rho_f = float(RHO_FLUID)

    t = _as_np_clean(d["t"])
    s = _as_np_clean(d["s"])
    if len(t) != len(s):
        raise ValueError(f"{name}: t y s deben tener igual longitud.")

    # Recorta ventana si corresponde
    mask = _window_mask(t, d.get("fit_window"))
    t, s = t[mask], s[mask]

    # Opcional: invierte signo de s para tener s decreciente -> vt positiva (magnitud)
    if d.get("invert_sign", False):
        s = -s

    # Ajuste
    (s0, vt, tau, t0), (ds0, dvt, dtau, dt0), R2 = fit_stokes(t, s)

    # Estimación de η por ambas vías
    eta_vt, eta_tau = viscosity_from_vt_tau(vt, tau, r_m, rho_s, rho_f)

    # Promedio simple y compatibilidad (diferencia relativa)
    eta_mean = 0.5*(eta_vt + eta_tau)
    compat_rel = abs(eta_vt - eta_tau) / max(eta_mean, 1e-12)

    rows.append({
        "name": name,
        "sphere_id": sid,
        "r_m": r_m,
        "rho_s": rho_s,
        "rho_f": rho_f,
        "s0": s0, "ds0": ds0,
        "t0": t0, "dt0": dt0,
        "vt": vt, "dvt": dvt,
        "tau": tau, "dtau": dtau,
        "eta_from_vt": eta_vt,
        "eta_from_tau": eta_tau,
        "eta_mean": eta_mean,
        "compat_rel": compat_rel,
        "R2": R2
    })

df = pd.DataFrame(rows)

# — Promedios por esfera (si hay repetidas del mismo sphere_id)
grp = df.groupby("sphere_id", as_index=False).agg(
    r_m=("r_m", "first"),
    rho_s=("rho_s", "first"),
    rho_f=("rho_f", "first"),
    n=("name", "count"),
    vt_avg=("vt", "mean"),
    vt_sd =("vt", "std"),
    tau_avg=("tau", "mean"),
    tau_sd =("tau", "std"),
    eta_vt_avg=("eta_from_vt", "mean"),
    eta_vt_sd =("eta_from_vt", "std"),
    eta_tau_avg=("eta_from_tau", "mean"),
    eta_tau_sd =("eta_from_tau", "std"),
    eta_mean_avg=("eta_mean", "mean"),
    eta_mean_sd =("eta_mean", "std"),
    R2_avg=("R2", "mean")
)

# Exporta CSV
df.to_csv("exp5_corridas_stokes.csv", index=False)
grp.to_csv("exp5_resumen_por_esfera.csv", index=False)

# ——— Salida de consola ———
pd.set_option("display.width", 140)
pd.set_option("display.max_columns", 999)

print("=== Resumen por corrida (ajuste Stokes) ===")
if len(df):
    cols = ["name","sphere_id","r_m","rho_s","rho_f","vt","dvt","tau","dtau",
            "eta_from_vt","eta_from_tau","eta_mean","compat_rel","R2"]
    print(df[cols].to_string(index=False, justify="right",
                             formatters={
                                 "r_m": "{:.6f}".format,
                                 "vt": "{:.6f}".format, "dvt": "{:.6f}".format,
                                 "tau":"{:.6f}".format, "dtau":"{:.6f}".format,
                                 "eta_from_vt":"{:.4f}".format, "eta_from_tau":"{:.4f}".format,
                                 "eta_mean":"{:.4f}".format, "compat_rel":"{:.3f}".format,
                                 "R2":"{:.4f}".format
                             }))
else:
    print("(Aún no hay corridas en DATASETS.)")

print("\n=== Resumen por esfera (promedios) ===")
if len(grp):
    cols = ["sphere_id","n","r_m","vt_avg","vt_sd","tau_avg","tau_sd",
            "eta_vt_avg","eta_vt_sd","eta_tau_avg","eta_tau_sd","eta_mean_avg","eta_mean_sd","R2_avg"]
    print(grp[cols].to_string(index=False, justify="right",
                              formatters={
                                  "r_m":"{:.6f}".format,
                                  "vt_avg":"{:.6f}".format, "vt_sd":"{:.6f}".format,
                                  "tau_avg":"{:.6f}".format, "tau_sd":"{:.6f}".format,
                                  "eta_vt_avg":"{:.4f}".format, "eta_vt_sd":"{:.4f}".format,
                                  "eta_tau_avg":"{:.4f}".format, "eta_tau_sd":"{:.4f}".format,
                                  "eta_mean_avg":"{:.4f}".format, "eta_mean_sd":"{:.4f}".format,
                                  "R2_avg":"{:.4f}".format
                              }))
else:
    print("(Un solo ensayo por esfera o sin datos.)")

print("\nNotas:")
print(" - 'compat_rel' mide la diferencia relativa entre η_(v_t) y η_(τ): valores << 1 indican buena consistencia.")
print(" - Si tu 's(t)' decrece (eje hacia abajo en Tracker), usa 'invert_sign': True para reportar v_t > 0.")
print(" - Define RHO_FLUID según la medición con densitómetro y los radios reales de tus esferas.")
