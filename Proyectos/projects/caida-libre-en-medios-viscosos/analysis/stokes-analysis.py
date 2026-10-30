# -*- coding: utf-8 -*-
"""
Análisis de Caída en Medio Viscoso - Ley de Stokes
Método: VELOCIDAD TERMINAL (equilibrio de fuerzas)

Este código:
1) Identifica la región de movimiento uniforme (velocidad terminal)
2) Calcula v_terminal mediante ajuste lineal: s = s₀ + v_t·t
3) Determina la viscosidad: η = 2r²g(ρₛ - ρf)/(9v_t)
4) Genera gráficos de posición, velocidad y comparación entre esferas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Constantes
g = 9.80665  # m/s²

# ============ CONFIGURACIÓN EXPERIMENTAL ============
RHO_FLUID = 1247.0   # kg/m³ - MEDIR CON DENSITÓMETRO

# Esferas (MEDIR radios con calibre)
SPHERES = [
    {"sphere_id": "esfera1", "radio_mm": 3.00, "rho_s": 7850.0},
    {"sphere_id": "esfera2", "radio_mm": 4.00, "rho_s": 7850.0},
    {"sphere_id": "esfera3", "radio_mm": 5.00, "rho_s": 7850.0},
]

SPH_MAP = {s["sphere_id"]: s for s in SPHERES}

# ============ DATOS EXPERIMENTALES ============
DATASETS = [
    {
        "name": "esfera1_corrida_1",
        "sphere_id": "esfera1",
        "invert_sign": True,
        "terminal_window": [0.00, 2.384000],  # Región de velocidad constante
        "t": [
        0.000000E0, 1.598333E-2, 3.213333E-2, 4.810000E-2, 6.398333E-2, 8.010000E-2,
        9.590000E-2, 1.119833E-1, 1.279000E-1, 1.439500E-1, 1.599667E-1, 1.759667E-1,
        1.920833E-1, 2.079500E-1, 2.239667E-1, 2.399000E-1, 2.560500E-1, 2.719333E-1,
        3.039167E-1, 3.199000E-1, 3.359167E-1, 3.518333E-1, 3.678333E-1, 3.838833E-1,
        3.998333E-1, 4.158333E-1, 4.319500E-1, 4.479500E-1, 4.639000E-1, 4.799000E-1,
        4.958667E-1, 5.119000E-1, 5.279167E-1, 5.440667E-1, 5.600500E-1, 5.762000E-1,
        5.920500E-1, 6.080667E-1, 6.241333E-1, 6.400167E-1, 6.560500E-1, 6.720500E-1,
        7.040500E-1, 7.198667E-1, 7.359333E-1, 7.518667E-1, 7.678667E-1, 7.838667E-1,
        7.999667E-1, 8.158833E-1, 8.318833E-1, 8.479333E-1, 8.639000E-1, 8.798500E-1,
        8.959500E-1, 9.120000E-1, 9.279667E-1, 9.438667E-1, 9.600167E-1, 9.759333E-1,
        9.920167E-1, 1.007950E0, 1.024000E0, 1.039900E0, 1.056083E0, 1.072000E0,
        1.103817E0, 1.119933E0, 1.135950E0, 1.151933E0, 1.167967E0, 1.183983E0,
        1.199883E0, 1.216033E0, 1.232017E0, 1.248033E0, 1.263950E0, 1.279967E0,
        1.296033E0, 1.311967E0, 1.328000E0, 1.343967E0, 1.359867E0, 1.376050E0,
        1.392067E0, 1.408083E0, 1.423983E0, 1.440050E0, 1.456167E0, 1.472283E0,
        1.503950E0, 1.519950E0, 1.535933E0, 1.551983E0, 1.568017E0, 1.584017E0,
        1.600017E0, 1.615950E0, 1.631933E0, 1.648033E0, 1.663950E0, 1.679950E0,
        1.696017E0, 1.712133E0, 1.728300E0, 1.744250E0, 1.760050E0, 1.776017E0,
        1.792317E0, 1.808100E0, 1.824017E0, 1.840050E0, 1.856300E0, 1.872283E0,
        1.903950E0, 1.919983E0, 1.936017E0, 1.952000E0, 1.968033E0, 1.984017E0,
        2.000067E0, 2.016067E0, 2.032067E0, 2.048033E0, 2.064200E0, 2.080383E0,
        2.096033E0, 2.111983E0, 2.128000E0, 2.144050E0, 2.160200E0, 2.176000E0,
        2.192133E0, 2.208167E0, 2.224133E0, 2.240283E0, 2.256167E0, 2.272250E0,
        2.303950E0, 2.319933E0, 2.335967E0, 2.351967E0, 2.368000E0, 2.384000E0
    ],
    "s": [
        5.105176E-1, 5.011446E-1, 4.928202E-1, 4.844958E-1, 4.769970E-1, 4.678253E-1,
        4.615620E-1, 4.555476E-1, 4.471324E-1, 4.414746E-1, 4.356670E-1, 4.285593E-1,
        4.231558E-1, 4.185217E-1, 4.131331E-1, 4.081239E-1, 4.035469E-1, 3.984264E-1,
        3.934247E-1, 3.891953E-1, 3.843461E-1, 3.795203E-1, 3.752892E-1, 3.706558E-1,
        3.656249E-1, 3.614718E-1, 3.568508E-1, 3.521560E-1, 3.479294E-1, 3.434012E-1,
        3.388012E-1, 3.347737E-1, 3.307572E-1, 3.260437E-1, 3.218747E-1, 3.181898E-1,
        3.137056E-1, 3.098022E-1, 3.061205E-1, 3.016691E-1, 2.974417E-1, 2.936899E-1,
        2.894116E-1, 2.856303E-1, 2.823693E-1, 2.780953E-1, 2.739300E-1, 2.703855E-1,
        2.665860E-1, 2.629008E-1, 2.592135E-1, 2.555521E-1, 2.516075E-1, 2.483881E-1,
        2.447093E-1, 2.409962E-1, 2.376348E-1, 2.340488E-1, 2.303501E-1, 2.272972E-1,
        2.239876E-1, 2.199111E-1, 2.172490E-1, 2.142481E-1, 2.102589E-1, 2.077176E-1,
        2.046465E-1, 2.008090E-1, 1.979125E-1, 1.952441E-1, 1.917278E-1, 1.890744E-1,
        1.862918E-1, 1.824383E-1, 1.800446E-1, 1.771112E-1, 1.736603E-1, 1.708303E-1,
        1.680069E-1, 1.646468E-1, 1.618968E-1, 1.590048E-1, 1.563349E-1, 1.531621E-1,
        1.506151E-1, 1.479173E-1, 1.445791E-1, 1.427122E-1, 1.395303E-1, 1.367070E-1,
        1.347262E-1, 1.320346E-1, 1.295722E-1, 1.271035E-1, 1.246385E-1, 1.228244E-1,
        1.196987E-1, 1.169780E-1, 1.148566E-1, 1.119661E-1, 1.101106E-1, 1.074526E-1,
        1.056357E-1, 1.031228E-1, 1.007753E-1, 9.869926E-2, 9.642449E-2, 9.384571E-2,
        9.204479E-2, 8.981515E-2, 8.731589E-2, 8.245081E-2, 8.046526E-2, 7.814279E-2,
        7.589826E-2, 7.388530E-2, 7.127179E-2, 6.944972E-2, 6.691301E-2, 6.491276E-2,
        6.226451E-2, 6.058581E-2, 5.781405E-2, 5.591250E-2, 5.357715E-2, 5.161283E-2,
        4.929902E-2, 4.733031E-2, 4.509376E-2, 4.269360E-2, 4.076105E-2, 3.857917E-2,
        3.616972E-2, 3.419726E-2, 3.205877E-2, 2.996479E-2, 2.804895E-2, 2.621410E-2,
        2.436291E-2, 2.310184E-2, 2.197352E-2, 2.115577E-2, 2.094988E-2, 2.065470E-2
    ],
    },
    {
        "name": "corrida1_esfera2",
        "sphere_id": "esfera2",
        "invert_sign": True,
        "terminal_window": [0.0, 1.353066],
        "t": [
        0.000000E0, 1.618922E-2, 3.227777E-2, 4.841666E-2, 6.452199E-2, 8.079509E-2,
        9.691720E-2, 1.130058E-1, 1.291279E-1, 1.451158E-1, 1.614224E-1, 1.775781E-1,
        1.934989E-1, 2.253908E-1, 2.415129E-1, 2.576686E-1, 2.738242E-1, 2.899631E-1,
        3.062194E-1, 3.222073E-1, 3.383127E-1, 3.544515E-1, 3.703723E-1, 3.865783E-1,
        4.026837E-1, 4.189064E-1, 4.349447E-1, 4.510500E-1, 4.672728E-1, 4.833278E-1,
        4.994331E-1, 5.155888E-1, 5.316941E-1, 5.478162E-1, 5.640893E-1, 5.799430E-1,
        5.961993E-1, 6.281080E-1, 6.442301E-1, 6.604025E-1, 6.764911E-1, 6.928984E-1,
        7.086514E-1, 7.248071E-1, 7.409628E-1, 7.570178E-1, 7.731399E-1, 7.893962E-1,
        8.055351E-1, 8.217075E-1, 8.376619E-1, 8.538511E-1, 8.699900E-1, 8.860953E-1,
        9.021839E-1, 9.181718E-1, 9.343106E-1, 9.506005E-1, 9.664207E-1, 9.827273E-1,
        9.988326E-1, 1.030758E0, 1.046897E0, 1.062935E0, 1.079225E0, 1.095196E0,
        1.111318E0, 1.127407E0, 1.143563E0, 1.159819E0, 1.175840E0, 1.191929E0,
        1.207950E0, 1.224089E0, 1.240278E0, 1.256484E0, 1.272489E0, 1.288678E0,
        1.304918E0, 1.320855E0, 1.337129E0, 1.353066E0
    ],
    "s": [
        4.934051E-1, 4.805911E-1, 4.691340E-1, 4.588897E-1, 4.477403E-1, 4.381337E-1,
        4.294085E-1, 4.197480E-1, 4.112086E-1, 4.032590E-1, 3.953805E-1, 3.877547E-1,
        3.805263E-1, 3.733766E-1, 3.657095E-1, 3.593786E-1, 3.520752E-1, 3.452437E-1,
        3.386751E-1, 3.323708E-1, 3.252708E-1, 3.191180E-1, 3.131670E-1, 3.062379E-1,
        3.005335E-1, 2.944229E-1, 2.878234E-1, 2.820744E-1, 2.783223E-1, 2.720451E-1,
        2.665052E-1, 2.606681E-1, 2.532627E-1, 2.459536E-1, 2.407516E-1, 2.365808E-1,
        2.315554E-1, 2.266954E-1, 2.214134E-1, 2.135726E-1, 2.086397E-1, 2.037193E-1,
        1.987171E-1, 1.940431E-1, 1.895275E-1, 1.843700E-1, 1.795456E-1, 1.747261E-1,
        1.698860E-1, 1.652312E-1, 1.606445E-1, 1.558349E-1, 1.512253E-1, 1.468606E-1,
        1.421525E-1, 1.378221E-1, 1.333801E-1, 1.289545E-1, 1.248426E-1, 1.213690E-1,
        1.171425E-1, 1.131842E-1, 1.092360E-1, 1.051483E-1, 1.012538E-1, 9.810312E-2,
        9.403553E-2, 9.018853E-2, 8.651903E-2, 8.284967E-2, 7.929962E-2, 7.564720E-2,
        7.162887E-2, 6.734225E-2, 6.401885E-2, 6.035173E-2, 5.652027E-2, 5.298823E-2,
        4.958797E-2, 4.566716E-2, 4.231930E-2, 3.866153E-2
    ]
},
    {
        "name": "corrida1_esfera3",
        "sphere_id": "esfera3",
        "invert_sign": True,
        "terminal_window": [0, 0.95],
        "t": [0.000000E0, 1.600000E-2, 3.208333E-2, 4.780000E-2, 6.405000E-2, 8.003333E-2,
        9.601667E-2, 1.120667E-1, 1.280833E-1, 1.469500E-1, 1.600000E-1, 1.917333E-1,
        2.076667E-1, 2.237833E-1, 2.397667E-1, 2.558000E-1, 2.718833E-1, 2.879000E-1,
        3.036167E-1, 3.196167E-1, 3.356000E-1, 3.516167E-1, 3.677167E-1, 3.837333E-1,
        3.997167E-1, 4.158500E-1, 4.318667E-1, 4.478500E-1, 4.638000E-1, 4.800167E-1,
        4.960000E-1, 5.119833E-1, 5.277500E-1, 5.438333E-1, 5.599667E-1, 5.915833E-1,
        6.075833E-1, 6.235833E-1, 6.396000E-1, 6.556000E-1, 6.716833E-1, 6.876667E-1,
        7.036333E-1, 7.196167E-1, 7.357667E-1, 7.517000E-1, 7.678833E-1, 7.838833E-1,
        7.999833E-1, 8.159500E-1, 8.318667E-1, 8.479000E-1, 8.639333E-1, 8.801167E-1,
        8.958833E-1, 9.119000E-1, 9.280167E-1, 9.438167E-1
    ],
    "s": [
        5.083670E-1, 4.989778E-1, 4.874783E-1, 4.728729E-1, 4.593106E-1, 4.439535E-1,
        4.304206E-1, 4.196995E-1, 4.102899E-1, 3.992746E-1, 3.895115E-1, 3.791753E-1,
        3.681621E-1, 3.584271E-1, 3.491783E-1, 3.398293E-1, 3.302579E-1, 3.215289E-1,
        3.124646E-1, 3.026959E-1, 2.959988E-1, 2.875908E-1, 2.785804E-1, 2.702634E-1,
        2.633671E-1, 2.543924E-1, 2.463630E-1, 2.381304E-1, 2.305767E-1, 2.229840E-1,
        2.154818E-1, 2.084141E-1, 2.003260E-1, 1.944204E-1, 1.866809E-1, 1.798299E-1,
        1.731926E-1, 1.668494E-1, 1.598418E-1, 1.534996E-1, 1.472586E-1, 1.403698E-1,
        1.343147E-1, 1.285279E-1, 1.227424E-1, 1.159915E-1, 1.110198E-1, 1.053042E-1,
        9.962835E-2, 9.456428E-2, 8.880522E-2, 8.310832E-2, 7.837979E-2, 7.285267E-2,
        6.731978E-2, 6.244255E-2, 5.736306E-2, 5.167513E-2
    ]
}
]

# ============ FUNCIONES ============

def clean_data(x):
    """Convierte a array y elimina NaN/Inf"""
    a = np.asarray(x, dtype=float)
    return a[np.isfinite(a)]

def window_mask(t, window):
    """Máscara para seleccionar ventana temporal"""
    if window is None:
        return np.ones_like(t, dtype=bool)
    tmin, tmax = window
    return (t >= tmin) & (t <= tmax)

def fit_terminal_velocity(t, s):
    """
    Ajuste lineal: s = s₀ + v_t·t
    Retorna: (s0, vt), (ds0, dvt), R2
    """
    # Ajuste por mínimos cuadrados
    A = np.column_stack([np.ones_like(t), t])
    params, residuals, rank, singular = np.linalg.lstsq(A, s, rcond=None)
    s0, vt = params
    
    # Errores
    s_pred = s0 + vt * t
    dof = max(len(t) - 2, 1)
    sigma2 = np.sum((s - s_pred)**2) / dof
    cov = sigma2 * np.linalg.pinv(A.T @ A)
    ds0, dvt = np.sqrt(np.diag(cov))
    
    # R²
    ss_tot = np.sum((s - np.mean(s))**2)
    ss_res = np.sum((s - s_pred)**2)
    R2 = 1 - ss_res / ss_tot
    
    return (s0, vt), (ds0, dvt), R2, s_pred

def calculate_viscosity(vt, r_m, rho_s, rho_f):
    """
    Viscosidad desde velocidad terminal:
    η = 2r²g(ρₛ - ρf)/(9v_t)
    """
    eta = (2 * r_m**2 * g * (rho_s - rho_f)) / (9 * abs(vt))
    return eta

def propagate_error_eta(vt, dvt, r, dr, rho_s, rho_f, drho_f):
    """Propagación de errores para η"""
    # η ∝ r²/vt
    rel_error_r = 2 * dr / r if r > 0 else 0
    rel_error_vt = dvt / abs(vt) if vt != 0 else 0
    rel_error_rho = drho_f / abs(rho_s - rho_f) if abs(rho_s - rho_f) > 0 else 0
    
    eta = calculate_viscosity(vt, r, rho_s, rho_f)
    deta = eta * np.sqrt(rel_error_r**2 + rel_error_vt**2 + rel_error_rho**2)
    
    return deta

# ============ PROCESAMIENTO ============

print("="*70)
print("ANÁLISIS DE CAÍDA EN MEDIO VISCOSO - LEY DE STOKES")
print("Método: Velocidad Terminal (Equilibrio de Fuerzas)")
print("="*70)

results = []
plot_data = []

for d in DATASETS:
    name = d["name"]
    sid = d["sphere_id"]
    sph = SPH_MAP[sid]
    
    r_mm = sph["radio_mm"]
    r_m = r_mm * 1e-3
    rho_s = sph["rho_s"]
    rho_f = RHO_FLUID
    
    print(f"\n{'─'*70}")
    print(f"Procesando: {name}")
    print(f"Esfera: {sid}, r = {r_mm:.2f} mm")
    print(f"{'─'*70}")
    
    # Preparar datos
    t = clean_data(d["t"])
    s = clean_data(d["s"])
    
    if d.get("invert_sign", False):
        s = -s
        s_all = -clean_data(d["s"])
    else:
        s_all = s.copy()
    
    # Seleccionar región terminal
    mask_term = window_mask(t, d.get("terminal_window"))
    t_term = t[mask_term]
    s_term = s[mask_term]
    
    print(f"Puntos totales: {len(t)}")
    print(f"Puntos en región terminal: {len(t_term)}")
    print(f"Ventana temporal: [{t_term[0]:.2f}, {t_term[-1]:.2f}] s")
    
    # Ajuste lineal
    (s0, vt), (ds0, dvt), R2, s_pred = fit_terminal_velocity(t_term, s_term)
    
    print(f"\nAjuste lineal: s = s₀ + v_t·t")
    print(f"  s₀ = {s0:.6f} ± {ds0:.6f} m")
    print(f"  v_terminal = {abs(vt):.6f} ± {dvt:.6f} m/s")
    print(f"  R² = {R2:.6f}")
    
    if R2 < 0.95:
        print(f"  ⚠️  ADVERTENCIA: R² < 0.95, ajustar ventana terminal_window")
    
    # Calcular viscosidad
    eta = calculate_viscosity(vt, r_m, rho_s, rho_f)
    deta = propagate_error_eta(vt, dvt, r_m, 0.01e-3, rho_s, rho_f, 5.0)
    
    print(f"\nViscosidad calculada:")
    print(f"  η = {eta:.4f} ± {deta:.4f} Pa·s")
    
    # Número de Reynolds
    Re = (rho_f * abs(vt) * 2 * r_m) / eta
    print(f"  Número de Reynolds: Re = {Re:.6f}")
    if Re < 1:
        print(f"  ✓ Régimen laminar (Re << 1)")
    else:
        print(f"  ⚠️  Re ≥ 1: Stokes puede no ser válido")
    
    # Guardar resultados
    results.append({
        "name": name,
        "sphere_id": sid,
        "r_mm": r_mm,
        "r_m": r_m,
        "rho_s": rho_s,
        "rho_f": rho_f,
        "s0": s0,
        "ds0": ds0,
        "vt": abs(vt),
        "dvt": dvt,
        "eta": eta,
        "deta": deta,
        "Re": Re,
        "R2": R2,
        "n_points_terminal": len(t_term)
    })
    
    plot_data.append({
        "name": name,
        "sphere_id": sid,
        "t_all": t,
        "s_all": s_all,
        "t_term": t_term,
        "s_term": s_term,
        "s_pred": s_pred,
        "vt": abs(vt),
        "R2": R2
    })

# DataFrame de resultados
df = pd.DataFrame(results)

# Promedios por esfera
grp = df.groupby("sphere_id", as_index=False).agg(
    r_mm=("r_mm", "first"),
    r_m=("r_m", "first"),
    n=("name", "count"),
    vt_avg=("vt", "mean"),
    vt_std=("vt", "std"),
    eta_avg=("eta", "mean"),
    eta_std=("eta", "std"),
    Re_avg=("Re", "mean"),
    R2_avg=("R2", "mean")
)

# Exportar
df.to_csv("resultados_velocidad_terminal.csv", index=False)
grp.to_csv("resumen_por_esfera.csv", index=False)

# ============ RESUMEN ============

print("\n" + "="*70)
print("RESUMEN DE RESULTADOS")
print("="*70)

print(f"\n{'Corrida':<20} {'Esfera':<10} {'v_t(m/s)':<12} {'η(Pa·s)':<10} {'Re':<10} {'R²':<8}")
print("─"*70)
for _, row in df.iterrows():
    print(f"{row['name']:<20} {row['sphere_id']:<10} {row['vt']:<12.6f} "
          f"{row['eta']:<10.4f} {row['Re']:<10.6f} {row['R2']:<8.4f}")

print("\n" + "="*70)
print("PROMEDIOS POR ESFERA")
print("="*70)
print(f"\n{'Esfera':<10} {'n':<5} {'r(mm)':<8} {'v_t(m/s)':<15} {'η(Pa·s)':<15} {'Re':<10}")
print("─"*70)
for _, row in grp.iterrows():
    vt_str = f"{row['vt_avg']:.6f}" if row['n'] == 1 else f"{row['vt_avg']:.6f}±{row['vt_std']:.6f}"
    eta_str = f"{row['eta_avg']:.4f}" if row['n'] == 1 else f"{row['eta_avg']:.4f}±{row['eta_std']:.4f}"
    print(f"{row['sphere_id']:<10} {int(row['n']):<5} {row['r_mm']:<8.2f} {vt_str:<15} {eta_str:<15} {row['Re_avg']:<10.6f}")

# Viscosidad promedio global
eta_global = grp['eta_avg'].mean()
eta_global_std = grp['eta_avg'].std() if len(grp) > 1 else 0

print(f"\n{'='*70}")
print(f"VISCOSIDAD PROMEDIO GLOBAL: η = {eta_global:.4f} ± {eta_global_std:.4f} Pa·s")
print(f"Valor referencia (glicerina pura a 20°C): 1.41 Pa·s")
print(f"{'='*70}")

# Verificar v_t ∝ r²
if len(grp) > 1:
    print(f"\nVERIFICACIÓN DE LEY DE STOKES (v_t ∝ r²):")
    print("─"*70)
    grp_sorted = grp.sort_values('r_m')
    for i in range(len(grp_sorted) - 1):
        r1 = grp_sorted.iloc[i]['r_m']
        r2 = grp_sorted.iloc[i+1]['r_m']
        vt1 = grp_sorted.iloc[i]['vt_avg']
        vt2 = grp_sorted.iloc[i+1]['vt_avg']
        ratio_r2 = (r2/r1)**2
        ratio_vt = vt2/vt1
        diff = abs(ratio_r2 - ratio_vt) / ratio_r2 * 100
        print(f"  {grp_sorted.iloc[i]['sphere_id']} → {grp_sorted.iloc[i+1]['sphere_id']}: "
              f"(r₂/r₁)² = {ratio_r2:.3f}, v_t₂/v_t₁ = {ratio_vt:.3f} "
              f"(diferencia: {diff:.1f}%)")

# ============ GRÁFICOS ============

print("\n" + "="*70)
print("GENERANDO GRÁFICOS...")
print("="*70)

colors = ['#e74c3c', '#3498db', '#2ecc71']
n_datasets = len(plot_data)

# Figura principal
fig = plt.figure(figsize=(15, 5*n_datasets))
gs = GridSpec(n_datasets, 2, figure=fig, hspace=0.3, wspace=0.3)

for i, pdata in enumerate(plot_data):
    color = colors[i % len(colors)]
    
    # Panel 1: Posición vs Tiempo
    ax1 = fig.add_subplot(gs[i, 0])
    ax1.scatter(pdata['t_all'], pdata['s_all'], alpha=0.4, s=25, 
                color=color, label='Datos experimentales')
    ax1.plot(pdata['t_term'], pdata['s_pred'], 'k-', linewidth=2.5, 
             label=f"Ajuste lineal (R²={pdata['R2']:.4f})")
    ax1.axvline(pdata['t_term'][0], color='gray', linestyle='--', 
                alpha=0.6, linewidth=1.5, label='Región terminal')
    ax1.axvline(pdata['t_term'][-1], color='gray', linestyle='--', 
                alpha=0.6, linewidth=1.5)
    ax1.set_xlabel('Tiempo (s)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Posición (m)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{pdata["name"]} - {pdata["sphere_id"]}\nv_terminal = {pdata["vt"]*100:.4f} cm/s', 
                  fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Velocidad instantánea (aproximada por diferencias finitas)
    ax2 = fig.add_subplot(gs[i, 1])
    t_all = pdata['t_all']
    s_all = pdata['s_all']
    
    # Calcular velocidad por diferencias centrales
    v_inst = np.zeros_like(t_all)
    for j in range(1, len(t_all)-1):
        dt = (t_all[j+1] - t_all[j-1])
        ds = (s_all[j+1] - s_all[j-1])
        v_inst[j] = abs(ds / dt) if dt > 0 else 0
    v_inst[0] = v_inst[1]
    v_inst[-1] = v_inst[-2]
    
    # Suavizado simple (promedio móvil)
    window_size = 5
    v_smooth = np.convolve(v_inst, np.ones(window_size)/window_size, mode='same')
    
    ax2.plot(t_all, v_smooth*100, color=color, linewidth=2, alpha=0.7, 
             label='Velocidad instantánea (suavizada)')
    ax2.axhline(pdata['vt']*100, color='k', linestyle='--', linewidth=2, 
                label=f'v_terminal = {pdata["vt"]*100:.4f} cm/s')
    ax2.axvspan(pdata['t_term'][0], pdata['t_term'][-1], alpha=0.2, 
                color='gray', label='Región de ajuste')
    ax2.set_xlabel('Tiempo (s)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Velocidad (cm/s)', fontsize=12, fontweight='bold')
    ax2.set_title('Velocidad vs Tiempo', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(bottom=0)

plt.suptitle('Análisis de Velocidad Terminal - Ley de Stokes', 
             fontsize=15, fontweight='bold', y=0.995)
plt.savefig('analisis_velocidad_terminal.png', dpi=150, bbox_inches='tight')
print("✓ Gráfico guardado: analisis_velocidad_terminal.png")

# Gráfico comparativo
if len(grp) > 1:
    fig2, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Gráfico 1: v_t vs r²
    ax1 = axes[0]
    r_mm = grp['r_mm'].values
    vt_avg = grp['vt_avg'].values * 100  # a cm/s
    
    ax1.scatter(r_mm**2, vt_avg, s=200, color='#e74c3c', 
                edgecolor='black', linewidth=2, zorder=3)
    
    # Ajuste lineal
    if len(r_mm) > 1:
        coeffs = np.polyfit(r_mm**2, vt_avg, 1)
        r2_fit = np.linspace(r_mm.min()**2, r_mm.max()**2, 100)
        ax1.plot(r2_fit, coeffs[0]*r2_fit + coeffs[1], 'k--', linewidth=2.5, 
                 alpha=0.7, label=f'v_t = {coeffs[0]:.3f}·r² + {coeffs[1]:.3f}')
    
    for i, row in grp.iterrows():
        ax1.annotate(row['sphere_id'], 
                    (row['r_mm']**2, row['vt_avg']*100),
                    xytext=(8, 8), textcoords='offset points', 
                    fontsize=11, fontweight='bold')
    
    ax1.set_xlabel('Radio² (mm²)', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Velocidad terminal (cm/s)', fontsize=13, fontweight='bold')
    ax1.set_title('Verificación: v_t ∝ r²', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Gráfico 2: Viscosidad por esfera
    ax2 = axes[1]
    bar_colors = colors[:len(grp)]
    x_pos = np.arange(len(grp))
    
    ax2.bar(x_pos, grp['eta_avg'], 
            yerr=grp['eta_std'] if len(grp) > 1 and grp['eta_std'].notna().any() else None,
            color=bar_colors, edgecolor='black', linewidth=2, alpha=0.8, capsize=5)
    ax2.axhline(eta_global, color='k', linestyle='--', linewidth=2.5, 
                label=f'Promedio: {eta_global:.3f} Pa·s')
    ax2.axhline(1.41, color='gray', linestyle=':', linewidth=2, alpha=0.7, 
                label='Glicerina pura (ref.)')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(grp['sphere_id'], fontsize=12, fontweight='bold')
    ax2.set_ylabel('Viscosidad η (Pa·s)', fontsize=13, fontweight='bold')
    ax2.set_title('Viscosidad por Esfera', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Gráfico 3: Número de Reynolds
    ax3 = axes[2]
    ax3.bar(x_pos, grp['Re_avg'], color=bar_colors, 
            edgecolor='black', linewidth=2, alpha=0.8)
    ax3.axhline(1.0, color='r', linestyle='--', linewidth=2, 
                label='Re = 1 (límite régimen laminar)')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(grp['sphere_id'], fontsize=12, fontweight='bold')
    ax3.set_ylabel('Número de Reynolds', fontsize=13, fontweight='bold')
    ax3.set_title('Verificación Régimen Laminar', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig('comparacion_esferas.png', dpi=150, bbox_inches='tight')
    print("✓ Gráfico guardado: comparacion_esferas.png")

print("\n" + "="*70)
print("ANÁLISIS COMPLETADO")
print("="*70)
print("\nArchivos generados:")
print("  • resultados_velocidad_terminal.csv")
print("  • resumen_por_esfera.csv")
print("  • analisis_velocidad_terminal.png")
if len(grp) > 1:
    print("  • comparacion_esferas.png")

print("\n📋 INSTRUCCIONES PARA MEJORAR RESULTADOS:")
print("  1. Si R² < 0.95: Ajusta 'terminal_window' en DATASETS")
print("  2. Región terminal = zona donde velocidad es constante")
print("  3. Excluye inicio (aceleración) y final (fondo de probeta)")
print("  4. MIDE los radios reales con calibre")
print("  5. VERIFICA densidad del fluido con densitómetro")
print("  6. Temperatura afecta viscosidad (~2-3% por °C)")

print("\n✅ CRITERIOS DE VALIDACIÓN:")
print("  • R² > 0.95 en todos los ajustes")
print("  • Re < 1 para todas las esferas (régimen laminar)")
print("  • η entre 0.8-1.5 Pa·s (glicerina comercial)")
print("  • v_t₂/v_t₁ ≈ (r₂/r₁)² con diferencia < 15%")

print("\n" + "="*70)