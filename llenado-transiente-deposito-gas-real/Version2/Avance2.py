import numpy as np
import matplotlib.pyplot as plt
import shutil

from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from matplotlib.animation import (
    FuncAnimation,
    PillowWriter,
    FFMpegWriter,
)
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D


if "__file__" in globals():
    DIRECTORIO_SALIDA = Path(__file__).resolve().parent
else:
    DIRECTORIO_SALIDA = Path.cwd()

OUTPUT_DIR = DIRECTORIO_SALIDA

MP4_PATH = OUTPUT_DIR / "flujovariable_video_sim.mp4"
PNG_PATH = OUTPUT_DIR / "flujovariable_foto_sim.png"
TXT_PATH = OUTPUT_DIR / "flujovariable_reporte.txt"
CSV_PATH = OUTPUT_DIR / "flujovariable_datos.csv"
GIF_PATH = OUTPUT_DIR / "ICF130_Tema7_flujo_variable.gif"


FPS = 30
N_FRAMES = 240

SAVE_MP4 = True
SAVE_GIF = False
SAVE_PNG = True
SAVE_CSV = True
SAVE_TXT = True


V = 0.050

P0 = 1.0e5
T0 = 300.0

P_supply = 200.0e5
T_supply = 300.0

P_target = 190.0e5

mdot_reference = 0.020
Cd = 0.80


if P_target >= P_supply:
    raise ValueError(
        "P_target debe ser menor que P_supply. "
        "Esta simulación se detiene antes de igualar "
        "la presión del suministro."
    )

if P0 >= P_target:
    raise ValueError(
        "P0 debe ser menor que P_target."
    )

if V <= 0:
    raise ValueError(
        "El volumen debe ser positivo."
    )


M = 28.0134e-3
R_u = 8.314462618
R = R_u / M

SH_A = 28.98641
SH_B = 1.853978
SH_C = -9.647459
SH_D = 16.63537
SH_E = 0.000117


def cp_n2_shomate(T):
    tau = T / 1000.0

    cp_molar = (
        SH_A
        + SH_B * tau
        + SH_C * tau**2
        + SH_D * tau**3
        + SH_E / tau**2
    )

    return cp_molar / M


cp = cp_n2_shomate(T_supply)
cv = cp - R
gamma = cp / cv

Tc = 126.19
Pc = 33.978e5

a = (27.0 / 64.0) * R**2 * Tc**2 / Pc
b = R * Tc / (8.0 * Pc)


def pressure_ideal(m, T):
    return m * R * T / V


def pressure_vdw_from_v(T, v):
    return R * T / (v - b) - a / v**2


def pressure_vdw(m, T):
    v = V / m

    if np.any(v <= b):
        raise ValueError(
            "Estado Van der Waals no físico: v <= b."
        )

    return pressure_vdw_from_v(T, v)


def specific_volume_vdw(P, T):
    def residual(v):
        return pressure_vdw_from_v(T, v) - P

    v_min = b * (1.0 + 1e-10)
    v_ideal = R * T / P

    v_max = max(
        100.0 * v_ideal,
        1.0e4 * b,
    )

    while residual(v_max) > 0:
        v_max *= 2.0

    return brentq(
        residual,
        v_min,
        v_max,
        xtol=1e-12,
        rtol=1e-12,
    )


def u_ideal(T):
    return cv * T


def u_vdw(T, v):
    return cv * T - a / v


v_supply_ideal = R * T_supply / P_supply
h0_ideal = cp * T_supply

v_supply_vdw = specific_volume_vdw(
    P_supply,
    T_supply,
)

u_supply_vdw = u_vdw(
    T_supply,
    v_supply_vdw,
)

h0_vdw = (
    u_supply_vdw
    + P_supply * v_supply_vdw
)


critical_pressure_ratio = (
    2.0 / (gamma + 1.0)
)**(
    gamma / (gamma - 1.0)
)

P_critical_tank = (
    critical_pressure_ratio * P_supply
)

choked_factor = (
    np.sqrt(
        gamma / (R * T_supply)
    )
    * (
        2.0 / (gamma + 1.0)
    )**(
        (gamma + 1.0)
        / (2.0 * (gamma - 1.0))
    )
)

A_orifice = (
    mdot_reference
    / (Cd * P_supply * choked_factor)
)

d_orifice = np.sqrt(
    4.0 * A_orifice / np.pi
)


def mass_flow_rate(P_downstream):
    if P_downstream >= P_supply:
        return 0.0

    ratio = max(
        P_downstream / P_supply,
        0.0,
    )

    if ratio <= critical_pressure_ratio:
        mdot = (
            Cd
            * A_orifice
            * P_supply
            * np.sqrt(
                gamma / (R * T_supply)
            )
            * (
                2.0 / (gamma + 1.0)
            )**(
                (gamma + 1.0)
                / (2.0 * (gamma - 1.0))
            )
        )

    else:
        pressure_term = (
            ratio**(2.0 / gamma)
            - ratio**(
                (gamma + 1.0) / gamma
            )
        )

        pressure_term = max(
            pressure_term,
            0.0,
        )

        mdot = (
            Cd
            * A_orifice
            * P_supply
            / np.sqrt(R * T_supply)
            * np.sqrt(
                (
                    2.0 * gamma
                    / (gamma - 1.0)
                )
                * pressure_term
            )
        )

    return mdot


def flow_regime(P_downstream):
    if P_downstream >= P_supply:
        return "SIN FLUJO"

    ratio = P_downstream / P_supply

    if ratio <= critical_pressure_ratio:
        return "ESTRANGULADO"

    return "SUBSONICO"


m0_ideal = P0 * V / (R * T0)

v0_vdw = specific_volume_vdw(P0, T0)
m0_vdw = V / v0_vdw


def rhs_ideal(t, y):
    m, T = y

    P_tank = pressure_ideal(m, T)
    mdot_now = mass_flow_rate(P_tank)

    dm_dt = mdot_now

    dT_dt = (
        mdot_now / m
        * (gamma * T_supply - T)
    )

    return [dm_dt, dT_dt]


def rhs_vdw(t, y):
    m, T = y

    v = V / m

    P_tank = pressure_vdw(m, T)
    mdot_now = mass_flow_rate(P_tank)

    u = u_vdw(T, v)
    dP_dT_v = R / (v - b)

    thermodynamic_term = (
        T * dP_dT_v - P_tank
    )

    dm_dt = mdot_now

    dT_dt = (
        mdot_now / (m * cv)
        * (
            h0_vdw
            - u
            + v * thermodynamic_term
        )
    )

    return [dm_dt, dT_dt]


def event_ideal(t, y):
    m, T = y
    return pressure_ideal(m, T) - P_target


event_ideal.terminal = True
event_ideal.direction = 1


def event_vdw(t, y):
    m, T = y
    return pressure_vdw(m, T) - P_target


event_vdw.terminal = True
event_vdw.direction = 1


T_MAX = 5000.0

print()
print("Integrando modelo de gas ideal...")

sol_ideal = solve_ivp(
    rhs_ideal,
    (0.0, T_MAX),
    (m0_ideal, T0),
    method="RK45",
    events=event_ideal,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.5,
)

print("Integrando modelo Van der Waals...")

sol_vdw = solve_ivp(
    rhs_vdw,
    (0.0, T_MAX),
    (m0_vdw, T0),
    method="RK45",
    events=event_vdw,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.5,
)

if len(sol_ideal.t_events[0]) == 0:
    raise RuntimeError(
        "El modelo ideal no alcanzó P_target."
    )

if len(sol_vdw.t_events[0]) == 0:
    raise RuntimeError(
        "Van der Waals no alcanzó P_target."
    )

tf_ideal = sol_ideal.t_events[0][0]
tf_vdw = sol_vdw.t_events[0][0]


def transition_time(
    solution,
    pressure_function,
    tf,
):
    if P_target <= P_critical_tank:
        return None

    def residual(t):
        m, T = solution.sol(t)

        return (
            pressure_function(m, T)
            - P_critical_tank
        )

    return brentq(
        residual,
        0.0,
        tf,
    )


tcrit_ideal = transition_time(
    sol_ideal,
    pressure_ideal,
    tf_ideal,
)

tcrit_vdw = transition_time(
    sol_vdw,
    pressure_vdw,
    tf_vdw,
)


tf_global = max(
    tf_ideal,
    tf_vdw,
)

time = np.linspace(
    0,
    tf_global,
    N_FRAMES,
)

time_I = np.minimum(
    time,
    tf_ideal,
)

time_V = np.minimum(
    time,
    tf_vdw,
)

mI, TI = sol_ideal.sol(time_I)
mV, TV = sol_vdw.sol(time_V)

PI = pressure_ideal(mI, TI)

PV = np.array([
    pressure_vdw(m, T)
    for m, T in zip(mV, TV)
])

mdot_I = np.array([
    mass_flow_rate(P)
    for P in PI
])

mdot_V = np.array([
    mass_flow_rate(P)
    for P in PV
])

TI_C = TI - 273.15
TV_C = TV - 273.15

PI_bar = PI / 1e5
PV_bar = PV / 1e5

vI = V / mI
vV = V / mV

ZI = PI * vI / (R * TI)
ZV = PV * vV / (R * TV)


mI_final, TI_final = sol_ideal.sol(tf_ideal)
mV_final, TV_final = sol_vdw.sol(tf_vdw)

PI_final = pressure_ideal(
    mI_final,
    TI_final,
)

PV_final = pressure_vdw(
    mV_final,
    TV_final,
)

mdot_I_final = mass_flow_rate(PI_final)
mdot_V_final = mass_flow_rate(PV_final)

ZI_final = (
    PI_final
    * (V / mI_final)
    / (R * TI_final)
)

ZV_final = (
    PV_final
    * (V / mV_final)
    / (R * TV_final)
)


T_analytical = (
    gamma * T_supply
    + (m0_ideal / mI)
    * (T0 - gamma * T_supply)
)

error_T_ideal = np.max(
    np.abs(TI - T_analytical)
)


def energy_validation_ideal():
    t_val = np.linspace(
        0,
        tf_ideal,
        1200,
    )

    m, T = sol_ideal.sol(t_val)

    U = m * u_ideal(T)
    U0 = m0_ideal * u_ideal(T0)

    Ein = (
        (m - m0_ideal)
        * h0_ideal
    )

    residual = U - U0 - Ein

    relative_error = np.max(
        np.abs(residual)
        / np.maximum(np.abs(U), 1.0)
    )

    return relative_error


def energy_validation_vdw():
    t_val = np.linspace(
        0,
        tf_vdw,
        1200,
    )

    m, T = sol_vdw.sol(t_val)

    v = V / m

    U = m * u_vdw(T, v)

    U0 = (
        m0_vdw
        * u_vdw(T0, v0_vdw)
    )

    Ein = (
        (m - m0_vdw)
        * h0_vdw
    )

    residual = U - U0 - Ein

    relative_error = np.max(
        np.abs(residual)
        / np.maximum(np.abs(U), 1.0)
    )

    return relative_error


energy_error_I = energy_validation_ideal()
energy_error_V = energy_validation_vdw()

mdot_avg_I = (
    (mI_final - m0_ideal)
    / tf_ideal
)

mdot_avg_V = (
    (mV_final - m0_vdw)
    / tf_vdw
)


if SAVE_CSV:
    csv_data = np.column_stack([
        time,
        PI_bar,
        PV_bar,
        TI_C,
        TV_C,
        mI,
        mV,
        mdot_I,
        mdot_V,
        ZI,
        ZV,
    ])

    csv_header = (
        "t_s,"
        "P_ideal_bar,"
        "P_vdw_bar,"
        "T_ideal_C,"
        "T_vdw_C,"
        "m_ideal_kg,"
        "m_vdw_kg,"
        "mdot_ideal_kg_s,"
        "mdot_vdw_kg_s,"
        "Z_ideal,"
        "Z_vdw"
    )

    np.savetxt(
        CSV_PATH,
        csv_data,
        delimiter=",",
        header=csv_header,
        comments="",
        fmt="%.10e",
    )

    print(f"CSV generado: {CSV_PATH}")


def fmt_time(value):
    if value is None:
        return "No se alcanza en esta simulación"

    return f"{value:.3f} s"


def create_report():
    report = f"""ICF130 — PROYECTO COMPUTACIONAL
VERSIÓN 2 — LLENADO TRANSIENTE CON CAUDAL VARIABLE
GAS IDEAL VS VAN DER WAALS
======================================================================

1. FENÓMENO FÍSICO Y OBJETIVO
----------------------------------------------------------------------

Se estudia el llenado transiente de un depósito rígido inicialmente
presurizado con nitrógeno.

Entran masa y energía, aumentando la masa almacenada y modificándose
la temperatura promedio y la presión.

Se comparan dos cierres termodinámicos:

    1) Gas ideal.
    2) Van der Waals.

El cambio principal respecto de la Versión 1 es que el caudal másico
ya no se impone constante durante todo el llenado.

Ahora se calcula mediante una ley de flujo compresible a través de un
orificio, distinguiendo flujo estrangulado y subsónico.

Esta versión continúa siendo adiabática y concentrada, de tipo 0D.
No incorpora transferencia de calor al ambiente.


2. CONDICIONES Y PARÁMETROS DEL PROBLEMA
----------------------------------------------------------------------

Gas: nitrógeno.

Volumen:
    V = {V:.12g} m³ = {V * 1000:.12g} L

Estado inicial:
    P0 = {P0:.12g} Pa = {P0 / 1e5:.12g} bar
    T0 = {T0:.12g} K

Suministro:
    P_supply = {P_supply:.12g} Pa
             = {P_supply / 1e5:.12g} bar
    T_supply = {T_supply:.12g} K

Presión objetivo:
    P_target = {P_target:.12g} Pa
             = {P_target / 1e5:.12g} bar

Caudal de referencia para calibración:
    mdot_reference = {mdot_reference:.12g} kg/s

Coeficiente de descarga:
    Cd = {Cd:.12g}

Las ecuaciones de estado y las razones de presión requieren
presiones absolutas.


3. PROPIEDADES DEL NITRÓGENO
----------------------------------------------------------------------

Masa molar:
    M = {M:.12g} kg/mol

Constante universal:
    R_u = {R_u:.12g} J/(mol K)

Constante específica:
    R = R_u/M = {R:.12g} J/(kg K)

Parámetros críticos utilizados:
    Tc = {Tc:.12g} K
    Pc = {Pc:.12g} Pa

Parámetros específicos de Van der Waals:

    a = (27/64)*R²*Tc²/Pc
      = {a:.12g} Pa·m⁶/kg²

    b = R*Tc/(8*Pc)
      = {b:.12g} m³/kg

El volumen específico se expresa en m³/kg. Por tanto, a y b son
parámetros específicos por unidad de masa, no parámetros molares.


4. CAPACIDADES CALORÍFICAS Y CORRELACIÓN SHOMATE
----------------------------------------------------------------------

cp se obtiene evaluando una sola vez la correlación:

    tau = T/1000

    cp_molar = A + B*tau + C*tau² + D*tau³ + E/tau²

    cp = cp_molar/M

Coeficientes:

    A = {SH_A}
    B = {SH_B}
    C = {SH_C}
    D = {SH_D}
    E = {SH_E}

El rango indicado en el código original para estos coeficientes
es de 100 a 500 K.

La evaluación se realiza en T_supply:

    cp = {cp:.12g} J/(kg K)

Luego se define:

    cv = cp - R = {cv:.12g} J/(kg K)

    gamma = cp/cv = {gamma:.12g}

Aunque existe una función cp_n2_shomate(T), el programa no vuelve
a evaluar cp con la temperatura instantánea del depósito.

Esta versión sigue usando propiedades caloríficas constantes.

En Van der Waals se mantiene cv constante y la entalpía se obtiene
mediante h = u + Pv. El cp anterior no se impone como capacidad
calorífica real constante del modelo Van der Waals.


5. HIPÓTESIS Y SIMPLIFICACIONES
----------------------------------------------------------------------

5.1. Depósito

Volumen de control fijo y rígido.
Una entrada y ninguna salida durante la etapa de llenado.

5.2. Intercambios energéticos

Se supone:

    - Qdot = 0.
    - Trabajo de eje nulo.
    - Trabajo por cambio de volumen nulo.
    - Energía cinética y potencial almacenadas despreciables.

La energía de la corriente se contabiliza mediante la entalpía total
h0 del suministro, supuesto prácticamente en reposo.

El trabajo de flujo está incluido en h0 y no se añade nuevamente.

5.3. Mezcla uniforme

Se considera una sola temperatura T(t), una sola presión P(t)
y una masa total m(t).

No se resuelven velocidades, turbulencia, estratificación ni
gradientes espaciales de temperatura.

5.4. Pared y ambiente

No se calcula temperatura de pared ni almacenamiento térmico en ella.
No existe transferencia de calor con el ambiente.

5.5. Suministro

P_supply y T_supply permanecen constantes.
No se modela agotamiento ni enfriamiento del reservorio.

5.6. Orificio

Área y Cd constantes.
Sin apertura dinámica, inercia de válvula ni pérdidas detalladas
en tuberías o accesorios.

No se modela flujo inverso: si P_tank >= P_supply se devuelve mdot=0.

5.7. Fases

No se consideran cambios de fase ni equilibrio líquido-vapor.

5.8. Modelo híbrido

La ley de caudal de gas ideal se utiliza en ambos modelos.

El depósito Van der Waals no utiliza una solución de flujo compresible
de gas real a través del orificio.


6. RAZÓN CRÍTICA Y REGÍMENES DE FLUJO
----------------------------------------------------------------------

Se define:

    r = P_tank/P_supply

La razón crítica utilizada es:

    r_crit = (2/(gamma+1))^(gamma/(gamma-1))
           = {critical_pressure_ratio:.12g}

La presión crítica del tanque es:

    P_critical_tank = r_crit*P_supply
                    = {P_critical_tank / 1e5:.12g} bar

Clasificación:

    P_tank >= P_supply:
        SIN FLUJO

    r <= r_crit y P_tank < P_supply:
        ESTRANGULADO

    r_crit < r < 1:
        SUBSÓNICO

En régimen estrangulado, con suministro y área constantes,
la fórmula produce un caudal constante.

En régimen subsónico, el caudal disminuye conforme aumenta
la contrapresión del depósito.

Estas relaciones corresponden a una aproximación de gas ideal
calóricamente perfecto para flujo compresible cuasiestacionario,
multiplicada por el coeficiente de descarga Cd.

No se calculan perfiles de velocidad ni número de Mach.
No se introduce Z ni una velocidad del sonido real en el orificio.


7. ECUACIONES DE CAUDAL
----------------------------------------------------------------------

FLUJO ESTRANGULADO:

    mdot =
        Cd*A_orifice*P_supply
        *sqrt(gamma/(R*T_supply))
        *(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

FLUJO SUBSÓNICO:

    mdot =
        Cd*A_orifice*P_supply/sqrt(R*T_supply)
        *sqrt(
            2*gamma/(gamma-1)
            *(
                r^(2/gamma)
                - r^((gamma+1)/gamma)
            )
        )

SIN DIFERENCIA FAVORABLE:

    Si P_tank >= P_supply:
        mdot = 0

Protecciones numéricas:

    r = max(P_tank/P_supply, 0)

    El término dentro de la raíz de la rama subsónica se limita
    inferiormente a cero.


8. CALIBRACIÓN DEL ORIFICIO
----------------------------------------------------------------------

Se define:

    F = sqrt(gamma/(R*T_supply))
        *(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

Luego:

    A_orifice = mdot_reference/(Cd*P_supply*F)

    d_orifice = sqrt(4*A_orifice/pi)

Resultados:

    A_orifice = {A_orifice:.12g} m²

    d_orifice = {d_orifice * 1000:.12g} mm

A_orifice es el área geométrica equivalente de la fórmula.
Cd*A_orifice es el área ponderada por el coeficiente de descarga.

La calibración fija el caudal de la rama estrangulada.
Coincide con el caudal inicial deseado si el estado inicial
pertenece a esa rama.

Estado inicial de esta ejecución:

    Régimen = {flow_regime(P0)}

    Caudal calculado = {mass_flow_rate(P0):.12g} kg/s

Ambos modelos comparten la misma función mdot(P_tank).
A igual presión predicen el mismo caudal.

A igual tiempo pueden tener caudales distintos porque sus presiones
pueden evolucionar de manera diferente.


9. BALANCES DE MASA Y ENERGÍA
----------------------------------------------------------------------

Durante el llenado activo:

    dm/dt = mdot(P_tank)

    dU/dt = mdot(P_tank)*h0_supply

    U = m*u

La masa ya no crece linealmente durante todo el llenado.
En la rama estrangulada sí crece linealmente con suministro fijo.

Como h0_supply es constante para cada modelo:

    U - U0 = (m - m0)*h0_supply

Esta igualdad sigue siendo válida aunque el caudal cambie con el tiempo.

El aislamiento térmico no obliga a mantener temperatura constante:
la corriente entrante transporta energía mediante su entalpía.


10. MODELO DE GAS IDEAL
----------------------------------------------------------------------

Ecuación de estado:

    P = m*R*T/V

    v = V/m

Energía interna:

    u = cv*T

Estado de suministro:

    v_supply_ideal = {v_supply_ideal:.12g} m³/kg

    h0_ideal = cp*T_supply
             = {h0_ideal:.12g} J/kg

Masa inicial:

    m0_ideal = P0*V/(R*T0)
             = {m0_ideal:.12g} kg

Ecuaciones diferenciales:

    dm/dt = mdot(P)

    dT/dt = (mdot(P)/m)*(gamma*T_supply - T)

Mientras existe entrada puede dividirse por dm/dt:

    dT/dm = (gamma*T_supply - T)/m

Solución analítica:

    T(m) = gamma*T_supply
           + (m0_ideal/m)*(T0 - gamma*T_supply)

Esta es una solución de T en función de m, no una expresión explícita
de T(t) para caudal variable. El programa obtiene m(t) numéricamente.

El límite formal para masa ilimitada sería:

    gamma*T_supply = {gamma * T_supply:.12g} K

No implica que esa temperatura se alcance en este llenado con
suministro de presión finita y detención por evento.


11. MODELO DE VAN DER WAALS
----------------------------------------------------------------------

Ecuación de estado:

    P = R*T/(v-b) - a/v²

    v = V/m

Dominio admitido:

    v > b

Energía interna:

    u = cv*T - a/v

    U = m*cv*T - a*m²/V

Identidad termodinámica:

    du = cv*dT + [T*(dP/dT)_v - P]*dv

Para Van der Waals:

    (dP/dT)_v = R/(v-b)

    T*(dP/dT)_v - P = a/v²

Con depósito rígido:

    dv/dt = -(v/m)*mdot

Ecuación de temperatura implementada:

    dT/dt =
        mdot/(m*cv)
        *[h0_vdw - u + v*(T*(dP/dT)_v - P)]

Forma equivalente:

    dT/dt =
        mdot/(m*cv)
        *[h0_vdw - cv*T + 2*a/v]

Condiciones iniciales:

    v0_vdw = {v0_vdw:.12g} m³/kg

    m0_vdw = V/v0_vdw = {m0_vdw:.12g} kg

Estado del suministro:

    v_supply_vdw = {v_supply_vdw:.12g} m³/kg

    u_supply_vdw = {u_supply_vdw:.12g} J/kg

    h0_vdw = u_supply_vdw + P_supply*v_supply_vdw
           = {h0_vdw:.12g} J/kg

No se sustituye h0_vdw por cp*T_supply.

El depósito y su entalpía de suministro usan Van der Waals,
pero la válvula conserva la ley de caudal ideal.


12. DIFERENCIAS ENTRE MODELOS Y FACTOR Z
----------------------------------------------------------------------

b representa volumen excluido y a atracción intermolecular.

Cambian:

    - La relación P-v-T.
    - La masa inicial a iguales P0, T0 y V.
    - La energía interna.
    - La entalpía del suministro.
    - La evolución de temperatura y presión.
    - El tiempo de cambio de régimen.
    - El tiempo y masa necesarios para alcanzar el objetivo.

Factor de compresibilidad:

    Z = P*v/(R*T)

Ideal:

    Z = 1

Van der Waals:

    Z = v/(v-b) - a/(R*T*v)

Z mide desviación respecto a la relación ideal.
No cuantifica por sí solo la exactitud frente a datos experimentales.

Van der Waals es una aproximación, no una referencia exacta
del comportamiento del nitrógeno.


13. MÉTODO NUMÉRICO Y EVENTOS
----------------------------------------------------------------------

Se integran dos sistemas independientes con:

    y = [m, T]

Configuración de solve_ivp:

    Método: RK45 adaptativo
    Intervalo: 0 a {T_MAX:.12g} s
    rtol: 1e-9
    atol: 1e-11
    max_step: 0.5 s
    dense_output: True

Los tiempos de los cuadros no son los pasos adaptativos del integrador.

Cada modelo se detiene al resolver:

    P - P_target = 0

con:

    terminal = True
    direction = 1

Si no aparece el evento se lanza RuntimeError.

Al inicio se comprueba:

    P0 < P_target < P_supply
    V > 0

Para obtener v a presión y temperatura conocidas, brentq resuelve:

    P_vdw(T,v) - P = 0

Límites:

    v_min = b*(1 + 1e-10)

    v_max = max(100*R*T/P, 1e4*b)

El límite superior se duplica mientras su residuo sea positivo.

Tolerancias:

    xtol = 1e-12
    rtol = 1e-12

pressure_vdw rechaza v<=b. Este control no verifica por sí solo
estabilidad termodinámica ni selecciona fases para otros parámetros.

La transición a subsónico se busca con brentq en la salida densa:

    P(t) - P_critical_tank = 0

Se busca si P_target > P_critical_tank.

El procedimiento presupone que la presión crítica se cruza dentro
del intervalo. No es un detector general para cualquier estado inicial.


14. RELOJ COMÚN Y CONGELAMIENTO DE LAS CURVAS
----------------------------------------------------------------------

Se usan {N_FRAMES} muestras hasta:

    tf_global = max(tf_ideal, tf_vdw)

Para evaluar cada modelo:

    time_I = minimum(time, tf_ideal)

    time_V = minimum(time, tf_vdw)

Después del evento se repite exactamente el último estado.

La meseta posterior es un congelamiento de la presentación.
No se integra una etapa posterior ni el cierre físico de la válvula.

El programa calcula los caudales mostrados evaluando mass_flow_rate
en las presiones muestreadas, incluso después de congelar cada estado.

Como P_target < P_supply, el caudal de esa ley al evento puede
seguir siendo positivo.

Después del evento ese valor es un diagnóstico de la válvula
a la presión final; no representa masa que siga entrando
en una simulación físicamente integrada.

flow_regime también clasifica la presión final repetida, no
el estado mecánico de una válvula cerrada.

Esta advertencia se aplica tanto a los gráficos como al CSV.


15. VALIDACIONES NUMÉRICAS
----------------------------------------------------------------------

15.1. Solución analítica del gas ideal

Se compara TI con:

    T(mI) = gamma*T_supply
            + (m0_ideal/mI)*(T0 - gamma*T_supply)

Máximo error absoluto en las muestras del reloj común:

    {error_T_ideal:.8e} K

15.2. Balance de energía

Para cada modelo se evalúan 1200 tiempos entre 0 y su propio tf.

Residuo:

    U - U0 - (m-m0)*h0_supply

Indicador:

    max(|residuo| / max(|U|, 1.0))

U se evalúa en joules y el piso del denominador es 1 J.

Resultados:

    Ideal = {energy_error_I:.8e}

    VdW = {energy_error_V:.8e}

Estas comprobaciones verifican consistencia numérica y de balances.
No constituyen validación experimental ni estudio de convergencia.

El programa reporta indicadores sin imponer un umbral automático
de aceptación.

No se comprueba de forma independiente la integral temporal del caudal:
la energía entrante se calcula a partir de la masa integrada.

Estado del integrador ideal:

    success = {sol_ideal.success}
    message = {sol_ideal.message}

Estado del integrador Van der Waals:

    success = {sol_vdw.success}
    message = {sol_vdw.message}


16. RESULTADOS DE ESTA EJECUCIÓN
----------------------------------------------------------------------

GAS IDEAL

    Tiempo final = {tf_ideal:.9g} s

    Transición a subsónico = {fmt_time(tcrit_ideal)}

    Temperatura final = {TI_final:.9g} K
                      = {TI_final - 273.15:.9g} °C

    Presión final = {PI_final / 1e5:.9g} bar

    Masa final = {mI_final:.9g} kg

    Z final = {ZI_final:.9g}

    Caudal de la ley al evento = {mdot_I_final:.9g} kg/s

    Caudal promedio = {mdot_avg_I:.9g} kg/s

VAN DER WAALS

    Tiempo final = {tf_vdw:.9g} s

    Transición a subsónico = {fmt_time(tcrit_vdw)}

    Temperatura final = {TV_final:.9g} K
                      = {TV_final - 273.15:.9g} °C

    Presión final = {PV_final / 1e5:.9g} bar

    Masa final = {mV_final:.9g} kg

    Z final = {ZV_final:.9g}

    Caudal de la ley al evento = {mdot_V_final:.9g} kg/s

    Caudal promedio = {mdot_avg_V:.9g} kg/s

El caudal promedio se calcula como:

    (m_final-m0)/tf


17. COMPARACIÓN E INTERPRETACIÓN
----------------------------------------------------------------------

Diferencias Van der Waals menos ideal en sus respectivos eventos:

    Delta tiempo = {tf_vdw - tf_ideal:+.9g} s

    Delta temperatura = {TV_final - TI_final:+.9g} K

    Delta masa = {mV_final - mI_final:+.9g} kg

Se comparan estados a igual presión objetivo, no necesariamente
a igual tiempo.

Para comparar simultáneamente los modelos debe comprobarse
si alguno ya está congelado.

No deben atribuirse todas las diferencias solo a Z:
también cambian la energía interna y la entalpía del suministro.

Durante la rama estrangulada el caudal permanece constante.
Después de la transición disminuye con la contrapresión.

El calentamiento se obtiene del balance de energía, no de una
transferencia de calor externa.


18. ACERCAMIENTO A LA PRESIÓN DEL SUMINISTRO
----------------------------------------------------------------------

Al acercarse r a 1 desde abajo, la ley subsónica da:

    mdot -> 0

El programa termina antes, en P_target < P_supply.

No simula ni demuestra el tiempo de igualación exacta de presiones.

La condición mdot->0 por sí sola no demuestra que el tiempo
de igualación sea infinito.

Esta ley se comporta como sqrt(1-r) cerca de r=1.
Por tanto, no se afirma un acercamiento necesariamente asintótico.


19. MAPA TÉRMICO 2D: CONSTRUCCIÓN Y LIMITACIONES
----------------------------------------------------------------------

El modelo físico calcula únicamente T(t).
El mapa de colores es ilustrativo, no una solución física T(x,y,t).

Malla:

    150 columnas
    80 filas

El campo parte de T0 expresada en °C.

Se aplican cinco actualizaciones visuales por cuadro con:

    Difusión: r_num = 0.075
    Advección: 0.045
    Relajación global: 0.055

Dos gaussianas representan entrada y pluma.

El objetivo local incorpora contribuciones con pesos 0.20 y 0.07
respecto al calentamiento positivo sobre T0.

Las relajaciones locales tienen pesos 0.10 y 0.03.

La posición vertical oscila según:

    0.055*NY*sin(2*pi*k/70)

Los bordes copian celdas vecinas como aislamiento visual.

No hay tamaños físicos de celda, paso físico de transporte,
velocidad ni difusividad asociados a esos coeficientes.

No se resuelven una EDP física de calor ni Navier-Stokes.

El mapa no utiliza mdot para ajustar la intensidad del chorro:
se construye a partir de T promedio y del índice de cuadro.

Por ello, no representa físicamente la reducción del caudal.

La media espacial no está obligada a coincidir exactamente con T(t).
Puede haber retraso y sobretemperaturas artificiales.

El campo puede seguir cambiando después de congelarse T promedio.
Su aspecto depende de N_FRAMES.

Los mapas comparten escala de colores, pero no predicen temperaturas
de pared, puntos calientes, flujos de calor ni tiempos de mezcla.

El campo no retroalimenta las ecuaciones de masa y energía.


20. ALCANCE DE LA VERSIÓN 2
----------------------------------------------------------------------

Se amplía la Versión 1 incorporando:

    - Caudal variable según la presión.
    - Rama estrangulada y subsónica.
    - Calibración del orificio.
    - Cálculo del tiempo de transición.
    - Exportación de resultados en CSV.

Se mantienen:

    - Depósito adiabático 0D.
    - Suministro constante.
    - cv constante.
    - Van der Waals como aproximación de gas real.
    - Visualización térmica ilustrativa.

No se incorporan:

    - Calor a pared o ambiente.
    - Capacidad térmica de pared.
    - cp(T) dinámico.
    - Modelo de válvula para gas real.
    - Agotamiento del suministro.
    - Flujo inverso.
    - CFD.
    - Cambios de fase.
    - Validación experimental.

Es una comparación didáctica de modelos.
No constituye certificación de diseño de equipos de alta presión.

Posibles ampliaciones: pared y ambiente, capacidades variables
y flujo de gas real, verificando nuevamente balances y dominio.


21. ARCHIVOS Y DIRECTORIO DE SALIDA
----------------------------------------------------------------------

Todos los archivos se guardan en:

    {OUTPUT_DIR}

Reporte TXT:
    {TXT_PATH}

Datos CSV:
    {CSV_PATH}

Imagen PNG:
    {PNG_PATH}

Video MP4:
    {MP4_PATH}

GIF opcional:
    {GIF_PATH}

Configuración:

    SAVE_TXT = {SAVE_TXT}
    SAVE_CSV = {SAVE_CSV}
    SAVE_PNG = {SAVE_PNG}
    SAVE_MP4 = {SAVE_MP4}
    SAVE_GIF = {SAVE_GIF}

Animación:

    N_FRAMES = {N_FRAMES}
    FPS = {FPS}

Los FPS son velocidad de reproducción, no el paso del integrador.
El GIF utiliza 18 FPS.

En un script se utiliza su propia carpeta.
En un notebook se utiliza la carpeta de trabajo.

Los archivos de igual nombre se reemplazan al ejecutar nuevamente.

TXT y CSV se generan antes del mapa térmico y no requieren FFmpeg.
El MP4 requiere FFmpeg.

Las rutas anteriores indican destinos configurados; la generación
de cada formato depende de su opción y de que termine la exportación.

El CSV contiene {N_FRAMES} muestras del reloj común con:

    tiempo,
    presiones,
    temperaturas en °C,
    masas,
    caudales de la ley de válvula,
    factores Z.

Incluye estados repetidos después del evento de cada modelo.
No contiene todos los pasos internos de RK45 ni una etapa física
posterior al cierre.

Los parámetros y resultados numéricos de este reporte proceden
directamente de la ejecución actual.

======================================================================
FIN DEL REPORTE — VERSIÓN 2
======================================================================
"""

    TXT_PATH.write_text(
        report,
        encoding="utf-8",
    )

    print(f"TXT generado: {TXT_PATH}")


if SAVE_TXT:
    create_report()


NX = 150
NY = 80


def create_thermal_field(T_average_C):
    field = np.full(
        (NY, NX),
        T0 - 273.15,
        dtype=np.float32,
    )

    Y, X = np.mgrid[0:NY, 0:NX]

    result = np.zeros(
        (len(time), NY, NX),
        dtype=np.float32,
    )

    for k in range(len(time)):
        target = T_average_C[k]

        delta = max(
            target - (T0 - 273.15),
            0.0,
        )

        center_y = (
            NY / 2
            + 0.055 * NY
            * np.sin(
                2.0 * np.pi * k / 70.0
            )
        )

        source = np.exp(
            -((X - 8) / 18.0)**2
            - ((Y - center_y) / 9.0)**2
        )

        plume = np.exp(
            -((X - 32) / 48.0)**2
            - ((Y - center_y) / 18.0)**2
        )

        local_target = (
            target
            + 0.20 * delta * source
            + 0.07 * delta * plume
        )

        for _ in range(5):
            new = field.copy()

            r_num = 0.075

            new[1:-1, 1:-1] = (
                (1 - 4 * r_num)
                * field[1:-1, 1:-1]
                + r_num * (
                    field[1:-1, 2:]
                    + field[1:-1, :-2]
                    + field[2:, 1:-1]
                    + field[:-2, 1:-1]
                )
            )

            advection = 0.045

            new[:, 1:] = (
                (1 - advection)
                * new[:, 1:]
                + advection
                * new[:, :-1]
            )

            new += (
                0.055 * (target - new)
            )

            new += (
                0.10 * source
                * (local_target - new)
            )

            new += (
                0.03 * plume
                * (local_target - new)
            )

            new[:, 0] = new[:, 1]
            new[:, -1] = new[:, -2]
            new[0, :] = new[1, :]
            new[-1, :] = new[-2, :]

            field = new

        result[k] = field

    return result


print()
print("Generando campo térmico ideal...")

field_I = create_thermal_field(TI_C)

print("Generando campo térmico Van der Waals...")

field_V = create_thermal_field(TV_C)

T_color_min = min(
    field_I.min(),
    field_V.min(),
)

T_color_max = max(
    field_I.max(),
    field_V.max(),
)


BG = "#070A0F"
PANEL = "#10151D"
GRID = "#293442"

IDEAL_COLOR = "#47D7FF"
VDW_COLOR = "#FF9F43"

WHITE = "#F4F7FB"
MUTED = "#AAB4C3"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titleweight": "bold",
    "axes.labelcolor": "#D6DCE5",
    "xtick.color": "#AEB8C6",
    "ytick.color": "#AEB8C6",
    "text.color": WHITE,
})


fig = plt.figure(
    figsize=(16, 9),
    facecolor=BG,
)

gs = fig.add_gridspec(
    4,
    4,
    left=0.045,
    right=0.955,
    bottom=0.080,
    top=0.955,
    wspace=0.34,
    hspace=0.48,
    height_ratios=[
        0.13,
        0.24,
        1.07,
        0.82,
    ],
)


ax_header = fig.add_subplot(gs[0, :])
ax_header.axis("off")

ax_header.text(
    0.00,
    0.76,
    "LLENADO TRANSIENTE DE UN DEPÓSITO",
    fontsize=22,
    fontweight="bold",
    color=WHITE,
    ha="left",
    va="center",
)

ax_header.text(
    0.00,
    0.12,
    (
        "Gas ideal vs. Van der Waals"
        "  |  flujo compresible variable"
        "  |  ICF130"
    ),
    fontsize=10.5,
    color=MUTED,
    ha="left",
    va="center",
)

ax_header.text(
    1.00,
    0.12,
    "Mapa térmico: representación visual, no CFD",
    fontsize=8.5,
    color=MUTED,
    ha="right",
    va="center",
)


def create_kpi(axis, title):
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")

    card = FancyBboxPatch(
        (0.02, 0.08),
        0.96,
        0.84,
        boxstyle=(
            "round,pad=0.015,"
            "rounding_size=0.04"
        ),
        facecolor=PANEL,
        edgecolor="#293340",
        linewidth=1.0,
    )

    axis.add_patch(card)

    axis.text(
        0.08,
        0.72,
        title,
        fontsize=8.5,
        color=MUTED,
        ha="left",
        va="center",
    )

    value = axis.text(
        0.08,
        0.35,
        "",
        fontsize=15,
        fontweight="bold",
        color=WHITE,
        ha="left",
        va="center",
    )

    return value


kpi_time = create_kpi(
    fig.add_subplot(gs[1, 0]),
    "TIEMPO",
)

kpi_pressure = create_kpi(
    fig.add_subplot(gs[1, 1]),
    "PRESIÓN I / VdW",
)

kpi_flow = create_kpi(
    fig.add_subplot(gs[1, 2]),
    "CAUDAL I / VdW",
)

kpi_regime = create_kpi(
    fig.add_subplot(gs[1, 3]),
    "RÉGIMEN IDEAL / VdW",
)


ax_I = fig.add_subplot(gs[2, 0:2])
ax_V = fig.add_subplot(gs[2, 2:4])

for ax in [ax_I, ax_V]:
    ax.set_facecolor(PANEL)
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_color("#303B48")

im_I = ax_I.imshow(
    field_I[0],
    origin="lower",
    aspect="auto",
    cmap="inferno",
    vmin=T_color_min,
    vmax=T_color_max,
)

im_V = ax_V.imshow(
    field_V[0],
    origin="lower",
    aspect="auto",
    cmap="inferno",
    vmin=T_color_min,
    vmax=T_color_max,
)

ax_I.set_title(
    "GAS IDEAL",
    color=IDEAL_COLOR,
    fontsize=14,
    pad=10,
)

ax_V.set_title(
    "VAN DER WAALS",
    color=VDW_COLOR,
    fontsize=14,
    pad=10,
)


for ax in [ax_I, ax_V]:
    ax.annotate(
        "",
        xy=(5, NY / 2),
        xytext=(-13, NY / 2),
        arrowprops=dict(
            arrowstyle="-|>",
            linewidth=3,
            color=WHITE,
        ),
    )


status_I = ax_I.text(
    0.025,
    0.04,
    "",
    transform=ax_I.transAxes,
    fontsize=9,
    color=WHITE,
    ha="left",
    va="bottom",
    bbox=dict(
        boxstyle="round,pad=0.4",
        facecolor=(0, 0, 0, 0.50),
        edgecolor="none",
    ),
)

status_V = ax_V.text(
    0.025,
    0.04,
    "",
    transform=ax_V.transAxes,
    fontsize=9,
    color=WHITE,
    ha="left",
    va="bottom",
    bbox=dict(
        boxstyle="round,pad=0.4",
        facecolor=(0, 0, 0, 0.50),
        edgecolor="none",
    ),
)


cbar = fig.colorbar(
    im_I,
    ax=[ax_I, ax_V],
    orientation="vertical",
    fraction=0.020,
    pad=0.015,
)

cbar.set_label(
    "Temperatura local ilustrativa [°C]",
    color=MUTED,
)

cbar.ax.tick_params(colors=MUTED)


ax_T = fig.add_subplot(gs[3, 0])
ax_P = fig.add_subplot(gs[3, 1])
ax_F = fig.add_subplot(gs[3, 2])
ax_Z = fig.add_subplot(gs[3, 3])

charts = [
    ax_T,
    ax_P,
    ax_F,
    ax_Z,
]

for ax in charts:
    ax.set_facecolor(PANEL)

    ax.grid(
        color=GRID,
        alpha=0.30,
        linewidth=0.7,
    )

    for spine in ax.spines.values():
        spine.set_color("#303B48")

    ax.tick_params(labelsize=8)


ax_T.set_title("Temperatura promedio")
ax_T.set_xlabel("t [s]")
ax_T.set_ylabel("T [°C]")

ax_T.plot(
    time,
    TI_C,
    color=IDEAL_COLOR,
    alpha=0.18,
)

ax_T.plot(
    time,
    TV_C,
    color=VDW_COLOR,
    alpha=0.18,
)

line_TI, = ax_T.plot(
    [],
    [],
    color=IDEAL_COLOR,
    linewidth=2.4,
)

line_TV, = ax_T.plot(
    [],
    [],
    color=VDW_COLOR,
    linewidth=2.4,
)


ax_P.set_title("Presión")
ax_P.set_xlabel("t [s]")
ax_P.set_ylabel("P [bar]")

ax_P.plot(
    time,
    PI_bar,
    color=IDEAL_COLOR,
    alpha=0.18,
)

ax_P.plot(
    time,
    PV_bar,
    color=VDW_COLOR,
    alpha=0.18,
)

line_PI, = ax_P.plot(
    [],
    [],
    color=IDEAL_COLOR,
    linewidth=2.4,
)

line_PV, = ax_P.plot(
    [],
    [],
    color=VDW_COLOR,
    linewidth=2.4,
)

ax_P.axhline(
    P_critical_tank / 1e5,
    color="#E9C46A",
    linestyle="--",
    linewidth=1.0,
    alpha=0.55,
    label="P crítica",
)

ax_P.axhline(
    P_target / 1e5,
    color=WHITE,
    linestyle=":",
    linewidth=1.0,
    alpha=0.35,
)


ax_F.set_title("Caudal másico")
ax_F.set_xlabel("t [s]")
ax_F.set_ylabel("ṁ [kg/s]")

ax_F.plot(
    time,
    mdot_I,
    color=IDEAL_COLOR,
    alpha=0.18,
)

ax_F.plot(
    time,
    mdot_V,
    color=VDW_COLOR,
    alpha=0.18,
)

line_FI, = ax_F.plot(
    [],
    [],
    color=IDEAL_COLOR,
    linewidth=2.4,
)

line_FV, = ax_F.plot(
    [],
    [],
    color=VDW_COLOR,
    linewidth=2.4,
)

if tcrit_ideal is not None:
    ax_F.axvline(
        tcrit_ideal,
        color=IDEAL_COLOR,
        linestyle="--",
        linewidth=1,
        alpha=0.35,
    )

if tcrit_vdw is not None:
    ax_F.axvline(
        tcrit_vdw,
        color=VDW_COLOR,
        linestyle="--",
        linewidth=1,
        alpha=0.35,
    )


ax_Z.set_title("Factor de compresibilidad")
ax_Z.set_xlabel("t [s]")
ax_Z.set_ylabel("Z")

ax_Z.plot(
    time,
    ZI,
    color=IDEAL_COLOR,
    alpha=0.18,
)

ax_Z.plot(
    time,
    ZV,
    color=VDW_COLOR,
    alpha=0.18,
)

line_ZI, = ax_Z.plot(
    [],
    [],
    color=IDEAL_COLOR,
    linewidth=2.4,
)

line_ZV, = ax_Z.plot(
    [],
    [],
    color=VDW_COLOR,
    linewidth=2.4,
)

ax_Z.axhline(
    1.0,
    color=WHITE,
    linestyle="--",
    linewidth=1,
    alpha=0.30,
)


for ax in charts:
    ax.set_xlim(
        0,
        tf_global,
    )

ax_T.set_ylim(
    min(TI_C.min(), TV_C.min()) - 5,
    max(TI_C.max(), TV_C.max()) + 8,
)

ax_P.set_ylim(
    0,
    P_target / 1e5 + 10,
)

ax_F.set_ylim(
    0,
    mdot_reference * 1.15,
)

ax_Z.set_ylim(
    min(ZI.min(), ZV.min()) - 0.04,
    max(ZI.max(), ZV.max()) + 0.04,
)


legend_elements = [
    Line2D(
        [0],
        [0],
        color=IDEAL_COLOR,
        lw=3,
        label="Gas ideal",
    ),
    Line2D(
        [0],
        [0],
        color=VDW_COLOR,
        lw=3,
        label="Van der Waals",
    ),
]

fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=2,
    frameon=False,
    bbox_to_anchor=(0.5, 0.017),
)


def update(k):
    im_I.set_data(field_I[k])
    im_V.set_data(field_V[k])

    sl = slice(0, k + 1)

    line_TI.set_data(
        time[sl],
        TI_C[sl],
    )

    line_TV.set_data(
        time[sl],
        TV_C[sl],
    )

    line_PI.set_data(
        time[sl],
        PI_bar[sl],
    )

    line_PV.set_data(
        time[sl],
        PV_bar[sl],
    )

    line_FI.set_data(
        time[sl],
        mdot_I[sl],
    )

    line_FV.set_data(
        time[sl],
        mdot_V[sl],
    )

    line_ZI.set_data(
        time[sl],
        ZI[sl],
    )

    line_ZV.set_data(
        time[sl],
        ZV[sl],
    )

    kpi_time.set_text(
        f"{time[k]:.1f} s"
    )

    kpi_pressure.set_text(
        f"{PI_bar[k]:.1f} / "
        f"{PV_bar[k]:.1f} bar"
    )

    kpi_flow.set_text(
        f"{mdot_I[k]:.4f} / "
        f"{mdot_V[k]:.4f}"
    )

    regime_I = flow_regime(PI[k])
    regime_V = flow_regime(PV[k])

    kpi_regime.set_text(
        f"{regime_I}\n{regime_V}"
    )

    text_I = (
        f"T = {TI_C[k]:6.1f} °C\n"
        f"P = {PI_bar[k]:6.1f} bar\n"
        f"m = {mI[k]:6.3f} kg\n"
        f"ṁ = {mdot_I[k]:.4f} kg/s\n"
        f"{regime_I}"
    )

    text_V = (
        f"T = {TV_C[k]:6.1f} °C\n"
        f"P = {PV_bar[k]:6.1f} bar\n"
        f"m = {mV[k]:6.3f} kg\n"
        f"ṁ = {mdot_V[k]:.4f} kg/s\n"
        f"{regime_V}"
    )

    if time[k] >= tf_ideal:
        text_I += "\nOBJETIVO ALCANZADO"

    if time[k] >= tf_vdw:
        text_V += "\nOBJETIVO ALCANZADO"

    status_I.set_text(text_I)
    status_V.set_text(text_V)

    return (
        im_I,
        im_V,
        line_TI,
        line_TV,
        line_PI,
        line_PV,
        line_FI,
        line_FV,
        line_ZI,
        line_ZV,
        kpi_time,
        kpi_pressure,
        kpi_flow,
        kpi_regime,
        status_I,
        status_V,
    )


animation = FuncAnimation(
    fig,
    update,
    frames=N_FRAMES,
    interval=1000 / FPS,
    blit=False,
    repeat=True,
)


print()
print("=" * 78)
print("ICF130 — VERSIÓN 2 — RESULTADOS")
print("=" * 78)

print()
print("PROPIEDADES N2")
print("-" * 78)
print(f"R                       = {R:.6f} J/(kg K)")
print(f"cp                      = {cp:.6f} J/(kg K)")
print(f"cv                      = {cv:.6f} J/(kg K)")
print(f"gamma                   = {gamma:.8f}")

print()
print("FLUJO COMPRESIBLE")
print("-" * 78)
print(
    f"Razón crítica           = "
    f"{critical_pressure_ratio:.8f}"
)
print(
    f"P crítica tanque        = "
    f"{P_critical_tank / 1e5:.4f} bar"
)
print(
    f"Área orificio           = "
    f"{A_orifice:.8e} m²"
)
print(
    f"Diámetro equivalente    = "
    f"{d_orifice * 1000:.5f} mm"
)

print()
print("GAS IDEAL")
print("-" * 78)
print(f"Tiempo final            = {tf_ideal:.4f} s")
print(
    f"Transición flujo        = "
    f"{fmt_time(tcrit_ideal)}"
)
print(
    f"T final                 = "
    f"{TI_final - 273.15:.4f} °C"
)
print(f"Masa final              = {mI_final:.6f} kg")
print(
    f"mdot de la ley al evento= "
    f"{mdot_I_final:.6f} kg/s"
)

print()
print("VAN DER WAALS")
print("-" * 78)
print(f"Tiempo final            = {tf_vdw:.4f} s")
print(
    f"Transición flujo        = "
    f"{fmt_time(tcrit_vdw)}"
)
print(
    f"T final                 = "
    f"{TV_final - 273.15:.4f} °C"
)
print(f"Masa final              = {mV_final:.6f} kg")
print(
    f"mdot de la ley al evento= "
    f"{mdot_V_final:.6f} kg/s"
)
print(f"Z final                 = {ZV_final:.6f}")

print()
print("VALIDACIONES")
print("-" * 78)
print(
    f"Error T analítico ideal = "
    f"{error_T_ideal:.4e} K"
)
print(
    f"Error energía ideal     = "
    f"{energy_error_I:.4e}"
)
print(
    f"Error energía VdW       = "
    f"{energy_error_V:.4e}"
)
print("=" * 78)


if SAVE_PNG:
    update(N_FRAMES - 1)

    fig.savefig(
        PNG_PATH,
        dpi=200,
        facecolor=BG,
    )

    print()
    print(f"PNG generado:\n{PNG_PATH}")


ffmpeg_path = shutil.which("ffmpeg")

if ffmpeg_path is not None:
    plt.rcParams[
        "animation.ffmpeg_path"
    ] = ffmpeg_path


mp4_generado = False
gif_generado = False

if SAVE_MP4:
    if ffmpeg_path is None:
        print()
        print("FFmpeg no fue encontrado por Python.")
        print("Se omite la exportación MP4.")

    else:
        try:
            print()
            print("Generando video MP4...")

            writer = FFMpegWriter(
                fps=FPS,
                bitrate=5000,
                metadata={
                    "title": "ICF130 Tema 7",
                    "artist": "Proyecto Computacional",
                },
            )

            animation.save(
                str(MP4_PATH),
                writer=writer,
                dpi=150,
            )

            mp4_generado = True

            print(f"MP4 generado:\n{MP4_PATH}")

        except Exception as error:
            print()
            print("No fue posible generar el MP4.")
            print(error)


if SAVE_GIF:
    try:
        animation.save(
            str(GIF_PATH),
            writer=PillowWriter(fps=18),
            dpi=105,
        )

        gif_generado = True

        print(f"GIF generado:\n{GIF_PATH}")

    except Exception as error:
        print()
        print("No fue posible generar el GIF.")
        print(error)


print()
print("=" * 78)
print("ARCHIVOS GENERADOS EN ESTA EJECUCIÓN")
print("=" * 78)

if SAVE_TXT:
    print(f"Reporte TXT : {TXT_PATH}")

if SAVE_CSV:
    print(f"Datos CSV   : {CSV_PATH}")

if SAVE_PNG:
    print(f"Imagen PNG  : {PNG_PATH}")

if mp4_generado:
    print(f"Video MP4   : {MP4_PATH}")

if gif_generado:
    print(f"Animación GIF: {GIF_PATH}")

print()
print(f"Directorio de salida: {OUTPUT_DIR}")
print("=" * 78)

plt.show()