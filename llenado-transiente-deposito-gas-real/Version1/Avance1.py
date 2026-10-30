import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

if "__file__" in globals():
    DIRECTORIO_SALIDA = Path(__file__).resolve().parent
else:
    DIRECTORIO_SALIDA = Path.cwd()
    
FPS = 30
N_FRAMES = 240

SAVE_MP4 = True
SAVE_GIF = False
SAVE_STATIC_PNG = True

BG = "#070A0F"
PANEL = "#10151D"
GRID = "#2A3340"

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


V = 0.050

P0 = 1.0e5
T0 = 300.0

P_in = 200.0e5
T_in = 300.0

mdot = 0.020
P_target = 100.0e5


M = 28.0134e-3
R_u = 8.314462618
R = R_u / M

cp = 1039.0
cv = cp - R
gamma = cp / cv

Tc = 126.19
Pc = 33.534 * 101325.0

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
            "Estado no físico en Van der Waals: v <= b."
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


v_in_ideal = R * T_in / P_in
h_in_ideal = cp * T_in

v_in_vdw = specific_volume_vdw(P_in, T_in)
u_in_vdw = u_vdw(T_in, v_in_vdw)
h_in_vdw = u_in_vdw + P_in * v_in_vdw

m0_ideal = P0 * V / (R * T0)

v0_vdw = specific_volume_vdw(P0, T0)
m0_vdw = V / v0_vdw


def rhs_ideal(t, y):
    m, T = y

    dm_dt = mdot

    dT_dt = (
        mdot / m
        * (gamma * T_in - T)
    )

    return [dm_dt, dT_dt]


def rhs_vdw(t, y):
    m, T = y

    v = V / m
    P = pressure_vdw(m, T)
    u = u_vdw(T, v)

    dP_dT_v = R / (v - b)

    thermodynamic_term = T * dP_dT_v - P

    dm_dt = mdot

    dT_dt = (
        mdot / (m * cv)
        * (
            h_in_vdw
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


sol_ideal = solve_ivp(
    rhs_ideal,
    (0.0, 1000.0),
    (m0_ideal, T0),
    method="RK45",
    events=event_ideal,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.25,
)

sol_vdw = solve_ivp(
    rhs_vdw,
    (0.0, 1000.0),
    (m0_vdw, T0),
    method="RK45",
    events=event_vdw,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.25,
)

if len(sol_ideal.t_events[0]) == 0:
    raise RuntimeError(
        "Gas ideal no alcanzó la presión objetivo."
    )

if len(sol_vdw.t_events[0]) == 0:
    raise RuntimeError(
        "Van der Waals no alcanzó la presión objetivo."
    )

tf_ideal = sol_ideal.t_events[0][0]
tf_vdw = sol_vdw.t_events[0][0]

tf_global = max(tf_ideal, tf_vdw)

time = np.linspace(
    0,
    tf_global,
    N_FRAMES,
)

time_I = np.minimum(time, tf_ideal)
time_V = np.minimum(time, tf_vdw)

mI, TI = sol_ideal.sol(time_I)
mV, TV = sol_vdw.sol(time_V)

PI = pressure_ideal(mI, TI)

PV = np.array([
    pressure_vdw(m, T)
    for m, T in zip(mV, TV)
])

TI_C = TI - 273.15
TV_C = TV - 273.15

PI_bar = PI / 1e5
PV_bar = PV / 1e5

vI = V / mI
vV = V / mV

ZI = PI * vI / (R * TI)
ZV = PV * vV / (R * TV)

T_analytical = (
    gamma * T_in
    + (m0_ideal / mI)
    * (T0 - gamma * T_in)
)

error_T_ideal = np.max(
    np.abs(TI - T_analytical)
)

U_I = mI * u_ideal(TI)
U0_I = m0_ideal * u_ideal(T0)
E_in_I = (mI - m0_ideal) * h_in_ideal

energy_residual_I = U_I - U0_I - E_in_I

energy_error_I = np.max(
    np.abs(energy_residual_I)
    / np.maximum(np.abs(U_I), 1.0)
)

U_V = mV * u_vdw(TV, vV)
U0_V = m0_vdw * u_vdw(T0, v0_vdw)
E_in_V = (mV - m0_vdw) * h_in_vdw

energy_residual_V = U_V - U0_V - E_in_V

energy_error_V = np.max(
    np.abs(energy_residual_V)
    / np.maximum(np.abs(U_V), 1.0)
)


def generar_informe_version_1():
    ruta = DIRECTORIO_SALIDA / "ideal_vs_vdw_reporte.txt"
    mass_error_ideal = np.max(
        np.abs(mI - m0_ideal - mdot * time_I)
    )

    mass_error_vdw = np.max(
        np.abs(mV - m0_vdw - mdot * time_V)
    )

    texto = f"""ICF130 — PROYECTO COMPUTACIONAL
VERSIÓN 1 — LLENADO TRANSIENTE DE UN DEPÓSITO
COMPARACIÓN: GAS IDEAL VS VAN DER WAALS
======================================================================

1. FENÓMENO FÍSICO Y OBJETIVO
----------------------------------------------------------------------

Se estudia el llenado transiente de un depósito rígido con nitrógeno.
El depósito contiene inicialmente gas y recibe una corriente de entrada.

Durante el llenado se acumulan masa y energía interna. Como consecuencia,
la temperatura promedio y la presión cambian con el tiempo.

El objetivo es comparar las predicciones de dos modelos termodinámicos:

    1) Gas ideal.
    2) Gas descrito mediante la ecuación de Van der Waals.

Ambos emplean balances de masa y energía para un volumen de control
abierto, pero difieren en sus relaciones de presión, energía interna
y entalpía de entrada.

El proceso es adiabático: no existe transferencia de calor a través
de la frontera del depósito. Esto no implica temperatura constante.

La corriente entrante transporta entalpía. Esta incluye la energía
interna del gas y el trabajo de flujo necesario para introducirlo.


2. CONDICIONES DEL PROBLEMA
----------------------------------------------------------------------

Gas: nitrógeno.

Volumen del depósito:
    V = {V:.12g} m³ = {V * 1000:.12g} L

Estado inicial:
    P0 = {P0:.12g} Pa = {P0 / 1e5:.12g} bar
    T0 = {T0:.12g} K = {T0 - 273.15:.12g} °C

Condiciones del suministro:
    P_in = {P_in:.12g} Pa = {P_in / 1e5:.12g} bar
    T_in = {T_in:.12g} K = {T_in - 273.15:.12g} °C

Caudal másico impuesto:
    mdot = {mdot:.12g} kg/s

Presión objetivo:
    P_target = {P_target:.12g} Pa
             = {P_target / 1e5:.12g} bar

Las presiones se interpretan como absolutas en las ecuaciones de estado.


3. PROPIEDADES Y PARÁMETROS UTILIZADOS
----------------------------------------------------------------------

Masa molar:
    M = {M:.12g} kg/mol

Constante universal:
    R_u = {R_u:.12g} J/(mol K)

Constante específica:
    R = R_u/M = {R:.12g} J/(kg K)

Referencia calorífica ideal:
    cp = {cp:.12g} J/(kg K)
    cv = cp - R = {cv:.12g} J/(kg K)
    gamma = cp/cv = {gamma:.12g}

Parámetros críticos usados:
    Tc = {Tc:.12g} K
    Pc = 33.534 * 101325 = {Pc:.12g} Pa

Parámetros de Van der Waals:

    a = (27/64) * R² * Tc² / Pc
      = {a:.12g} Pa·m⁶/kg²

    b = R * Tc / (8 * Pc)
      = {b:.12g} m³/kg

El volumen específico v se expresa en m³/kg. Por ello, a y b son
parámetros específicos por unidad de masa, no parámetros molares.


4. HIPÓTESIS Y SIMPLIFICACIONES
----------------------------------------------------------------------

4.1. Volumen de control

El depósito se trata como un volumen de control fijo y rígido.
Existe una entrada y ninguna salida durante la etapa integrada.

4.2. Trabajo y calor

Se supone:
    - Transferencia de calor nula: Qdot = 0.
    - Trabajo de eje nulo.
    - Trabajo de frontera por variación de volumen nulo.
    - Energías cinética y potencial despreciables.

El trabajo de flujo de entrada ya está incluido en h_in.

4.3. Modelo concentrado 0D

Se supone que las propiedades son uniformes dentro del depósito.
En cada instante se calcula una sola temperatura promedio T(t),
una sola presión P(t) y una masa total m(t).

No se resuelven gradientes espaciales, estratificación, turbulencia,
campos de velocidad ni geometría física del chorro.

4.4. Suministro y caudal

P_in, T_in y mdot se mantienen constantes.

No se modelan:
    - Agotamiento o enfriamiento del suministro.
    - Dinámica de una válvula.
    - Pérdidas de presión en tuberías.
    - Variación del caudal con la diferencia de presión.
    - Flujo compresible a través de una restricción.

P_in se usa para determinar el estado de entrada. No existe una ley
dinámica que calcule mdot a partir de P_in y la presión del depósito.

4.5. Propiedades caloríficas

Se usa cv constante.

cp = cv + R se utiliza como referencia calorífica del modelo ideal.
En Van der Waals no se impone ese cp como capacidad calorífica real:
la entalpía de entrada se calcula mediante h = u + Pv.

4.6. Pared y fases

No se calcula temperatura de pared ni su capacidad térmica.
No se modelan cambios de fase o equilibrio líquido-vapor.


5. BALANCE DE MASA
----------------------------------------------------------------------

Balance general:

    dm/dt = mdot_in - mdot_out

Como no hay salida y el caudal de entrada es constante:

    dm/dt = mdot

Mientras el llenado permanece activo:

    m(t) = m0 + mdot*t

La pendiente de la curva de masa es mdot. Las masas iniciales de los
dos modelos pueden diferir porque cada uno obtiene la densidad inicial
con su propia ecuación de estado.


6. BALANCE DE ENERGÍA
----------------------------------------------------------------------

Para las hipótesis adoptadas:

    dU/dt = mdot*h_in

donde:

    U = m*u

La entalpía de entrada es constante para cada modelo. Por tanto:

    U(t) - U0 = mdot*h_in*t

y, usando el balance de masa:

    U(t) - U0 = (m(t) - m0)*h_in

Esta última igualdad se usa para verificar conservación de energía.


7. MODELO DE GAS IDEAL
----------------------------------------------------------------------

Ecuación de estado:

    P*v = R*T
    v = V/m
    P = m*R*T/V

Energía interna específica:

    u = cv*T

Entalpía específica de entrada:

    h_in_ideal = cp*T_in
               = {h_in_ideal:.12g} J/kg

Volumen específico de entrada:

    v_in_ideal = R*T_in/P_in
               = {v_in_ideal:.12g} m³/kg

Masa inicial:

    m0_ideal = P0*V/(R*T0)
             = {m0_ideal:.12g} kg

Expandiendo el balance de energía:

    cv*(m*dT/dt + T*dm/dt) = mdot*cp*T_in

Con dm/dt = mdot:

    dT/dt = (mdot/m)*(gamma*T_in - T)

Solución analítica usada para validar:

    T(t) = gamma*T_in
           + (m0_ideal/m(t))*(T0 - gamma*T_in)

Si el llenado ideal continuara indefinidamente con estas hipótesis,
la temperatura tendería formalmente a:

    T_limite = gamma*T_in = {gamma * T_in:.12g} K

Este límite no implica que se alcance exactamente antes de llegar
a la presión objetivo.


8. MODELO DE VAN DER WAALS
----------------------------------------------------------------------

Ecuación de estado:

    P = R*T/(v-b) - a/v²

Dominio admitido:

    v > b

El parámetro b representa el efecto de volumen excluido y el término
a/v² introduce una corrección asociada a atracciones intermoleculares.

Relación calorífica:

    du = cv*dT + [T*(dP/dT)_v - P]*dv

Para Van der Waals:

    (dP/dT)_v = R/(v-b)

    T*(dP/dT)_v - P = a/v²

Por integración con cv constante:

    u(T,v) = cv*T - a/v

Se adopta una referencia energética compatible con el modelo ideal.

La energía interna total es:

    U = m*cv*T - a*m²/V

Estado inicial:

    v0_vdw = {v0_vdw:.12g} m³/kg
    m0_vdw = V/v0_vdw = {m0_vdw:.12g} kg

El volumen específico inicial se determina resolviendo:

    P_vdw(T0,v0_vdw) = P0

Estado de entrada:

    v_in_vdw = {v_in_vdw:.12g} m³/kg

    u_in_vdw = cv*T_in - a/v_in_vdw
             = {u_in_vdw:.12g} J/kg

    h_in_vdw = u_in_vdw + P_in*v_in_vdw
             = {h_in_vdw:.12g} J/kg

No se sustituye h_in_vdw por cp*T_in.

Como v = V/m y V es constante:

    dv/dt = -(v/m)*mdot

La ecuación diferencial implementada es:

    dT/dt = mdot/(m*cv)
            * [h_in_vdw - u + v*(T*(dP/dT)_v - P)]

Para estas relaciones resulta equivalente a:

    dT/dt = mdot/(m*cv)
            * [h_in_vdw - cv*T + 2*a/v]


9. DIFERENCIAS ENTRE LOS MODELOS
----------------------------------------------------------------------

Los dos modelos usan los mismos valores de V, P0, T0, P_in, T_in,
mdot y presión objetivo.

Sin embargo, cambian:
    - La ecuación de estado.
    - La masa inicial asociada a P0, T0 y V.
    - El volumen específico del suministro.
    - La energía interna específica.
    - La entalpía de entrada.
    - La evolución de temperatura y presión.
    - La masa y el tiempo necesarios para alcanzar el objetivo.

No toda diferencia puede atribuirse únicamente a la ecuación de presión:
también cambia el cierre energético.

Factor de compresibilidad:

    Z = P*v/(R*T)

Para el gas ideal:

    Z = 1

Para Van der Waals:

    Z = v/(v-b) - a/(R*T*v)

Z expresa desviación respecto a la relación ideal. No cuantifica por
sí solo el error frente a datos experimentales de nitrógeno.

Van der Waals es una aproximación de gas real, no una referencia exacta.


10. MÉTODO NUMÉRICO
----------------------------------------------------------------------

Se integran dos sistemas independientes con variables:

    y = [m, T]

Se utiliza scipy.integrate.solve_ivp con:

    Método: RK45 adaptativo
    Intervalo máximo: 0 a 1000 s
    rtol: 1e-9
    atol: 1e-11
    max_step: 0.25 s
    dense_output: True

La salida densa permite evaluar la solución en los tiempos requeridos
por la animación sin imponer esos tiempos como pasos del integrador.

Los volúmenes específicos de Van der Waals se obtienen con brentq:

    Residuo: P_vdw(T,v) - P
    v_min = b*(1 + 1e-10)
    v_max = max(100*R*T/P, 1e4*b)

El límite superior se duplica mientras su residuo sea positivo.

Tolerancias de brentq:

    xtol = 1e-12
    rtol = 1e-12

El control v>b evita evaluar el estado fuera del dominio admitido.
No constituye una verificación completa de estabilidad termodinámica
ni una selección de fases para otros parámetros.

Las tolerancias numéricas no representan incertidumbre física.


11. DETENCIÓN Y RELOJ COMÚN
----------------------------------------------------------------------

Cada modelo tiene su propio evento:

    P(t) - P_target = 0

El evento es terminal y solo detecta cruces crecientes:

    terminal = True
    direction = 1

Si no se alcanza el objetivo dentro del intervalo, el programa genera
un RuntimeError.

Para mostrar ambos resultados se usa un reloj común hasta:

    tf_global = max(tf_ideal, tf_vdw)

La evaluación de cada solución se limita mediante:

    time_I = minimum(time, tf_ideal)
    time_V = minimum(time, tf_vdw)

Cuando un modelo llega primero al objetivo, su estado final se repite
mientras el otro continúa.


12. EXPLICACIÓN DE LA RECTA HORIZONTAL FINAL
----------------------------------------------------------------------

La meseta exacta posterior al evento aparece porque el programa congela
deliberadamente masa, temperatura, presión y Z del modelo detenido.

No representa un llenado que sigue recibiendo masa con todas esas
variables constantes.

Tampoco se integra explícitamente una etapa posterior al cierre de una
válvula ni una relajación térmica posterior.

Es necesario distinguir dos efectos:

    a) Antes del evento:
       La temperatura ideal se aproxima progresivamente a gamma*T_in.

    b) Después del evento:
       Se repite exactamente el último estado por el uso de minimum.

Por ello, la línea constante final no debe interpretarse solamente
como una asíntota física de temperatura.


13. VALIDACIONES IMPLEMENTADAS
----------------------------------------------------------------------

13.1. Solución analítica ideal

Se compara TI con:

    T_analitica = gamma*T_in
                  + (m0_ideal/mI)*(T0 - gamma*T_in)

Error máximo absoluto:

    max(|TI - T_analitica|) = {error_T_ideal:.6e} K

13.2. Conservación de energía

Se calcula para cada modelo:

    residuo = U - U0 - (m-m0)*h_in

El indicador reportado es:

    max(|residuo| / max(|U|, 1.0))

U se evalúa numéricamente en joules; el piso del denominador es 1 J.

Resultados:

    Error energético ideal = {energy_error_I:.6e}
    Error energético VdW   = {energy_error_V:.6e}

Los indicadores se evalúan sobre las muestras de la animación, incluidos
los estados finales repetidos. No comprueban todos los puntos del
continuo temporal.

El programa reporta estos errores sin imponer un umbral automático
de aceptación.

Son verificaciones de coherencia numérica y de balances. No sustituyen
una validación experimental ni un estudio de convergencia.

13.3. Diagnósticos adicionales incluidos en este TXT

Máximo residuo de masa:

    Ideal = {mass_error_ideal:.6e} kg
    VdW   = {mass_error_vdw:.6e} kg

Desviación de presión final respecto al objetivo:

    Ideal = {PI[-1] - P_target:+.6e} Pa
    VdW   = {PV[-1] - P_target:+.6e} Pa

Margen mínimo de dominio en las muestras VdW:

    min(v-b) = {np.min(vV - b):.6e} m³/kg

Estado del integrador ideal:

    success = {sol_ideal.success}
    message = {sol_ideal.message}

Estado del integrador VdW:

    success = {sol_vdw.success}
    message = {sol_vdw.message}


14. RESULTADOS DE ESTA EJECUCIÓN
----------------------------------------------------------------------

GAS IDEAL

    Tiempo final       = {tf_ideal:.9g} s
    Temperatura final  = {TI[-1]:.9g} K
                       = {TI_C[-1]:.9g} °C
    Presión final      = {PI_bar[-1]:.9g} bar
    Masa final         = {mI[-1]:.9g} kg
    Z final            = {ZI[-1]:.9g}

VAN DER WAALS

    Tiempo final       = {tf_vdw:.9g} s
    Temperatura final  = {TV[-1]:.9g} K
                       = {TV_C[-1]:.9g} °C
    Presión final      = {PV_bar[-1]:.9g} bar
    Masa final         = {mV[-1]:.9g} kg
    Z final            = {ZV[-1]:.9g}

DIFERENCIAS: VAN DER WAALS MENOS IDEAL

    Delta tiempo       = {tf_vdw - tf_ideal:+.9g} s
    Delta temperatura  = {TV[-1] - TI[-1]:+.9g} K
    Delta masa         = {mV[-1] - mI[-1]:+.9g} kg

Estas diferencias comparan los estados al alcanzar la misma presión
objetivo, aunque ocurran en tiempos distintos.


15. INTERPRETACIÓN DE LOS RESULTADOS
----------------------------------------------------------------------

La masa aumenta linealmente durante la etapa de llenado activo porque
el caudal está impuesto como constante.

La evolución de temperatura resulta del balance de energía y del modelo
de propiedades. El aislamiento térmico no impide el calentamiento.

La presión depende simultáneamente de la masa, temperatura y ecuación
de estado. La curva de presión no debe interpretarse solo a partir
de la curva de temperatura.

Para comparar los dos modelos en un mismo tiempo, debe comprobarse si
alguno ya alcanzó su evento y está congelado en la representación.

La comparación final es a igual presión objetivo; no corresponde
necesariamente a igual tiempo físico.


16. MAPA TÉRMICO 2D: CONSTRUCCIÓN Y LIMITACIONES
----------------------------------------------------------------------

El modelo termodinámico es 0D. El mapa 2D es una representación
ilustrativa y no calcula un campo físico T(x,y,t).

La malla visual tiene:

    150 columnas
    80 filas

El campo parte de T0 expresada en °C.

Para cada cuadro se aplican cinco actualizaciones visuales con:

    Difusión visual: r = 0.075
    Advección visual: adv = 0.045
    Relajación hacia la T promedio: 0.055

Se usan dos distribuciones gaussianas para representar entrada y pluma.

El objetivo local añade contribuciones con pesos 0.20 y 0.07 respecto
al calentamiento positivo sobre T0.

También se aplican relajaciones locales con pesos 0.10 y 0.03.

La posición vertical de la entrada oscila según el índice del cuadro.
Los bordes copian sus celdas vecinas como aislamiento visual.

Estas operaciones no tienen asociados:
    - Tamaños físicos de celda.
    - Un paso temporal físico de transporte.
    - Una difusividad térmica real.
    - Una velocidad física del gas.

Por ello, no constituyen una solución física de transporte de calor
o dinámica de fluidos computacional.

La media espacial del campo no está obligada a ser exactamente T(t).
Puede presentar retraso o sobretemperaturas artificiales.

El campo puede seguir cambiando después de congelarse la temperatura
promedio, porque sus operaciones visuales siguen ejecutándose.

Cambiar N_FRAMES puede modificar el aspecto del mapa.

Los dos mapas comparten escala de colores, pero sus valores locales
no predicen temperaturas de pared, puntos calientes reales, flujos
de calor ni tiempos físicos de mezcla.

El mapa no retroalimenta las ecuaciones diferenciales: no modifica
masa, temperatura promedio, presión ni energía interna del modelo 0D.


17. ALCANCE DE ESTA VERSIÓN 1
----------------------------------------------------------------------

Esta versión estudia el llenado transiente adiabático de un depósito
rígido mediante balances concentrados de masa y energía.

Permite:
    - Comparar gas ideal y Van der Waals.
    - Analizar T(t), P(t), m(t) y Z(t).
    - Comparar tiempos y estados a igual presión objetivo.
    - Verificar la solución ideal contra una expresión analítica.
    - Comprobar el cierre del balance de energía.
    - Mostrar una animación ilustrativa del proceso.

Quedan fuera:
    - Transferencia de calor con pared y ambiente.
    - Capacidad térmica de la pared.
    - Propiedades caloríficas variables.
    - Dinámica de válvulas y suministro.
    - Resolución espacial real y CFD.
    - Cambios de fase.
    - Validación experimental.

Van der Waals no se adopta como descripción exacta del nitrógeno.
El programa tiene alcance didáctico y de comparación de modelos;
no constituye una certificación de diseño de equipos a alta presión.

Una versión posterior puede incorporar pared, ambiente, una ley de
caudal y propiedades más precisas, verificando nuevamente balances
y dominio de validez.


18. ARCHIVOS Y CONFIGURACIÓN DE SALIDA
----------------------------------------------------------------------

Animación:

    N_FRAMES = {N_FRAMES}
    FPS = {FPS}

Los FPS determinan la reproducción, no el paso del integrador.

Opciones:

    SAVE_MP4 = {SAVE_MP4}
    SAVE_GIF = {SAVE_GIF}
    SAVE_STATIC_PNG = {SAVE_STATIC_PNG}

El TXT se genera antes del mapa térmico y de las exportaciones.
Su generación no requiere FFmpeg.

Se guarda junto al script. En un notebook sin __file__, se guarda
en la carpeta de trabajo.

El archivo se reemplaza al ejecutar de nuevo el programa.

Los parámetros y resultados numéricos de este documento proceden
directamente de las variables y soluciones de esta ejecución.

======================================================================
FIN DEL INFORME — VERSIÓN 1
======================================================================
"""

    ruta.write_text(texto, encoding="utf-8")

    print()
    print("TXT explicativo generado:")
    print(ruta)

    return ruta


ruta_informe = generar_informe_version_1()


NX = 150
NY = 80


def create_thermal_field(T_average):
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
        target = T_average[k]

        delta_T = max(
            target - (T0 - 273.15),
            0.0,
        )

        center_y = (
            NY / 2
            + 0.05 * NY
            * np.sin(2 * np.pi * k / 70)
        )

        source = np.exp(
            -((X - 8) / 18)**2
            - ((Y - center_y) / 9)**2
        )

        plume = np.exp(
            -((X - 30) / 45)**2
            - ((Y - center_y) / 17)**2
        )

        local_target = (
            target
            + 0.20 * delta_T * source
            + 0.07 * delta_T * plume
        )

        for _ in range(5):
            new = field.copy()

            r = 0.075

            new[1:-1, 1:-1] = (
                (1 - 4 * r) * field[1:-1, 1:-1]
                + r * (
                    field[1:-1, 2:]
                    + field[1:-1, :-2]
                    + field[2:, 1:-1]
                    + field[:-2, 1:-1]
                )
            )

            adv = 0.045

            new[:, 1:] = (
                (1 - adv) * new[:, 1:]
                + adv * new[:, :-1]
            )

            new += 0.055 * (target - new)

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


print("Generando campo térmico del gas ideal...")
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


fig = plt.figure(
    figsize=(16, 9),
    facecolor=BG,
)

gs = fig.add_gridspec(
    nrows=4,
    ncols=4,
    left=0.045,
    right=0.955,
    bottom=0.075,
    top=0.955,
    wspace=0.34,
    hspace=0.46,
    height_ratios=[0.13, 0.24, 1.05, 0.82],
)

ax_header = fig.add_subplot(gs[0, :])
ax_header.axis("off")

ax_header.text(
    0.00,
    0.75,
    "LLENADO TRANSIENTE DE UN DEPÓSITO",
    fontsize=22,
    fontweight="bold",
    color=WHITE,
    ha="left",
    va="center",
)

ax_header.text(
    0.00,
    0.15,
    "Gas ideal vs. Van der Waals  |  Balance transiente de masa y energía",
    fontsize=10.5,
    color=MUTED,
    ha="left",
    va="center",
)

ax_header.text(
    1.00,
    0.15,
    "Campo térmico 2D: visualización ilustrativa",
    fontsize=9,
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
        boxstyle="round,pad=0.015,rounding_size=0.04",
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
        0.36,
        "",
        fontsize=16,
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
    "PRESIÓN IDEAL / REAL",
)

kpi_temperature = create_kpi(
    fig.add_subplot(gs[1, 2]),
    "TEMPERATURA IDEAL / REAL",
)

kpi_Z = create_kpi(
    fig.add_subplot(gs[1, 3]),
    "FACTOR Z REAL",
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
        xy=(4, NY / 2),
        xytext=(-14, NY / 2),
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
        boxstyle="round,pad=0.35",
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
        boxstyle="round,pad=0.35",
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
ax_M = fig.add_subplot(gs[3, 2])
ax_Z = fig.add_subplot(gs[3, 3])

charts = [ax_T, ax_P, ax_M, ax_Z]

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
ax_T.set_xlabel("Tiempo [s]")
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
ax_P.set_xlabel("Tiempo [s]")
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
    P_target / 1e5,
    linestyle="--",
    color=WHITE,
    alpha=0.25,
)


ax_M.set_title("Masa acumulada")
ax_M.set_xlabel("Tiempo [s]")
ax_M.set_ylabel("m [kg]")

ax_M.plot(
    time,
    mI,
    color=IDEAL_COLOR,
    alpha=0.18,
)

ax_M.plot(
    time,
    mV,
    color=VDW_COLOR,
    alpha=0.18,
)

line_MI, = ax_M.plot(
    [],
    [],
    color=IDEAL_COLOR,
    linewidth=2.4,
)

line_MV, = ax_M.plot(
    [],
    [],
    color=VDW_COLOR,
    linewidth=2.4,
)


ax_Z.set_title("Factor de compresibilidad")
ax_Z.set_xlabel("Tiempo [s]")
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
    linestyle="--",
    color=WHITE,
    alpha=0.25,
)


for ax in charts:
    ax.set_xlim(0, tf_global)

ax_T.set_ylim(
    min(TI_C.min(), TV_C.min()) - 5,
    max(TI_C.max(), TV_C.max()) + 8,
)

ax_P.set_ylim(0, 107)

ax_M.set_ylim(
    min(mI.min(), mV.min()) * 0.90,
    max(mI.max(), mV.max()) * 1.06,
)

ax_Z.set_ylim(
    min(ZI.min(), ZV.min()) - 0.05,
    max(ZI.max(), ZV.max()) + 0.05,
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
    bbox_to_anchor=(0.5, 0.018),
)


def update(k):
    im_I.set_data(field_I[k])
    im_V.set_data(field_V[k])

    sl = slice(0, k + 1)

    line_TI.set_data(time[sl], TI_C[sl])
    line_TV.set_data(time[sl], TV_C[sl])

    line_PI.set_data(time[sl], PI_bar[sl])
    line_PV.set_data(time[sl], PV_bar[sl])

    line_MI.set_data(time[sl], mI[sl])
    line_MV.set_data(time[sl], mV[sl])

    line_ZI.set_data(time[sl], ZI[sl])
    line_ZV.set_data(time[sl], ZV[sl])

    kpi_time.set_text(
        f"{time[k]:.1f} s"
    )

    kpi_pressure.set_text(
        f"{PI_bar[k]:.1f} / {PV_bar[k]:.1f} bar"
    )

    kpi_temperature.set_text(
        f"{TI_C[k]:.1f} / {TV_C[k]:.1f} °C"
    )

    kpi_Z.set_text(
        f"{ZV[k]:.4f}"
    )

    text_I = (
        f"T = {TI_C[k]:6.1f} °C\n"
        f"P = {PI_bar[k]:6.1f} bar\n"
        f"m = {mI[k]:6.3f} kg"
    )

    if time[k] >= tf_ideal:
        text_I += "\n100 bar alcanzados"

    text_V = (
        f"T = {TV_C[k]:6.1f} °C\n"
        f"P = {PV_bar[k]:6.1f} bar\n"
        f"m = {mV[k]:6.3f} kg"
    )

    if time[k] >= tf_vdw:
        text_V += "\n100 bar alcanzados"

    status_I.set_text(text_I)
    status_V.set_text(text_V)

    return (
        im_I,
        im_V,
        line_TI,
        line_TV,
        line_PI,
        line_PV,
        line_MI,
        line_MV,
        line_ZI,
        line_ZV,
        kpi_time,
        kpi_pressure,
        kpi_temperature,
        kpi_Z,
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
print("=" * 72)
print("ICF130 — TEMA 7: LLENADO TRANSIENTE")
print("=" * 72)

print()
print("CONDICIONES")
print("-" * 72)
print(f"Volumen tanque       = {V * 1000:.1f} L")
print(f"Presión inicial      = {P0 / 1e5:.2f} bar")
print(f"Temperatura inicial  = {T0:.2f} K")
print(f"Presión suministro   = {P_in / 1e5:.2f} bar")
print(f"Caudal               = {mdot:.4f} kg/s")
print(f"Presión objetivo     = {P_target / 1e5:.2f} bar")

print()
print("GAS IDEAL")
print("-" * 72)
print(f"Tiempo final         = {tf_ideal:.3f} s")
print(f"Temperatura final    = {TI_C[-1]:.3f} °C")
print(f"Masa final           = {mI[-1]:.5f} kg")
print(f"Factor Z final       = {ZI[-1]:.6f}")

print()
print("VAN DER WAALS")
print("-" * 72)
print(f"Tiempo final         = {tf_vdw:.3f} s")
print(f"Temperatura final    = {TV_C[-1]:.3f} °C")
print(f"Masa final           = {mV[-1]:.5f} kg")
print(f"Factor Z final       = {ZV[-1]:.6f}")

print()
print("DIFERENCIAS")
print("-" * 72)
print(f"Δ tiempo             = {tf_vdw - tf_ideal:+.3f} s")
print(
    f"Δ temperatura        = "
    f"{TV_C[-1] - TI_C[-1]:+.3f} °C"
)
print(f"Δ masa               = {mV[-1] - mI[-1]:+.5f} kg")

print()
print("VALIDACIÓN")
print("-" * 72)
print(f"Error analítico T ideal = {error_T_ideal:.3e} K")
print(f"Error energía ideal     = {energy_error_I:.3e}")
print(f"Error energía VdW       = {energy_error_V:.3e}")

print()
print(f"Informe TXT: {ruta_informe}")
print("=" * 72)



if SAVE_STATIC_PNG:
    update(N_FRAMES - 1)

    ruta_png = DIRECTORIO_SALIDA / "ideal_vs_vdw_sim_final.png"

    plt.savefig(
        ruta_png,
        dpi=200,
        facecolor=BG,
        bbox_inches="tight",
    )

    print(f"\nPNG generado: {ruta_png}")


if SAVE_MP4:
    try:
        from matplotlib.animation import FFMpegWriter

        ruta_mp4 = (
            DIRECTORIO_SALIDA
            / "ideal_vs_vdw_sim.mp4"
        )

        writer = FFMpegWriter(
            fps=FPS,
            bitrate=5000,
            metadata={
                "title": "ICF130 Tema 7 - Llenado transiente",
                "artist": "Proyecto Computacional",
            },
        )

        animation.save(
            str(ruta_mp4),
            writer=writer,
            dpi=150,
        )

        print(f"\nMP4 generado: {ruta_mp4}")

    except Exception as error:
        print("\nNo fue posible crear el MP4.")
        print(error)


if SAVE_GIF:
    ruta_gif = (
        DIRECTORIO_SALIDA
        / "ideal_vs_vdw.gif"
    )

    animation.save(
        str(ruta_gif),
        writer=PillowWriter(fps=18),
        dpi=105,
    )

    print(f"\nGIF generado: {ruta_gif}")


print(f"\nDirectorio de salida: {DIRECTORIO_SALIDA}")

plt.show()