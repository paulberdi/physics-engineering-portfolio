# ======================================================================
# ICF130 — PROYECTO COMPUTACIONAL
# TEMA 7: LLENADO TRANSIENTE DE UN DEPÓSITO CON GAS REAL
#
# VERSIÓN 3
# CAUDAL COMPRESIBLE VARIABLE + PARED TÉRMICA
# ======================================================================
#
# PROGRESIÓN DEL MODELO
#
# V1
#   - caudal constante
#   - llenado adiabático
#   - gas ideal vs Van der Waals
#
# V2
#   - caudal variable
#   - flujo estrangulado / subsónico
#   - llenado adiabático
#
# V3
#   - conserva el modelo de caudal de V2
#   - incorpora una pared con capacidad térmica
#   - transferencia de calor gas -> pared
#   - transferencia de calor pared -> ambiente
#
# VARIABLES DINÁMICAS:
#
#   m(t)   masa dentro del tanque
#   Tg(t)  temperatura promedio del gas
#   Tw(t)  temperatura promedio de la pared
#
# MODELO:
#
#   dm/dt = mdot
#
#   dUg/dt =
#       mdot*h_in - Qdot_gw
#
#   Cw*dTw/dt =
#       Qdot_gw - Qdot_wa
#
# donde:
#
#   Qdot_gw = UA_gw*(Tg-Tw)
#
#   Qdot_wa = UA_wa*(Tw-Tamb)
#
# IMPORTANTE:
#
# El modelo termodinámico es 0D.
# El mapa 2D mostrado en la animación es únicamente ilustrativo.
# No corresponde a CFD ni retroalimenta las EDO.
#
# ======================================================================


# ======================================================================
# 0. LIBRERÍAS
# ======================================================================

from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from matplotlib.animation import (
    FuncAnimation,
    FFMpegWriter,
    PillowWriter
)


# ======================================================================
# 1. DIRECTORIO Y OPCIONES
# ======================================================================

DIRECTORIO_SALIDA = (
    Path(__file__).resolve().parent
    if "__file__" in globals()
    else Path.cwd()
)


SAVE_TXT = True
SAVE_CSV = True

SAVE_PNG = True
SAVE_MP4 = True
SAVE_GIF = False

# Comparación V2 vs V3 en una figura aparte.
# No aparece en el dashboard principal.

SAVE_COMPARACION_V2 = True


FPS = 30

N_FRAMES = 240


# ======================================================================
# 2. DATOS FÍSICOS
# ======================================================================

# ----------------------------------------------------------------------
# TANQUE
# ----------------------------------------------------------------------

V = 0.050
# m^3
# 50 litros


# ----------------------------------------------------------------------
# ESTADO INICIAL
# ----------------------------------------------------------------------

P0 = 1.0e5
# Pa = 1 bar

T0 = 300.0
# K


# ----------------------------------------------------------------------
# SUMINISTRO
# ----------------------------------------------------------------------

P_SUPPLY = 200.0e5
# Pa = 200 bar

T_SUPPLY = 300.0
# K


# ----------------------------------------------------------------------
# PRESIÓN OBJETIVO
# ----------------------------------------------------------------------

P_TARGET = 190.0e5
# Pa = 190 bar


# ======================================================================
# 3. MODELO DE CAUDAL
# ======================================================================

MDOT_REFERENCIA = 0.020
# kg/s

CD = 0.80


# ======================================================================
# 4. PARED Y AMBIENTE — NUEVO EN V3
# ======================================================================

T_AMBIENTE = 300.0
# K

T_PARED_INICIAL = 300.0
# K


# Masa efectiva de pared

MASA_PARED = 8.0
# kg


# Capacidad calorífica específica efectiva

CP_PARED = 500.0
# J/(kg K)


# Conductancia térmica efectiva gas -> pared

UA_GAS_PARED = 8.0
# W/K


# Conductancia térmica efectiva pared -> ambiente

UA_PARED_AMBIENTE = 3.0
# W/K


# Capacidad térmica total de la pared

CAPACIDAD_PARED = (
    MASA_PARED
    * CP_PARED
)
# J/K


# ======================================================================
# 5. PROPIEDADES DEL NITRÓGENO
# ======================================================================

M = 28.0134e-3
# kg/mol


R_U = 8.314462618
# J/(mol K)


R = R_U / M
# J/(kg K)


# ----------------------------------------------------------------------
# Propiedades caloríficas constantes
#
# Se mantienen como en V1/V2 para que la única nueva física de V3
# sea el intercambio térmico con la pared.
# ----------------------------------------------------------------------

CP = 1039.0
# J/(kg K)


CV = CP - R


GAMMA = CP / CV


# ======================================================================
# 6. VAN DER WAALS
# ======================================================================

TC = 126.19
# K


PC = 33.978e5
# Pa


A_VDW = (
    (27.0 / 64.0)
    * R**2
    * TC**2
    / PC
)


B_VDW = (
    R
    * TC
    / (8.0 * PC)
)


# ======================================================================
# 7. FLUJO COMPRESIBLE
# ======================================================================

# Razón crítica de presiones

RAZON_CRITICA = (
    2.0
    / (GAMMA + 1.0)
)**(
    GAMMA
    / (GAMMA - 1.0)
)


P_CRITICA = (
    RAZON_CRITICA
    * P_SUPPLY
)


# Factor del flujo estrangulado

FACTOR_CHOKED = (
    np.sqrt(
        GAMMA
        / (R * T_SUPPLY)
    )
    *
    (
        2.0
        / (GAMMA + 1.0)
    )**(
        (GAMMA + 1.0)
        /
        (
            2.0
            * (GAMMA - 1.0)
        )
    )
)


# Área efectiva calibrada para que:
#
# mdot_choked = 0.020 kg/s

AREA_ORIFICIO = (
    MDOT_REFERENCIA
    /
    (
        CD
        * P_SUPPLY
        * FACTOR_CHOKED
    )
)


DIAMETRO_ORIFICIO = np.sqrt(
    4.0
    * AREA_ORIFICIO
    / np.pi
)


# ======================================================================
# 8. ESTILO — MISMO LENGUAJE VISUAL DE V1 / V2
# ======================================================================

BG = "#070A0F"

PANEL = "#10151D"

GRID = "#2A3340"


IDEAL_COLOR = "#47D7FF"

VDW_COLOR = "#FF9F43"


WHITE = "#F4F7FB"

MUTED = "#AAB4C3"


# Malla del mapa visual

NX = 150

NY = 80


plt.rcParams.update({

    "font.family":
    "DejaVu Sans",

    "font.size":
    10,

    "figure.facecolor":
    BG,

    "axes.facecolor":
    PANEL,

    "axes.edgecolor":
    "#303B48",

    "axes.labelcolor":
    "#D6DCE5",

    "axes.titlecolor":
    WHITE,

    "xtick.color":
    "#AEB8C6",

    "ytick.color":
    "#AEB8C6",

    "text.color":
    WHITE
})


# ======================================================================
# 9. VALIDACIÓN DE PARÁMETROS
# ======================================================================

def validar_parametros():

    if not (
        0.0
        < P0
        < P_TARGET
        < P_SUPPLY
    ):

        raise ValueError(
            "Debe cumplirse "
            "0 < P0 < P_TARGET < P_SUPPLY."
        )


    if V <= 0.0:

        raise ValueError(
            "El volumen debe ser positivo."
        )


    if min(
        T0,
        T_SUPPLY,
        T_AMBIENTE,
        T_PARED_INICIAL
    ) <= 0.0:

        raise ValueError(
            "Las temperaturas absolutas "
            "deben ser positivas."
        )


    if min(
        M,
        R,
        CP,
        CV,
        MASA_PARED,
        CP_PARED,
        MDOT_REFERENCIA,
        CD
    ) <= 0.0:

        raise ValueError(
            "Los parámetros físicos "
            "deben ser positivos."
        )


    if min(
        UA_GAS_PARED,
        UA_PARED_AMBIENTE
    ) < 0.0:

        raise ValueError(
            "Las conductancias térmicas "
            "no pueden ser negativas."
        )


# ======================================================================
# 10. ECUACIONES DE ESTADO
# ======================================================================

def presion(
    modelo,
    masa,
    temperatura
):

    # ------------------------------------------------------------------
    # GAS IDEAL
    # ------------------------------------------------------------------

    if modelo == "Ideal":

        return (
            masa
            * R
            * temperatura
            / V
        )


    # ------------------------------------------------------------------
    # VAN DER WAALS
    # ------------------------------------------------------------------

    volumen = (
        V / masa
    )


    if np.any(
        volumen <= B_VDW
    ):

        raise ValueError(
            "Van der Waals fuera de dominio: "
            "v <= b."
        )


    return (
        R
        * temperatura
        / (
            volumen
            - B_VDW
        )
        -
        A_VDW
        / volumen**2
    )


# ======================================================================
# 11. VOLUMEN ESPECÍFICO VAN DER WAALS
# ======================================================================

def volumen_vdw(
    P,
    T
):

    def residuo(v):

        return (
            R
            * T
            / (
                v
                - B_VDW
            )
            -
            A_VDW
            / v**2
            -
            P
        )


    v_min = (
        B_VDW
        * (
            1.0
            + 1e-10
        )
    )


    v_max = max(

        100.0
        * R
        * T
        / P,

        1.0e4
        * B_VDW
    )


    while (
        residuo(
            v_max
        )
        > 0.0
    ):

        v_max *= 2.0


    return brentq(

        residuo,

        v_min,

        v_max,

        xtol=1e-12,

        rtol=1e-12
    )


# ======================================================================
# 12. ENERGÍA INTERNA TOTAL DEL GAS
# ======================================================================

def energia_gas(
    modelo,
    masa,
    temperatura
):

    # ------------------------------------------------------------------
    # IDEAL
    #
    # U = m cv T
    # ------------------------------------------------------------------

    if modelo == "Ideal":

        return (
            masa
            * CV
            * temperatura
        )


    # ------------------------------------------------------------------
    # VAN DER WAALS
    #
    # u = cv T - a/v
    #
    # U = m cv T - a m^2/V
    # ------------------------------------------------------------------

    return (
        masa
        * CV
        * temperatura
        -
        A_VDW
        * masa**2
        / V
    )


# ======================================================================
# 13. CAUDAL COMPRESIBLE
# ======================================================================

def caudal(
    P_tanque
):

    # No se modela flujo inverso

    if P_tanque >= P_SUPPLY:

        return 0.0


    razon = np.clip(

        P_tanque
        / P_SUPPLY,

        0.0,

        1.0
    )


    # ------------------------------------------------------------------
    # FLUJO ESTRANGULADO
    # ------------------------------------------------------------------

    if razon <= RAZON_CRITICA:

        return (
            CD
            * AREA_ORIFICIO
            * P_SUPPLY
            * FACTOR_CHOKED
        )


    # ------------------------------------------------------------------
    # FLUJO SUBSÓNICO
    # ------------------------------------------------------------------

    termino = (

        razon**(
            2.0
            / GAMMA
        )

        -

        razon**(
            (GAMMA + 1.0)
            / GAMMA
        )
    )


    termino = max(
        termino,
        0.0
    )


    return (

        CD

        * AREA_ORIFICIO

        * P_SUPPLY

        / np.sqrt(
            R
            * T_SUPPLY
        )

        * np.sqrt(

            2.0
            * GAMMA
            / (
                GAMMA
                - 1.0
            )

            * termino
        )
    )


# ======================================================================
# 14. RÉGIMEN DEL FLUJO
# ======================================================================

def regimen(
    P_tanque
):

    if P_tanque >= P_SUPPLY:

        return "SIN FLUJO"


    if P_tanque <= P_CRITICA:

        return "ESTRANGULADO"


    return "SUBSONICO"


# ======================================================================
# 15. SIMULACIÓN
# ======================================================================

def simular(
    modelo,
    UA_interior
):

    # ==================================================================
    # ESTADO INICIAL
    # ==================================================================

    if modelo == "Ideal":

        masa_inicial = (
            P0
            * V
            / (
                R
                * T0
            )
        )


        h_entrada = (
            CP
            * T_SUPPLY
        )


    else:

        v0 = volumen_vdw(
            P0,
            T0
        )


        masa_inicial = (
            V / v0
        )


        v_supply = volumen_vdw(
            P_SUPPLY,
            T_SUPPLY
        )


        u_supply = (
            CV
            * T_SUPPLY
            -
            A_VDW
            / v_supply
        )


        h_entrada = (
            u_supply
            +
            P_SUPPLY
            * v_supply
        )


    # ==================================================================
    # SISTEMA DE EDO
    #
    # estado =
    #
    # [
    #   masa,
    #   Tg,
    #   Tw,
    #   Ein,
    #   Qgw,
    #   Qwa
    # ]
    #
    # Convención:
    #
    # Qgw > 0:
    #   calor gas -> pared
    #
    # Qwa > 0:
    #   calor pared -> ambiente
    # ==================================================================

    def derivadas(
        tiempo,
        estado
    ):

        (
            masa,
            Tg,
            Tw,
            Ein,
            Qgw,
            Qwa
        ) = estado


        P = presion(
            modelo,
            masa,
            Tg
        )


        mdot = caudal(
            P
        )


        # --------------------------------------------------------------
        # GAS -> PARED
        # --------------------------------------------------------------

        qdot_gw = (
            UA_interior
            * (
                Tg
                - Tw
            )
        )


        # --------------------------------------------------------------
        # PARED -> AMBIENTE
        # --------------------------------------------------------------

        qdot_wa = (
            UA_PARED_AMBIENTE
            * (
                Tw
                - T_AMBIENTE
            )
        )


        # --------------------------------------------------------------
        # BALANCE DE MASA
        # --------------------------------------------------------------

        dm_dt = mdot


        # --------------------------------------------------------------
        # BALANCE DE ENERGÍA DEL GAS
        # --------------------------------------------------------------

        if modelo == "Ideal":

            # U = m cv Tg
            #
            # dU/dt =
            #
            # mdot*h_in
            # -
            # Qdot_gw
            #
            #
            # m cv dT/dt =
            #
            # mdot(
            #     h_in
            #     -
            #     cv*Tg
            # )
            #
            # -
            #
            # Qdot_gw

            dTg_dt = (

                mdot
                * (
                    h_entrada
                    -
                    CV
                    * Tg
                )

                -

                qdot_gw

            ) / (

                masa
                * CV
            )


        else:

            # Van der Waals:
            #
            # U =
            #
            # m cv Tg
            #
            # -
            #
            # a m^2/V
            #
            #
            # Derivando:
            #
            # m cv dT/dt =
            #
            # mdot(
            #   h_in
            #   -
            #   cv Tg
            #   +
            #   2 a m/V
            # )
            #
            # -
            # Qdot_gw

            dTg_dt = (

                mdot

                * (

                    h_entrada

                    -

                    CV
                    * Tg

                    +

                    2.0
                    * A_VDW
                    * masa
                    / V
                )

                -

                qdot_gw

            ) / (

                masa
                * CV
            )


        # --------------------------------------------------------------
        # BALANCE DE ENERGÍA DE LA PARED
        #
        # Cw dTw/dt =
        #
        # Qdot_gw
        # -
        # Qdot_wa
        # --------------------------------------------------------------

        dTw_dt = (

            qdot_gw
            -
            qdot_wa

        ) / CAPACIDAD_PARED


        # --------------------------------------------------------------
        # ACUMULADORES
        # --------------------------------------------------------------

        dEin_dt = (
            mdot
            * h_entrada
        )


        dQgw_dt = (
            qdot_gw
        )


        dQwa_dt = (
            qdot_wa
        )


        return [

            dm_dt,

            dTg_dt,

            dTw_dt,

            dEin_dt,

            dQgw_dt,

            dQwa_dt
        ]


    # ==================================================================
    # EVENTO
    #
    # Detener cuando:
    #
    # P = P_TARGET
    # ==================================================================

    def evento_objetivo(
        tiempo,
        estado
    ):

        masa = estado[0]

        Tg = estado[1]


        return (

            presion(
                modelo,
                masa,
                Tg
            )

            -

            P_TARGET
        )


    evento_objetivo.terminal = True

    evento_objetivo.direction = 1


    # ==================================================================
    # ESTADO INICIAL
    # ==================================================================

    estado_inicial = [

        masa_inicial,

        T0,

        T_PARED_INICIAL,

        0.0,

        0.0,

        0.0
    ]


    # ==================================================================
    # INTEGRACIÓN
    # ==================================================================

    solucion = solve_ivp(

        derivadas,

        (
            0.0,
            5000.0
        ),

        estado_inicial,

        method="RK45",

        events=evento_objetivo,

        dense_output=True,

        rtol=1e-9,

        atol=1e-11,

        max_step=0.5
    )


    if not solucion.success:

        raise RuntimeError(

            f"{modelo}: "
            f"{solucion.message}"
        )


    if len(
        solucion.t_events[0]
    ) == 0:

        raise RuntimeError(

            f"{modelo}: "
            "no alcanzó la presión objetivo."
        )


    tiempo_final = (
        solucion.t_events[0][0]
    )


    # ==================================================================
    # TRANSICIÓN ESTRANGULADO -> SUBSÓNICO
    # ==================================================================

    tiempo_transicion = None


    if (
        P0
        <
        P_CRITICA
        <
        P_TARGET
    ):

        def condicion_transicion(
            tiempo
        ):

            estado = solucion.sol(
                tiempo
            )


            return (

                presion(

                    modelo,

                    estado[0],

                    estado[1]

                )

                -

                P_CRITICA
            )


        tiempo_transicion = brentq(

            condicion_transicion,

            0.0,

            tiempo_final
        )


    return {

        "modelo":
        modelo,

        "sol":
        solucion,

        "tf":
        tiempo_final,

        "m0":
        masa_inicial,

        "h_in":
        h_entrada,

        "UA":
        UA_interior,

        "transicion":
        tiempo_transicion
    }


# ======================================================================
# 16. EVALUAR RESULTADO
# ======================================================================

def evaluar(
    resultado,
    tiempos
):

    tiempos = np.asarray(
        tiempos,
        dtype=float
    )


    tiempos_locales = np.minimum(

        tiempos,

        resultado["tf"]
    )


    valores = (
        resultado[
            "sol"
        ].sol(
            tiempos_locales
        )
    )


    masa = valores[0]

    Tg = valores[1]

    Tw = valores[2]

    Ein = valores[3]

    Qgw = valores[4]

    Qwa = valores[5]


    P = presion(

        resultado["modelo"],

        masa,

        Tg
    )


    mdot = np.array([

        caudal(
            Pi
        )

        for Pi
        in np.atleast_1d(P)
    ])


    qdot_gw = (

        resultado["UA"]

        * (
            Tg
            -
            Tw
        )
    )


    qdot_wa = (

        UA_PARED_AMBIENTE

        * (
            Tw
            -
            T_AMBIENTE
        )
    )


    Z = (

        P
        * V

        /

        (
            masa
            * R
            * Tg
        )
    )


    activo = (

        tiempos
        <=
        resultado["tf"]
        + 1e-12
    )


    return {

        "t_local":
        tiempos_locales,

        "m":
        masa,

        "Tg":
        Tg,

        "Tw":
        Tw,

        "P":
        P,

        "mdot":
        mdot,

        "Z":
        Z,

        "Ein":
        Ein,

        "Qgw":
        Qgw,

        "Qwa":
        Qwa,

        "qdot_gw":
        qdot_gw,

        "qdot_wa":
        qdot_wa,

        "activo":
        activo
    }


# ======================================================================
# 17. VALIDACIÓN DE ENERGÍA
# ======================================================================

def validar_balance(
    resultado
):

    tiempos = np.linspace(

        0.0,

        resultado["tf"],

        2000
    )


    datos = evaluar(
        resultado,
        tiempos
    )


    # ------------------------------------------------------------------
    # CAMBIO DE ENERGÍA DEL GAS
    # ------------------------------------------------------------------

    Ug = energia_gas(

        resultado["modelo"],

        datos["m"],

        datos["Tg"]
    )


    Ug0 = energia_gas(

        resultado["modelo"],

        resultado["m0"],

        T0
    )


    delta_Ug = (
        Ug
        -
        Ug0
    )


    # ------------------------------------------------------------------
    # CAMBIO DE ENERGÍA DE LA PARED
    # ------------------------------------------------------------------

    delta_Uw = (

        CAPACIDAD_PARED

        * (
            datos["Tw"]
            -
            T_PARED_INICIAL
        )
    )


    # ------------------------------------------------------------------
    # BALANCE DEL GAS
    #
    # Delta Ug =
    #
    # Ein
    # -
    # Qgw
    # ------------------------------------------------------------------

    residuo_gas = (

        delta_Ug

        -

        datos["Ein"]

        +

        datos["Qgw"]
    )


    # ------------------------------------------------------------------
    # BALANCE DE LA PARED
    #
    # Delta Uw =
    #
    # Qgw
    # -
    # Qwa
    # ------------------------------------------------------------------

    residuo_pared = (

        delta_Uw

        -

        datos["Qgw"]

        +

        datos["Qwa"]
    )


    # ------------------------------------------------------------------
    # BALANCE TOTAL
    #
    # Delta(Ug+Uw) =
    #
    # Ein
    # -
    # Qwa
    # ------------------------------------------------------------------

    residuo_total = (

        delta_Ug

        +

        delta_Uw

        -

        datos["Ein"]

        +

        datos["Qwa"]
    )


    escala = np.maximum(

        np.abs(
            datos["Ein"]
        ),

        1.0
    )


    return {

        "gas":

        np.max(

            np.abs(
                residuo_gas
            )

            / escala
        ),


        "pared":

        np.max(

            np.abs(
                residuo_pared
            )

            / escala
        ),


        "total":

        np.max(

            np.abs(
                residuo_total
            )

            / escala
        ),


        "total_J":

        np.max(

            np.abs(
                residuo_total
            )
        )
    }


# ======================================================================
# 18. VALIDACIÓN ANALÍTICA DEL CASO IDEAL ADIABÁTICO
# ======================================================================

def error_analitico_ideal(
    referencia
):

    tiempos = np.linspace(

        0.0,

        referencia["tf"],

        2000
    )


    estado = (
        referencia[
            "sol"
        ].sol(
            tiempos
        )
    )


    masa = estado[0]

    T_numerica = estado[1]


    # Incluso con mdot variable:
    #
    # dT/dm =
    #
    # (gamma*T_supply - T)/m
    #
    # por lo que la relación T(m) sigue siendo válida.

    T_analitica = (

        GAMMA
        * T_SUPPLY

        +

        referencia["m0"]
        / masa

        * (

            T0

            -

            GAMMA
            * T_SUPPLY
        )
    )


    return np.max(

        np.abs(

            T_numerica

            -

            T_analitica
        )
    )


# ======================================================================
# 19. CAMPO TÉRMICO ILUSTRATIVO
# ======================================================================
#
# Solo se usa para conservar el lenguaje visual de V1/V2.
#
# NO es CFD.
#
# NO se resuelven:
#
#   - Navier-Stokes
#   - campo real de velocidades
#   - campo físico T(x,y,t)
#
# La termodinámica real utilizada por el modelo es Tg(t).
#
# ======================================================================

def crear_campos_visuales(
    temperaturas,
    activos
):

    temperaturas = np.asarray(
        temperaturas
    )


    activos = np.asarray(
        activos
    )


    campo = np.full(

        (
            NY,
            NX
        ),

        T0
        -
        273.15,

        dtype=np.float32
    )


    yy, xx = np.mgrid[
        0:NY,
        0:NX
    ]


    campos = np.zeros(

        (
            len(
                temperaturas
            ),
            NY,
            NX
        ),

        dtype=np.float32
    )


    for k, objetivo in enumerate(
        temperaturas
    ):

        # --------------------------------------------------------------
        # Una vez que el modelo llega al objetivo,
        # congelamos también la visualización.
        # --------------------------------------------------------------

        if (
            k > 0
            and
            not activos[k]
        ):

            campos[k] = (
                campos[
                    k - 1
                ]
            )

            continue


        delta = max(

            objetivo

            -

            (
                T0
                -
                273.15
            ),

            0.0
        )


        centro_y = (

            NY / 2

            +

            0.05
            * NY

            * np.sin(

                2.0
                * np.pi
                * k
                / 70.0
            )
        )


        fuente = np.exp(

            -

            (
                (
                    xx
                    - 8.0
                )
                / 18.0
            )**2

            -

            (
                (
                    yy
                    - centro_y
                )
                / 9.0
            )**2
        )


        pluma = np.exp(

            -

            (
                (
                    xx
                    - 32.0
                )
                / 48.0
            )**2

            -

            (
                (
                    yy
                    - centro_y
                )
                / 18.0
            )**2
        )


        objetivo_local = (

            objetivo

            +

            0.20
            * delta
            * fuente

            +

            0.07
            * delta
            * pluma
        )


        # Cinco pasos internos puramente visuales

        for _ in range(5):

            nuevo = (
                campo.copy()
            )


            # Difusión visual

            nuevo[
                1:-1,
                1:-1
            ] = (

                0.70

                * campo[
                    1:-1,
                    1:-1
                ]

                +

                0.075

                * (

                    campo[
                        1:-1,
                        2:
                    ]

                    +

                    campo[
                        1:-1,
                        :-2
                    ]

                    +

                    campo[
                        2:,
                        1:-1
                    ]

                    +

                    campo[
                        :-2,
                        1:-1
                    ]
                )
            )


            # Advección visual hacia la derecha

            nuevo[
                :,
                1:
            ] = (

                0.955

                * nuevo[
                    :,
                    1:
                ]

                +

                0.045

                * nuevo[
                    :,
                    :-1
                ]
            )


            # Llevar el promedio hacia Tg

            nuevo += (

                0.055

                * (
                    objetivo
                    -
                    nuevo
                )
            )


            # Efecto local de entrada

            nuevo += (

                0.10

                * fuente

                * (
                    objetivo_local
                    -
                    nuevo
                )
            )


            nuevo += (

                0.03

                * pluma

                * (
                    objetivo_local
                    -
                    nuevo
                )
            )


            # Bordes

            nuevo[:, 0] = (
                nuevo[:, 1]
            )

            nuevo[:, -1] = (
                nuevo[:, -2]
            )

            nuevo[0, :] = (
                nuevo[1, :]
            )

            nuevo[-1, :] = (
                nuevo[-2, :]
            )


            campo = nuevo


        campos[k] = campo


    return campos


# ======================================================================
# 20. DASHBOARD PRINCIPAL
# ======================================================================
# ======================================================================
# 20. DASHBOARD PRINCIPAL — LAYOUT CORREGIDO
# ======================================================================

def crear_dashboard(
    resultados
):

    tiempo_maximo = max(
        resultado["tf"]
        for resultado in resultados
    )

    tiempos = np.linspace(
        0.0,
        tiempo_maximo,
        N_FRAMES
    )

    datos = [
        evaluar(
            resultado,
            tiempos
        )
        for resultado in resultados
    ]

    colores = [
        IDEAL_COLOR,
        VDW_COLOR
    ]


    # ==================================================================
    # MAPAS VISUALES
    # ==================================================================

    campos = [
        crear_campos_visuales(
            serie["Tg"] - 273.15,
            serie["activo"]
        )
        for serie in datos
    ]

    Tmin = min(
        np.min(campo)
        for campo in campos
    )

    Tmax = max(
        np.max(campo)
        for campo in campos
    )


    # ==================================================================
    # FIGURA
    #
    # Más alta que antes para evitar superposiciones.
    # ==================================================================

    figura = plt.figure(
        figsize=(16, 12.3),
        facecolor=BG
    )


    # IMPORTANTE:
    #
    # top = 0.76
    #
    # Todo lo que está sobre 0.76 queda reservado
    # exclusivamente para el encabezado.

    rejilla = figura.add_gridspec(
        3,
        12,

        left=0.055,
        right=0.935,

        bottom=0.105,
        top=0.795,

        hspace=0.62,
        wspace=0.95,

        height_ratios=[
            1.75,
            0.9,
            0.9
        ]
    )


    # ==================================================================
    # ENCABEZADO
    # ==================================================================

    figura.text(
        0.055,
        0.970,

        "VERSIÓN 3",

        fontsize=11,
        color=MUTED
    )


    figura.text(
        0.055,
        0.940,

        "LLENADO TRANSIENTE DE UN DEPÓSITO",

        fontsize=21,
        weight="bold",
        color=WHITE
    )


    figura.text(
        0.055,
        0.908,

        "Gas ideal vs. Van der Waals"
        "   |   flujo compresible variable"
        "   |   intercambio térmico con la pared",

        fontsize=11,
        color=MUTED
    )


    # Línea de condiciones más corta.
    # Ya no colocamos aquí el texto "no CFD".

    figura.text(
        0.055,
        0.875,

        f"Depósito {V*1000:.0f} L"
        f"   |   Inicio {P0/1e5:.0f} bar"
        f"   |   Suministro {P_SUPPLY/1e5:.0f} bar"
        f"   |   Objetivo {P_TARGET/1e5:.0f} bar"
        f"   |   Ambiente {T_AMBIENTE-273.15:.1f} °C",

        fontsize=9.5,
        color=MUTED
    )


    reloj = figura.text(
        0.935,
        0.962,

        "",

        ha="right",

        fontsize=14,
        weight="bold",
        color=WHITE
    )


    # ==================================================================
    # TANQUES
    # ==================================================================

    ejes_tanque = [

        figura.add_subplot(
            rejilla[
                0,
                0:6
            ]
        ),

        figura.add_subplot(
            rejilla[
                0,
                6:12
            ]
        )
    ]


    imagenes = []
    textos_estado = []
    titulos_tanque = []


    for i, eje in enumerate(
        ejes_tanque
    ):

        eje.set_facecolor(
            PANEL
        )

        imagen = eje.imshow(
            campos[i][0],

            origin="lower",
            aspect="auto",

            cmap="inferno",

            vmin=Tmin,
            vmax=Tmax
        )

        imagenes.append(
            imagen
        )

        eje.set_xticks([])
        eje.set_yticks([])

        for borde in eje.spines.values():

            borde.set_color(
                "#303B48"
            )

            borde.set_linewidth(
                1.1
            )


        # El título ahora queda claramente separado
        # del encabezado general.

        titulo = eje.set_title(
            "",

            fontsize=14,
            weight="bold",

            color=colores[i],

            pad=10
        )

        titulos_tanque.append(
            titulo
        )


        # Cuadro de información ligeramente más pequeño.

        texto = eje.text(
            0.020,
            0.055,

            "",

            transform=eje.transAxes,

            fontsize=10.2,

            color=WHITE,

            verticalalignment="bottom",

            bbox=dict(
                facecolor=(
                    0,
                    0,
                    0,
                    0.72
                ),

                edgecolor="none",

                pad=6
            )
        )

        textos_estado.append(
            texto
        )


    # ==================================================================
    # BARRA DE COLOR
    # ==================================================================

    posicion = (
        ejes_tanque[1]
        .get_position()
    )


    eje_color = figura.add_axes([
        0.946,
        posicion.y0,
        0.008,
        posicion.height
    ])


    barra = figura.colorbar(
        imagenes[0],
        cax=eje_color
    )


    barra.set_label(
        "T ilustrativa [°C]",

        fontsize=8.5,
        color=WHITE,

        labelpad=8
    )


    barra.ax.tick_params(
        colors=MUTED,
        labelsize=7.5
    )


    # ==================================================================
    # GRÁFICOS
    # ==================================================================

    ax_T = figura.add_subplot(
        rejilla[
            1,
            0:6
        ]
    )


    ax_P = figura.add_subplot(
        rejilla[
            1,
            6:12
        ]
    )


    ax_mdot = figura.add_subplot(
        rejilla[
            2,
            0:6
        ]
    )


    ax_Q = figura.add_subplot(
        rejilla[
            2,
            6:12
        ]
    )


    graficos = [
        ax_T,
        ax_P,
        ax_mdot,
        ax_Q
    ]


    titulos = [
        "Temperatura",
        "Presión del depósito",
        "Caudal másico",
        "Transferencia de calor"
    ]


    unidades = [
        "T [°C]",
        "P [bar]",
        "ṁ [kg/s]",
        "Potencia [W]"
    ]


    for eje, titulo, unidad in zip(
        graficos,
        titulos,
        unidades
    ):

        eje.set_facecolor(
            PANEL
        )

        eje.set_title(
            titulo,

            fontsize=11,
            pad=11,

            color=WHITE
        )

        eje.set_xlabel(
            "Tiempo [s]",
            fontsize=9,

            labelpad=5
        )

        eje.set_ylabel(
            unidad,
            fontsize=9,

            labelpad=5
        )

        eje.set_xlim(
            0.0,
            tiempo_maximo
        )

        eje.grid(
            color=GRID,
            alpha=0.42,
            linewidth=0.8
        )

        eje.tick_params(
            labelsize=8
        )

        for borde in eje.spines.values():

            borde.set_color(
                "#303B48"
            )


    # ==================================================================
    # REFERENCIAS
    # ==================================================================

    # Temperatura ambiente

    ax_T.axhline(
        T_AMBIENTE - 273.15,

        color="#ADB8C6",

        linestyle=":",

        linewidth=1.1,

        alpha=0.60
    )


    # Presión objetivo

    ax_P.axhline(
        P_TARGET / 1e5,

        color=WHITE,

        linestyle=":",

        linewidth=1.2,

        alpha=0.55
    )


    # Presión crítica

    ax_P.axhline(
        P_CRITICA / 1e5,

        color="#E9C46A",

        linestyle="--",

        linewidth=1.1,

        alpha=0.60
    )


    ax_mdot.set_ylim(
        0.0,
        MDOT_REFERENCIA * 1.18
    )


    ax_P.set_ylim(
        0.0,
        P_TARGET / 1e5 * 1.08
    )


    # ==================================================================
    # CURVAS
    # ==================================================================

    artistas = []


    for i, (
        resultado,
        serie
    ) in enumerate(
        zip(
            resultados,
            datos
        )
    ):

        color = colores[i]

        modelo = resultado[
            "modelo"
        ]


        # ==============================================================
        # TEMPERATURA GAS
        # ==============================================================

        linea_Tg, = ax_T.plot(
            [],
            [],

            color=color,

            linewidth=2.3,

            label=f"Gas {modelo}"
        )


        # ==============================================================
        # TEMPERATURA PARED
        # ==============================================================

        linea_Tw, = ax_T.plot(
            [],
            [],

            color=color,

            linestyle="--",

            linewidth=1.7,

            alpha=0.85,

            label=f"Pared {modelo}"
        )


        # ==============================================================
        # PRESIÓN
        # ==============================================================

        linea_P, = ax_P.plot(
            [],
            [],

            color=color,

            linewidth=2.3,

            label=modelo
        )


        # ==============================================================
        # CAUDAL
        # ==============================================================

        linea_mdot, = ax_mdot.plot(
            [],
            [],

            color=color,

            linewidth=2.3,

            label=modelo
        )


        # ==============================================================
        # CALOR GAS -> PARED
        # ==============================================================

        linea_Q, = ax_Q.plot(
            [],
            [],

            color=color,

            linewidth=2.3,

            label=modelo
        )


        artistas.append({

            "Tg":
            (
                linea_Tg,
                serie["Tg"] - 273.15
            ),

            "Tw":
            (
                linea_Tw,
                serie["Tw"] - 273.15
            ),

            "P":
            (
                linea_P,
                serie["P"] / 1e5
            ),

            "mdot":
            (
                linea_mdot,
                serie["mdot"]
            ),

            "Q":
            (
                linea_Q,
                serie["qdot_gw"]
            )
        })


    # ==================================================================
    # LEYENDAS
    # ==================================================================

    # La de temperatura tiene cuatro entradas.
    # La dejamos arriba a la izquierda en dos columnas.

    ax_T.legend(
        loc="upper left",

        ncol=2,

        fontsize=7.3,

        frameon=True,

        framealpha=0.86,

        facecolor=PANEL,

        edgecolor="#303B48",

        borderpad=0.5,

        handlelength=2.0
    )


    # El resto tiene solo dos entradas.

    for eje in [
        ax_P,
        ax_mdot,
        ax_Q
    ]:

        eje.legend(
            loc="lower right",

            fontsize=8,

            frameon=True,

            framealpha=0.82,

            facecolor=PANEL,

            edgecolor="#303B48"
        )


    # ==================================================================
    # ESCALAS DE TEMPERATURA Y CALOR
    # ==================================================================

    temperaturas = []

    potencias = []


    for serie in datos:

        temperaturas.extend(
            serie["Tg"][
                serie["activo"]
            ]
            - 273.15
        )

        temperaturas.extend(
            serie["Tw"][
                serie["activo"]
            ]
            - 273.15
        )

        potencias.extend(
            serie["qdot_gw"][
                serie["activo"]
            ]
        )


    Tmin_grafico = min(
        temperaturas
    )

    Tmax_grafico = max(
        temperaturas
    )


    margen_T = max(
        5.0,
        0.08
        * (
            Tmax_grafico
            -
            Tmin_grafico
        )
    )


    ax_T.set_ylim(
        Tmin_grafico - margen_T,
        Tmax_grafico + margen_T
    )


    Qmin = min(
        potencias
    )

    Qmax = max(
        potencias
    )


    margen_Q = max(
        25.0,
        0.10
        * (
            Qmax
            -
            Qmin
        )
    )


    ax_Q.set_ylim(
        min(
            0.0,
            Qmin - margen_Q
        ),
        Qmax + margen_Q
    )


    # ==================================================================
    # PIE INFERIOR
    #
    # Ahora queda fuera de la zona de los ejes.
    # ==================================================================

    figura.text(
        0.055,
        0.060,

        "Temperatura: línea continua = gas"
        "   ·   línea discontinua = pared"
        "   ·   presión: línea amarilla = condición crítica",

        fontsize=8.7,

        color=MUTED
    )


    figura.text(
        0.055,
        0.032,

        "Mapa térmico ilustrativo (no CFD). "
        "El modelo termodinámico es concentrado (0D).",

        fontsize=8.7,

        color=MUTED
    )


    # ==================================================================
    # ACTUALIZACIÓN DE LA ANIMACIÓN
    # ==================================================================

    def actualizar(
        k
    ):

        tiempo_actual = tiempos[k]


        reloj.set_text(
            f"t = {tiempo_actual:6.1f}"
            f" / {tiempo_maximo:.1f} s"
        )


        # ==============================================================
        # TANQUES
        # ==============================================================

        for i, (
            resultado,
            serie
        ) in enumerate(
            zip(
                resultados,
                datos
            )
        ):

            imagenes[i].set_data(
                campos[i][k]
            )


            P_actual = (
                serie["P"][k]
                / 1e5
            )


            Tg_actual = (
                serie["Tg"][k]
                - 273.15
            )


            Tw_actual = (
                serie["Tw"][k]
                - 273.15
            )


            masa_actual = (
                serie["m"][k]
            )


            qdot_actual = (
                serie["qdot_gw"][k]
            )


            if (
                tiempo_actual
                <
                resultado["tf"]
            ):

                estado = "LLENANDO"

                mdot_texto = (
                    f"ṁ = "
                    f"{serie['mdot'][k]:.4f} kg/s"
                )

            else:

                estado = (
                    "OBJETIVO ALCANZADO"
                )

                mdot_texto = (
                    "llenado finalizado"
                )


            regimen_actual = regimen(
                serie["P"][k]
            )


            titulos_tanque[i].set_text(

                f"{resultado['modelo'].upper()}"
                f"  ·  {estado}"
                f"  ·  {regimen_actual}"
            )


            # Tres líneas limpias.
            # Evitamos llenar demasiado el tanque.

            textos_estado[i].set_text(

                f"T gas = {Tg_actual:5.1f} °C"
                f"      "
                f"T pared = {Tw_actual:5.1f} °C\n"

                f"P = {P_actual:5.1f} bar"
                f"      "
                f"m = {masa_actual:5.3f} kg\n"

                f"{mdot_texto}"
                f"      "
                f"Q̇g→w = {qdot_actual:5.1f} W"
            )


        # ==============================================================
        # CURVAS
        # ==============================================================

        for i, resultado in enumerate(
            resultados
        ):

            indice_final = (
                k + 1
            )


            mascara = (

                tiempos[
                    :indice_final
                ]

                <=

                resultado["tf"]
                + 1e-12
            )


            t_plot = (

                tiempos[
                    :indice_final
                ][
                    mascara
                ]
            )


            for (
                linea,
                valores
            ) in artistas[i].values():

                linea.set_data(

                    t_plot,

                    valores[
                        :indice_final
                    ][
                        mascara
                    ]
                )


        return tuple(

            imagenes

            +

            textos_estado

            +

            titulos_tanque

            +

            [reloj]

            +

            [
                elemento[0]

                for grupo
                in artistas

                for elemento
                in grupo.values()
            ]
        )


    # ==================================================================
    # ANIMACIÓN
    # ==================================================================

    animacion = FuncAnimation(
        figura,

        actualizar,

        frames=N_FRAMES,

        interval=1000 / FPS,

        blit=False,

        repeat=True
    )


    # Frame final para PNG

    actualizar(
        N_FRAMES - 1
    )


    return (
        figura,
        animacion,
        tiempos,
        datos
    )
# ======================================================================
# 21. COMPARACIÓN SIMPLE V2 VS V3
# ======================================================================

def crear_comparacion_v2_v3(
    resultados,
    referencias
):

    figura, eje = plt.subplots(

        figsize=(
            10,
            6
        ),

        facecolor=BG
    )


    eje.set_facecolor(
        PANEL
    )


    colores = {

        "Ideal":
        IDEAL_COLOR,

        "VdW":
        VDW_COLOR
    }


    for resultado, referencia in zip(

        resultados,

        referencias
    ):

        modelo = (
            resultado["modelo"]
        )


        color = (
            colores[modelo]
        )


        tiempos_v3 = np.linspace(

            0.0,

            resultado["tf"],

            800
        )


        tiempos_v2 = np.linspace(

            0.0,

            referencia["tf"],

            800
        )


        datos_v3 = evaluar(

            resultado,

            tiempos_v3
        )


        datos_v2 = evaluar(

            referencia,

            tiempos_v2
        )


        eje.plot(

            tiempos_v2,

            datos_v2["Tg"]
            -
            273.15,

            color=color,

            linestyle=":",

            linewidth=2,

            label=f"{modelo} — V2 adiabática"
        )


        eje.plot(

            tiempos_v3,

            datos_v3["Tg"]
            -
            273.15,

            color=color,

            linewidth=2.4,

            label=f"{modelo} — V3"
        )


    eje.axhline(

        T_AMBIENTE
        -
        273.15,

        color=MUTED,

        linestyle="--",

        alpha=0.55
    )


    eje.set_title(

        "EFECTO DE LA PARED TÉRMICA: V2 VS. V3",

        fontsize=14,

        weight="bold"
    )


    eje.set_xlabel(
        "Tiempo [s]"
    )


    eje.set_ylabel(
        "Temperatura del gas [°C]"
    )


    eje.grid(

        color=GRID,

        alpha=0.40
    )


    eje.legend(
        frameon=False
    )


    for borde in eje.spines.values():

        borde.set_color(
            "#303B48"
        )


    plt.tight_layout()


    return figura


# ======================================================================
# 22. REPORTE TXT
# ======================================================================

def crear_reporte(
    resultados,
    referencias,
    validaciones,
    error_analitico
):

    lineas = []


    lineas.append(
        "ICF130 — VERSIÓN 3"
    )

    lineas.append(
        "LLENADO TRANSIENTE CON PARED TÉRMICA"
    )

    lineas.append(
        "=" * 78
    )


    lineas.append("")

    lineas.append(
        "1. MODELO FÍSICO"
    )

    lineas.append("")

    lineas.append(
        "dm/dt = mdot(P)"
    )

    lineas.append(
        "dUg/dt = mdot*h_in - Qdot_gw"
    )

    lineas.append(
        "Cw*dTw/dt = Qdot_gw - Qdot_wa"
    )

    lineas.append(
        "Qdot_gw = UA_gw*(Tg-Tw)"
    )

    lineas.append(
        "Qdot_wa = UA_wa*(Tw-Tamb)"
    )


    lineas.append("")

    lineas.append(
        "2. CONDICIONES"
    )

    lineas.append("")


    lineas.append(
        f"Volumen = {V:.6g} m3"
    )

    lineas.append(
        f"P inicial = {P0/1e5:.3f} bar"
    )

    lineas.append(
        f"T inicial = {T0:.3f} K"
    )

    lineas.append(
        f"P suministro = {P_SUPPLY/1e5:.3f} bar"
    )

    lineas.append(
        f"T suministro = {T_SUPPLY:.3f} K"
    )

    lineas.append(
        f"P objetivo = {P_TARGET/1e5:.3f} bar"
    )

    lineas.append(
        f"P crítica = {P_CRITICA/1e5:.6f} bar"
    )

    lineas.append(
        f"Caudal estrangulado de referencia = "
        f"{MDOT_REFERENCIA:.6f} kg/s"
    )

    lineas.append(
        f"Cd = {CD:.6f}"
    )

    lineas.append(
        f"Área efectiva = {AREA_ORIFICIO:.12e} m2"
    )

    lineas.append(
        f"Diámetro efectivo = "
        f"{DIAMETRO_ORIFICIO*1000:.6f} mm"
    )


    lineas.append("")

    lineas.append(
        "3. PARED TÉRMICA"
    )

    lineas.append("")


    lineas.append(
        f"T ambiente = {T_AMBIENTE:.3f} K"
    )

    lineas.append(
        f"T pared inicial = {T_PARED_INICIAL:.3f} K"
    )

    lineas.append(
        f"Masa pared = {MASA_PARED:.6f} kg"
    )

    lineas.append(
        f"cp pared = {CP_PARED:.6f} J/(kg K)"
    )

    lineas.append(
        f"Capacidad térmica pared = "
        f"{CAPACIDAD_PARED:.6f} J/K"
    )

    lineas.append(
        f"UA gas-pared = {UA_GAS_PARED:.6f} W/K"
    )

    lineas.append(
        f"UA pared-ambiente = "
        f"{UA_PARED_AMBIENTE:.6f} W/K"
    )


    lineas.append("")

    lineas.append(
        "Los parámetros térmicos de pared son valores de estudio, "
        "no datos experimentales medidos."
    )


    for resultado, referencia in zip(

        resultados,

        referencias
    ):

        modelo = (
            resultado["modelo"]
        )


        final = evaluar(

            resultado,

            np.array([
                resultado["tf"]
            ])
        )


        final_v2 = evaluar(

            referencia,

            np.array([
                referencia["tf"]
            ])
        )


        balance = (
            validaciones[modelo]
        )


        lineas.append("")

        lineas.append(
            "=" * 78
        )

        lineas.append(
            f"RESULTADOS — {modelo.upper()}"
        )

        lineas.append(
            "=" * 78
        )


        lineas.append(
            f"Tiempo final V3 = "
            f"{resultado['tf']:.9f} s"
        )

        lineas.append(
            f"Tiempo final V2 = "
            f"{referencia['tf']:.9f} s"
        )

        lineas.append(
            f"Temperatura gas final V3 = "
            f"{final['Tg'][0]:.9f} K"
        )

        lineas.append(
            f"Temperatura gas final V3 = "
            f"{final['Tg'][0]-273.15:.9f} C"
        )

        lineas.append(
            f"Temperatura pared final = "
            f"{final['Tw'][0]-273.15:.9f} C"
        )

        lineas.append(
            f"Temperatura gas final V2 = "
            f"{final_v2['Tg'][0]-273.15:.9f} C"
        )

        lineas.append(
            f"Masa final = "
            f"{final['m'][0]:.9f} kg"
        )

        lineas.append(
            f"Presión final = "
            f"{final['P'][0]/1e5:.9f} bar"
        )

        lineas.append(
            f"Factor Z final = "
            f"{final['Z'][0]:.9f}"
        )

        lineas.append(
            f"Caudal de la ley al evento = "
            f"{final['mdot'][0]:.9f} kg/s"
        )

        lineas.append(
            f"Energía de entrada = "
            f"{final['Ein'][0]/1000:.9f} kJ"
        )

        lineas.append(
            f"Calor gas->pared = "
            f"{final['Qgw'][0]/1000:.9f} kJ"
        )

        lineas.append(
            f"Calor pared->ambiente = "
            f"{final['Qwa'][0]/1000:.9f} kJ"
        )


        if (
            resultado["transicion"]
            is not None
        ):

            lineas.append(
                "Tiempo de transición a subsónico = "
                f"{resultado['transicion']:.9f} s"
            )


        lineas.append(
            f"Error relativo balance gas = "
            f"{balance['gas']:.6e}"
        )

        lineas.append(
            f"Error relativo balance pared = "
            f"{balance['pared']:.6e}"
        )

        lineas.append(
            f"Error relativo balance total = "
            f"{balance['total']:.6e}"
        )

        lineas.append(
            f"Residuo energético total máximo = "
            f"{balance['total_J']:.6e} J"
        )


    lineas.append("")

    lineas.append(
        "=" * 78
    )

    lineas.append(
        "VALIDACIÓN ADICIONAL"
    )

    lineas.append(
        "=" * 78
    )


    lineas.append(
        "Error máximo solución analítica ideal adiabática = "
        f"{error_analitico:.6e} K"
    )


    lineas.append("")

    lineas.append(
        "INTERPRETACIÓN DE V3"
    )

    lineas.append("")


    lineas.append(
        "Respecto de V2, parte de la energía que entra con la corriente "
        "de gas ya no queda almacenada únicamente en el gas."
    )


    lineas.append(
        "El gas transfiere energía a la pared, la pared aumenta su "
        "temperatura y simultáneamente intercambia calor con el ambiente."
    )


    lineas.append(
        "Como consecuencia, la temperatura del gas es menor que en "
        "la referencia adiabática."
    )


    lineas.append(
        "A una misma presión objetivo, una menor temperatura implica "
        "que debe acumularse una cantidad diferente de masa y cambia "
        "el tiempo necesario para completar el llenado."
    )


    ruta = (
        DIRECTORIO_SALIDA
        /
        "V3_reporte.txt"
    )


    ruta.write_text(

        "\n".join(
            lineas
        ),

        encoding="utf-8"
    )


    print(
        f"TXT: {ruta}"
    )


# ======================================================================
# 23. EXPORTAR CSV
# ======================================================================

def exportar_csv(
    resultados
):

    tiempo_maximo = max(

        resultado["tf"]

        for resultado
        in resultados
    )


    tiempos = np.linspace(

        0.0,

        tiempo_maximo,

        1500
    )


    columnas = [
        tiempos
    ]


    nombres = [
        "tiempo_global_s"
    ]


    for resultado in resultados:

        datos = evaluar(
            resultado,
            tiempos
        )


        modelo = (
            resultado["modelo"]
            .replace(
                " ",
                "_"
            )
        )


        variables = [

            (
                "tiempo_local_s",
                datos["t_local"]
            ),

            (
                "masa_kg",
                datos["m"]
            ),

            (
                "Tgas_K",
                datos["Tg"]
            ),

            (
                "Tpared_K",
                datos["Tw"]
            ),

            (
                "presion_Pa",
                datos["P"]
            ),

            (
                "mdot_kg_s",
                datos["mdot"]
            ),

            (
                "Z",
                datos["Z"]
            ),

            (
                "Ein_J",
                datos["Ein"]
            ),

            (
                "Q_gas_pared_J",
                datos["Qgw"]
            ),

            (
                "Q_pared_ambiente_J",
                datos["Qwa"]
            ),

            (
                "Qdot_gas_pared_W",
                datos["qdot_gw"]
            ),

            (
                "Qdot_pared_ambiente_W",
                datos["qdot_wa"]
            ),

            (
                "activo",
                datos["activo"].astype(
                    float
                )
            )
        ]


        for nombre, valores in variables:

            nombres.append(
                f"{modelo}_{nombre}"
            )


            columnas.append(
                valores
            )


    matriz = np.column_stack(
        columnas
    )


    ruta = (
        DIRECTORIO_SALIDA
        /
        "V3_datos.csv"
    )


    np.savetxt(

        ruta,

        matriz,

        delimiter=",",

        header=",".join(
            nombres
        ),

        comments="",

        fmt="%.10e"
    )


    print(
        f"CSV: {ruta}"
    )


# ======================================================================
# 24. RESULTADOS EN CONSOLA
# ======================================================================

def imprimir_resultados(
    resultados,
    referencias,
    validaciones,
    error_analitico
):

    print()

    print(
        "=" * 78
    )

    print(
        "ICF130 — VERSIÓN 3"
    )

    print(
        "CAUDAL VARIABLE + PARED TÉRMICA"
    )

    print(
        "=" * 78
    )


    print()

    print(
        "CONDICIONES"
    )

    print(
        "-" * 78
    )


    print(
        f"Volumen tanque          = "
        f"{V*1000:.2f} L"
    )


    print(
        f"Presión inicial         = "
        f"{P0/1e5:.3f} bar"
    )


    print(
        f"Presión suministro      = "
        f"{P_SUPPLY/1e5:.3f} bar"
    )


    print(
        f"Presión objetivo        = "
        f"{P_TARGET/1e5:.3f} bar"
    )


    print(
        f"Presión crítica         = "
        f"{P_CRITICA/1e5:.3f} bar"
    )


    print(
        f"Caudal estrangulado     = "
        f"{MDOT_REFERENCIA:.5f} kg/s"
    )


    print(
        f"Diámetro efectivo       = "
        f"{DIAMETRO_ORIFICIO*1000:.4f} mm"
    )


    print(
        f"UA gas-pared            = "
        f"{UA_GAS_PARED:.3f} W/K"
    )


    print(
        f"UA pared-ambiente       = "
        f"{UA_PARED_AMBIENTE:.3f} W/K"
    )


    print(
        f"Capacidad pared         = "
        f"{CAPACIDAD_PARED:.2f} J/K"
    )


    for resultado, referencia in zip(

        resultados,

        referencias
    ):

        modelo = (
            resultado["modelo"]
        )


        final = evaluar(

            resultado,

            np.array([
                resultado["tf"]
            ])
        )


        final_v2 = evaluar(

            referencia,

            np.array([
                referencia["tf"]
            ])
        )


        print()

        print(
            modelo.upper()
        )

        print(
            "-" * 78
        )


        print(
            f"Tiempo final V3         = "
            f"{resultado['tf']:.3f} s"
        )


        print(
            f"Tiempo final V2         = "
            f"{referencia['tf']:.3f} s"
        )


        print(
            f"T gas final V3          = "
            f"{final['Tg'][0]-273.15:.3f} °C"
        )


        print(
            f"T pared final           = "
            f"{final['Tw'][0]-273.15:.3f} °C"
        )


        print(
            f"T gas final V2          = "
            f"{final_v2['Tg'][0]-273.15:.3f} °C"
        )


        print(
            f"Masa final              = "
            f"{final['m'][0]:.6f} kg"
        )


        print(
            f"Factor Z final          = "
            f"{final['Z'][0]:.6f}"
        )


        print(
            f"Caudal ley al evento    = "
            f"{final['mdot'][0]:.6f} kg/s"
        )


        print(
            f"Q gas -> pared          = "
            f"{final['Qgw'][0]/1000:.3f} kJ"
        )


        print(
            f"Q pared -> ambiente     = "
            f"{final['Qwa'][0]/1000:.3f} kJ"
        )


        if (
            resultado["transicion"]
            is not None
        ):

            print(
                "Transición subsónica     = "
                f"{resultado['transicion']:.3f} s"
            )


        print(
            "Error balance total      = "
            f"{validaciones[modelo]['total']:.3e}"
        )


    print()

    print(
        "VALIDACIÓN"
    )

    print(
        "-" * 78
    )


    print(
        "Error analítico ideal V2 = "
        f"{error_analitico:.3e} K"
    )


    print(
        "=" * 78
    )


# ======================================================================
# 25. MAIN
# ======================================================================

def main():

    validar_parametros()


    print(
        "Integrando casos V3..."
    )


    # ==================================================================
    # V3 — CON PARED TÉRMICA
    # ==================================================================

    resultados = [

        simular(
            "Ideal",
            UA_GAS_PARED
        ),

        simular(
            "VdW",
            UA_GAS_PARED
        )
    ]


    # ==================================================================
    # V2 — REFERENCIA ADIABÁTICA
    #
    # UA_gas-pared = 0
    #
    # La pared queda desacoplada del gas.
    # ==================================================================

    referencias = [

        simular(
            "Ideal",
            0.0
        ),

        simular(
            "VdW",
            0.0
        )
    ]


    # ==================================================================
    # VALIDACIONES
    # ==================================================================

    validaciones = {

        resultado["modelo"]:

        validar_balance(
            resultado
        )

        for resultado
        in resultados
    }


    error_analitico = (
        error_analitico_ideal(
            referencias[0]
        )
    )


    # ==================================================================
    # CONSOLA
    # ==================================================================

    imprimir_resultados(

        resultados,

        referencias,

        validaciones,

        error_analitico
    )


    # ==================================================================
    # TXT
    # ==================================================================

    if SAVE_TXT:

        crear_reporte(

            resultados,

            referencias,

            validaciones,

            error_analitico
        )


    # ==================================================================
    # CSV
    # ==================================================================

    if SAVE_CSV:

        exportar_csv(
            resultados
        )


    # ==================================================================
    # DASHBOARD
    # ==================================================================

    print(
        "Generando dashboard..."
    )


    (
        figura,
        animacion,
        tiempos,
        datos
    ) = crear_dashboard(
        resultados
    )


    # ==================================================================
    # PNG
    # ==================================================================

    if SAVE_PNG:

        ruta_png = (
            DIRECTORIO_SALIDA
            /
            "V3_resultado_final.png"
        )


        figura.savefig(

            ruta_png,

            dpi=220,

            facecolor=BG,

            bbox_inches="tight"
        )


        print(
            f"PNG: {ruta_png}"
        )


    # ==================================================================
    # COMPARACIÓN V2 VS V3
    # ==================================================================

    if SAVE_COMPARACION_V2:

        figura_comparacion = (
            crear_comparacion_v2_v3(

                resultados,

                referencias
            )
        )


        ruta_comparacion = (

            DIRECTORIO_SALIDA

            /

            "V3_comparacion_V2.png"
        )


        figura_comparacion.savefig(

            ruta_comparacion,

            dpi=220,

            facecolor=BG,

            bbox_inches="tight"
        )


        print(
            f"Comparación V2-V3: "
            f"{ruta_comparacion}"
        )


    # ==================================================================
    # MP4
    # ==================================================================

    if SAVE_MP4:

        ffmpeg = shutil.which(
            "ffmpeg"
        )


        if ffmpeg is None:

            try:

                import imageio_ffmpeg


                ffmpeg = (
                    imageio_ffmpeg
                    .get_ffmpeg_exe()
                )

            except Exception:

                ffmpeg = None


        if ffmpeg is not None:

            plt.rcParams[
                "animation.ffmpeg_path"
            ] = ffmpeg


            ruta_mp4 = (

                DIRECTORIO_SALIDA

                /

                "V3_simulacion.mp4"
            )


            try:

                writer = FFMpegWriter(

                    fps=FPS,

                    bitrate=6500,

                    metadata={

                        "title":
                        "ICF130 Tema 7 - Version 3",

                        "artist":
                        "Proyecto Computacional"
                    }
                )


                print(
                    "Exportando MP4..."
                )


                animacion.save(

                    str(
                        ruta_mp4
                    ),

                    writer=writer,

                    dpi=115
                )


                print(
                    f"MP4: {ruta_mp4}"
                )


            except Exception as error:

                print(
                    "No fue posible generar el MP4:"
                )

                print(
                    error
                )


        else:

            print(
                "FFmpeg no disponible. "
                "Se omite el MP4."
            )


    # ==================================================================
    # GIF OPCIONAL
    # ==================================================================

    if SAVE_GIF:

        ruta_gif = (

            DIRECTORIO_SALIDA

            /

            "V3_simulacion.gif"
        )


        try:

            animacion.save(

                str(
                    ruta_gif
                ),

                writer=PillowWriter(
                    fps=18
                ),

                dpi=100
            )


            print(
                f"GIF: {ruta_gif}"
            )


        except Exception as error:

            print(
                "No fue posible generar el GIF:"
            )

            print(
                error
            )


    # ==================================================================
    # FIN
    # ==================================================================

    print()

    print(
        "Archivos generados en:"
    )

    print(
        DIRECTORIO_SALIDA
    )


    plt.show()


    return (

        resultados,

        referencias,

        validaciones
    )


# ======================================================================
# 26. EJECUCIÓN
# ======================================================================

if __name__ == "__main__":

    (
        resultados,
        referencias,
        validaciones
    ) = main()