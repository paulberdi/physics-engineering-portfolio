from machine import ADC, Pin
import utime
import math

adc = ADC(26)

# Relés asociados a cada filtro
rele_graves = Pin(2, Pin.OUT)
rele_medios = Pin(3, Pin.OUT)
rele_agudos = Pin(4, Pin.OUT)

FS = 8000
N = 256
BIAS = 32768

FREQ_GRAVE = 120
FREQ_MEDIO = 1000
FREQ_AGUDO = 3000

ultimo_filtro = "medio"

# Cambia esto a True si tu módulo de relé se activa con LOW
RELE_ACTIVO_EN_LOW = False


def rele_on(rele):
    rele.value(0 if RELE_ACTIVO_EN_LOW else 1)


def rele_off(rele):
    rele.value(1 if RELE_ACTIVO_EN_LOW else 0)


def seleccionar_filtro(nombre):
    global ultimo_filtro

    # Apaga todos primero
    rele_off(rele_graves)
    rele_off(rele_medios)
    rele_off(rele_agudos)

    # Pequeña pausa para evitar dos caminos conectados a la vez
    utime.sleep_ms(20)

    if nombre == "grave":
        rele_on(rele_graves)
    elif nombre == "medio":
        rele_on(rele_medios)
    elif nombre == "agudo":
        rele_on(rele_agudos)

    ultimo_filtro = nombre


def goertzel(samples, freq, fs):
    k = int(0.5 + ((N * freq) / fs))
    omega = (2.0 * math.pi * k) / N
    coeff = 2.0 * math.cos(omega)

    q0 = 0
    q1 = 0
    q2 = 0

    for sample in samples:
        q0 = coeff * q1 - q2 + sample
        q2 = q1
        q1 = q0

    return q1 * q1 + q2 * q2 - q1 * q2 * coeff


def leer_bloque():
    samples = []
    periodo_us = int(1000000 / FS)

    for _ in range(N):
        valor = adc.read_u16()
        valor = valor - BIAS
        samples.append(valor)
        utime.sleep_us(periodo_us)

    return samples


def decidir_filtro(grave, medio, agudo):
    energias = {
        "grave": grave,
        "medio": medio,
        "agudo": agudo
    }

    candidato = max(energias, key=energias.get)

    energia_actual = energias[ultimo_filtro]
    energia_candidato = energias[candidato]

    # Evita cambios nerviosos entre relés
    if energia_candidato > energia_actual * 1.50:
        return candidato

    return ultimo_filtro


seleccionar_filtro("medio")

while True:
    muestras = leer_bloque()

    energia_grave = goertzel(muestras, FREQ_GRAVE, FS)
    energia_medio = goertzel(muestras, FREQ_MEDIO, FS)
    energia_agudo = goertzel(muestras, FREQ_AGUDO, FS)

    filtro = decidir_filtro(
        energia_grave,
        energia_medio,
        energia_agudo
    )

    seleccionar_filtro(filtro)

    print("G:", int(energia_grave),
          "M:", int(energia_medio),
          "A:", int(energia_agudo),
          "RELE:", filtro)