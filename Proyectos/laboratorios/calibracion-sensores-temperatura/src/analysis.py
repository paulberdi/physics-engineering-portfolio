import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Datos act 2 - PT100
T_pt100 = np.array([
    10, 14, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52,
    53, 53, 58, 60, 62, 64, 65, 65, 67, 66, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 76, 75, 78, 77, 80, 81, 82, 84, 88, 86, 85, 87, 90,
    89, 90, 91, 90, 91, 92, 94, 95, 96, 96, 98, 99, 99, 98, 100, 97
])

R_pt100 = np.array([
    103.9, 104.4, 110.4, 111.9, 112.3, 112.7, 113.4,
    114.3, 114.6, 115, 115.5, 115.8, 116.1, 116.7, 117,
    117.4, 118, 118.8, 119, 119.6, 120.4, 120.6, 121.3,
    121.6, 121.7, 121.8, 122.2, 122.4, 122.8, 123, 123.4,
    123.5, 123.7, 123.9, 124.6, 124.8, 125.1, 125.4, 125.7,
    125.9, 126, 127.1, 128.6, 128.7, 128.9, 129.5, 129.7,
    129.7, 129.9, 130.5, 130.5, 131, 131.1, 131.9, 132.3,
    132.5, 132.8, 133.4, 133.5, 134.3, 134.6, 134.7, 135.1
])

dT_real = 1
dR = 0.1

# Datos act 1 - Gas / columna de mercurio
T_gas = np.array([
    10, 14, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52,
    53, 53, 58, 60, 62, 64, 65, 65, 67, 66, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 76, 75, 78, 77, 80, 81, 82, 84, 88, 86, 85, 87, 90,
    89, 90, 91, 90, 91, 92, 94, 95, 96, 96, 98, 99, 99, 98, 100, 97
])

# Altura del mercurio - lado izquierdo
h_i = np.array([
    27.3, 26.6, 26.3, 26.2, 26.2, 26.1, 26.1,
    26.1, 26, 26, 25.9, 25.9, 25.8, 25.7, 25.7,
    25.7, 25.7, 25.4, 25.4, 25.3, 25.1, 25.1, 25,
    24.9, 24.9, 24.9, 24.8, 24.8, 24.8, 24.7, 24.7,
    24.6, 24.5, 24.6, 24.3, 24.3, 24.2, 24.2, 24.1,
    24.1, 24.1, 23.9, 23.7, 23.6, 23.6, 23.4, 23.4,
    23.4, 23.4, 23.3, 23.3, 23.3, 23.2, 22.9, 22.9,
    22.9, 22.9, 22.8, 22.8, 22.5, 22.5, 22.4, 22.3
])

# Altura del mercurio - lado derecho
h_d = np.array([
    25.1, 25.6, 26.1, 26.1, 26.1, 26.2, 26.2,
    26.3, 26.3, 26.4, 26.4, 26.5, 26.5, 26.7, 26.7,
    26.7, 26.8, 26.9, 26.9, 27, 27.1, 27.3, 27.5,
    27.5, 27.5, 27.5, 27.5, 27.5, 27.6, 27.6, 27.7,
    27.8, 27.8, 27.8, 28, 28.1, 28.2, 28.2, 28.3,
    28.4, 28.3, 28.6, 28.9, 28.9, 28.9, 29, 29,
    29, 29.1, 29.2, 29.2, 29.3, 29.4, 29.6, 29.7,
    29.7, 29.6, 29.7, 29.7, 30, 29.9, 30.1, 30.2
])

dh = 0.05
rho = 13595       # kg/m^3
g = 9.81           # m/s^2
Patm = 101325      # Pa

# Cálculos act 1
delta_h_cm = h_d - h_i
delta_h_m = delta_h_cm / 100.0
P_gas = Patm + rho * g * delta_h_m
error_delta_h = np.sqrt(dh**2 + dh**2)
error_P = rho * g * (error_delta_h / 100)


def calc_err_intercept(x, err_slope):
    n = len(x)
    return err_slope * np.sqrt(np.sum(x**2) / n)


# Ajuste Gráfico 1: PT100
m1, b1, r1, p1, err_m1 = linregress(T_pt100, R_pt100)
err_b1 = calc_err_intercept(T_pt100, err_m1)
R0 = b1
alpha = m1 / b1
sigma_alpha = alpha * np.sqrt((err_m1 / m1)**2 + (err_b1 / R0)**2)

# Ajuste Gráfico 2: Altura vs Temperatura
m2, b2, r2, p2, err_m2 = linregress(T_gas, delta_h_cm)
err_b2 = calc_err_intercept(T_gas, err_m2)

# Ajuste Gráfico 3: Presión vs Temperatura (kPa)
P_gas_kPa = P_gas / 1000
error_P_kPa = error_P / 1000
m3, b3, r3, p3, err_m3 = linregress(T_gas, P_gas_kPa)
err_b3 = calc_err_intercept(T_gas, err_m3)

# Modelo Gas Ideal para Gráfico 3
T_K = T_gas + 273.15
P_ideal_kPa = P_gas_kPa[0] * T_K / T_K[0]

print("==================================================")
print("RESULTADOS DE LOS AJUSTES LINEALES (y = mx + b)")
print("==================================================\n")

print("--- GRÁFICO 1: Resistencia PT100 vs Temperatura ---")
print(f"Ecuación   : R = {m1:.4f}*T + {b1:.4f}")
print(f"Pendiente  : {m1:.4f} ± {err_m1:.4f} Ω/°C")
print(f"Intercepto : {b1:.4f} ± {err_b1:.4f} Ω")
print(f"R²         : {r1**2:.5f}")
print(f"Extras     : α = {alpha:.5f} ± {sigma_alpha:.5f} 1/°C\n")

print("--- GRÁFICO 2: Δh Mercurio vs Temperatura ---")
print(f"Ecuación   : Δh = {m2:.4f}*T {b2:+.4f}")
print(f"Pendiente  : {m2:.4f} ± {err_m2:.4f} cm/°C")
print(f"Intercepto : {b2:.4f} ± {err_b2:.4f} cm")
print(f"R²         : {r2**2:.5f}")
print(f"Error inst.: lectura dh = {dh} cm, Δh propagado = {error_delta_h:.4f} cm\n")

print("--- GRÁFICO 3: Presión Absoluta vs Temperatura ---")
print(f"Ecuación   : P = {m3:.4f}*T + {b3:.4f}")
print(f"Pendiente  : {m3:.4f} ± {err_m3:.4f} kPa/°C")
print(f"Intercepto : {b3:.4f} ± {err_b3:.4f} kPa")
print(f"R²         : {r3**2:.5f}")
print(f"Error inst.: P propagado = {error_P:.2f} Pa ({error_P_kPa:.4f} kPa)\n")

props_texto = dict(boxstyle="round", facecolor="white", alpha=0.8)

# Gráfico 1: Resistencia vs Temperatura (PT100)
plt.figure(figsize=(6, 4))
plt.errorbar(
    T_pt100,
    R_pt100,
    yerr=dR,
    xerr=dT_real,
    fmt="o",
    markersize=3,
    capsize=4,
    elinewidth=1.5,
    label="Datos experimentales",
)
plt.plot(T_pt100, b1 + m1 * T_pt100, "r-", label="Ajuste lineal")
texto_eq1 = f"$R = {m1:.4f}T + {b1:.2f}$\n$R^2 = {r1**2:.4f}$"
plt.text(
    0.05,
    0.95,
    texto_eq1,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    bbox=props_texto,
)
plt.xlabel("Temperatura ± 1 (°C)")
plt.ylabel("Resistencia PT100 ± 0.1 (Ω)")
plt.legend(loc="lower right")
plt.grid(True)
plt.tight_layout()
plt.savefig("resistencia_vs_temperatura.png", dpi=300)
plt.show()

# Gráfico 2: Altura mercurio vs Temperatura
plt.figure(figsize=(6, 4))
plt.errorbar(
    T_gas,
    delta_h_cm,
    yerr=error_delta_h,
    xerr=dT_real,
    fmt="o",
    markersize=3,
    capsize=4,
    elinewidth=1.5,
    label="Datos experimentales",
)
plt.plot(T_gas, b2 + m2 * T_gas, "r-", label="Ajuste lineal")
texto_eq2 = f"$\\Delta h = {m2:.4f}T {b2:+.2f}$\n$R^2 = {r2**2:.4f}$"
plt.text(
    0.05,
    0.95,
    texto_eq2,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    bbox=props_texto,
)
plt.xlabel("Temperatura ± 1 (°C)")
plt.ylabel(r"$\Delta h$ mercurio ± 0.07 (cm)")
plt.legend(loc="lower right")
plt.grid(True)
plt.tight_layout()
plt.savefig("altura_mercurio_vs_temperatura.png", dpi=300)
plt.show()

# Gráfico 3: Presión del bulbo vs Temperatura y gas ideal
plt.figure(figsize=(6, 4))
plt.errorbar(
    T_gas,
    P_gas_kPa,
    yerr=error_P_kPa,
    xerr=dT_real,
    fmt="o",
    markersize=3,
    capsize=4,
    elinewidth=1.5,
    label="Presión experimental",
)
plt.plot(T_gas, b3 + m3 * T_gas, "r-", label="Ajuste lineal exp.")
plt.plot(T_gas, P_ideal_kPa, "g--", label="Modelo Gas Ideal")
texto_eq3 = f"$P = {m3:.4f}T + {b3:.2f}$\n$R^2 = {r3**2:.4f}$"
plt.text(
    0.05,
    0.95,
    texto_eq3,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    bbox=props_texto,
)
plt.xlabel("Temperatura ± 1 (°C)")
plt.ylabel("Presión Absoluta ± 0.09 (kPa)")
plt.legend(loc="lower right")
plt.grid(True)
plt.tight_layout()
plt.savefig("presion_vs_temperatura.png", dpi=300)
plt.show()
