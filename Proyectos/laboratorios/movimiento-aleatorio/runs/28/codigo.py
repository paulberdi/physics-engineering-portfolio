import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')  # Para VS Code
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
import glob
import os
import warnings
from datetime import datetime
warnings.filterwarnings('ignore')

# Configuración para gráficos científicos según estándares académicos
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.titlesize': 16,
    'font.family': 'serif',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.dpi': 100
})

class AnalizadorMovimientoAleatorio:
    """
    Analizador de movimiento aleatorio siguiendo la guía FIS 200 UTFSM.
    Implementa la metodología específica del laboratorio universitario.
    """
    
    def __init__(self):
        self.tiempo = None
        self.posicion_x = None
        self.posicion_y = None
        self.archivo_origen = None
        self.dt = None
        self.n_puntos = None
        self.resultados = {}
        
        # Parámetros específicos de la guía UTFSM
        self.dimension = 2  # Análisis en 2D según la guía
        self.regimenes = {
            'subdifusivo': (0, 1),
            'difusivo_clasico': (0.9, 1.1),  # Tolerancia para α ≈ 1
            'subbalítico': (1, 2),
            'balístico': (1.9, 2.1)  # Tolerancia para α ≈ 2
        }
    
    def cargar_datos_tracker(self):
        """Busca y carga datos de Tracker SOLO en la carpeta actual (no en subcarpetas)."""
        print("🔍 Buscando archivos de datos del experimento en la carpeta actual...")
        
        # Obtener la carpeta donde está el script
        carpeta_actual = os.path.dirname(os.path.abspath(__file__))
        print(f"📁 Buscando en: {carpeta_actual}")
        
        # Buscar archivos SOLO en la carpeta actual (no recursivo)
        archivos_datos = []
        extensiones = ['*.xlsx', '*.xls', '*.csv']
        
        for extension in extensiones:
            # Usar os.path.join para buscar solo en carpeta actual
            patron = os.path.join(carpeta_actual, extension)
            archivos_encontrados = glob.glob(patron)
            archivos_datos.extend(archivos_encontrados)
        
        if not archivos_datos:
            print("❌ No se encontraron archivos de datos en esta carpeta")
            print("💡 Coloque un archivo Excel/CSV con datos de Tracker en la MISMA carpeta que este script")
            print("   Formato esperado: t | x | y (tiempo en segundos, posiciones en metros)")
            print(f"   Carpeta actual: {carpeta_actual}")
            return False
        
        # Mostrar solo los nombres de archivo (sin ruta completa)
        nombres_archivos = [os.path.basename(archivo) for archivo in archivos_datos]
        print(f"📄 Archivos encontrados en la carpeta actual: {nombres_archivos}")
        
        # Usar el primer archivo encontrado
        archivo = archivos_datos[0]
        nombre_archivo = os.path.basename(archivo)
        print(f"📊 Procesando archivo: {nombre_archivo}")
        
        try:
            # Leer archivo (Excel o CSV)
            if archivo.endswith('.csv'):
                df = pd.read_csv(archivo)
            else:
                df = pd.read_excel(archivo)
            
            print(f"✅ Archivo cargado: {df.shape[0]} filas, {df.shape[1]} columnas")
            print(f"📋 Columnas disponibles: {list(df.columns)}")
            
            # Identificar columnas automáticamente
            columnas = df.columns.tolist()
            
            # Buscar columnas de tiempo, x, y (común en Tracker)
            col_tiempo = None
            col_x = None
            col_y = None
            
            for col in columnas:
                col_lower = str(col).lower()
                if any(x in col_lower for x in ['t', 'time', 'tiempo']) and col_tiempo is None:
                    col_tiempo = col
                elif any(x in col_lower for x in ['x', 'pos x', 'position x']) and col_x is None:
                    col_x = col
                elif any(x in col_lower for x in ['y', 'pos y', 'position y']) and col_y is None:
                    col_y = col
            
            # Si no encontró por nombre, usar las primeras 3 columnas
            if col_tiempo is None:
                col_tiempo = columnas[0]
            if col_x is None:
                col_x = columnas[1] if len(columnas) > 1 else None
            if col_y is None:
                col_y = columnas[2] if len(columnas) > 2 else None
            
            if col_x is None or col_y is None:
                print("❌ No se pudieron identificar columnas X e Y")
                return False
            
            print(f"📊 Columnas identificadas:")
            print(f"   • Tiempo: '{col_tiempo}'")
            print(f"   • Posición X: '{col_x}'")
            print(f"   • Posición Y: '{col_y}'")
            
            # Extraer y limpiar datos
            self.tiempo = df[col_tiempo].dropna().values
            self.posicion_x = df[col_x].dropna().values
            self.posicion_y = df[col_y].dropna().values
            
            # Asegurar mismo tamaño
            min_length = min(len(self.tiempo), len(self.posicion_x), len(self.posicion_y))
            self.tiempo = self.tiempo[:min_length]
            self.posicion_x = self.posicion_x[:min_length]
            self.posicion_y = self.posicion_y[:min_length]
            
            # Guardar solo el nombre del archivo (sin ruta)
            self.archivo_origen = nombre_archivo
            self.n_puntos = len(self.tiempo)
            
            # Calcular dt promedio
            if self.n_puntos > 1:
                self.dt = np.mean(np.diff(self.tiempo))
            else:
                self.dt = 1.0
            
            print(f"📈 Datos procesados: {self.n_puntos} puntos válidos")
            print(f"⏱️  Δt promedio: {self.dt:.4f} segundos")
            print(f"🕐 Duración total: {self.tiempo[-1] - self.tiempo[0]:.1f} segundos")
            print(f"📏 Rango X: [{np.min(self.posicion_x):.4f}, {np.max(self.posicion_x):.4f}] m")
            print(f"📏 Rango Y: [{np.min(self.posicion_y):.4f}, {np.max(self.posicion_y):.4f}] m")
            
            return True
            
        except Exception as e:
            print(f"❌ Error procesando {nombre_archivo}: {e}")
            return False
    
    def calcular_msd_segun_guia(self):
        """
        Calcula el MSD siguiendo la ecuación de la guía UTFSM:
        ⟨r²(t)⟩ = ⟨[x(t) − x(0)]² + [y(t) − y(0)]²⟩
        """
        print("\n🧮 Calculando Desplazamiento Cuadrático Medio (MSD)...")
        
        # Posición inicial
        x0 = self.posicion_x[0]
        y0 = self.posicion_y[0]
        
        # Calcular desplazamientos desde el origen
        desplazamiento_x = self.posicion_x - x0
        desplazamiento_y = self.posicion_y - y0
        
        # MSD según ecuación (1) de la guía
        msd_instantaneo = desplazamiento_x**2 + desplazamiento_y**2
        
        # Para análisis, también calculamos MSD para diferentes lags
        max_lag = min(self.n_puntos // 4, 1000)  # Máximo 25% de los datos
        lags = np.arange(1, max_lag + 1)
        msd_promedio = []
        
        for lag in lags:
            if lag < self.n_puntos:
                # Calcular MSD para este lag específico
                dx_lag = self.posicion_x[lag:] - self.posicion_x[:-lag]
                dy_lag = self.posicion_y[lag:] - self.posicion_y[:-lag]
                msd_lag = np.mean(dx_lag**2 + dy_lag**2)
                msd_promedio.append(msd_lag)
        
        # Tiempo correspondiente a cada lag
        tiempo_lags = lags * self.dt
        
        msd_data = {
            'tiempo_total': self.tiempo - self.tiempo[0],  # Tiempo desde t=0
            'msd_instantaneo': msd_instantaneo,
            'lags': lags,
            'tiempo_lags': tiempo_lags,
            'msd_promedio': np.array(msd_promedio),
            'posicion_inicial': (x0, y0),
            'desplazamiento_x': desplazamiento_x,
            'desplazamiento_y': desplazamiento_y
        }
        
        self.resultados['msd'] = msd_data
        print(f"✅ MSD calculado para {len(lags)} valores de lag")
        print(f"   Posición inicial: ({x0:.4f}, {y0:.4f}) m")
        print(f"   MSD final: {msd_instantaneo[-1]:.6f} m²")
        
        return msd_data
    
    def analisis_ley_escala(self):
        """
        Análisis de la ley de escala ⟨r²(t)⟩ = A·t^α mediante regresión log-log.
        Sigue el procedimiento de los pasos 12-14 de la guía.
        """
        print("\n📊 Realizando análisis de ley de escala (log-log)...")
        
        msd_data = self.resultados['msd']
        
        # Usar solo valores positivos para el análisis log-log
        mask_positivos = (msd_data['tiempo_lags'] > 0) & (msd_data['msd_promedio'] > 0)
        tiempo_log = msd_data['tiempo_lags'][mask_positivos]
        msd_log = msd_data['msd_promedio'][mask_positivos]
        
        if len(tiempo_log) < 3:
            print("❌ Insuficientes puntos para análisis log-log")
            return None
        
        # Regresión lineal en escala logarítmica: log(MSD) = log(A) + α·log(t)
        log_tiempo = np.log(tiempo_log)
        log_msd = np.log(msd_log)
        
        # Ajuste lineal
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_tiempo, log_msd)
        
        # Parámetros de la ley de escala
        alpha = slope  # Exponente de escalamiento
        A = np.exp(intercept)  # Constante pre-exponencial
        r_squared = r_value**2
        
        # Clasificar régimen según la guía
        regimen = self.clasificar_regimen(alpha)
        
        # Calcular coeficiente de difusión si es régimen difusivo
        coef_difusion = None
        if 0.9 <= alpha <= 1.1:  # Régimen difusivo (α ≈ 1)
            # Según la guía: ⟨r²(t)⟩ = 2dDt, donde d=2 para 2D
            # Por lo tanto: A = 2dD = 4D → D = A/4
            coef_difusion = A / 4
            print(f"🌊 Régimen difusivo detectado - Calculando coeficiente de difusión")
        
        escala_data = {
            'alpha': alpha,
            'A': A,
            'r_squared': r_squared,
            'p_value': p_value,
            'std_error': std_err,
            'regimen': regimen,
            'coeficiente_difusion': coef_difusion,
            'log_tiempo': log_tiempo,
            'log_msd': log_msd,
            'tiempo_fit': tiempo_log,
            'msd_fit': msd_log
        }
        
        self.resultados['ley_escala'] = escala_data
        
        print(f"✅ Análisis de ley de escala completado:")
        print(f"   • Exponente α = {alpha:.4f} ± {std_err:.4f}")
        print(f"   • Constante A = {A:.6f}")
        print(f"   • R² = {r_squared:.4f}")
        print(f"   • Régimen identificado: {regimen}")
        if coef_difusion:
            print(f"   • Coeficiente de difusión D = {coef_difusion:.6f} m²/s")
        
        return escala_data
    
    def clasificar_regimen(self, alpha):
        """Clasifica el régimen según el exponente α y los criterios de la guía."""
        if 0 < alpha < 1:
            return "Subdifusivo anómalo"
        elif 0.9 <= alpha <= 1.1:
            return "Difusivo (clásico browniano)"
        elif 1 < alpha < 2:
            return "Superdifusivo (subbalístico)"
        elif 1.9 <= alpha <= 2.1:
            return "Balístico"
        elif alpha > 2:
            return "Superdifusivo extremo"
        else:
            return "Régimen no clasificado"
    
    def analisis_estadistico_basico(self):
        """Análisis estadístico básico de las trayectorias."""
        print("\n📈 Realizando análisis estadístico básico...")
        
        # Estadísticas de posición
        stats_posicion = {
            'x': {
                'media': np.mean(self.posicion_x),
                'std': np.std(self.posicion_x),
                'rango': np.max(self.posicion_x) - np.min(self.posicion_x),
                'min': np.min(self.posicion_x),
                'max': np.max(self.posicion_x)
            },
            'y': {
                'media': np.mean(self.posicion_y),
                'std': np.std(self.posicion_y),
                'rango': np.max(self.posicion_y) - np.min(self.posicion_y),
                'min': np.min(self.posicion_y),
                'max': np.max(self.posicion_y)
            }
        }
        
        # Estadísticas de desplazamiento
        msd_data = self.resultados['msd']
        stats_desplazamiento = {
            'desplazamiento_total': np.sqrt(msd_data['msd_instantaneo'][-1]),
            'desplazamiento_maximo': np.sqrt(np.max(msd_data['msd_instantaneo'])),
            'velocidad_media_x': np.mean(np.diff(self.posicion_x) / self.dt),
            'velocidad_media_y': np.mean(np.diff(self.posicion_y) / self.dt),
            'distancia_recorrida': np.sum(np.sqrt(np.diff(self.posicion_x)**2 + np.diff(self.posicion_y)**2))
        }
        
        self.resultados['estadisticas'] = {
            'posicion': stats_posicion,
            'desplazamiento': stats_desplazamiento
        }
        
        return self.resultados['estadisticas']
    
    def crear_graficos_segun_guia(self):
        """Crea gráficos siguiendo la metodología específica de la guía UTFSM."""
        print("\n📊 Creando visualizaciones según metodología FIS 200...")
        
        # Configurar figura con 6 subgráficos
        fig = plt.figure(figsize=(18, 12))
        
        # Título principal
        nombre_archivo = os.path.splitext(self.archivo_origen)[0] if self.archivo_origen else "Datos"
        fig.suptitle(f'Análisis de Movimiento Aleatorio - Experimento {nombre_archivo}\n' +
                    f'Universidad Técnica Federico Santa María - FIS 200', 
                    fontsize=16, fontweight='bold', y=0.95)
        
        msd_data = self.resultados['msd']
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        # 1. Trayectoria 2D (como en Tracker)
        ax1 = plt.subplot(2, 3, 1)
        ax1.plot(self.posicion_x, self.posicion_y, 'b-', alpha=0.7, linewidth=1.0)
        ax1.plot(self.posicion_x[0], self.posicion_y[0], 'go', markersize=8, label='Inicio')
        ax1.plot(self.posicion_x[-1], self.posicion_y[-1], 'ro', markersize=8, label='Final')
        ax1.set_xlabel('Posición X [m]')
        ax1.set_ylabel('Posición Y [m]')
        ax1.set_title('Trayectoria de la Esfera', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.axis('equal')
        
        # Información del experimento
        info_exp = f'N = {self.n_puntos} puntos\nΔt = {self.dt:.4f} s\nT = {self.tiempo[-1]-self.tiempo[0]:.1f} s'
        ax1.text(0.02, 0.98, info_exp, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        # 2. Posición X vs Tiempo
        ax2 = plt.subplot(2, 3, 2)
        ax2.plot(self.tiempo, self.posicion_x, 'b-', alpha=0.8, linewidth=0.8)
        ax2.set_xlabel('Tiempo [s]')
        ax2.set_ylabel('Posición X [m]')
        ax2.set_title('Posición X vs Tiempo', fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Estadísticas
        stats_x = stats_data['posicion']['x']
        info_x = f'⟨x⟩ = {stats_x["media"]:.4f} m\nσₓ = {stats_x["std"]:.4f} m\nRango = {stats_x["rango"]:.4f} m'
        ax2.text(0.02, 0.98, info_x, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        # 3. Posición Y vs Tiempo
        ax3 = plt.subplot(2, 3, 3)
        ax3.plot(self.tiempo, self.posicion_y, 'r-', alpha=0.8, linewidth=0.8)
        ax3.set_xlabel('Tiempo [s]')
        ax3.set_ylabel('Posición Y [m]')
        ax3.set_title('Posición Y vs Tiempo', fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # Estadísticas
        stats_y = stats_data['posicion']['y']
        info_y = f'⟨y⟩ = {stats_y["media"]:.4f} m\nσᵧ = {stats_y["std"]:.4f} m\nRango = {stats_y["rango"]:.4f} m'
        ax3.text(0.02, 0.98, info_y, transform=ax3.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        
        # 4. MSD vs Tiempo (escala lineal) - Paso 12 de la guía
        ax4 = plt.subplot(2, 3, 4)
        ax4.plot(msd_data['tiempo_lags'], msd_data['msd_promedio'], 'ko-', markersize=3, linewidth=1.5)
        
        # Ajuste teórico si es difusivo
        if escala_data['coeficiente_difusion']:
            D = escala_data['coeficiente_difusion']
            msd_teorico = 4 * D * msd_data['tiempo_lags']  # ⟨r²(t)⟩ = 2dDt con d=2
            ax4.plot(msd_data['tiempo_lags'], msd_teorico, 'r--', linewidth=2, 
                    label=f'Teórico: 4Dt\nD = {D:.6f} m²/s')
            ax4.legend()
        
        ax4.set_xlabel('Tiempo [s]')
        ax4.set_ylabel('⟨r²(t)⟩ [m²]')
        ax4.set_title('Desplazamiento Cuadrático Medio', fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        # Ecuación MSD
        ecuacion_msd = r'$\langle r^2(t) \rangle = \langle [x(t)-x(0)]^2 + [y(t)-y(0)]^2 \rangle$'
        ax4.text(0.02, 0.98, ecuacion_msd, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # 5. Análisis Log-Log (Pasos 13-14 de la guía) - GRÁFICO PRINCIPAL
        ax5 = plt.subplot(2, 3, 5)
        ax5.loglog(escala_data['tiempo_fit'], escala_data['msd_fit'], 'bo', markersize=4, alpha=0.7, label='Datos experimentales')
        
        # Línea de ajuste
        msd_fit_teorico = escala_data['A'] * escala_data['tiempo_fit']**escala_data['alpha']
        ax5.loglog(escala_data['tiempo_fit'], msd_fit_teorico, 'r-', linewidth=2, 
                  label=f'Ajuste: At^α\nα = {escala_data["alpha"]:.3f}')
        
        # Líneas de referencia para diferentes regímenes
        t_ref = escala_data['tiempo_fit']
        msd_ref_base = escala_data['msd_fit'][len(escala_data['msd_fit'])//2]
        t_ref_base = escala_data['tiempo_fit'][len(escala_data['tiempo_fit'])//2]
        
        # Referencias teóricas
        if not (0.9 <= escala_data['alpha'] <= 1.1):  # Si no es difusivo, mostrar referencia
            msd_difusivo = msd_ref_base * (t_ref / t_ref_base)**1.0
            ax5.loglog(t_ref, msd_difusivo, 'k--', alpha=0.5, label='α = 1 (difusivo)')
        
        ax5.set_xlabel('Tiempo [s]')
        ax5.set_ylabel('⟨r²(t)⟩ [m²]')
        ax5.set_title('Análisis Log-Log: ⟨r²(t)⟩ ∝ t^α', fontweight='bold')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        
        # Información del ajuste
        info_ajuste = (f'⟨r²(t)⟩ = A·t^α\n'
                      f'A = {escala_data["A"]:.2e} m²/s^α\n'
                      f'α = {escala_data["alpha"]:.4f} ± {escala_data["std_error"]:.4f}\n'
                      f'R² = {escala_data["r_squared"]:.4f}\n'
                      f'Régimen: {escala_data["regimen"]}')
        ax5.text(0.02, 0.98, info_ajuste, transform=ax5.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
        
        # 6. Interpretación de Regímenes
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')  # Sin ejes, solo texto
        
        # Tabla de regímenes
        regimenes_texto = """CLASIFICACIÓN DE REGÍMENES (Según Guía FIS 200):

• Subdifusivo anómalo: 0 < α < 1
  Movimiento restringido por barreras físicas
  
• Difusivo clásico: α = 1
  Movimiento browniano ideal
  
• Superdifusivo subbalístico: 1 < α < 2  
  Movimiento más persistente que difusivo
  
• Balístico: α = 2
  Movimiento con dirección definida

RESULTADO EXPERIMENTAL:"""
        
        resultado_texto = (f'\nα = {escala_data["alpha"]:.4f}\n'
                          f'Régimen: {escala_data["regimen"]}\n\n')
        
        if escala_data['coeficiente_difusion']:
            resultado_texto += f'Coeficiente de difusión:\nD = {escala_data["coeficiente_difusion"]:.2e} m²/s\n\n'
        
        # Evaluación crítica
        if 0.9 <= escala_data['alpha'] <= 1.1:
            evaluacion = "✅ Comportamiento consistente con\nmovimiento browniano clásico"
        else:
            evaluacion = f"⚠️ Comportamiento anómalo\n(α ≠ 1)"
        
        ax6.text(0.05, 0.95, regimenes_texto + resultado_texto + evaluacion,
                transform=ax6.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightsteelblue', alpha=0.8),
                family='monospace')
        
        plt.tight_layout()
        
        # Guardar figura
        nombre_figura = f"analisis_movimiento_aleatorio_{os.path.splitext(self.archivo_origen)[0]}.png"
        plt.savefig(nombre_figura, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"📊 Gráficos guardados en: {nombre_figura}")
        
        # Mostrar gráficos
        try:
            plt.show(block=True)
            print("✅ Gráficos mostrados correctamente")
        except Exception as e:
            print(f"⚠️ Error mostrando gráficos: {e}")
            print("💡 Los gráficos se guardaron correctamente en el archivo PNG")
        
        return fig
    
    def generar_reporte_paper(self):
        """Genera un reporte estilo paper científico según solicita la guía."""
        print("\n📄 Generando paper científico...")
        
        fecha_actual = datetime.now().strftime("%Y-%m-%d")
        nombre_archivo = os.path.splitext(self.archivo_origen)[0] if self.archivo_origen else "experimento"
        nombre_reporte = f"paper_movimiento_aleatorio_{nombre_archivo}.txt"
        
        msd_data = self.resultados['msd']
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        with open(nombre_reporte, 'w', encoding='utf-8') as f:
            # Encabezado estilo paper
            f.write("="*80 + "\n")
            f.write("ANÁLISIS EXPERIMENTAL DE MOVIMIENTO ALEATORIO\n")
            f.write("Universidad Técnica Federico Santa María\n")
            f.write("Física Experimental - FIS 200\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"Fecha: {fecha_actual}\n")
            f.write(f"Archivo de datos: {self.archivo_origen}\n")
            f.write(f"Duración del experimento: {self.tiempo[-1]-self.tiempo[0]:.1f} segundos\n")
            f.write(f"Número de mediciones: {self.n_puntos}\n\n")
            
            # RESUMEN
            f.write("RESUMEN\n")
            f.write("-"*40 + "\n")
            f.write("Se analizó el movimiento aleatorio de una esfera en un recipiente vibrante\n")
            f.write("utilizando análisis de video mediante software Tracker. Se determinó el\n")
            f.write(f"régimen dinámico como '{escala_data['regimen']}' con exponente de\n")
            f.write(f"escalamiento α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}.\n")
            
            if escala_data['coeficiente_difusion']:
                f.write(f"El coeficiente de difusión estimado es D = {escala_data['coeficiente_difusion']:.2e} m²/s.\n")
            f.write("\n")
            
            # INTRODUCCIÓN
            f.write("1. INTRODUCCIÓN\n")
            f.write("-"*40 + "\n")
            f.write("El movimiento browniano describe el desplazamiento aleatorio de partículas\n")
            f.write("microscópicas suspendidas en un fluido. En este experimento se simula este\n")
            f.write("comportamiento mediante esferas en un recipiente vibrante.\n\n")
            
            f.write("El desplazamiento cuadrático medio (MSD) se define como:\n")
            f.write("⟨r²(t)⟩ = ⟨[x(t) − x(0)]² + [y(t) − y(0)]²⟩\n\n")
            
            f.write("La evolución temporal sigue la ley de escala:\n")
            f.write("⟨r²(t)⟩ = A·t^α\n\n")
            
            f.write("donde α caracteriza el régimen dinámico:\n")
            f.write("• α < 1: Subdifusivo anómalo\n")
            f.write("• α = 1: Difusivo clásico (browniano)\n")
            f.write("• α > 1: Superdifusivo\n\n")
            
            # METODOLOGÍA
            f.write("2. METODOLOGÍA EXPERIMENTAL\n")
            f.write("-"*40 + "\n")
            f.write("• Se utilizó una esfera en recipiente vibrante\n")
            f.write("• Grabación de video durante 120 segundos\n")
            f.write("• Análisis de trayectoria con software Tracker\n")
            f.write("• Calibración dimensional del sistema\n")
            f.write(f"• Frecuencia de muestreo: {1/self.dt:.1f} Hz\n\n")
            
            # RESULTADOS
            f.write("3. RESULTADOS\n")
            f.write("-"*40 + "\n")
            
            f.write("3.1 Características del movimiento:\n")
            stats_desp = stats_data['desplazamiento']
            f.write(f"• Desplazamiento total: {stats_desp['desplazamiento_total']:.4f} m\n")
            f.write(f"• Desplazamiento máximo: {stats_desp['desplazamiento_maximo']:.4f} m\n")
            f.write(f"• Distancia recorrida: {stats_desp['distancia_recorrida']:.4f} m\n")
            f.write(f"• Velocidad media X: {stats_desp['velocidad_media_x']:.6f} m/s\n")
            f.write(f"• Velocidad media Y: {stats_desp['velocidad_media_y']:.6f} m/s\n\n")
            
            f.write("3.2 Análisis de ley de escala:\n")
            f.write(f"• Exponente α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}\n")
            f.write(f"• Constante A = {escala_data['A']:.2e} m²/s^α\n")
            f.write(f"• Coeficiente de correlación R² = {escala_data['r_squared']:.4f}\n")
            f.write(f"• p-valor = {escala_data['p_value']:.2e}\n")
            f.write(f"• Régimen identificado: {escala_data['regimen']}\n\n")
            
            if escala_data['coeficiente_difusion']:
                f.write("3.3 Coeficiente de difusión:\n")
                f.write("Para régimen difusivo (α ≈ 1), se aplica:\n")
                f.write("⟨r²(t)⟩ = 2dDt, donde d = 2 (2D)\n")
                f.write(f"D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
            
            # DISCUSIÓN
            f.write("4. DISCUSIÓN\n")
            f.write("-"*40 + "\n")
            
            # Evaluación del comportamiento
            if 0.9 <= escala_data['alpha'] <= 1.1:
                f.write("4.1 Comportamiento difusivo:\n")
                f.write("El exponente α ≈ 1 confirma un comportamiento difusivo clásico,\n")
                f.write("consistente con movimiento browniano. Esto indica que:\n")
                f.write("• Las colisiones son aleatorias e independientes\n")
                f.write("• No hay efectos de memoria significativos\n")
                f.write("• El sistema ha alcanzado equilibrio estadístico\n\n")
            else:
                f.write("4.1 Comportamiento anómalo:\n")
                f.write(f"El exponente α = {escala_data['alpha']:.4f} indica desviación del\n")
                f.write("comportamiento difusivo ideal. Posibles causas:\n")
                if escala_data['alpha'] < 1:
                    f.write("• Restricciones geométricas del recipiente\n")
                    f.write("• Efectos de fricción o atrapamiento\n")
                    f.write("• Heterogeneidad en la vibración\n")
                else:
                    f.write("• Movimiento balístico a cortos tiempos\n")
                    f.write("• Efectos de inercia de la esfera\n")
                    f.write("• Correlaciones en la vibración\n")
                f.write("\n")
            
            f.write("4.2 Incertidumbres experimentales:\n")
            f.write("• Resolución espacial del video\n")
            f.write("• Precisión del software de tracking\n")
            f.write("• Efectos de bordes del recipiente\n")
            f.write("• Variaciones en la amplitud de vibración\n\n")
            
            # CONCLUSIONES
            f.write("5. CONCLUSIONES\n")
            f.write("-"*40 + "\n")
            f.write(f"1. Se identificó un régimen {escala_data['regimen'].lower()}\n")
            f.write(f"   con exponente de escalamiento α = {escala_data['alpha']:.4f}\n\n")
            
            if escala_data['coeficiente_difusion']:
                f.write(f"2. El coeficiente de difusión efectivo es D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
            
            f.write("3. El análisis log-log permite discriminar efectivamente\n")
            f.write("   entre diferentes regímenes de transporte\n\n")
            
            if 0.9 <= escala_data['alpha'] <= 1.1:
                f.write("4. Los resultados son consistentes con la teoría\n")
                f.write("   del movimiento browniano clásico\n")
            else:
                f.write("4. Se observan desviaciones del comportamiento\n")
                f.write("   browniano ideal que requieren mayor investigación\n")
            
            f.write("\n")
            
            # REFERENCIAS
            f.write("REFERENCIAS\n")
            f.write("-"*40 + "\n")
            f.write("[1] Einstein, A. (1905). Über die von der molekularkinetischen Theorie\n")
            f.write("    der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten\n")
            f.write("    suspendierten Teilchen. Ann. Phys. 17, 549-560.\n\n")
            f.write("[2] Metzler, R., & Klafter, J. (2000). The random walk's guide to\n")
            f.write("    anomalous diffusion. Physics Reports, 339(1), 1-77.\n\n")
            f.write("[3] Guía de Laboratorio FIS 200 - Universidad Técnica Federico Santa María\n")
            
            f.write("\n" + "="*80 + "\n")
            
    def generar_reporte_paper(self):
        """Genera un reporte estilo paper científico según solicita la guía."""
        print("\n📄 Generando paper científico...")
        
        fecha_actual = datetime.now().strftime("%Y-%m-%d")
        nombre_archivo = os.path.splitext(self.archivo_origen)[0] if self.archivo_origen else "experimento"
        nombre_reporte = f"paper_movimiento_aleatorio_{nombre_archivo}.txt"
        
    def generar_reporte_paper(self):
        """Genera un reporte estilo paper científico según solicita la guía."""
        print("\n📄 Generando paper científico...")
        
        fecha_actual = datetime.now().strftime("%Y-%m-%d")
        nombre_archivo = os.path.splitext(self.archivo_origen)[0] if self.archivo_origen else "experimento"
        nombre_reporte = f"paper_movimiento_aleatorio_{nombre_archivo}.txt"
        
    def generar_reporte_paper(self):
        """Genera un reporte estilo paper científico según solicita la guía."""
        print("\n📄 Generando paper científico...")
        
        # OBTENER CARPETA EXACTA DEL SCRIPT (mismo método que para PNG)
        import sys
        if hasattr(sys, '_getframe'):
            script_dir = os.path.dirname(os.path.abspath(sys._getframe().f_code.co_filename))
        else:
            script_dir = os.path.dirname(os.path.realpath(__file__))
        
        fecha_actual = datetime.now().strftime("%Y-%m-%d")
        nombre_archivo = os.path.splitext(self.archivo_origen)[0] if self.archivo_origen else "experimento"
        nombre_reporte = f"paper_movimiento_aleatorio_{nombre_archivo}.txt"
        
        # CREAR RUTA EN LA MISMA CARPETA QUE EL SCRIPT
        ruta_completa_reporte = os.path.join(script_dir, nombre_reporte)
        
        # Debug para mostrar dónde se guardará
        print(f"📁 Carpeta del script: {script_dir}")
        print(f"📁 Directorio de trabajo: {os.getcwd()}")
        print(f"📄 Archivo de reporte: {nombre_reporte}")
        print(f"🔗 Ruta completa del reporte: {ruta_completa_reporte}")
        
        msd_data = self.resultados['msd']
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        try:
            with open(ruta_completa_reporte, 'w', encoding='utf-8') as f:
                # Encabezado estilo paper
                f.write("="*80 + "\n")
                f.write("ANÁLISIS EXPERIMENTAL DE MOVIMIENTO ALEATORIO\n")
                f.write("Universidad Técnica Federico Santa María\n")
                f.write("Física Experimental - FIS 200\n")
                f.write("="*80 + "\n\n")
                
                f.write(f"Fecha: {fecha_actual}\n")
                f.write(f"Archivo de datos: {self.archivo_origen}\n")
                f.write(f"Duración del experimento: {self.tiempo[-1]-self.tiempo[0]:.1f} segundos\n")
                f.write(f"Número de mediciones: {self.n_puntos}\n\n")
                
                # RESUMEN
                f.write("RESUMEN\n")
                f.write("-"*40 + "\n")
                f.write("Se analizó el movimiento aleatorio de una esfera en un recipiente vibrante\n")
                f.write("utilizando análisis de video mediante software Tracker. Se determinó el\n")
                f.write(f"régimen dinámico como '{escala_data['regimen']}' con exponente de\n")
                f.write(f"escalamiento α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}.\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write(f"El coeficiente de difusión estimado es D = {escala_data['coeficiente_difusion']:.2e} m²/s.\n")
                f.write("\n")
                
                # INTRODUCCIÓN
                f.write("1. INTRODUCCIÓN\n")
                f.write("-"*40 + "\n")
                f.write("El movimiento browniano describe el desplazamiento aleatorio de partículas\n")
                f.write("microscópicas suspendidas en un fluido. En este experimento se simula este\n")
                f.write("comportamiento mediante esferas en un recipiente vibrante.\n\n")
                
                f.write("El desplazamiento cuadrático medio (MSD) se define como:\n")
                f.write("⟨r²(t)⟩ = ⟨[x(t) − x(0)]² + [y(t) − y(0)]²⟩\n\n")
                
                f.write("La evolución temporal sigue la ley de escala:\n")
                f.write("⟨r²(t)⟩ = A·t^α\n\n")
                
                f.write("donde α caracteriza el régimen dinámico:\n")
                f.write("• α < 1: Subdifusivo anómalo\n")
                f.write("• α = 1: Difusivo clásico (browniano)\n")
                f.write("• α > 1: Superdifusivo\n\n")
                
                # METODOLOGÍA
                f.write("2. METODOLOGÍA EXPERIMENTAL\n")
                f.write("-"*40 + "\n")
                f.write("• Se utilizó una esfera en recipiente vibrante\n")
                f.write("• Grabación de video durante 120 segundos\n")
                f.write("• Análisis de trayectoria con software Tracker\n")
                f.write("• Calibración dimensional del sistema\n")
                f.write(f"• Frecuencia de muestreo: {1/self.dt:.1f} Hz\n\n")
                
                # RESULTADOS
                f.write("3. RESULTADOS\n")
                f.write("-"*40 + "\n")
                
                f.write("3.1 Características del movimiento:\n")
                stats_desp = stats_data['desplazamiento']
                f.write(f"• Desplazamiento total: {stats_desp['desplazamiento_total']:.4f} m\n")
                f.write(f"• Desplazamiento máximo: {stats_desp['desplazamiento_maximo']:.4f} m\n")
                f.write(f"• Distancia recorrida: {stats_desp['distancia_recorrida']:.4f} m\n")
                f.write(f"• Velocidad media X: {stats_desp['velocidad_media_x']:.6f} m/s\n")
                f.write(f"• Velocidad media Y: {stats_desp['velocidad_media_y']:.6f} m/s\n\n")
                
                f.write("3.2 Análisis de ley de escala:\n")
                f.write(f"• Exponente α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}\n")
                f.write(f"• Constante A = {escala_data['A']:.2e} m²/s^α\n")
                f.write(f"• Coeficiente de correlación R² = {escala_data['r_squared']:.4f}\n")
                f.write(f"• p-valor = {escala_data['p_value']:.2e}\n")
                f.write(f"• Régimen identificado: {escala_data['regimen']}\n\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write("3.3 Coeficiente de difusión:\n")
                    f.write("Para régimen difusivo (α ≈ 1), se aplica:\n")
                    f.write("⟨r²(t)⟩ = 2dDt, donde d = 2 (2D)\n")
                    f.write(f"D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
                
                # DISCUSIÓN
                f.write("4. DISCUSIÓN\n")
                f.write("-"*40 + "\n")
                
                # Evaluación del comportamiento
                if 0.9 <= escala_data['alpha'] <= 1.1:
                    f.write("4.1 Comportamiento difusivo:\n")
                    f.write("El exponente α ≈ 1 confirma un comportamiento difusivo clásico,\n")
                    f.write("consistente con movimiento browniano. Esto indica que:\n")
                    f.write("• Las colisiones son aleatorias e independientes\n")
                    f.write("• No hay efectos de memoria significativos\n")
                    f.write("• El sistema ha alcanzado equilibrio estadístico\n\n")
                else:
                    f.write("4.1 Comportamiento anómalo:\n")
                    f.write(f"El exponente α = {escala_data['alpha']:.4f} indica desviación del\n")
                    f.write("comportamiento difusivo ideal. Posibles causas:\n")
                    if escala_data['alpha'] < 1:
                        f.write("• Restricciones geométricas del recipiente\n")
                        f.write("• Efectos de fricción o atrapamiento\n")
                        f.write("• Heterogeneidad en la vibración\n")
                    else:
                        f.write("• Movimiento balístico a cortos tiempos\n")
                        f.write("• Efectos de inercia de la esfera\n")
                        f.write("• Correlaciones en la vibración\n")
                    f.write("\n")
                
                f.write("4.2 Incertidumbres experimentales:\n")
                f.write("• Resolución espacial del video\n")
                f.write("• Precisión del software de tracking\n")
                f.write("• Efectos de bordes del recipiente\n")
                f.write("• Variaciones en la amplitud de vibración\n\n")
                
                # CONCLUSIONES
                f.write("5. CONCLUSIONES\n")
                f.write("-"*40 + "\n")
                f.write(f"1. Se identificó un régimen {escala_data['regimen'].lower()}\n")
                f.write(f"   con exponente de escalamiento α = {escala_data['alpha']:.4f}\n\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write(f"2. El coeficiente de difusión efectivo es D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
                
                f.write("3. El análisis log-log permite discriminar efectivamente\n")
                f.write("   entre diferentes regímenes de transporte\n\n")
                
                if 0.9 <= escala_data['alpha'] <= 1.1:
                    f.write("4. Los resultados son consistentes con la teoría\n")
                    f.write("   del movimiento browniano clásico\n")
                else:
                    f.write("4. Se observan desviaciones del comportamiento\n")
                    f.write("   browniano ideal que requieren mayor investigación\n")
                
                f.write("\n")
                
                # REFERENCIAS
                f.write("REFERENCIAS\n")
                f.write("-"*40 + "\n")
                f.write("[1] Einstein, A. (1905). Über die von der molekularkinetischen Theorie\n")
                f.write("    der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten\n")
                f.write("    suspendierten Teilchen. Ann. Phys. 17, 549-560.\n\n")
                f.write("[2] Metzler, R., & Klafter, J. (2000). The random walk's guide to\n")
                f.write("    anomalous diffusion. Physics Reports, 339(1), 1-77.\n\n")
                f.write("[3] Guía de Laboratorio FIS 200 - Universidad Técnica Federico Santa María\n")
                
                f.write("\n" + "="*80 + "\n")
            
            # Verificar que se guardó
            if os.path.exists(ruta_completa_reporte):
                print(f"✅ Paper científico guardado en: {ruta_completa_reporte}")
            else:
                print(f"❌ Error: El reporte no se guardó correctamente")
                
        except Exception as e:
            print(f"❌ Error generando paper: {e}")
            return None
            
        return ruta_completa_reporte
        
        msd_data = self.resultados['msd']
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        try:
            with open(ruta_completa_reporte, 'w', encoding='utf-8') as f:
                # Encabezado estilo paper
                f.write("="*80 + "\n")
                f.write("ANÁLISIS EXPERIMENTAL DE MOVIMIENTO ALEATORIO\n")
                f.write("Universidad Técnica Federico Santa María\n")
                f.write("Física Experimental - FIS 200\n")
                f.write("="*80 + "\n\n")
                
                f.write(f"Fecha: {fecha_actual}\n")
                f.write(f"Archivo de datos: {self.archivo_origen}\n")
                f.write(f"Duración del experimento: {self.tiempo[-1]-self.tiempo[0]:.1f} segundos\n")
                f.write(f"Número de mediciones: {self.n_puntos}\n\n")
                
                # RESUMEN
                f.write("RESUMEN\n")
                f.write("-"*40 + "\n")
                f.write("Se analizó el movimiento aleatorio de una esfera en un recipiente vibrante\n")
                f.write("utilizando análisis de video mediante software Tracker. Se determinó el\n")
                f.write(f"régimen dinámico como '{escala_data['regimen']}' con exponente de\n")
                f.write(f"escalamiento α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}.\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write(f"El coeficiente de difusión estimado es D = {escala_data['coeficiente_difusion']:.2e} m²/s.\n")
                f.write("\n")
                
                # INTRODUCCIÓN
                f.write("1. INTRODUCCIÓN\n")
                f.write("-"*40 + "\n")
                f.write("El movimiento browniano describe el desplazamiento aleatorio de partículas\n")
                f.write("microscópicas suspendidas en un fluido. En este experimento se simula este\n")
                f.write("comportamiento mediante esferas en un recipiente vibrante.\n\n")
                
                f.write("El desplazamiento cuadrático medio (MSD) se define como:\n")
                f.write("⟨r²(t)⟩ = ⟨[x(t) − x(0)]² + [y(t) − y(0)]²⟩\n\n")
                
                f.write("La evolución temporal sigue la ley de escala:\n")
                f.write("⟨r²(t)⟩ = A·t^α\n\n")
                
                f.write("donde α caracteriza el régimen dinámico:\n")
                f.write("• α < 1: Subdifusivo anómalo\n")
                f.write("• α = 1: Difusivo clásico (browniano)\n")
                f.write("• α > 1: Superdifusivo\n\n")
                
                # METODOLOGÍA
                f.write("2. METODOLOGÍA EXPERIMENTAL\n")
                f.write("-"*40 + "\n")
                f.write("• Se utilizó una esfera en recipiente vibrante\n")
                f.write("• Grabación de video durante 120 segundos\n")
                f.write("• Análisis de trayectoria con software Tracker\n")
                f.write("• Calibración dimensional del sistema\n")
                f.write(f"• Frecuencia de muestreo: {1/self.dt:.1f} Hz\n\n")
                
                # RESULTADOS
                f.write("3. RESULTADOS\n")
                f.write("-"*40 + "\n")
                
                f.write("3.1 Características del movimiento:\n")
                stats_desp = stats_data['desplazamiento']
                f.write(f"• Desplazamiento total: {stats_desp['desplazamiento_total']:.4f} m\n")
                f.write(f"• Desplazamiento máximo: {stats_desp['desplazamiento_maximo']:.4f} m\n")
                f.write(f"• Distancia recorrida: {stats_desp['distancia_recorrida']:.4f} m\n")
                f.write(f"• Velocidad media X: {stats_desp['velocidad_media_x']:.6f} m/s\n")
                f.write(f"• Velocidad media Y: {stats_desp['velocidad_media_y']:.6f} m/s\n\n")
                
                f.write("3.2 Análisis de ley de escala:\n")
                f.write(f"• Exponente α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}\n")
                f.write(f"• Constante A = {escala_data['A']:.2e} m²/s^α\n")
                f.write(f"• Coeficiente de correlación R² = {escala_data['r_squared']:.4f}\n")
                f.write(f"• p-valor = {escala_data['p_value']:.2e}\n")
                f.write(f"• Régimen identificado: {escala_data['regimen']}\n\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write("3.3 Coeficiente de difusión:\n")
                    f.write("Para régimen difusivo (α ≈ 1), se aplica:\n")
                    f.write("⟨r²(t)⟩ = 2dDt, donde d = 2 (2D)\n")
                    f.write(f"D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
                
                # DISCUSIÓN
                f.write("4. DISCUSIÓN\n")
                f.write("-"*40 + "\n")
                
                # Evaluación del comportamiento
                if 0.9 <= escala_data['alpha'] <= 1.1:
                    f.write("4.1 Comportamiento difusivo:\n")
                    f.write("El exponente α ≈ 1 confirma un comportamiento difusivo clásico,\n")
                    f.write("consistente con movimiento browniano. Esto indica que:\n")
                    f.write("• Las colisiones son aleatorias e independientes\n")
                    f.write("• No hay efectos de memoria significativos\n")
                    f.write("• El sistema ha alcanzado equilibrio estadístico\n\n")
                else:
                    f.write("4.1 Comportamiento anómalo:\n")
                    f.write(f"El exponente α = {escala_data['alpha']:.4f} indica desviación del\n")
                    f.write("comportamiento difusivo ideal. Posibles causas:\n")
                    if escala_data['alpha'] < 1:
                        f.write("• Restricciones geométricas del recipiente\n")
                        f.write("• Efectos de fricción o atrapamiento\n")
                        f.write("• Heterogeneidad en la vibración\n")
                    else:
                        f.write("• Movimiento balístico a cortos tiempos\n")
                        f.write("• Efectos de inercia de la esfera\n")
                        f.write("• Correlaciones en la vibración\n")
                    f.write("\n")
                
                f.write("4.2 Incertidumbres experimentales:\n")
                f.write("• Resolución espacial del video\n")
                f.write("• Precisión del software de tracking\n")
                f.write("• Efectos de bordes del recipiente\n")
                f.write("• Variaciones en la amplitud de vibración\n\n")
                
                # CONCLUSIONES
                f.write("5. CONCLUSIONES\n")
                f.write("-"*40 + "\n")
                f.write(f"1. Se identificó un régimen {escala_data['regimen'].lower()}\n")
                f.write(f"   con exponente de escalamiento α = {escala_data['alpha']:.4f}\n\n")
                
                if escala_data['coeficiente_difusion']:
                    f.write(f"2. El coeficiente de difusión efectivo es D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
                
                f.write("3. El análisis log-log permite discriminar efectivamente\n")
                f.write("   entre diferentes regímenes de transporte\n\n")
                
                if 0.9 <= escala_data['alpha'] <= 1.1:
                    f.write("4. Los resultados son consistentes con la teoría\n")
                    f.write("   del movimiento browniano clásico\n")
                else:
                    f.write("4. Se observan desviaciones del comportamiento\n")
                    f.write("   browniano ideal que requieren mayor investigación\n")
                
                f.write("\n")
                
                # REFERENCIAS
                f.write("REFERENCIAS\n")
                f.write("-"*40 + "\n")
                f.write("[1] Einstein, A. (1905). Über die von der molekularkinetischen Theorie\n")
                f.write("    der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten\n")
                f.write("    suspendierten Teilchen. Ann. Phys. 17, 549-560.\n\n")
                f.write("[2] Metzler, R., & Klafter, J. (2000). The random walk's guide to\n")
                f.write("    anomalous diffusion. Physics Reports, 339(1), 1-77.\n\n")
                f.write("[3] Guía de Laboratorio FIS 200 - Universidad Técnica Federico Santa María\n")
                
                f.write("\n" + "="*80 + "\n")
            
            # Verificar que se guardó
            if os.path.exists(ruta_completa_reporte):
                print(f"✅ Paper científico guardado en: {ruta_completa_reporte}")
            else:
                print(f"❌ Error: El reporte no se guardó correctamente")
                
        except Exception as e:
            print(f"❌ Error generando paper: {e}")
            return None
            
        return ruta_completa_reporte
        
        msd_data = self.resultados['msd']
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        with open(ruta_completa_reporte, 'w', encoding='utf-8') as f:
            # Encabezado estilo paper
            f.write("="*80 + "\n")
            f.write("ANÁLISIS EXPERIMENTAL DE MOVIMIENTO ALEATORIO\n")
            f.write("Universidad Técnica Federico Santa María\n")
            f.write("Física Experimental - FIS 200\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"Fecha: {fecha_actual}\n")
            f.write(f"Archivo de datos: {self.archivo_origen}\n")
            f.write(f"Duración del experimento: {self.tiempo[-1]-self.tiempo[0]:.1f} segundos\n")
            f.write(f"Número de mediciones: {self.n_puntos}\n\n")
            
            # RESUMEN
            f.write("RESUMEN\n")
            f.write("-"*40 + "\n")
            f.write("Se analizó el movimiento aleatorio de una esfera en un recipiente vibrante\n")
            f.write("utilizando análisis de video mediante software Tracker. Se determinó el\n")
            f.write(f"régimen dinámico como '{escala_data['regimen']}' con exponente de\n")
            f.write(f"escalamiento α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}.\n")
            
            if escala_data['coeficiente_difusion']:
                f.write(f"El coeficiente de difusión estimado es D = {escala_data['coeficiente_difusion']:.2e} m²/s.\n")
            f.write("\n")
            
            # INTRODUCCIÓN
            f.write("1. INTRODUCCIÓN\n")
            f.write("-"*40 + "\n")
            f.write("El movimiento browniano describe el desplazamiento aleatorio de partículas\n")
            f.write("microscópicas suspendidas en un fluido. En este experimento se simula este\n")
            f.write("comportamiento mediante esferas en un recipiente vibrante.\n\n")
            
            f.write("El desplazamiento cuadrático medio (MSD) se define como:\n")
            f.write("⟨r²(t)⟩ = ⟨[x(t) − x(0)]² + [y(t) − y(0)]²⟩\n\n")
            
            f.write("La evolución temporal sigue la ley de escala:\n")
            f.write("⟨r²(t)⟩ = A·t^α\n\n")
            
            f.write("donde α caracteriza el régimen dinámico:\n")
            f.write("• α < 1: Subdifusivo anómalo\n")
            f.write("• α = 1: Difusivo clásico (browniano)\n")
            f.write("• α > 1: Superdifusivo\n\n")
            
            # METODOLOGÍA
            f.write("2. METODOLOGÍA EXPERIMENTAL\n")
            f.write("-"*40 + "\n")
            f.write("• Se utilizó una esfera en recipiente vibrante\n")
            f.write("• Grabación de video durante 120 segundos\n")
            f.write("• Análisis de trayectoria con software Tracker\n")
            f.write("• Calibración dimensional del sistema\n")
            f.write(f"• Frecuencia de muestreo: {1/self.dt:.1f} Hz\n\n")
            
            # RESULTADOS
            f.write("3. RESULTADOS\n")
            f.write("-"*40 + "\n")
            
            f.write("3.1 Características del movimiento:\n")
            stats_desp = stats_data['desplazamiento']
            f.write(f"• Desplazamiento total: {stats_desp['desplazamiento_total']:.4f} m\n")
            f.write(f"• Desplazamiento máximo: {stats_desp['desplazamiento_maximo']:.4f} m\n")
            f.write(f"• Distancia recorrida: {stats_desp['distancia_recorrida']:.4f} m\n")
            f.write(f"• Velocidad media X: {stats_desp['velocidad_media_x']:.6f} m/s\n")
            f.write(f"• Velocidad media Y: {stats_desp['velocidad_media_y']:.6f} m/s\n\n")
            
            f.write("3.2 Análisis de ley de escala:\n")
            f.write(f"• Exponente α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}\n")
            f.write(f"• Constante A = {escala_data['A']:.2e} m²/s^α\n")
            f.write(f"• Coeficiente de correlación R² = {escala_data['r_squared']:.4f}\n")
            f.write(f"• p-valor = {escala_data['p_value']:.2e}\n")
            f.write(f"• Régimen identificado: {escala_data['regimen']}\n\n")
            
            if escala_data['coeficiente_difusion']:
                f.write("3.3 Coeficiente de difusión:\n")
                f.write("Para régimen difusivo (α ≈ 1), se aplica:\n")
                f.write("⟨r²(t)⟩ = 2dDt, donde d = 2 (2D)\n")
                f.write(f"D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
            
            # DISCUSIÓN
            f.write("4. DISCUSIÓN\n")
            f.write("-"*40 + "\n")
            
            # Evaluación del comportamiento
            if 0.9 <= escala_data['alpha'] <= 1.1:
                f.write("4.1 Comportamiento difusivo:\n")
                f.write("El exponente α ≈ 1 confirma un comportamiento difusivo clásico,\n")
                f.write("consistente con movimiento browniano. Esto indica que:\n")
                f.write("• Las colisiones son aleatorias e independientes\n")
                f.write("• No hay efectos de memoria significativos\n")
                f.write("• El sistema ha alcanzado equilibrio estadístico\n\n")
            else:
                f.write("4.1 Comportamiento anómalo:\n")
                f.write(f"El exponente α = {escala_data['alpha']:.4f} indica desviación del\n")
                f.write("comportamiento difusivo ideal. Posibles causas:\n")
                if escala_data['alpha'] < 1:
                    f.write("• Restricciones geométricas del recipiente\n")
                    f.write("• Efectos de fricción o atrapamiento\n")
                    f.write("• Heterogeneidad en la vibración\n")
                else:
                    f.write("• Movimiento balístico a cortos tiempos\n")
                    f.write("• Efectos de inercia de la esfera\n")
                    f.write("• Correlaciones en la vibración\n")
                f.write("\n")
            
            f.write("4.2 Incertidumbres experimentales:\n")
            f.write("• Resolución espacial del video\n")
            f.write("• Precisión del software de tracking\n")
            f.write("• Efectos de bordes del recipiente\n")
            f.write("• Variaciones en la amplitud de vibración\n\n")
            
            # CONCLUSIONES
            f.write("5. CONCLUSIONES\n")
            f.write("-"*40 + "\n")
            f.write(f"1. Se identificó un régimen {escala_data['regimen'].lower()}\n")
            f.write(f"   con exponente de escalamiento α = {escala_data['alpha']:.4f}\n\n")
            
            if escala_data['coeficiente_difusion']:
                f.write(f"2. El coeficiente de difusión efectivo es D = {escala_data['coeficiente_difusion']:.2e} m²/s\n\n")
            
            f.write("3. El análisis log-log permite discriminar efectivamente\n")
            f.write("   entre diferentes regímenes de transporte\n\n")
            
            if 0.9 <= escala_data['alpha'] <= 1.1:
                f.write("4. Los resultados son consistentes con la teoría\n")
                f.write("   del movimiento browniano clásico\n")
            else:
                f.write("4. Se observan desviaciones del comportamiento\n")
                f.write("   browniano ideal que requieren mayor investigación\n")
            
            f.write("\n")
            
            # REFERENCIAS
            f.write("REFERENCIAS\n")
            f.write("-"*40 + "\n")
            f.write("[1] Einstein, A. (1905). Über die von der molekularkinetischen Theorie\n")
            f.write("    der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten\n")
            f.write("    suspendierten Teilchen. Ann. Phys. 17, 549-560.\n\n")
            f.write("[2] Metzler, R., & Klafter, J. (2000). The random walk's guide to\n")
            f.write("    anomalous diffusion. Physics Reports, 339(1), 1-77.\n\n")
            f.write("[3] Guía de Laboratorio FIS 200 - Universidad Técnica Federico Santa María\n")
            
            f.write("\n" + "="*80 + "\n")
            
        print(f"📄 Paper científico generado: {ruta_completa_reporte}")
        return ruta_completa_reporte
    
    def mostrar_resumen_resultados(self):
        """Muestra un resumen de resultados en pantalla."""
        print("\n" + "="*70)
        print("📊 RESUMEN DE RESULTADOS - MOVIMIENTO ALEATORIO FIS 200")
        print("="*70)
        
        escala_data = self.resultados['ley_escala']
        stats_data = self.resultados['estadisticas']
        
        print(f"\n🎯 EXPERIMENTO:")
        print(f"   • Archivo: {self.archivo_origen}")
        print(f"   • Duración: {self.tiempo[-1]-self.tiempo[0]:.1f} segundos")
        print(f"   • Puntos: {self.n_puntos}")
        print(f"   • Δt: {self.dt:.4f} s")
        
        print(f"\n📐 ANÁLISIS DE LEY DE ESCALA:")
        print(f"   • Exponente α = {escala_data['alpha']:.4f} ± {escala_data['std_error']:.4f}")
        print(f"   • R² = {escala_data['r_squared']:.4f}")
        print(f"   • Régimen: {escala_data['regimen']}")
        
        if escala_data['coeficiente_difusion']:
            print(f"\n🌊 COEFICIENTE DE DIFUSIÓN:")
            print(f"   • D = {escala_data['coeficiente_difusion']:.2e} m²/s")
        
        stats_desp = stats_data['desplazamiento']
        print(f"\n📏 CARACTERÍSTICAS DEL MOVIMIENTO:")
        print(f"   • Desplazamiento final: {stats_desp['desplazamiento_total']:.4f} m")
        print(f"   • Distancia recorrida: {stats_desp['distancia_recorrida']:.4f} m")
        
        print(f"\n🔬 EVALUACIÓN TEÓRICA:")
        if 0.9 <= escala_data['alpha'] <= 1.1:
            print("   ✅ Comportamiento difusivo clásico (browniano)")
        elif escala_data['alpha'] < 1:
            print("   ⚠️ Comportamiento subdifusivo (movimiento restringido)")
        elif 1 < escala_data['alpha'] < 2:
            print("   ⚠️ Comportamiento superdifusivo (movimiento persistente)")
        elif escala_data['alpha'] >= 2:
            print("   ⚠️ Comportamiento balístico")
        
        print("="*70)
    
    def ejecutar_analisis_completo(self):
        """Ejecuta el análisis completo siguiendo la guía FIS 200."""
        print("🔬 ANÁLISIS DE MOVIMIENTO ALEATORIO - GUÍA FIS 200 UTFSM")
        print("="*60)
        
        # 1. Cargar datos de Tracker
        if not self.cargar_datos_tracker():
            return False
        
        # 2. Cálculo de MSD
        self.calcular_msd_segun_guia()
        
        # 3. Análisis estadístico básico
        self.analisis_estadistico_basico()
        
        # 4. Análisis de ley de escala (pasos 12-14 de la guía)
        self.analisis_ley_escala()
        
        # 5. Crear gráficos según metodología
        self.crear_graficos_segun_guia()
        
        # 6. Generar paper científico
        self.generar_reporte_paper()
        
        # 7. Mostrar resumen
        self.mostrar_resumen_resultados()
        
        print(f"\n✅ ANÁLISIS COMPLETADO SEGÚN GUÍA FIS 200")
        print("="*60)
        nombre_base = os.path.splitext(self.archivo_origen)[0]
        print("📁 Archivos generados:")
        print(f"   • Gráficos: analisis_movimiento_aleatorio_{nombre_base}.png")
        print(f"   • Paper: paper_movimiento_aleatorio_{nombre_base}.txt")
        print("\n💡 Utilice estos archivos para su informe científico")
        
        return True

def main():
    """Función principal del programa."""
    print("🧬 ANALIZADOR DE MOVIMIENTO ALEATORIO")
    print("Universidad Técnica Federico Santa María - FIS 200")
    print("="*50)
    print("📋 Sigue la metodología específica de la guía de laboratorio")
    print("📊 Análisis de MSD y determinación de regímenes dinámicos")
    print("-"*50)
    
    # Crear analizador
    analizador = AnalizadorMovimientoAleatorio()
    
    # Ejecutar análisis completo
    exito = analizador.ejecutar_analisis_completo()
    
    if exito:
        print("\n🎉 Análisis completado exitosamente")
        print("📄 Revise los archivos generados para su informe")
        print("📊 Los gráficos incluyen las ecuaciones teóricas")
        print("📝 El paper sigue el formato científico solicitado")
    else:
        print("\n❌ No se pudo completar el análisis")
        print("💡 Verifique que tiene datos de Tracker en formato Excel/CSV")
        print("📋 Formato esperado: columnas de tiempo, posición X, posición Y")

if __name__ == "__main__":
    main()