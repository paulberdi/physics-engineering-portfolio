# Verificación numérica

Esta carpeta documenta resultados **previos** de pruebas automáticas, diagnósticos y controles de consistencia.

- **Pruebas/** contiene el registro histórico de pytest.
- **Diagnosticos/** contiene diagnósticos numéricos, audiovisuales y del entorno.
- **Registro_Validacion.md** resume las comprobaciones originales.

Para verificar el estado actual del código en tu computador:

~~~bash
cd codigo
python -m pip install -r requirements.txt
python -m pytest
~~~

Estos archivos no implican que las pruebas se hayan reejecutado tras la reorganización. **Tampoco equivalen a validación experimental del sistema físico.**
