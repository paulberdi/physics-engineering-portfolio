# Entrega y presentación

Las presentaciones finales se conservaron completas, pero se guardan como partes binarias debido a su tamaño:

- **Presentacion_Final.pdf.parts/** — PDF de diapositivas.
- **Presentacion_Final.pptx.parts/** — presentación editable PowerPoint.

Para reconstruirlas, desde la carpeta principal del proyecto ejecuta:

~~~bash
python reconstruir_multimedia.py
~~~

La reconstrucción verifica el tamaño y SHA-256 de cada archivo. No hay un PDF o PPTX directamente descargable hasta ejecutar este comando sobre una copia local del repositorio.

El programa ejecutable está en [codigo/](../codigo/); la animación se recupera en [resultados/Animaciones_MP4/](../resultados/Animaciones_MP4/). Consulta las instrucciones generales en el [README principal](../README.md).
