**Interpolación Espacial de Información Geotécnico-Geológica del Subsuelo

Polígono Providencia–Country, Guadalajara, Jalisco**

Este repositorio contiene los scripts utilizados para el procesamiento, análisis geoestadístico e interpolación espacial de la información geotécnica y geológica recopilada dentro del polígono Providencia–Country en Guadalajara, Jalisco, México.
El objetivo principal es generar productos geoespaciales que permitan visualizar, interpretar y comprender la variabilidad del subsuelo, mediante la integración de sondeos geotécnicos, ensayos geofísicos y mediciones hidrológicas realizadas en la zona durante los últimos 12 años.

📌 Contenido del repositorio

El repositorio incluye:

Scripts en R para:

Limpieza y preparación de datos espaciales

Generación de variogramas y ajustes de modelos

Interpolaciones mediante kriging ordinario

Conversión de resultados a formatos ráster (GeoTIFF)

Visualizaciones preliminares con ggplot2

Archivos de apoyo (CSV) con coordenadas y valores medidos en campo y laboratorio

Modelos espaciales utilizados para la construcción de mapas geotécnicos y geofísicos

📍 Información procesada

Los análisis están basados en:

69 sondeos SPT

9 ensayos geosísmicos Down-Hole

11 piezómetros

21 pozos de monitoreo

Toda la información proviene del acervo técnico generado por la empresa a lo largo de diversos proyectos desarrollados en la zona de estudio.

🗺️ Productos generados

Mediante los scripts incluidos, se generan:

Mapas de profundidad a roca basal

Mapas de variación de Vs30 y detección de estrato denso (N>50)

Mapas de resistencia UCS de la roca

Archivos ráster compatibles con QGIS para maquetación y análisis espacial avanzado

⚙️ Dependencias principales

Los scripts fueron desarrollados principalmente con las siguientes librerías de R:

gstat – Geoestadística y kriging

sp / sf – Manejo de objetos espaciales

raster – Creación y exportación de modelos ráster

ggplot2 – Visualización preliminar de resultados

viridis – Escalas de color perceptualmente uniformes

📄 Objetivo del repositorio

Este repositorio sirve como respaldo técnico y científico para la investigación de maestría sobre la caracterización del subsuelo en el área de Providencia–Country, y facilita:

La revisión transparente del proceso analítico

La reproducibilidad de los resultados

La futura actualización del modelo conforme se integren nuevos sondeos o ensayos en la zona

🤝 Contribuciones

Este repositorio está abierto a mejoras, sugerencias y ampliaciones futuras conforme se genere nueva información geotécnica en el polígono de estudio.

📧 Contacto

Para preguntas, dudas o colaboración técnica:
Iván Mariano Nieto Cárdenas — Ingeniero Civil / Geotecnia & Geociencias Aplicadas
ivannieto_95@hotmail.com
