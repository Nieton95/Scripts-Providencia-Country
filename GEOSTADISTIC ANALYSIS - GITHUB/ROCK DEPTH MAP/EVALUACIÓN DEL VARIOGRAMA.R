# =======================================================================================
# ANÁLISIS DEL VARIOGRAMA PARA EL MAPA DE LA PROFUNDIDAD DE LA ROCA BASAL
# POLÍGONO PROVIDENCIA-COUNTRY, GUADALAJARA, JALISCO
# =======================================================================================

# --- Importar la Base de Datos ---
Rock_Depth <- read.csv("C:/Users/GProyectos/Desktop/Tésis Iván/GEOSTADISTIC ANALYSIS/ROCK DEPTH MAP/BD Bedrock Depth.csv")
View(Rock_Depth)

# --- Cargar librerías necesarias ---
install.packages(c("ggplot2", "sf", "gstat", "sp", "raster", "osmdata", "cowplot"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)
library(osmdata)
library(cowplot)

# =======================================================================================
# CONVERSIÓN A COORDENADAS MÉTRICAS (UTM ZONA 13N)
# =======================================================================================

# Crear objeto sf con las columnas de coordenadas
Rock_Depth_sf <- st_as_sf(Rock_Depth, coords = c("Long", "Lat"), crs = 4326)

# Transformar a sistema de coordenadas proyectadas (UTM Zona 13N - EPSG:6372)
Rock_Depth_sf <- st_transform(Rock_Depth_sf, 6372)

# Convertir a objeto SpatialPointsDataFrame para usar con gstat
Rock_Depth <- as(Rock_Depth_sf, "Spatial")

# =======================================================================================
# CÁLCULO Y AJUSTE DEL VARIOGRAMA
# =======================================================================================

# Cálculo del variograma empírico
variograma <- variogram(Prof_Inicial ~ 1, Rock_Depth)

# Ajuste con tres modelos distintos
modelo_esferico <- fit.variogram(variograma, vgm(model = "Sph"))
modelo_exponencial <- fit.variogram(variograma, vgm(model = "Exp"))
modelo_gaussiano <- fit.variogram(variograma, vgm(model = "Gau"))

# Graficar los tres modelos juntos para comparar
plot(variograma, modelo_esferico, main = "Variograma - Modelo Esférico")
plot(variograma, modelo_exponencial, main = "Variograma - Modelo Exponencial")
plot(variograma, modelo_gaussiano, main = "Variograma - Modelo Gaussiano")

# Comparación visual en un solo gráfico (opcional)
plot(variograma, modelo_esferico, col = "blue", main = "Comparación de Modelos")
plot(variograma, modelo_exponencial, add = TRUE, col = "red")
plot(variograma, modelo_gaussiano, add = TRUE, col = "green")
legend("bottomright",
       legend = c("Esférico", "Exponencial", "Gaussiano"),
       col = c("blue", "red", "green"),
       lty = 1, cex = 0.9)


# Validación cruzada para cada modelo
cv_esf <- krige.cv(Prof_Inicial ~ 1, Rock_Depth, modelo_esferico)
cv_exp <- krige.cv(Prof_Inicial ~ 1, Rock_Depth, modelo_exponencial)
cv_gau <- krige.cv(Prof_Inicial ~ 1, Rock_Depth, modelo_gaussiano)

# Comparar el error cuadrático medio (RMSE). 
# El modelo con el menor RMSE suele ser el que mejor representa la estructura espacial de tus datos.
data.frame(
  Modelo = c("Esférico", "Exponencial", "Gaussiano"),
  RMSE = c(
    sqrt(mean(cv_esf$residual^2, na.rm = TRUE)),
    sqrt(mean(cv_exp$residual^2, na.rm = TRUE)),
    sqrt(mean(cv_gau$residual^2, na.rm = TRUE))
  )
)


# # Calcular variograma empírico
# variograma <- variogram(Prof_Inicial ~ 1, Rock_Depth)
# 
# # Graficar variograma empírico
# plot(variograma, main = "Variograma Empírico - Profundidad de Roca",
#      xlab = "Distancia (m)", ylab = "Semivarianza")
# 
# # Ajustar modelo teórico (esférico)
# modelo_variograma <- fit.variogram(variograma, vgm(model = "Sph"))
# 
# # Graficar variograma con modelo ajustado
# plot(variograma, modelo_variograma, main = "Variograma Ajustado (Modelo Esférico)",
#      xlab = "Distancia (m)", ylab = "Semivarianza")
# 
# # Mostrar los parámetros del modelo ajustado
# print(modelo_variograma)
# 
# # =======================================================================================
# # INTERPRETACIÓN DEL ALCANCE (RANGE)
# # =======================================================================================
# 
# # Extraer el valor del "range" del modelo ajustado
# alcance <- modelo_variograma$range[2]
# 
# # Imprimir interpretación automática
# cat("\n-----------------------------------------------\n")
# cat("📏 Distancia máxima de confianza (alcance del variograma):", round(alcance, 2), "metros\n")
# cat("🔎 Esto indica que los valores de profundidad de roca presentan correlación espacial\n")
# cat("     significativa hasta aproximadamente", round(alcance, 0), "m de distancia.\n")
# cat("     Más allá de esa distancia, la interpolación kriging pierde fiabilidad.\n")
# cat("-----------------------------------------------\n")











