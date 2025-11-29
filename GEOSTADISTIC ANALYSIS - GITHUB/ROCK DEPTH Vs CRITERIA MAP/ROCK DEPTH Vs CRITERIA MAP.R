#MAPA DE PROFUNDIDAD DEL BASAMENTO ROCOSO (Vs > 720 m/s) DENTRO DEL POLIGONO PROVIDENCIA-COUNTRY, EN GUADALAJARA, JALISCO, CON BASE
#EN LA UBICACIÓN DE LOS SONDEOS GEOSÍSMICOS TIPO DOWN-HOLE REALIZADOS

#Importar la Base de Datos
Rock_Vs_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/ROCK DEPTH Vs CRITERIA MAP/Vs Rock Depth.csv")
View(Rock_Vs_Depth)

#----------------------------------------------------------------------------------------------------------

#MODELADO GEOESTADÍSTICO

install.packages(c("ggplot2", "sf", "gstat", "sp", "raster"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)


# Crear un objeto SpatialPointsDataFrame a partir de los datos
coordinates(Rock_Vs_Depth) <- ~Long + Lat

# Definir el modelo de variograma
variograma <- variogram(Prof_rock ~ 1, Rock_Vs_Depth)

# Ajustar el modelo de variograma
modelo_variograma <- fit.variogram(
  variograma, 
  vgm(psill = 1, model = "Sph", range = 100, nugget = 0.1)
)

# Expandir los límites de la malla de interpolación
expand_factor <- 0.001  # Ajusta este valor para definir cuánto se extiende la interpolación más allá de los sondeos
bbox <- st_bbox(Rock_Depth)  # Obtener los límites actuales de los datos
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)  # Expansión del área de interpolación

# Generar una nueva malla de puntos en el área expandida
grid_expanded <- spsample(as_Spatial(st_as_sfc(bbox_expanded)), type = "regular", n = 1000000)  # Aumentar el número de puntos en la malla para mejor cobertura
gridded(grid_expanded) <- TRUE  # Convertir la malla en un objeto gridded adecuado para interpolación

# Realizar la interpolación mediante kriging
kriging_result_expanded <- krige(Prof_rock ~ 1, Rock_Vs_Depth, grid_expanded, model = modelo_variograma, nmax=1)

# Convertir los resultados a un data frame para ggplot2
kriging_df_expanded <- as.data.frame(kriging_result_expanded)
colnames(kriging_df_expanded) <- c("Long", "Lat", "Prof_rock")  # Renombrar columnas para ggplot2
kriging_df_expanded <- na.omit(kriging_df_expanded)

# Redondear valores para mejorar la visualización
kriging_df_expanded$Prof_rock <- round(kriging_df_expanded$Prof_rock, 2)

# Asegurarse de que Rock_Depth es un data frame
Rock_Vs_Depth_df <- as.data.frame(Rock_Vs_Depth)


#----------------------------------------------------------------------------------------------------------

#EXPORTAR LA INTERPOLACIÓN A QGIS

# Convertir los resultados de kriging a un raster
kriging_raster <- rasterFromXYZ(kriging_df_expanded[, c("Long", "Lat", "Prof_rock")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(kriging_raster, filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_output.tif", format = "GTiff", overwrite = TRUE)

