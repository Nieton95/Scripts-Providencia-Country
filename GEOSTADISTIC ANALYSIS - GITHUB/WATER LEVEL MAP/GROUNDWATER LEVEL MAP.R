#MAPA DEL TIRANTE DE AGUA SUBTERRÁNEA DENTRO DEL POLIGONO PROVIDENCIA-COUNTRY EN GUADALAJARA, JALISCO, CON BASE
#EN LA UBICACIÓN DE LOS SONDEOS REALIZADOS

#Importar la Base de Datos
Water_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/WATER LEVEL MAP/DB Sondeo Agua Subterránea.csv")
View(Water_Depth)

#----------------------------------------------------------------------------------------------------------

#MODELADO GEOESTADÍSTICO

install.packages(c("ggplot2", "sf", "gstat", "sp", "raster"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)


# Crear un objeto SpatialPointsDataFrame a partir de los datos
coordinates(Water_Depth) <- ~Long + Lat

# Definir el modelo de variograma
variograma <- variogram(Prof_agua ~ 1, Water_Depth)

# Ajustar el modelo de variograma 
modelo_variograma <- fit.variogram(
  variograma, 
  vgm(psill = 1, model = "Sph", range = 100, nugget = 0.1)
)

# Expandir los límites de la malla de interpolación
expand_factor <- 0.001  # Ajusta este valor para definir cuánto se extiende la interpolación más allá de los sondeos
bbox <- st_bbox(Water_Depth)  # Obtener los límites actuales de los datos
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)  # Expansión del área de interpolación

# Generar una nueva malla de puntos en el área expandida
grid_expanded <- spsample(as_Spatial(st_as_sfc(bbox_expanded)), type = "regular", n = 1000000)  # Aumentar el número de puntos en la malla para mejor cobertura
gridded(grid_expanded) <- TRUE  # Convertir la malla en un objeto gridded adecuado para interpolación

# Realizar la interpolación mediante kriging
kriging_result_expanded <- krige(Prof_agua ~ 1, Water_Depth, grid_expanded, model = modelo_variograma, nmax = 1)  # Aplicar kriging con nmax=10 para una extrapolación más suave

# Convertir los resultados a un data frame para ggplot2
kriging_df_expanded <- as.data.frame(kriging_result_expanded)
colnames(kriging_df_expanded) <- c("Long", "Lat", "Prof_agua")  # Renombrar columnas para ggplot2
kriging_df_expanded <- na.omit(kriging_df_expanded)

# Redondear valores para mejorar la visualización
kriging_df_expanded$Prof_agua <- round(kriging_df_expanded$Prof_agua, 2)

Water_Depth_df <- as.data.frame(Water_Depth)


#----------------------------------------------------------------------------------------------------------

#EXPORTAR LA INTERPOLACIÓN A QGIS

# Convertir los resultados de kriging a un raster
kriging_raster <- rasterFromXYZ(kriging_df_expanded[, c("Long", "Lat", "Prof_agua")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(kriging_raster, filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_output.tif", format = "GTiff", overwrite = TRUE)





