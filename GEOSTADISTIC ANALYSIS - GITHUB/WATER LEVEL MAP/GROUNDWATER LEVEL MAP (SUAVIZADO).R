#MAPA DEL TIRANTE DE AGUA SUBTERRÁNEA DENTRO DEL POLIGONO PROVIDENCIA-COUNTRY EN GUADALAJARA, JALISCO, CON BASE
#EN LA UBICACIÓN DE LOS SONDEOS REALIZADOS

# Cargar librerías necesarias
install.packages(c("ggplot2", "sf", "gstat", "sp", "raster"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)

#Importar la Base de Datos
Water_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/WATER LEVEL MAP/DB Sondeo Agua Subterránea.csv")

# Crear un objeto SpatialPointsDataFrame a partir de los datos
coordinates(Water_Depth) <- ~Long + Lat

# Definir el modelo de variograma
variograma <- variogram(Prof_agua ~ 1, Water_Depth)

# Ajustar el modelo de variograma 
modelo_variograma <- fit.variogram(
  variograma, 
  vgm(model = "Sph")
)

# Expandir los límites de la malla de interpolación
expand_factor <- 0.001  # Ajusta este valor para definir cuánto se extiende la interpolación más allá de los sondeos
bbox <- st_bbox(Water_Depth)  # Obtener los límites actuales de los datos
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)  # Expansión del área de interpolación

# Definir resolución de la malla (ajusta según la densidad deseada)
res <- 0.0001  # ~10 m de resolución en GDL aprox.
x_seq <- seq(bbox_expanded["xmin"], bbox_expanded["xmax"], by = res)
y_seq <- seq(bbox_expanded["ymin"], bbox_expanded["ymax"], by = res)

grid_df <- expand.grid(Long = x_seq, Lat = y_seq)
coordinates(grid_df) <- ~Long + Lat
gridded(grid_df) <- TRUE

# Realizar la interpolación mediante kriging
kriging_result <- krige(Prof_agua ~ 1, Water_Depth, grid_df, model = modelo_variograma, nmax = 10)  # Aplicar kriging con nmax=10 para una extrapolación más suave

# Convertir los resultados a un data frame para ggplot2
kriging_df <- as.data.frame(kriging_result)
colnames(kriging_df) [1:2] <- c("Long", "Lat")  # Renombrar columnas para ggplot2

# VISUALIZACIÓN EN R
# ------------------------------------------------------------------------------
ggplot() +
  geom_tile(data = kriging_df, aes(x = Long, y = Lat, fill = var1.pred)) +
  scale_fill_viridis_c(option = "C", name = "Profundidad (m)", direction = -1) +
  labs(
    title = "Mapa de profundidad de agua subterránea (Kriging)",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#----------------------------------------------------------------------------------------------------------

#EXPORTAR LA INTERPOLACIÓN A QGIS

# Convertir los resultados de kriging a un raster
kriging_raster <- rasterFromXYZ(kriging_df[, c("Long", "Lat", "var1.pred")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(kriging_raster, filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_suavizado.tif", format = "GTiff", overwrite = TRUE)





