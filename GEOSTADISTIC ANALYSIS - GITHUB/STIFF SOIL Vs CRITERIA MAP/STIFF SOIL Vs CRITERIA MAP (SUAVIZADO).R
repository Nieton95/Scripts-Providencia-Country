#MAPA DE PROFUNDIDAD DEL ESTRATO DE SUELO DENSO (N>50) DENTRO DEL POLIGONO PROVIDENCIA-COUNTRY EN GUADALAJARA, JALISCO, CON BASE
#EN LA UBICACIÓN DE LOS SONDEOS GEOSÍSMICOS TIPO DOWN-HOLE REALIZADOS

#-----------------------------------------------------------------------------------------------------------
# CARGAR LAS LIBRERÍAS NECESARIAS
install.packages(c("ggplot2", "sf", "gstat", "sp", "raster"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)

#-----------------------------------------------------------------------------------------------------------
# IMPORTAR LA BASE DE DATOS
Soil_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/STIFF SOIL Vs CRITERIA MAP/Vs Stiff Soil Depth.csv")

# Crear un objeto SpatialPointsDataFrame a partir de los datos
coordinates(Soil_Depth) <- ~Long + Lat

#-----------------------------------------------------------------------------------------------------------
#KRIGING
# Definir el modelo de variograma
variograma <- variogram(Prof_stiff ~ 1, Soil_Depth)
plot(variograma, main = "Variograma empírico - Profundidad suelo denso")

# Definir modelo inicial estable
modelo_inicial <- vgm(
  psill = var(Soil_Depth$Prof_stiff, na.rm = TRUE) / 2,
  model = "Sph",
  range = 0.005,
  nugget = var(Soil_Depth$Prof_stiff, na.rm = TRUE) / 10
)

# Ajustar modelo con control de errores
modelo_variograma <- tryCatch({
  fit.variogram(variograma, modelo_inicial)
}, error = function(e) {
  message("Ajuste automático del variograma falló. Usando modelo inicial por defecto.")
  modelo_inicial
})

# Mostrar y graficar el modelo ajustado
print(modelo_variograma)
plot(variograma, modelo_variograma, main = "Ajuste del variograma")

# Expandir los límites de la malla de interpolación
expand_factor <- 0.005  # Ajusta este valor para definir cuánto se extiende la interpolación más allá de los sondeos
bbox <- st_bbox(Soil_Depth)  # Obtener los límites actuales de los datos
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)  # Expansión del área de interpolación

# Definir resolución de la malla (ajusta según la densidad deseada)
res <- 0.0001  # ~10 m de resolución en GDL aprox.
x_seq <- seq(bbox_expanded["xmin"], bbox_expanded["xmax"], by = res)
y_seq <- seq(bbox_expanded["ymin"], bbox_expanded["ymax"], by = res)

grid_df <- expand.grid(Long = x_seq, Lat = y_seq)
coordinates(grid_df) <- ~Long + Lat
gridded(grid_df) <- TRUE

# Realizar la interpolación mediante kriging
kriging_result <- krige(Prof_stiff ~ 1, Soil_Depth, grid_df, model = modelo_variograma, nmax=3)

# Convertir los resultados a un data frame para ggplot2
kriging_df <- as.data.frame(kriging_result)
colnames(kriging_df) [1:2] <- c("Long", "Lat")  # Renombrar columnas para ggplot2

#----------------------------------------------------------------------------------------------------------
# VISUALIZACIÓN EN R
ggplot() +
  geom_tile(data = kriging_df, aes(x = Long, y = Lat, fill = var1.pred)) +
  scale_fill_viridis_c(option = "C", name = "Profundidad (m)", direction = -1) +
  labs(
    title = "Mapa de profundidad a estrato de suelo denso (Kriging) con base en Vs",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#-----------------------------------------------------------------------------------------------------------
#EXPORTAR LA INTERPOLACIÓN A QGIS
# Convertir los resultados de kriging a un raster (usa var1.pred)
kriging_raster <- rasterFromXYZ(kriging_df[, c("Long", "Lat", "var1.pred")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(
  kriging_raster, 
  filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_output_SUAVIZADO.tif", 
  format = "GTiff", overwrite = TRUE
  )


