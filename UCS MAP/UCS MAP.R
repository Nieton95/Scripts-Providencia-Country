# MAPA DE LA ZONIFICACIÓN GEOLÓGICA CON BASE EN LOS VALORES DEL UCS DE LA ROCA BASAL DEL POLÍGONO PROVIDENCIA-COUNTRY EN 
#GUADALAJARA, JALISCO CON BASE EN EL CHART DE RQD INDICADO EN LA PÁGINA 246 DEL LIBRO DE GONZÁLEZ DE VALLEJO.

#-----------------------------------------------------------------------------------------------------------
# CARGAR LAS LIBRERÍAS NECESARIAS
install.packages(c("ggplot2", "sf", "gstat", "sp", "raster", "rgdal"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)
library(rgdal)

#-----------------------------------------------------------------------------------------------------------
# IMPORTAR LA BASE DE DATOS
UCS_var <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/UCS MAP/UCS Map.csv")

# Convertir a objeto espacial
coordinates(UCS_var) <- ~Long + Lat

#-----------------------------------------------------------------------------------------------------------
#KRIGING
# Definir el modelo de variograma
variograma <- variogram(UCS ~ 1, UCS_var)

# Ajustar el modelo de variograma 
modelo_variograma <- fit.variogram(
  variograma,
  vgm(psill = 400, model = "Sph", range = 0.01, nugget = 50)
)

# Crear una cuadrícula regular para la interpolación
expand_factor <- 0.001
bbox <- st_bbox(st_as_sf(UCS_var))
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)

# Definir resolución de la malla (ajusta según la densidad deseada)
res <- 0.0001  # ~10 m de resolución en GDL aprox.
x_seq <- seq(bbox_expanded["xmin"], bbox_expanded["xmax"], by = res)
y_seq <- seq(bbox_expanded["ymin"], bbox_expanded["ymax"], by = res)

grid_df <- expand.grid(Long = x_seq, Lat = y_seq)
coordinates(grid_df) <- ~Long + Lat
gridded(grid_df) <- TRUE

# Realizar interpolación kriging
kriging_result <- krige(UCS ~ 1, UCS_var, grid_df, model = modelo_variograma, nmax = 10)
kriging_df <- as.data.frame(kriging_result)
colnames(kriging_df)[1:2] <- c("Long", "Lat")
kriging_df$Zona <- sapply(kriging_df$var1.pred, define_zone_numeric)

#-----------------------------------------------------------------------------------------------------------
# VISUALIZACIÓN EN R
ggplot() +
  geom_tile(data = kriging_df, aes(x = Long, y = Lat, fill = var1.pred)) +
  scale_fill_viridis_c(option = "C", name = "UCS (kg/cm2)", direction = -1) +
  labs(
    title = "Mapa de distribución del valor de UCS",
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

# Definir el sistema de coordenadas (EPSG:4326)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar el GeoTIFF con valores interpolados
writeRaster(
  kriging_raster,
  filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_completo.tif",
  format = "GTiff",
  overwrite = TRUE
)
















