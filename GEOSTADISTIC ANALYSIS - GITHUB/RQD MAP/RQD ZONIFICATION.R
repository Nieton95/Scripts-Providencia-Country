#MAPA DE LA ZONIFICACIÓN GEOLÓGICA CON BASE EN LOS VALORES DEL RQD DE LA ROCA BASAL DEL POLIGONO PROVIDENCIA-COUNTRY EN GUADALAJARA, JALISCO, 
#CON BASE EN EL CHART DE RQD INDICADO EN LA PÁGINA 326 DEL LIBRO DE GONZALEZ DE VALLEJO.

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
RQD_var <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/RQD MAP/RQD Map.csv")

# Convertir a objeto espacial
coordinates(RQD_var) <- ~Long + Lat

# Definir el modelo de variograma
variograma <- variogram(RQD ~ 1, RQD_var)

# Ajustar el modelo de variograma 
modelo_variograma <- fit.variogram(
  variograma, 
  vgm(psill = 400, model = "Sph", range = 0.01, nugget = 50)
  )

# Expandir los límites de la malla de interpolación
expand_factor <- 0.001
bbox <- st_bbox(st_as_sf(RQD_var))
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)

# Definir resolución de la malla (ajusta según la densidad deseada)
res <- 0.0001  # ~10 m de resolución en GDL aprox.
x_seq <- seq(bbox_expanded["xmin"], bbox_expanded["xmax"], by = res)
y_seq <- seq(bbox_expanded["ymin"], bbox_expanded["ymax"], by = res)

grid_df <- expand.grid(Long = x_seq, Lat = y_seq)
coordinates(grid_df) <- ~Long + Lat
gridded(grid_df) <- TRUE

# Realizar interpolación kriging
kriging_result <- krige(RQD ~ 1, RQD_var, grid_df, model = modelo_variograma, nmax = 10)
kriging_df <- as.data.frame(kriging_result)
colnames(kriging_df)[1:2] <- c("Long", "Lat")

# Clasificar en zonas geológicas
define_zone_numeric <- function(rqd) {
  if (rqd >= 90) {
    return(1)  # Zona 1 - Muy Buena Calidad
  } else if (rqd >= 75) {
    return(2)  # Zona 2 - Buena Calidad
  } else if (rqd >= 50) {
    return(3)  # Zona 3 - Calidad Media
  } else if (rqd >= 25) {
    return(4)  # Zona 4 - Mala Calidad
  } else {
    return(5)  # Zona 5 - Muy Mala Calidad
  }
}


kriging_df$Zona <- sapply(kriging_df$var1.pred, define_zone_numeric)
kriging_df$Zona <- factor(kriging_df$Zona, levels = c(1, 2, 3, 4, 5))

#-----------------------------------------------------------------------------------------------------------
# VISUALIZACIÓN EN R
ggplot() +
  geom_tile(data = kriging_df, aes(x = Long, y = Lat, fill = Zona)) +
  scale_fill_manual(
    values = c("#2E8B57", "#66C2A5", "#FDD835", "#E67E22", "#D73027"),
    breaks = c(1, 2, 3, 4, 5),
    labels = c("Muy Buena Calidad", "Buena Calidad", "Calidad Media", "Mala Calidad", "Muy Mala Calidad")
  ) +
  labs(
    title = "Zonificación Geológica según RQD", 
    fill = "Calidad de Roca"
  ) +
  theme_minimal()

#-----------------------------------------------------------------------------------------------------------
#EXPORTAR LA INTERPOLACIÓN A QGIS
# Convertir los resultados de kriging a un raster
kriging_df$Zona_num <- as.numeric(kriging_df$Zona)
zona_raster <- rasterFromXYZ(kriging_df[, c("Long", "Lat", "Zona_num")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(zona_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(
  zona_raster,
  filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_output_zonification_RQD_SUAVIZADO.tif",
  format = "GTiff",
  overwrite = TRUE
)








