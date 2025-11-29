#MAPA DE LA PROFUNDIDAD DE LA ROCA BASAL DEL POLIGONO PROVIDENCIA-COUNTRY EN GUADALAJARA, JALISCO, CON BASE
#EN LA UBICACIÓN DE LOS SONDEOS REALIZADOS

# Importar la Base de Datos
Rock_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/ROCK DEPTH MAP/BD Bedrock Depth.csv")
View(Rock_Depth)

# Cargar librerías necesarias
install.packages(c("ggplot2", "sf", "gstat", "sp", "raster", "osmdata", "cowplot"))
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)
library(osmdata)
library(cowplot)

# Crear un objeto SpatialPointsDataFrame a partir de los datos
coordinates(Rock_Depth) <- ~Long + Lat

# Definir el modelo de variograma
variograma <- variogram(Prof_Inicial ~ 1, Rock_Depth)
# modelo_variograma <- fit.variogram(variograma, vgm(psill = 1, model = "Sph", range = 0.01, nugget = 0.1))

modelo_variograma <- fit.variogram(variograma, vgm(model = "Sph"))

# Expandir los límites de la malla de interpolación
expand_factor <- 0.001  # Ajusta este valor para definir cuánto se extiende la interpolación más allá de los sondeos
bbox <- st_bbox(Rock_Depth)  # Obtener los límites actuales de los datos
bbox_expanded <- bbox + c(-expand_factor, -expand_factor, expand_factor, expand_factor)  # Expansión del área de interpolación

# Generar una nueva malla de puntos en el área expandida
grid_expanded <- spsample(as_Spatial(st_as_sfc(bbox_expanded)), type = "regular", n = 1000000)  # Aumentar el número de puntos en la malla para mejor cobertura
gridded(grid_expanded) <- TRUE  # Convertir la malla en un objeto gridded adecuado para interpolación

# Realizar la interpolación mediante kriging
kriging_result_expanded <- krige(Prof_Inicial ~ 1, Rock_Depth, grid_expanded, model = modelo_variograma, nmax = 1)  # Aplicar kriging con nmax=10 para una extrapolación más suave

# Convertir los resultados a un data frame para ggplot2
kriging_df_expanded <- as.data.frame(kriging_result_expanded)
colnames(kriging_df_expanded) <- c("Long", "Lat", "Prof_Inicial")  # Renombrar columnas para ggplot2
kriging_df_expanded <- na.omit(kriging_df_expanded)


# Redondear valores para mejorar la visualización
kriging_df_expanded$Prof_Inicial <- round(kriging_df_expanded$Prof_Inicial, 2)

Rock_Depth <- as.data.frame(Rock_Depth)


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# Obtener datos de OpenStreetMap
placebb <- getbb("Guadalajara")
calles <- placebb %>% opq() %>% add_osm_feature(key = "highway") %>% osmdata_sf()
kriging_sf <- st_as_sf(kriging_df_expanded, coords = c("Long", "Lat"), crs = 4326)
kriging_extent <- st_bbox(kriging_sf)
kriging_polygon <- st_as_sfc(kriging_extent)
calles_recortadas <- st_intersection(calles$osm_lines, kriging_polygon)

# Crear el mapa
mapa_base <- ggplot() +
  geom_tile(data = kriging_df_expanded, aes(x = Long, y = Lat, fill = Prof_Inicial), alpha = 0.8) +
  geom_sf(data = calles_recortadas, color = "black", size = 0.5) +
  geom_point(data = Rock_Depth, aes(x = Long, y = Lat), color = "red", size = 3, shape = 21, fill = "white") +
  scale_fill_gradientn(name = "Profundidad de la Roca (m)",
                       colors = c("lightgrey", "#C0C0C0", "#696969", "darkslategray"),
                       breaks = round(seq(min(kriging_df_expanded$Prof_Inicial), max(kriging_df_expanded$Prof_Inicial), length.out = 10), 2)
  ) +
  coord_sf() +
  labs(title = "Profundidad de la Roca Basal",
       subtitle = "Polígono Providencia-Country, Guadalajara",
       x = "Longitud", y = "Latitud") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Crear anotaciones
tabla_anotaciones <- ggplot() +
  annotate("text", x = 1.5, y = 0.95, label = "Ubicación de los sondeos exploratorios:", color = "black", size = 4, hjust = 0) +
  annotate("point", x = 1.45, y = 0.95, color = "red", size = 4, shape = 21, fill = "white") +
  xlim(0, 2) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Combinar mapa y anotaciones
mapa_final <- plot_grid(mapa_base, tabla_anotaciones, ncol = 1, rel_heights = c(0.90, 0.10))
print(mapa_final)


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#EXPORTAR LA INTERPOLACIÓN A QGIS

# Convertir los resultados de kriging a un raster
kriging_raster <- rasterFromXYZ(kriging_df_expanded[, c("Long", "Lat", "Prof_Inicial")])

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(kriging_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(kriging_raster, filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/kriging_output.tif", format = "GTiff", overwrite = TRUE)


