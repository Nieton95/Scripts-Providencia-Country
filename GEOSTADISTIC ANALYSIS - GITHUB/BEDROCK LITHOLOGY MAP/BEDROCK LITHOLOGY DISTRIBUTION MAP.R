#MAPA DE LA DISTRIBUCIÓN DE LA LITOLOGÍA IDENTIFICADA DENTRO DEL POLÍGONO DE ESTUDIO, MEDIANTE LA TÉCNICA IDW

#El método IDW (Inverse Distance Weighting) es una técnica de interpolación geoespacial que asigna un valor estimado a un punto no muestreado 
#en función de los valores de los puntos cercanos. La idea principal es que los puntos más cercanos al punto de estimación tienen más influencia 
#en el valor final, y esta influencia disminuye a medida que aumenta la distancia. Los valores de los puntos cercanos se ponderan inversamente 
#a su distancia, lo que significa que los puntos más cercanos tienen un peso mayor en la estimación. Es un método simple y rápido, pero no 
#considera la estructura espacial ni la autocorrelación de los datos.

#-----------------------------------------------------------------------------------------------------------
# CARGAR LIBRERÍAS NECESARIAS
install.packages(c("sf", "gstat", "sp", "ggplot2", "viridis", "osmdata", "cowplot", "raster"))
library(sf)
library(gstat)
library(sp)
library(ggplot2)
library(viridis)
library(osmdata)
library(cowplot)
library(raster)

#-----------------------------------------------------------------------------------------------------------
# IMPORTAR LA BASE DE DATOS
df <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/BEDROCK LITHOLOGY MAP/BD Bedrock Type.csv", 
               header = TRUE, stringsAsFactors = FALSE)

# Convertir a objeto espacial (CRS 4326)
df_sf <- st_as_sf(df, coords = c("Long", "Lat"), crs = 4326)

#-----------------------------------------------------------------------------------------------------------
#KRIGING
# Crear variables binarias para la interpolación
df_sf$Basalto <- ifelse(df$Tipo_roca == "Basalto", 1, 0)
df_sf$Ignimbrita <- ifelse(df$Tipo_roca == "Ignimbrita", 1, 0)

# Definir un margen de extensión (por ejemplo, 10% del rango)
margen_x <- (max(df$Long) - min(df$Long)) * 0.1
margen_y <- (max(df$Lat) - min(df$Lat)) * 0.1

# Crear una malla de interpolación extendida
x_range <- seq(min(df$Long) - margen_x, max(df$Long) + margen_x, length.out = 1500)
y_range <- seq(min(df$Lat) - margen_y, max(df$Lat) + margen_y, length.out = 1500)
grid <- expand.grid(Long = x_range, Lat = y_range)
grid_sf <- st_as_sf(grid, coords = c("Long", "Lat"), crs = 4326)

# Convertir a objeros spatial para la interpolación
df_sp <- as(df_sf, "Spatial")
grid_sp <- as(grid_sf, "Spatial")

# Interpolación mediante IDW para basalto e ignimbrita
idw_basalto <- idw(Basalto ~ 1, locations = df_sp, newdata = grid_sp, idp = 2)
idw_ignimbrita <- idw(Ignimbrita ~ 1, locations = df_sp, newdata = grid_sp, idp = 2)

# Combinar resultados y determinar litología dominante
kriged_df <- data.frame(Long = coordinates(grid_sp)[,1], 
                        Lat = coordinates(grid_sp)[,2], 
                        Basalto = idw_basalto$var1.pred, 
                        Ignimbrita = idw_ignimbrita$var1.pred)
kriged_df$Tipo_roca <- ifelse(kriged_df$Basalto > kriged_df$Ignimbrita, "Basalto", "Ignimbrita")

# Convertir a SF para limitar la extensión
kriging_sf <- st_as_sf(kriged_df, coords = c("Long", "Lat"), crs = 4326)

# Definir extensión del mapa
kriging_extent <- st_bbox(kriging_sf)  # Esto crea un 'bounding box' alrededor de los datos de kriging
kriging_polygon <- st_as_sfc(kriging_extent)  # Convertir el bounding box a un polígono sf


# Definir los límites espaciales basados en la extensión de kriging_sf
x_lim <- range(st_coordinates(grid_sf)[,1])  # Extensión en X
y_lim <- range(st_coordinates(grid_sf)[,2])  # Extensión en Y


#-----------------------------------------------------------------------------------------------------------
#EXPORTAR LA INTERPOLACIÓN A QGIS
# Convertir los resultados de IDW en un data.frame con coordenadas explícitas
kriging_df <- cbind(st_coordinates(kriging_sf), as.data.frame(kriging_sf))

# Asignar valores numéricos a cada litología
kriging_df$Tipo_roca_num <- ifelse(kriging_df$Tipo_roca == "Basalto", 1, 2)

# Crear el raster desde el data.frame (asegurando que las coordenadas sean correctas)
idw_raster <- rasterFromXYZ(kriging_df[, c("X", "Y", "Tipo_roca_num")])  # X = Longitud, Y = Latitud

# Definir el sistema de coordenadas (EPSG:4326 para coordenadas geográficas)
crs(idw_raster) <- CRS("+init=epsg:4326")

# Exportar como GeoTIFF
writeRaster(idw_raster, filename = "C:/Users/ivann/OneDrive/Escritorio/Thesis/idw_litologia_suavizada.tif", 
            format = "GTiff", datatype = "INT1U", overwrite = TRUE)










