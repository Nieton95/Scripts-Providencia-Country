# VALIDACIÓN GEOESTADÍSTICA – PROFUNDIDAD DE ROCA BASAL
# Subsuelo del polígono Providencia–Country, Guadalajara, Jalisco
#---------------------------------------------------------------------------------

# CARGAR LIBRERÍAS NECESARIAS (instala sólo si no están)
pkgs <- c("ggplot2","sf","gstat","sp","raster","dplyr","viridis")
for(p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)
library(dplyr)
library(viridis)

#---------------------------------------------------------------------------------
# 1) IMPORTACIÓN Y PREPARACIÓN DE DATOS
# Ajusta la ruta si es necesario
Rock_Depth <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/ROCK DEPTH MAP/BD Bedrock Depth.csv",
                       stringsAsFactors = FALSE)

# Convertir a sf con CRS WGS84 (long/lat)
Rock_Depth_sf <- st_as_sf(Rock_Depth, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM (zona 13N, WGS84). EPSG:32613 es UTM zone 13N (metros).
# Ajusta el EPSG si tu zona corresponde a otra UTM.
Rock_Depth_utm <- st_transform(Rock_Depth_sf, crs = 32613)

# Convertir a Spatial (objeto requerido por gstat)
Rock_Depth_utm_sp <- as(Rock_Depth_utm, "Spatial")

# Chequeo rápido
cat("Número de puntos importados:", nrow(Rock_Depth_utm), "\n")
print(st_crs(Rock_Depth_utm))

#---------------------------------------------------------------------------------
# 2. VARIOGRAMA EXPERIMENTAL
variograma_emp <- variogram(Prof_Inicial ~ 1, Rock_Depth_utm_sp)
plot(variograma_emp, main = "Variograma experimental – Profundidad a roca basal")

#---------------------------------------------------------------------------------
# 3. EVALUACIÓN DE ANISOTROPÍA
variograma_dir <- variogram(Prof_Inicial ~ 1, Rock_Depth_utm_sp, alpha = c(0, 45, 90, 135))
plot(variograma_dir, main = "Variogramas direccionales - Profundidad a roca basal")

#---------------------------------------------------------------------------------
# 4. AJUSTE DEL MODELO ISOTRÓPICO (como punto de partida)
# Inicializamos con parámetros aproximados en unidades de metros
modelo_inicial <- vgm(psill = var(Rock_Depth_utm_sp$Prof_Inicial, na.rm = TRUE)/2,
                      model = "Sph",
                      range = as.numeric(diff(range(st_coordinates(Rock_Depth_utm)[,1])))/4, # heurístico
                      nugget = var(Rock_Depth_utm_sp$Prof_Inicial, na.rm = TRUE)/10)

modelo_variograma <- fit.variogram(variograma_emp, modelo_inicial)

plot(variograma_emp, modelo_variograma,
     main = "Ajuste del variograma – Modelo esférico isotrópico")
print(modelo_variograma)

#---------------------------------------------------------------------------------
# 5. GENERACIÓN DE MALLA EN UTM (metros)
bb <- st_bbox(Rock_Depth_utm)

# Convertir a numéricos explícitamente
xmin <- as.numeric(bb["xmin"])
ymin <- as.numeric(bb["ymin"])
xmax <- as.numeric(bb["xmax"])
ymax <- as.numeric(bb["ymax"])

expand_factor <- 100  # metros

bbox_exp <- c(
  xmin = xmin - expand_factor,
  ymin = ymin - expand_factor,
  xmax = xmax + expand_factor,
  ymax = ymax + expand_factor
)

# Definir resolución
res_m <- 10

x_seq <- seq(bbox_exp["xmin"], bbox_exp["xmax"], by = res_m)
y_seq <- seq(bbox_exp["ymin"], bbox_exp["ymax"], by = res_m)


grid_pts <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_pts) <- ~x + y
gridded(grid_pts) <- TRUE
proj4string(grid_pts) <- CRS("+init=epsg:32613")

#---------------------------------------------------------------------------------
# 6. KRIGING – INCERTIDUMBRE (MODELO ISOTRÓPICO)
kriging_uncert_iso <- krige(Prof_Inicial ~ 1,
                            Rock_Depth_utm_sp,
                            grid_pts,
                            model = modelo_variograma,
                            nmax = 5)

uncert_iso_df <- as.data.frame(kriging_uncert_iso)
colnames(uncert_iso_df)[1:2] <- c("X","Y")

# Visualización preliminar (varianza)
ggplot(uncert_iso_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(title = "Mapa de incertidumbre (varianza) – Modelo isotrópico (UTM)",
       x = "Easting (m)", y = "Northing (m)") +
  coord_equal() +
  theme_minimal()

#---------------------------------------------------------------------------------
# 7) VALIDACIÓN CRUZADA (LEAVE-ONE-OUT) - MODELO ISOTRÓPICO
cv_iso <- krige.cv(Prof_Inicial ~ 1, Rock_Depth_utm_sp, model = modelo_variograma, nfold = nrow(Rock_Depth_utm_sp))
cv_iso$residual <- cv_iso$observed - cv_iso$var1.pred

ME_iso  <- mean(cv_iso$residual, na.rm = TRUE)
MSE_iso <- mean(cv_iso$residual^2, na.rm = TRUE)
MAE_iso <- mean(abs(cv_iso$residual), na.rm = TRUE)
COR_iso <- cor(cv_iso$observed, cv_iso$var1.pred, use = "complete.obs")

cat("\n--- VALIDACIÓN (MODELO ISOTRÓPICO) ---\n")
cat("ME: ", ME_iso, "\n")
cat("MSE:", MSE_iso, "\n")
cat("MAE:", MAE_iso, "\n")
cat("COR:", COR_iso, "\n")

#--------------------------------------------------------------------------------
# 8. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#Convertir objeto de validación cruzada a data frame
res_df <- as.data.frame(cv_iso)

#Asegurar columnas clave
res_df$residual <- res_df$observed - res_df$var1.pred

#Estadísticos básicos de residuos
res_mean <- mean(res_df$residual, na.rm = TRUE)
res_sd <- sd(res_df$residual, na.rm = TRUE)
res_min <- min(res_df$residual, na.rm = TRUE)
res_max <- max(res_df$residual, na.rm = TRUE)

cat("\n--- ESTADÍSTICOS DE RESIDUOS ---\n")
cat("Media:", round(res_mean, 3), "\n")
cat("Desv. estándar:", round(res_sd, 3), "\n")
cat("Mínimo:", round(res_min, 3), "\n")
cat("Máximo:", round(res_max, 3), "\n")

#Histograma de residuos
p_hist_res <- ggplot(res_df, aes(x = residual)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Distribución de residuos del modelo de kriging",
    x = "Residuo (m)",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p_hist_res)

#QQ-plot de residuos (normalidad aproximada)
p_qq_res <- ggplot(res_df, aes(sample = residual)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red", linewidth = 0.8) +
  labs(
    title = "QQ-plot de residuos del modelo de kriging",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p_qq_res)

#Residuos vs valores estimados (heterocedasticidad)
p_res_pred <- ggplot(res_df, aes(x = var1.pred, y = residual)) +
  geom_point(color = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuos vs valores estimados",
    x = "Profundidad estimada (m)",
    y = "Residuo (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p_res_pred)

#Mapa espacial de residuos (detección de patrones)
#Extraer coordenadas UTM
coords_res <- coordinates(cv_iso)

res_map_df <- data.frame(
  X = coords_res[,1],
  Y = coords_res[,2],
  residual = res_df$residual
)

p_map_res <- ggplot(res_map_df) +
  geom_point(aes(x = X, y = Y, color = residual), size = 3) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Residuo (m)"
  ) +
  labs(
    title = "Distribución espacial de los residuos",
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p_map_res)

#Variograma de residuos (correlación espacial residual)
residual_sp <- cv_iso
residual_sp$residual <- res_df$residual

variograma_res <- variogram(residual ~ 1, residual_sp,
                            cutoff = cutoff_val,
                            width = width_val)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")


