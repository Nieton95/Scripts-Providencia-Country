# VALIDACIÓN GEOESTADÍSTICA – RQD (Zonificación y Mapa RQD)
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
RQD <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/RQD MAP/RQD Map.csv",
                stringsAsFactors = FALSE)

# Convertir a sf (long/lat WGS84)
RQD_sf <- st_as_sf(RQD, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM zona 13N (EPSG:32613)
RQD_utm <- st_transform(RQD_sf, crs = 32613)

# Convertir a Spatial
RQD_utm_sp <- as(RQD_utm, "Spatial")

cat("Puntos cargados:", nrow(RQD_utm), "\n")
print(st_crs(RQD_utm))

#---------------------------------------------------------------------------------
# 2. VARIOGRAMA EXPERIMENTAL
variograma_emp_RQD <- variogram(RQD ~ 1, RQD_utm_sp)
plot(variograma_emp_RQD, main = "Variograma experimental – RQD")

#---------------------------------------------------------------------------------
# 3. EVALUACIÓN DE ANISOTROPÍA
variograma_dir_RQD <- variogram(RQD ~ 1, RQD_utm_sp,
                                alpha = c(0, 45, 90, 135))
plot(variograma_dir_RQD,
     main = "Variogramas direccionales – RQD")

#---------------------------------------------------------------------------------
# 4. AJUSTE DEL MODELO ISOTRÓPICO
modelo_inicial_RQD <- vgm(
  psill  = var(RQD_utm_sp$RQD, na.rm = TRUE)/2,
  model  = "Sph",
  range  = diff(range(st_coordinates(RQD_utm)[,1]))/4,
  nugget = var(RQD_utm_sp$RQD, na.rm = TRUE)/10
)

modelo_RQD <- fit.variogram(variograma_emp_RQD, modelo_inicial_RQD)

plot(variograma_emp_RQD, modelo_RQD,
     main = "Ajuste del variograma – Modelo esférico isotrópico")
print(modelo_RQD)

#---------------------------------------------------------------------------------
# 5. GENERACIÓN DE MALLA EN UTM (metros)
bb_RQD <- st_bbox(RQD_utm)

xmin <- as.numeric(bb_RQD["xmin"])
ymin <- as.numeric(bb_RQD["ymin"])
xmax <- as.numeric(bb_RQD["xmax"])
ymax <- as.numeric(bb_RQD["ymax"])

expand_factor <- 100  # metros
bbox_exp_RQD <- c(
  xmin = xmin - expand_factor,
  ymin = ymin - expand_factor,
  xmax = xmax + expand_factor,
  ymax = ymax + expand_factor
)

res_m <- 10  # resolución del grid

x_seq <- seq(bbox_exp_RQD["xmin"], bbox_exp_RQD["xmax"], by = res_m)
y_seq <- seq(bbox_exp_RQD["ymin"], bbox_exp_RQD["ymax"], by = res_m)

grid_RQD <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_RQD) <- ~x + y
gridded(grid_RQD) <- TRUE
proj4string(grid_RQD) <- CRS("+init=epsg:32613")

#---------------------------------------------------------------------------------
# 6. KRIGING – MAPA DE INCERTIDUMBRE
kriging_RQD <- krige(
  RQD ~ 1,
  RQD_utm_sp,
  grid_RQD,
  model = modelo_RQD,
  nmax = 5
)

RQD_uncert_df <- as.data.frame(kriging_RQD)
colnames(RQD_uncert_df)[1:2] <- c("X","Y")

# Visualización preliminar (varianza)
ggplot(RQD_uncert_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(title = "Varianza de kriging – RQD (modelo isotrópico)",
       x = "Easting (m)", y = "Northing (m)") +
  coord_equal() + theme_minimal()

#---------------------------------------------------------------------------------
# 7) VALIDACIÓN CRUZADA (LEAVE-ONE-OUT) - MODELO ISOTRÓPICO
cv_RQD <- krige.cv(
  RQD ~ 1,
  RQD_utm_sp,
  model = modelo_RQD,
  nfold = nrow(RQD_utm_sp)
)

cv_RQD$residual <- cv_RQD$observed - cv_RQD$var1.pred

ME_RQD  <- mean(cv_RQD$residual, na.rm = TRUE)
MSE_RQD <- mean(cv_RQD$residual^2, na.rm = TRUE)
MAE_RQD <- mean(abs(cv_RQD$residual), na.rm = TRUE)
COR_RQD <- cor(cv_RQD$observed, cv_RQD$var1.pred, use = "complete.obs")

cat("\n--- VALIDACIÓN RQD (MODELO ISOTRÓPICO) ---\n")
cat("ME:  ", ME_RQD,  "\n")
cat("MSE: ", MSE_RQD, "\n")
cat("MAE: ", MAE_RQD, "\n")
cat("COR: ", COR_RQD, "\n")
cat("-------------------------------------------\n\n")

#--------------------------------------------------------------------------------
# 8. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#Convertir objeto de validación cruzada a data frame
res_df <- as.data.frame(cv_RQD)

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
coords_res <- coordinates(cv_RQD)

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
residual_sp <- cv_RQD
residual_sp$residual <- res_df$residual

variograma_res <- variogram(residual ~ 1, residual_sp,
                            #cutoff = cutoff_val,
                            #width = width_val)
)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")



