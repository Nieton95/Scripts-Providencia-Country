# VALIDACIÓN GEOESTADÍSTICA – Rock Depth Vs Criteria Map
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
# 1. IMPORTACIÓN Y PREPARACIÓN DE DATOS
#---------------------------------------------------------------------------------
# Ajusta la ruta si es necesario
VsRock <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/ROCK DEPTH VS CRITERIA MAP/Vs Rock Depth.csv",
                   stringsAsFactors = FALSE)

# Convertir a sf (WGS84)
VsRock_sf <- st_as_sf(VsRock, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM zona 13N
VsRock_utm <- st_transform(VsRock_sf, crs = 32613)

# Convertir a Spatial
VsRock_utm_sp <- as(VsRock_utm, "Spatial")

cat("Puntos cargados:", nrow(VsRock_utm), "\n")
print(st_crs(VsRock_utm))

#---------------------------------------------------------------------------------
# 2. ANÁLISIS DE NORMALIDAD – PROFUNDIDAD A ROCA (Vs)
#---------------------------------------------------------------------------------
rock_vals <- VsRock_utm_sp$Prof_rock

# Histograma
p_hist_norm <- ggplot(data.frame(Prof_rock = rock_vals), aes(x = Prof_rock)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    title = "Histograma – Profundidad a roca (Vs)",
    x = "Profundidad a roca (m)",
    y = "Frecuencia"
  ) +
  theme_minimal()

print(p_hist_norm)

# QQ-plot
p_qq_norm <- ggplot(data.frame(Prof_rock = rock_vals), aes(sample = Prof_rock)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ-plot – Profundidad a roca (Vs)",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados"
  ) +
  theme_minimal()

print(p_qq_norm)

# Prueba de Shapiro-Wilk
shapiro_test <- shapiro.test(rock_vals)

cat("\n--- PRUEBA DE NORMALIDAD (Shapiro-Wilk) – PROFUNDIDAD A ROCA (Vs) ---\n")
print(shapiro_test)

#---------------------------------------------------------------------------------
# 3. VARIOGRAMA EXPERIMENTAL
#---------------------------------------------------------------------------------
variograma_emp_Vs <- variogram(Prof_rock ~ 1, VsRock_utm_sp)
plot(variograma_emp_Vs, main = "Variograma experimental – Profundidad a roca (Vs)")

#---------------------------------------------------------------------------------
# 4. EVALUACIÓN DE ANISOTROPÍA
#---------------------------------------------------------------------------------
variograma_dir_Vs <- variogram(
  Prof_rock ~ 1, VsRock_utm_sp,
  alpha = c(0, 45, 90, 135)
)

plot(variograma_dir_Vs,
     main = "Variogramas direccionales – Profundidad a roca (Vs)")

#---------------------------------------------------------------------------------
# 5. AJUSTE DEL MODELO ISOTRÓPICO
#---------------------------------------------------------------------------------
modelo_inicial_Vs <- vgm(
  psill  = var(VsRock_utm_sp$Prof_rock, na.rm = TRUE) / 2,
  model  = "Sph",
  range  = diff(range(st_coordinates(VsRock_utm)[,1])) / 4,
  nugget = var(VsRock_utm_sp$Prof_rock, na.rm = TRUE) / 10
)

modelo_Vs <- fit.variogram(variograma_emp_Vs, modelo_inicial_Vs)

plot(variograma_emp_Vs, modelo_Vs,
     main = "Ajuste del variograma – Modelo esférico isotrópico (Rock Depth)")
print(modelo_Vs)

#---------------------------------------------------------------------------------
# 6. GENERACIÓN DE MALLA EN UTM (metros)
#---------------------------------------------------------------------------------
bb_Vs <- st_bbox(VsRock_utm)

xmin <- as.numeric(bb_Vs["xmin"])
ymin <- as.numeric(bb_Vs["ymin"])
xmax <- as.numeric(bb_Vs["xmax"])
ymax <- as.numeric(bb_Vs["ymax"])

expand_factor <- 100
bbox_exp_Vs <- c(
  xmin = xmin - expand_factor,
  ymin = ymin - expand_factor,
  xmax = xmax + expand_factor,
  ymax = ymax + expand_factor
)

res_m <- 10

x_seq <- seq(bbox_exp_Vs["xmin"], bbox_exp_Vs["xmax"], by = res_m)
y_seq <- seq(bbox_exp_Vs["ymin"], bbox_exp_Vs["ymax"], by = res_m)

grid_Vs <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_Vs) <- ~x + y
gridded(grid_Vs) <- TRUE
proj4string(grid_Vs) <- CRS("+init=epsg:32613")

#---------------------------------------------------------------------------------
# 7. KRIGING – MAPA DE INCERTIDUMBRE
#---------------------------------------------------------------------------------
kriging_Vs <- krige(
  Prof_rock ~ 1,
  VsRock_utm_sp,
  grid_Vs,
  model = modelo_Vs,
  nmax = 5
)

Vs_uncert_df <- as.data.frame(kriging_Vs)
colnames(Vs_uncert_df)[1:2] <- c("X","Y")

ggplot(Vs_uncert_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(
    title = "Varianza de kriging – Profundidad a roca (Vs)",
    x = "Easting (m)", y = "Northing (m)"
  ) +
  coord_equal() + theme_minimal()

#---------------------------------------------------------------------------------
# 8. VALIDACIÓN CRUZADA (LEAVE-ONE-OUT)
#---------------------------------------------------------------------------------
cv_Vs <- krige.cv(
  Prof_rock ~ 1,
  VsRock_utm_sp,
  model = modelo_Vs,
  nfold = nrow(VsRock_utm_sp)
)

cv_Vs$residual <- cv_Vs$observed - cv_Vs$var1.pred

ME_Vs  <- mean(cv_Vs$residual, na.rm = TRUE)
MSE_Vs <- mean(cv_Vs$residual^2, na.rm = TRUE)
MAE_Vs <- mean(abs(cv_Vs$residual), na.rm = TRUE)
COR_Vs <- cor(cv_Vs$observed, cv_Vs$var1.pred, use = "complete.obs")

cat("\n--- VALIDACIÓN (ROCK DEPTH – MODELO ISOTRÓPICO) ---\n")
cat("ME:  ", ME_Vs,  "\n")
cat("MSE: ", MSE_Vs, "\n")
cat("MAE: ", MAE_Vs, "\n")
cat("COR: ", COR_Vs, "\n")
cat("-------------------------------------------\n\n")

#--------------------------------------------------------------------------------
# 9. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#--------------------------------------------------------------------------------
#Convertir objeto de validación cruzada a data frame
res_df <- as.data.frame(cv_Vs)

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
coords_res <- coordinates(cv_Vs)

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
residual_sp <- cv_Vs
residual_sp$residual <- res_df$residual

variograma_res <- variogram(residual ~ 1, residual_sp,
                            cutoff = cutoff_val,
                            width = width_val)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")

#--------------------------------------------------------------------------------
# 10. ANÁLISIS DE SENSIBILIDAD (variación de nmax) – PROFUNDIDAD A ROCA (Vs)
#--------------------------------------------------------------------------------

nmax_vals <- c(3, 5, 8, 12)

sens_results <- data.frame(
  nmax = integer(),
  ME = numeric(),
  MAE = numeric(),
  RMSE = numeric(),
  COR = numeric()
)

for (n in nmax_vals) {
  
  cv_tmp <- krige.cv(
    Prof_rock ~ 1,
    VsRock_utm_sp,
    model = modelo_Vs,
    nfold = nrow(VsRock_utm_sp),
    nmax = n
  )
  
  residuals_tmp <- cv_tmp$observed - cv_tmp$var1.pred
  
  sens_results <- rbind(
    sens_results,
    data.frame(
      nmax = n,
      ME = mean(residuals_tmp, na.rm = TRUE),
      MAE = mean(abs(residuals_tmp), na.rm = TRUE),
      RMSE = sqrt(mean(residuals_tmp^2, na.rm = TRUE)),
      COR = cor(cv_tmp$observed, cv_tmp$var1.pred, use = "complete.obs")
    )
  )
}

cat("\n--- ANÁLISIS DE SENSIBILIDAD (nmax) – PROFUNDIDAD A ROCA (Vs) ---\n")
print(sens_results)

# Visualización del efecto de nmax
p_sens <- ggplot(sens_results, aes(x = nmax)) +
  geom_line(aes(y = RMSE), linewidth = 1) +
  geom_point(aes(y = RMSE), size = 3) +
  labs(
    title = "Análisis de sensibilidad – RMSE vs nmax (Profundidad a roca Vs)",
    x = "Número máximo de vecinos (nmax)",
    y = "RMSE"
  ) +
  theme_minimal()

print(p_sens)

# AVISO CRÍTICO (PROFUNDIDAD AL BASAMENTO ROCOSO): La validación cruzada para la profundidad de detección del basamento rocoso con base en los valores
#de Vs arrojó un COR de -0.577, indicando una relación inversa y una falla severa del modelo. El alto MAE (7.41 metros) confirma un error de predicción
#críticamente elevado. Estos resultados demuestran que la topografía del basamento es altamente irregular o caótica (no estacionaria) y está 
#influenciada por estructuras geológicas no espaciales (paleo-relieve, fallas, o erosión diferencial). 
#El mapa de Kriging generado no debe utilizarse para la toma de decisiones. 
#Se sugiere aumentar la cantidad de datos.
