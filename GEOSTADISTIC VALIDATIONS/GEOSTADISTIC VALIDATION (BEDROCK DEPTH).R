# VALIDACIÓN GEOESTADÍSTICA – PROFUNDIDAD DE ROCA BASAL (SCRIPT CORREGIDO)
# Subsuelo del polígono Providencia–Country, Guadalajara, Jalisco
# --------------------------------------------------------------------------------

#Paquetes (instalación condicional)
pkgs <- c("ggplot2","sf","gstat","sp","raster","dplyr","viridis")
for(p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(sf)
library(gstat)
library(sp)
library(raster)
library(dplyr)
library(viridis)

# --------------------------------------------------------------------------------
# 1. IMPORTACIÓN Y PREPARACIÓN DE DATOS (CSV -> sf -> UTM)
# --------------------------------------------------------------------------------
Rock_Depth <- read.csv(
  "C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/ROCK DEPTH MAP/BD Bedrock Depth.csv",
  stringsAsFactors = FALSE
)

# Convertir a sf con CRS WGS84 (long/lat)
Rock_Depth_sf <- st_as_sf(Rock_Depth, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM (zona 13N — EPSG:32613). Ajusta EPSG si tu UTM es distinto.
Rock_Depth_utm <- st_transform(Rock_Depth_sf, crs = 32613)

# Convertir a Spatial (objeto requerido por gstat)
Rock_Depth_utm_sp <- as(Rock_Depth_utm, "Spatial")

cat("Puntos cargados:", nrow(Rock_Depth_utm), "\n")
print(st_crs(Rock_Depth_utm))


# --------------------------------------------------------------------------------
# 2. ANÁLISIS DE NORMALIDAD (VARIABLE ORIGINAL)
# --------------------------------------------------------------------------------

depth_vals <- Rock_Depth_utm_sp$Prof_Inicial

# Histograma
p_hist_norm <- ggplot(data.frame(depth = depth_vals), aes(x = depth)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    title = "Histograma – Profundidad a roca basal",
    x = "Profundidad (m)",
    y = "Frecuencia"
  ) +
  theme_minimal()

print(p_hist_norm)

# QQ-plot
p_qq_norm <- ggplot(data.frame(depth = depth_vals), aes(sample = depth)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ-plot – Profundidad a roca basal",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados"
  ) +
  theme_minimal()

print(p_qq_norm)

# Prueba Shapiro-Wilk
shapiro_test <- shapiro.test(depth_vals)

cat("\n--- PRUEBA DE NORMALIDAD (Shapiro-Wilk) ---\n")
print(shapiro_test)


# --------------------------------------------------------------------------------
# 3. VARIOGRAMA EXPERIMENTAL (EN UTM - metros) — bins controlados (cutoff / width)
# --------------------------------------------------------------------------------
# Calculamos una distancia de corte razonable: diagonal del bounding box / 2
bb_num <- as.numeric(st_bbox(Rock_Depth_utm))
names(bb_num) <- c("xmin","ymin","xmax","ymax")
diag_m <- sqrt((bb_num["xmax"]-bb_num["xmin"])^2 + (bb_num["ymax"]-bb_num["ymin"])^2)
cutoff_val <- diag_m / 2          # evaluar hasta la mitad de la diagonal del área
width_val  <- cutoff_val / 12     # dividir en ~12 bins (heurístico)

# Asegurar valores sensatos
if(!is.finite(cutoff_val) || cutoff_val <= 0) cutoff_val <- 1000
if(!is.finite(width_val)  || width_val <= 0)  width_val <- cutoff_val/12

variograma_emp <- variogram(Prof_Inicial ~ 1, Rock_Depth_utm_sp,
                            cutoff = cutoff_val, width = width_val)

plot(variograma_emp, main = "Variograma experimental – Profundidad a roca basal")
cat("cutoff (m):", round(cutoff_val,1), " width (m):", round(width_val,1), "\n")

# --------------------------------------------------------------------------------
# 4. EVALUACIÓN DE ANISOTROPÍA (direccionales en UTM)
# --------------------------------------------------------------------------------
variograma_dir <- variogram(Prof_Inicial ~ 1, Rock_Depth_utm_sp,
                            alpha = c(0, 45, 90, 135),
                            cutoff = cutoff_val, width = width_val)
plot(variograma_dir, main = "Variogramas direccionales - evaluación de anisotropía")

# --------------------------------------------------------------------------------
# 5. AJUSTE DEL MODELO ISOTRÓPICO (REAJUSTADO PARA UTM)
# --------------------------------------------------------------------------------
# Parámetros iniciales razonables en metros (heurísticos)
var_data <- var(Rock_Depth_utm_sp$Prof_Inicial, na.rm = TRUE)
psill_init  <- max(var_data/2, 1e-6)
nugget_init <- max(var_data/10, 0)
# rango inicial: 1/4 del mayor span del bbox
span_x <- abs(bb_num["xmax"] - bb_num["xmin"])
span_y <- abs(bb_num["ymax"] - bb_num["ymin"])
range_init  <- max(min(span_x, span_y)/4, cutoff_val/6, 50)  # m, heurístico

modelo_inicial <- vgm(psill = psill_init, model = "Sph", range = range_init, nugget = nugget_init)

# Intentamos ajuste con control de errores
modelo_variograma <- tryCatch({
  fit.variogram(variograma_emp, modelo_inicial)
}, error = function(e){
  message("fit.variogram falló — se usará modelo inicial como fallback.")
  return(modelo_inicial)
})

plot(variograma_emp, modelo_variograma, main = "Ajuste variograma")
print(modelo_variograma)

# --------------------------------------------------------------------------------
# 6. GENERACIÓN DE MALLA EN UTM (metros) — con conversión numérica segura
# --------------------------------------------------------------------------------
expand_factor <- 100    # expansión en metros (ajustable)
xmin <- as.numeric(bb_num["xmin"]) - expand_factor
xmax <- as.numeric(bb_num["xmax"]) + expand_factor
ymin <- as.numeric(bb_num["ymin"]) - expand_factor
ymax <- as.numeric(bb_num["ymax"]) + expand_factor

# Definir resolución en metros (ej: 10 m). Aumenta si la malla es muy grande.
res_m <- 10
x_seq <- seq(xmin, xmax, by = res_m)
y_seq <- seq(ymin, ymax, by = res_m)

grid_pts <- expand.grid(X = x_seq, Y = y_seq)
coordinates(grid_pts) <- ~X + Y
gridded(grid_pts) <- TRUE
proj4string(grid_pts) <- CRS("+init=epsg:32613")

cat("Malla generada:", length(x_seq) * length(y_seq), "celdas (puede tardar si es muy grande)\n")

# --------------------------------------------------------------------------------
# 7. KRIGING – INCERTIDUMBRE (MODELO ISOTRÓPICO)
# --------------------------------------------------------------------------------
kriging_uncert_iso <- krige(Prof_Inicial ~ 1,
                            Rock_Depth_utm_sp,
                            grid_pts,
                            model = modelo_variograma,
                            nmax = 5)

uncert_iso_df <- as.data.frame(kriging_uncert_iso)
colnames(uncert_iso_df)[1:2] <- c("X","Y")

# Visualización preliminar (UTM)
p_unc_iso <- ggplot(uncert_iso_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(
    title = "Mapa de incertidumbre (varianza)",
    x = "Easting (m)", 
    y = "Northing (m)"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p_unc_iso)


# --------------------------------------------------------------------------------
# 8. VALIDACIÓN CRUZADA (MODELO ISOTRÓPICO)
# --------------------------------------------------------------------------------
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
# 9. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#--------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------
# 10. ANÁLISIS DE SENSIBILIDAD (variación de nmax)
# --------------------------------------------------------------------------------

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
    Prof_Inicial ~ 1,
    Rock_Depth_utm_sp,
    model = modelo_variograma,
    nfold = nrow(Rock_Depth_utm_sp),
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

cat("\n--- ANÁLISIS DE SENSIBILIDAD (nmax) ---\n")
print(sens_results)

# Visualización del efecto de nmax
p_sens <- ggplot(sens_results, aes(x = nmax)) +
  geom_line(aes(y = RMSE), linewidth = 1) +
  geom_point(aes(y = RMSE), size = 3) +
  labs(
    title = "Análisis de sensibilidad – RMSE vs nmax",
    x = "Número máximo de vecinos (nmax)",
    y = "RMSE"
  ) +
  theme_minimal()

print(p_sens)

