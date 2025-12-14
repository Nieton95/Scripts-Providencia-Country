# VALIDACIÓN GEOESTADÍSTICA – PROFUNDIDAD DE AGUA SUBTERRÁNEA
# Polígono Providencia–Country, Guadalajara, Jalisco
# ---------------------------------------------------------------------------------

# CARGAR LIBRERÍAS NECESARIAS (instala sólo si no están)
pkgs <- c("ggplot2","sf","gstat","sp","raster","dplyr","viridis")
for(p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
lapply(pkgs, library, character.only = TRUE)

#---------------------------------------------------------------------------------
# 1) IMPORTACIÓN Y PREPARACIÓN DE DATOS
# Ajusta la ruta si es necesario
csv_path <- "C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/WATER LEVEL MAP/DB Sondeo Agua Subterránea.csv"

Water_Depth <- read.csv(csv_path, stringsAsFactors = FALSE)

# Verificar que la columna exista
if (!"Prof_agua" %in% names(Water_Depth)) stop("La columna 'Prof_agua' no se encontró en el CSV.")

# Convertir a sf con CRS WGS84 (long/lat)
Water_sf <- st_as_sf(Water_Depth, coords = c("Long","Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM (zona 13N, WGS84). Ajusta EPSG si tu zona UTM es otra.
utm_epsg <- 32613
Water_utm <- st_transform(Water_sf, crs = utm_epsg)

# Convertir a objeto Spatial (requerido por gstat)
Water_sp <- as(Water_utm, "Spatial")

cat("Puntos importados:", nrow(Water_utm), "\n")
print(st_crs(Water_utm))

#---------------------------------------------------------------------------------
# 2) VARIOGRAMA EXPERIMENTAL
variograma_emp <- variogram(Prof_agua ~ 1, Water_sp)
plot(variograma_emp, main = "Variograma experimental - Profundidad de agua")

#---------------------------------------------------------------------------------
# 3) EVALUACIÓN DE ANISOTROPÍA (variogramas direccionales)
variograma_dir <- variogram(Prof_agua ~ 1, Water_sp, alpha = c(0, 45, 90, 135))
plot(variograma_dir, main = "Variogramas direccionales – Profundidad de agua")

#---------------------------------------------------------------------------------
# 4) AJUSTE DEL VARIOGRAMA (modelo isotrópico inicial)
# Parámetros iniciales heurísticos en metros
var_value <- var(Water_sp$Prof_agua, na.rm = TRUE)
range_guess <- as.numeric(diff(range(st_coordinates(Water_utm)[,1]))) / 4
if (!is.finite(range_guess) || range_guess <= 0) range_guess <- 500  # fallback

modelo_inicial <- vgm(psill = var_value/2, model = "Sph", range = range_guess, nugget = var_value/10)

modelo_variograma <- tryCatch(
  fit.variogram(variograma_emp, modelo_inicial),
  error = function(e){
    message("fit.variogram falló, usando el modelo inicial sin ajuste: ", e$message)
    modelo_inicial
  }
)

plot(variograma_emp, modelo_variograma, main = "Ajuste del variograma – Modelo esférico isotrópico")
print(modelo_variograma)

#---------------------------------------------------------------------------------
# 5) CREACIÓN DE LA MALLA EN UTM (metros)
bb <- st_bbox(Water_utm)
xmin <- as.numeric(bb["xmin"]); ymin <- as.numeric(bb["ymin"])
xmax <- as.numeric(bb["xmax"]); ymax <- as.numeric(bb["ymax"])

expand_m <- 100  # buffer en metros alrededor de la caja (ajustable)
bbox_exp <- c(xmin = xmin - expand_m, ymin = ymin - expand_m, xmax = xmax + expand_m, ymax = ymax + expand_m)

# resolución en metros (ej: 10 m). Ajusta si la malla debe ser más fina/gruesa.
res_m <- 10

# secuencias
x_seq <- seq(bbox_exp["xmin"], bbox_exp["xmax"], by = res_m)
y_seq <- seq(bbox_exp["ymin"], bbox_exp["ymax"], by = res_m)

# proteger contra vectores vacíos
if (length(x_seq) < 2 || length(y_seq) < 2) stop("La malla generada es demasiado pequeña. Revisa bbox y res_m.")

grid_pts <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_pts) <- ~x + y
gridded(grid_pts) <- TRUE
proj4string(grid_pts) <- CRS(paste0("+init=epsg:", utm_epsg))

#---------------------------------------------------------------------------------
# 6. KRIGING – MAPA DE INCERTIDUMBRE
kriging_iso <- krige(Prof_agua ~ 1, Water_sp, grid_pts, model = modelo_variograma, nmax = 5)

krig_iso_df <- as.data.frame(kriging_iso)
colnames(krig_iso_df)[1:2] <- c("X", "Y")

# Visualización preliminar (varianza)
p_var_iso <- ggplot(krig_iso_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(title = "Varianza de kriging (isotrópico) - Profundidad de agua", x = "Easting (m)", y = "Northing (m)") +
  coord_equal() + theme_minimal()
print(p_var_iso)

#---------------------------------------------------------------------------------
# 7) VALIDACIÓN CRUZADA (LEAVE-ONE-OUT) - MODELO ISOTRÓPICO
cv_iso <- tryCatch(
  krige.cv(Prof_agua ~ 1, Water_sp, model = modelo_variograma, nfold = nrow(Water_sp)),
  error = function(e){
    message("krige.cv falló (isotrópico): ", e$message)
    return(NULL)
  }
)

if (!is.null(cv_iso)) {
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
} else {
  cat("Validación cruzada isotrópica no se ejecutó.\n")
}

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
                            #cutoff = cutoff_val,
                            #width = width_val)
)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")

