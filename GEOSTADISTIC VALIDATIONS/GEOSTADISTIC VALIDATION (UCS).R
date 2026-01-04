# VALIDACIÓN GEOESTADÍSTICA – UCS (Zonificación y Mapa UCS)
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
UCS <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/UCS MAP/UCS Map.csv",
                stringsAsFactors = FALSE)

# Convertir a sf (long/lat WGS84)
UCS_sf <- st_as_sf(UCS, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM zona 13N (EPSG:32613)
UCS_utm <- st_transform(UCS_sf, crs = 32613)

# Convertir a Spatial
UCS_utm_sp <- as(UCS_utm, "Spatial")

cat("Puntos cargados:", nrow(UCS_utm), "\n")
print(st_crs(UCS_utm))

#---------------------------------------------------------------------------------
# 2. ANÁLISIS DE NORMALIDAD (UCS)
#---------------------------------------------------------------------------------

ucs_vals <- UCS_utm_sp$UCS

# Histograma
p_hist_norm <- ggplot(data.frame(UCS = ucs_vals), aes(x = UCS)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    title = "Histograma – UCS",
    x = "Resistencia a compresión simple (UCS)",
    y = "Frecuencia"
  ) +
  theme_minimal()

print(p_hist_norm)

# QQ-plot
p_qq_norm <- ggplot(data.frame(UCS = ucs_vals), aes(sample = UCS)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ-plot – UCS",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados"
  ) +
  theme_minimal()

print(p_qq_norm)

# Prueba Shapiro-Wilk
shapiro_test <- shapiro.test(ucs_vals)

cat("\n--- PRUEBA DE NORMALIDAD (Shapiro-Wilk) – UCS ---\n")
print(shapiro_test)

#---------------------------------------------------------------------------------
# 3. VARIOGRAMA EXPERIMENTAL
#---------------------------------------------------------------------------------
variograma_emp_UCS <- variogram(UCS ~ 1, UCS_utm_sp)
plot(variograma_emp_UCS, main = "Variograma experimental – UCS")

#---------------------------------------------------------------------------------
# 4. EVALUACIÓN DE ANISOTROPÍA
#---------------------------------------------------------------------------------
variograma_dir_UCS <- variogram(UCS ~ 1, UCS_utm_sp,
                                alpha = c(0, 45, 90, 135))
plot(variograma_dir_UCS,
     main = "Variogramas direccionales – UCS")

#---------------------------------------------------------------------------------
# 5. AJUSTE DEL MODELO ISOTRÓPICO
#---------------------------------------------------------------------------------
modelo_inicial_UCS <- vgm(
  psill  = var(UCS_utm_sp$UCS, na.rm = TRUE)/2,
  model  = "Sph",
  range  = diff(range(st_coordinates(UCS_utm)[,1]))/4,
  nugget = var(UCS_utm_sp$UCS, na.rm = TRUE)/10
)

modelo_UCS <- fit.variogram(variograma_emp_UCS, modelo_inicial_UCS)

plot(variograma_emp_UCS, modelo_UCS,
     main = "Ajuste del variograma – Modelo esférico isotrópico (UCS)")
print(modelo_UCS)

#---------------------------------------------------------------------------------
# 6. GENERACIÓN DE MALLA EN UTM (metros)
#---------------------------------------------------------------------------------
bb_UCS <- st_bbox(UCS_utm)

xmin <- as.numeric(bb_UCS["xmin"])
ymin <- as.numeric(bb_UCS["ymin"])
xmax <- as.numeric(bb_UCS["xmax"])
ymax <- as.numeric(bb_UCS["ymax"])

expand_factor <- 100  # metros
bbox_exp_UCS <- c(
  xmin = xmin - expand_factor,
  ymin = ymin - expand_factor,
  xmax = xmax + expand_factor,
  ymax = ymax + expand_factor
)

res_m <- 10  # resolución del grid

x_seq <- seq(bbox_exp_UCS["xmin"], bbox_exp_UCS["xmax"], by = res_m)
y_seq <- seq(bbox_exp_UCS["ymin"], bbox_exp_UCS["ymax"], by = res_m)

grid_UCS <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_UCS) <- ~x + y
gridded(grid_UCS) <- TRUE
proj4string(grid_UCS) <- CRS("+init=epsg:32613")

#---------------------------------------------------------------------------------
# 7. KRIGING – MAPA DE INCERTIDUMBRE
#---------------------------------------------------------------------------------
kriging_UCS <- krige(
  UCS ~ 1,
  UCS_utm_sp,
  grid_UCS,
  model = modelo_UCS,
  nmax = 5
)

UCS_uncert_df <- as.data.frame(kriging_UCS)
colnames(UCS_uncert_df)[1:2] <- c("X","Y")

# Visualización preliminar
ggplot(UCS_uncert_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(title = "Varianza de kriging – UCS (modelo isotrópico)",
       x = "Easting (m)", y = "Northing (m)") +
  coord_equal() + theme_minimal()

#---------------------------------------------------------------------------------
# 8. VALIDACIÓN CRUZADA (LEAVE-ONE-OUT) - MODELO ISOTRÓPICO
#---------------------------------------------------------------------------------
cv_UCS <- krige.cv(
  UCS ~ 1,
  UCS_utm_sp,
  model = modelo_UCS,
  nfold = nrow(UCS_utm_sp)
)

cv_UCS$residual <- cv_UCS$observed - cv_UCS$var1.pred

ME_UCS  <- mean(cv_UCS$residual, na.rm = TRUE)
MSE_UCS <- mean(cv_UCS$residual^2, na.rm = TRUE)
MAE_UCS <- mean(abs(cv_UCS$residual), na.rm = TRUE)
COR_UCS <- cor(cv_UCS$observed, cv_UCS$var1.pred, use = "complete.obs")

cat("\n--- VALIDACIÓN UCS (MODELO ISOTRÓPICO) ---\n")
cat("ME:  ", ME_UCS,  "\n")
cat("MSE: ", MSE_UCS, "\n")
cat("MAE: ", MAE_UCS, "\n")
cat("COR: ", COR_UCS, "\n")
cat("-------------------------------------------\n\n")

#--------------------------------------------------------------------------------
# 9. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#---------------------------------------------------------------------------------
#Convertir objeto de validación cruzada a data frame
res_df <- as.data.frame(cv_UCS)

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
coords_res <- coordinates(cv_UCS)

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
residual_sp <- cv_UCS
residual_sp$residual <- res_df$residual

variograma_res <- variogram(residual ~ 1, residual_sp,
                            #cutoff = cutoff_val,
                            #width = width_val)
)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")

#--------------------------------------------------------------------------------
# 10. ANÁLISIS DE SENSIBILIDAD (variación de nmax) – UCS
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
    UCS ~ 1,
    UCS_utm_sp,
    model = modelo_UCS,
    nfold = nrow(UCS_utm_sp),
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

cat("\n--- ANÁLISIS DE SENSIBILIDAD (nmax) – UCS ---\n")
print(sens_results)

# Visualización del efecto de nmax
p_sens <- ggplot(sens_results, aes(x = nmax)) +
  geom_line(aes(y = RMSE), linewidth = 1) +
  geom_point(aes(y = RMSE), size = 3) +
  labs(
    title = "Análisis de sensibilidad – RMSE vs nmax (UCS)",
    x = "Número máximo de vecinos (nmax)",
    y = "RMSE"
  ) +
  theme_minimal()

print(p_sens)



# AVISO CRÍTICO (UCS): La validación cruzada para la variable Resistencia a Compresión Simple (UCS) confirma un rendimiento predictivo inaceptable.
#El bajo COR (0.138) indica una ausencia de correlación espacial efectiva, y el MAE (200.46 unidades) refleja un error promedio de predicción 
#críticamente elevado. Estos resultados son típicos en parámetros geotécnicos con alta heterogeneidad litológica y estructural, donde la varianza a
#corta distancia es dominante. 
#El mapa de Kriging generado es considerado de baja fiabilidad. 
#Se recomienda enfáticamente el uso de Cokriging (incorporando datos correlacionados como RQD o Litología) y la delimitación de dominios geomecánicos 
#para estabilizar la dependencia espacial y mejorar la predicción.