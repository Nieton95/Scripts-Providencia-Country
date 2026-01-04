# VALIDACIÓN GEOESTADÍSTICA – Stiff Soil Vs Criteria Map
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
VsStiff <- read.csv("C:/Users/ivann/OneDrive/Escritorio/Thesis/GEOSTADISTIC ANALYSIS/STIFF SOIL VS CRITERIA MAP/Vs Stiff Soil Depth.csv",
                    stringsAsFactors = FALSE)

# Convertir a sf (WGS84)
VsStiff_sf <- st_as_sf(VsStiff, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# Reproyectar a UTM zona 13N
VsStiff_utm <- st_transform(VsStiff_sf, crs = 32613)

# Convertir a Spatial
VsStiff_utm_sp <- as(VsStiff_utm, "Spatial")

cat("Puntos cargados:", nrow(VsStiff_utm), "\n")
print(st_crs(VsStiff_utm))

#---------------------------------------------------------------------------------
# 2. ANÁLISIS DE NORMALIDAD – PROFUNDIDAD DE SUELO DENSO (Vs)
#---------------------------------------------------------------------------------
stiff_vals <- VsStiff_utm_sp$Prof_stiff

# Histograma
p_hist_norm <- ggplot(data.frame(Prof_stiff = stiff_vals), aes(x = Prof_stiff)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    title = "Histograma – Profundidad de suelo denso (Vs)",
    x = "Profundidad suelo denso (m)",
    y = "Frecuencia"
  ) +
  theme_minimal()

print(p_hist_norm)

# QQ-plot
p_qq_norm <- ggplot(data.frame(Prof_stiff = stiff_vals), aes(sample = Prof_stiff)) +
  stat_qq(color = "steelblue") +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ-plot – Profundidad de suelo denso (Vs)",
    x = "Cuantiles teóricos",
    y = "Cuantiles observados"
  ) +
  theme_minimal()

print(p_qq_norm)

# Prueba de Shapiro-Wilk
shapiro_test <- shapiro.test(stiff_vals)

cat("\n--- PRUEBA DE NORMALIDAD (Shapiro-Wilk) – PROFUNDIDAD DE SUELO DENSO (Vs) ---\n")
print(shapiro_test)

#---------------------------------------------------------------------------------
# 3. VARIOGRAMA EXPERIMENTAL
#---------------------------------------------------------------------------------
variograma_emp_Stiff <- variogram(Prof_stiff ~ 1, VsStiff_utm_sp)
plot(variograma_emp_Stiff, main = "Variograma experimental – Profundidad suelo denso (Vs)")

#---------------------------------------------------------------------------------
# 4. EVALUACIÓN DE ANISOTROPÍA
#---------------------------------------------------------------------------------
variograma_dir_Stiff <- variogram(
  Prof_stiff ~ 1, VsStiff_utm_sp,
  alpha = c(0, 45, 90, 135)
)

plot(variograma_dir_Stiff,
     main = "Variogramas direccionales – Suelo denso (Vs)")

#---------------------------------------------------------------------------------
# 5. AJUSTE DEL MODELO ISOTRÓPICO
#---------------------------------------------------------------------------------
modelo_inicial_Stiff <- vgm(
  psill  = var(VsStiff_utm_sp$Prof_stiff, na.rm = TRUE) / 2,
  model  = "Sph",
  range  = diff(range(st_coordinates(VsStiff_utm)[,1])) / 4,
  nugget = var(VsStiff_utm_sp$Prof_stiff, na.rm = TRUE) / 10
)

modelo_Stiff <- fit.variogram(variograma_emp_Stiff, modelo_inicial_Stiff)

plot(variograma_emp_Stiff, modelo_Stiff,
     main = "Ajuste del variograma – Modelo esférico isotrópico (Stiff Soil Depth)")
print(modelo_Stiff)

#---------------------------------------------------------------------------------
# 6. GENERACIÓN DE MALLA EN UTM (metros)
#---------------------------------------------------------------------------------
bb_Stiff <- st_bbox(VsStiff_utm)

xmin <- as.numeric(bb_Stiff["xmin"])
ymin <- as.numeric(bb_Stiff["ymin"])
xmax <- as.numeric(bb_Stiff["xmax"])
ymax <- as.numeric(bb_Stiff["ymax"])

expand_factor <- 100
bbox_exp_Stiff <- c(
  xmin = xmin - expand_factor,
  ymin = ymin - expand_factor,
  xmax = xmax + expand_factor,
  ymax = ymax + expand_factor
)

res_m <- 10

x_seq <- seq(bbox_exp_Stiff["xmin"], bbox_exp_Stiff["xmax"], by = res_m)
y_seq <- seq(bbox_exp_Stiff["ymin"], bbox_exp_Stiff["ymax"], by = res_m)

grid_Stiff <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid_Stiff) <- ~x + y
gridded(grid_Stiff) <- TRUE
proj4string(grid_Stiff) <- CRS("+init=epsg:32613")

#---------------------------------------------------------------------------------
# 7. KRIGING – MAPA DE INCERTIDUMBRE
#---------------------------------------------------------------------------------
kriging_Stiff <- krige(
  Prof_stiff ~ 1,
  VsStiff_utm_sp,
  grid_Stiff,
  model = modelo_Stiff,
  nmax = 5
)

Stiff_uncert_df <- as.data.frame(kriging_Stiff)
colnames(Stiff_uncert_df)[1:2] <- c("X","Y")

ggplot(Stiff_uncert_df) +
  geom_raster(aes(x = X, y = Y, fill = var1.var)) +
  scale_fill_viridis_c(option = "A", name = "Varianza") +
  labs(
    title = "Varianza de kriging – Profundidad suelo denso (modelo isotrópico)",
    x = "Easting (m)", y = "Northing (m)"
  ) +
  coord_equal() + theme_minimal()

#---------------------------------------------------------------------------------
# 8. VALIDACIÓN CRUZADA (LEAVE-ONE-OUT)
#---------------------------------------------------------------------------------
cv_Stiff <- krige.cv(
  Prof_stiff ~ 1,
  VsStiff_utm_sp,
  model = modelo_Stiff,
  nfold = nrow(VsStiff_utm_sp)
)

cv_Stiff$residual <- cv_Stiff$observed - cv_Stiff$var1.pred

ME_Stiff  <- mean(cv_Stiff$residual, na.rm = TRUE)
MSE_Stiff <- mean(cv_Stiff$residual^2, na.rm = TRUE)
MAE_Stiff <- mean(abs(cv_Stiff$residual), na.rm = TRUE)
COR_Stiff <- cor(cv_Stiff$observed, cv_Stiff$var1.pred, use = "complete.obs")

cat("\n--- VALIDACIÓN STIFF SOIL (MODELO ISOTRÓPICO) ---\n")
cat("ME:  ", ME_Stiff,  "\n")
cat("MSE: ", MSE_Stiff, "\n")
cat("MAE: ", MAE_Stiff, "\n")
cat("COR: ", COR_Stiff, "\n")
cat("-------------------------------------------\n\n")

#--------------------------------------------------------------------------------
# 9. ANÁLISIS DE RESIDUOS (POST–VALIDACIÓN CRUZADA)
#--------------------------------------------------------------------------------
#Convertir objeto de validación cruzada a data frame
res_df <- as.data.frame(cv_Stiff)

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
coords_res <- coordinates(cv_Stiff)

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
residual_sp <- cv_Stiff
residual_sp$residual <- res_df$residual

variograma_res <- variogram(residual ~ 1, residual_sp,
                            #cutoff = cutoff_val,
                            #width = width_val)
)

plot(variograma_res,
     main = "Variograma de residuos (verificación de independencia espacial)")

#--------------------------------------------------------------------------------
# 10. ANÁLISIS DE SENSIBILIDAD (variación de nmax) – SUELO DENSO (Vs)
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
    Prof_stiff ~ 1,
    VsStiff_utm_sp,
    model = modelo_Stiff,
    nfold = nrow(VsStiff_utm_sp),
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

cat("\n--- ANÁLISIS DE SENSIBILIDAD (nmax) – SUELO DENSO (Vs) ---\n")
print(sens_results)

# Visualización del efecto de nmax
p_sens <- ggplot(sens_results, aes(x = nmax)) +
  geom_line(aes(y = RMSE), linewidth = 1) +
  geom_point(aes(y = RMSE), size = 3) +
  labs(
    title = "Análisis de sensibilidad – RMSE vs nmax (Suelo denso Vs)",
    x = "Número máximo de vecinos (nmax)",
    y = "RMSE"
  ) +
  theme_minimal()

print(p_sens)

# AVISO CRÍTICO: El análisis de Validación Cruzada para la variable STIFF SOIL (Suelo Rígido) arrojó un índice de Correlación (COR) de -0.604, 
#indicando un fallo severo en la capacidad predictiva del modelo de Kriging Ordinario y una relación inversa entre valores predichos y reales. 
#Este resultado, junto con la fuerte varianza a corta distancia (alto 'nugget' en el variograma), confirma la ausencia de correlación espacial útil 
#atribuible a la alta heterogeneidad del subsuelo a la escala de muestreo actual. 
#El mapa generado es meramente ilustrativo de los datos existentes. 
#Se recomienda incrementar la densidad de sondeos y aplicar zonificación geomecánica en futuros análisis para establecer una dependencia espacial
#robusta.
