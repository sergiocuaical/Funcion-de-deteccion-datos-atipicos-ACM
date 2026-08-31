# ============================================================
# deteccion_atipicos_acm()
# Detección de atípicos multivariados sobre coordenadas de un ACM
# usando distancia de Mahalanobis, con umbral chi-cuadrado.
# ============================================================
#
# Parámetros:
#   res         : objeto resultante de FactoMineR::MCA()
#   datos       : data.frame original usado para el ACM (misma cantidad
#                 y orden de filas que res$ind$coord). Se usa para
#                 devolver las bases separadas (limpia / atípicos).
#   inercia_min : porcentaje mínimo de inercia acumulada (0-100) para
#                 escoger H (número de ejes retenidos). Por defecto 60.
#   alpha       : nivel de significancia para el umbral chi-cuadrado.
#                 Por defecto 0.05.
#
# Retorna una lista con:
#   $H              : número de ejes retenidos
#   $umbral         : valor umbral (qchisq(1-alpha, df = H))
#   $D2             : vector de distancias de Mahalanobis al cuadrado
#   $es_atipico     : vector lógico (TRUE = atípico)
#   $datos_limpios  : data.frame sin las filas atípicas
#   $datos_atipicos : data.frame solo con las filas atípicas
#   $resumen        : data.frame con n_total, n_atipicos, pct_atipicos,
#                      H y umbral usados
#
# ------------------------------------------------------------

deteccion_atipicos_acm <- function(res, datos, inercia_min = 60, alpha = 0.05) {
  
  # ---- validaciones básicas ----
  if (!inherits(res, "MCA")) {
    stop("`res` debe ser un objeto devuelto por FactoMineR::MCA().")
  }
  if (nrow(datos) != nrow(res$ind$coord)) {
    stop("`datos` debe tener el mismo número de filas que res$ind$coord.")
  }
  if (inercia_min <= 0 || inercia_min >= 100) {
    stop("`inercia_min` debe estar entre 0 y 100.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("`alpha` debe estar entre 0 y 1.")
  }
  
  coord <- res$ind$coord
  eig   <- res$eig
  
  # ---- selección de H según inercia acumulada ----
  H <- which(eig[, 3] >= inercia_min)[1]
  
  if (is.na(H)) {
    H <- ncol(coord)
  } else {
    H <- min(H, ncol(coord))
  }
  
  coord_H <- as.matrix(coord[, 1:H, drop = FALSE])
  
  # ---- umbral chi-cuadrado ----
  umbral <- qchisq(1 - alpha, df = H)
  
  # ---- distancia de Mahalanobis ----
  centro <- colMeans(coord_H)
  S      <- cov(coord_H)
  
  D2 <- mahalanobis(coord_H, center = centro, cov = S)
  
  es_atipico <- D2 > umbral
  
  # ---- separación de bases ----
  datos_limpios  <- datos[!es_atipico, , drop = FALSE]
  datos_atipicos <- datos[ es_atipico, , drop = FALSE]
  
  # ---- resumen ----
  n_total     <- nrow(datos)
  n_atipicos  <- sum(es_atipico)
  pct_atipicos <- round(100 * n_atipicos / n_total, 2)
  
  resumen <- data.frame(
    n_total      = n_total,
    n_atipicos   = n_atipicos,
    pct_atipicos = pct_atipicos,
    H            = H,
    alpha        = alpha,
    umbral       = round(umbral, 4)
  )
  
  # ---- salida ----
  list(
    H              = H,
    umbral         = umbral,
    D2             = D2,
    es_atipico     = es_atipico,
    datos_limpios  = datos_limpios,
    datos_atipicos = datos_atipicos,
    resumen        = resumen
  )
}



 