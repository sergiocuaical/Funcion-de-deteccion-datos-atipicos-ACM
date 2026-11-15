#===============================================================================
# Autores: Sergio Alejandro Florez y Alejandro Agredo Cruz
#===============================================================================


# ==============================================================================
# coord_atipicos_sup_acm()
# Coordenadas de individuos suplementarios (p. ej. atípicos) en un ACM,
# calculadas con la fórmula de transición y promediando solo sobre las
# modalidades disponibles en el ACM activo.
# ==============================================================================
#
# Motivación:
#   MCA(..., ind.sup = ...) falla cuando, al retirar los atípicos, alguna
#   modalidad queda sin observaciones en la base activa (columna con suma 0).
#   Aquí cada individuo se proyecta con las modalidades que SÍ existen en el
#   ACM activo; las no soportadas se omiten del promedio.
#
# Fórmula:
#   F_is = ( sum_{k disponibles} G_ks ) / ( Q_i * sqrt(lambda_s) )
#   con G_ks = coordenada de la modalidad k en el eje s (res$var$coord),
#   lambda_s = autovalor del eje s (res$eig[, 1], sin corregir) y
#   Q_i = número de modalidades disponibles del individuo i.
#
# Parámetros:
#   res       : objeto resultante de FactoMineR::MCA() ajustado sobre la
#               base SIN atípicos (activos).
#   datos_sup : data.frame con los individuos a proyectar. Debe contener
#               únicamente las variables activas del ACM (mismos nombres
#               de columna que la base usada en MCA()).
#   ncp       : número de dimensiones a calcular. Por defecto 5.
#
# Retorna una lista con:
#   $coord       : matriz n_sup x ncp con las coordenadas (misma escala
#                  que res$ind$coord)
#   $diagnostico : data.frame por individuo con Qi (modalidades
#                  disponibles), prop_no_soportada y las modalidades
#                  no soportadas
#   $resumen     : data.frame con n_sup, n_con_no_soportadas y ncp
#
# ------------------------------------------------------------

coord_atipicos_sup_acm <- function(res, datos_sup, ncp = 5) {
  
  # ---- validaciones básicas ----
  if (!inherits(res, "MCA")) {
    stop("`res` debe ser un objeto devuelto por FactoMineR::MCA().")
  }
  if (!is.data.frame(datos_sup) || nrow(datos_sup) == 0) {
    stop("`datos_sup` debe ser un data.frame con al menos una fila.")
  }
  if (ncp < 1 || ncp > min(ncol(res$var$coord), nrow(res$eig))) {
    stop("`ncp` debe estar entre 1 y el número de ejes disponibles en `res`.")
  }
  
  Gk     <- res$var$coord[, 1:ncp, drop = FALSE]   # rownames: "Variable_modalidad"
  lambda <- res$eig[1:ncp, 1]
  Q      <- ncol(datos_sup)
  
  # ---- clave "Variable_modalidad" de cada individuo ----
  key_sup <- sapply(names(datos_sup), function(v) paste(v, as.character(datos_sup[[v]]), sep = "_"))
  if (is.null(dim(key_sup))) key_sup <- matrix(key_sup, nrow = nrow(datos_sup))
  dimnames(key_sup) <- list(rownames(datos_sup), names(datos_sup))
  
  # ---- correspondencia con las modalidades del ACM activo ----
  idx <- matrix(match(key_sup, rownames(Gk)), nrow = nrow(key_sup))   # NA = modalidad no soportada
  
  if (all(is.na(idx))) {
    stop("Ninguna modalidad de `datos_sup` coincide con rownames(res$var$coord). ",
         "Se espera el formato 'Variable_modalidad'.")
  }
  
  Qi <- rowSums(!is.na(idx))                                          # modalidades disponibles por individuo
  
  if (any(Qi == 0)) {
    warning("Hay individuos sin ninguna modalidad disponible; sus coordenadas serán NA.")
  }
  
  # ---- coordenadas: suma de G_k disponibles / (Qi * sqrt(lambda)) ----
  den <- ifelse(Qi == 0, NA, Qi)
  
  coord_sup <- sapply(1:ncp, function(s) {
    g <- matrix(Gk[as.vector(idx), s], nrow = nrow(idx))
    g[is.na(g)] <- 0
    rowSums(g) / (den * sqrt(lambda[s]))
  })
  coord_sup <- matrix(coord_sup, nrow = nrow(datos_sup), ncol = ncp)
  rownames(coord_sup) <- rownames(datos_sup)
  colnames(coord_sup) <- paste0("Dim.", 1:ncp)
  
  # ---- diagnóstico por individuo ----
  diagnostico <- data.frame(
    Qi                        = Qi,
    prop_no_soportada         = round(1 - Qi / Q, 3),
    modalidades_no_soportadas = sapply(seq_len(nrow(key_sup)), function(i)
      paste(key_sup[i, is.na(idx[i, ])], collapse = "; "))
  )
  rownames(diagnostico) <- rownames(datos_sup)
  
  # ---- resumen ----
  resumen <- data.frame(
    n_sup               = nrow(datos_sup),
    n_con_no_soportadas = sum(Qi < Q),
    ncp                 = ncp
  )
  
  # ---- salida ----
  list(
    coord       = coord_sup,
    diagnostico = diagnostico,
    resumen     = resumen
  )
}


