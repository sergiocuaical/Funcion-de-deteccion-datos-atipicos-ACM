#===============================================================================
# Autores: Sergio Alejandro Florez y Alejandro Agredo Cruz
#===============================================================================

#=======================================Descripción=============================
#La siguiente función permite la detección de datos atipicos en conjuntos de datos
# multivariados, recrea el ACM usando la libreria FactoMineR, detecta posibles atipicos
# Usando las coordenadas factoriales por componente
#===============================================================================

# ================================== Paquetes requeridos========================
#install.packages("FactoMineR")
#install.packages("discFA")
#install.packages("car")
#===============================================================================

# Cargar librerías
library(FactoMineR)
library(discFA)
library(car)

#' Detección de datos atípicos mediante ACM y el criterio del Boxplot
#'
#' @param data Base de datos (data.frame o tibble).
#' @param vars Vector con los nombres o índices de las columnas que SÍ se desean incluir (por defecto NULL, usa todas).
#' @param excluir_vars Vector con los nombres o índices de las columnas que se desean OMITIR.
#' @param h Número de dimensiones factoriales a evaluar (por defecto h = 5).
#' @param c Factor multiplicativo del IQR para los límites del boxplot (por defecto c = 1.5).
#' @param incluir_acm Lógico. Si es TRUE, incluye el objeto completo del ACM en la salida.
#' @param ... Argumentos adicionales para la función MCA de FactoMineR.

detectar_atipicos_boxplot_acm_directo <- function(data, vars = NULL, excluir_vars = NULL, h = 5, c = 1.5, incluir_acm = FALSE, ...) {
  
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  # 1. Filtrar columnas según la selección del usuario
  if (!is.null(vars)) {
    data <- data[, vars, drop = FALSE]
  }
  
  if (!is.null(excluir_vars)) {
    data <- data[, !(names(data) %in% excluir_vars), drop = FALSE]
  }
  
  # 2. Convertir únicamente las columnas seleccionadas a factores
  data[] <- lapply(data, function(x) {
    if (!is.factor(x)) {
      as.factor(x)
    } else {
      x
    }
  })
  
  # 3. Ejecutar internamente el ACM utilizando FactoMineR::MCA
  res_acm <- MCA(data, graph = FALSE, ...)
  
  # Extraer coordenadas factoriales de los individuos
  coords <- res_acm$ind$coord
  
  # Validar dimensiones
  n_total <- nrow(coords)
  if (ncol(coords) < h) {
    warning(sprintf("Se solicitaron h = %d dimensiones, pero solo hay %d disponibles. Evaluando con %d dimensiones.", 
                    h, ncol(coords), ncol(coords)))
    h <- ncol(coords)
  }
  
  coords_h <- coords[, 1:h, drop = FALSE]
  ids_individuos <- if (!is.null(rownames(coords_h))) rownames(coords_h) else 1:n_total
  
  # 4. Evaluación univariada con Boxplot por cada dimensión
  matriz_atipicos_dim <- matrix(FALSE, nrow = n_total, ncol = h)
  colnames(matriz_atipicos_dim) <- paste0("Dim_", 1:h)
  
  limites_calculados <- list()
  
  for (j in 1:h) {
    col_j <- coords_h[, j]
    q1 <- quantile(col_j, 0.25, na.rm = TRUE)
    q3 <- quantile(col_j, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    lim_inf <- q1 - c * iqr
    lim_sup <- q3 + c * iqr
    
    limites_calculados[[paste0("Dim_", j)]] <- c(L_inf = lim_inf, Q1 = q1, Q3 = q3, L_sup = lim_sup)
    matriz_atipicos_dim[, j] <- (col_j < lim_inf) | (col_j > lim_sup)
  }
  
  es_atipico_global <- apply(matriz_atipicos_dim, 1, any)
  ids_atipicos <- ids_individuos[es_atipico_global]
  total_atipicos <- length(ids_atipicos)
  porcentaje_atipicos <- (total_atipicos / n_total) * 100
  
  # 5. Reporte en consola
  cat("======================================================================\n")
  cat("      DETECCIÓN DE DATOS ATÍPICOS VÍA ACM (CRITERIO BOXPLOT DIRECTO)     \n")
  cat("======================================================================\n")
  cat(sprintf(" Muestra evaluada: %d individuos | Dimensiones evaluadas (h): %d\n", n_total, h))
  cat(sprintf(" Variables utilizadas en el ACM: %d\n", ncol(data)))
  cat(sprintf(" Total de individuos atípicos detectados: %d (%.2f%%)\n", total_atipicos, porcentaje_atipicos))
  cat("======================================================================\n\n")
  
  resultado <- list(
    resumen = data.frame(
      N_Total = n_total,
      N_Atipicos = total_atipicos,
      Porcentaje_Atipicos = porcentaje_atipicos,
      h_evaluadas = h,
      c_utilizado = c
    ),
    ids_atipicos = ids_atipicos,
    posiciones_atipicos = which(es_atipico_global),
    es_atipico_vector = es_atipico_global,
    detalle_por_dimension = matriz_atipicos_dim,
    limites_dimensiones = limites_calculados
  )
  
  if (incluir_acm) {
    resultado$objeto_mca <- res_acm
    resultado$coordenadas_individuos <- coords
  }
  
  return(invisible(resultado))
}

# 1. Cargar y explorar la estructura de los datos
# data(car_data) # (Descomenta esta línea si aún no la tienes cargada en tu entorno)
datos <- car_data
str(datos)

# =====================================================================
# EJEMPLOS DE USO CON TU BASE DE DATOS
# =====================================================================

# Opción 1: Seleccionar explícitamente solo las columnas interesan (vars)
# Ejemplo usando las variables de car_data:
res_seleccion <- detectar_atipicos_boxplot_acm_directo(
  data = datos, 
  vars = c("Price", "Safety", "Exterior_Looks", "Space_comfort", "Technology"), # Solo estas entran al ACM
  h = 2, 
  c = 1.5,
  incluir_acm = TRUE
)

# Opción 2: Omitir columnas que no se deseen analizar (excluir_vars)
res_exclusion <- detectar_atipicos_boxplot_acm_directo(
  data = datos, 
  excluir_vars = c("Color", "Testimonials"), # Se excluyen estas dos y el resto se procesa
  h = 3, 
  c = 1.5,
  incluir_acm = FALSE
)

# 3. Ver los automóviles atípicos detectados
print(res_seleccion$ids_atipicos)
