library(FactoClass)
library(car)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggrepel)

library(discFA)

data(car_data)
str(car_data)



#acm original
car_data <- car_data %>%
  mutate(across(everything(), as.factor))

res.MCA2 <- MCA(car_data, graph = FALSE)
fviz_mca_biplot(res.MCA2, title = "ACM original (sin depurar)")
res.MCA2$eig[1:5, ]


#aplicar la función sobre la bd

resultado2 <- deteccion_atipicos_acm(res.MCA2, car_data, inercia_min = 60, alpha = 0.05)

resultado2$resumen
rownames(car_data)[resultado2$es_atipico]     
resultado2$datos_atipicos 

# impacto de quitar atipicos sobre el acm

res.MCA.limpio2 <- MCA(resultado2$datos_limpios, graph = FALSE)

#comparación de la inercia
comparacion_inercia2 <- data.frame(
  Dim               = 1:5,
  Inercia_original  = round(res.MCA2$eig[1:5, 2], 2),
  Inercia_depurado  = round(res.MCA.limpio2$eig[1:5, 2], 2),
  Acum_original     = round(res.MCA2$eig[1:5, 3], 2),
  Acum_depurado     = round(res.MCA.limpio2$eig[1:5, 3], 2)
)
comparacion_inercia2

## comparacion plano ind
fviz_mca_ind(
  res.MCA2,
  col.ind = ifelse(resultado2$es_atipico, "Atípico", "Normal"),
  palette = c("red", "blue"),
  legend.title = "Individuo",
  title = "ACM | Car",
  mean.point = FALSE,
  alpha.ind = 0.8,
  geom.ind = c("point","text"), 
  pointsize = 1.4,
  ggtheme = theme_bw(base_size = 11)
) +
  labs(
    subtitle = "Representación de individuos identificando atípicos"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )




fviz_mca_ind(
  res.MCA.limpio2,
  title = "ACM | Car",
  mean.point = FALSE,
  alpha.ind = 0.8,
  geom.ind = c("point","text"), 
  pointsize = 1.4,
  ggtheme = theme_bw(base_size = 11)
) +
  labs(
    subtitle = "Representación de individuos sin los atípicos"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


#Variables

fviz_mca_var(
  res.MCA2,
  title = "ACM | Tecnoestrés",
  alpha.var = 0.8,
  geom.var = c("point", "text"),
  shape.var = 17, 
  pointsize = 1.5,
  labelsize = 3.8,
  ggtheme = theme_bw(base_size = 11)
) +
  labs(
    subtitle = "Representación inicial de las categorías de las variables"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


fviz_mca_var(
  res.MCA.limpio2,
  title = "ACM | Tecnoestrés",
  alpha.var = 0.8,
  geom.var = c("point", "text"),
  shape.var = 17, 
  pointsize = 1.5,
  labelsize = 3.8,
  ggtheme = theme_bw(base_size = 11)
) +
  labs(
    subtitle = "Representación de las categorías luego de quitar atípicos"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


#incluir atip como supl
res.MCA.consup2 <- MCA(car_data, ind.sup = which(resultado2$es_atipico), graph = FALSE)












#incluir atip como supl (calculo manual con la funcion, promedio sobre modalidades disponibles)
#
sup <- coord_atipicos_sup_acm(res.MCA.limpio2,
                               car_data[resultado2$es_atipico, , drop = FALSE])
sup$coord
sup$diagnostico

# verificación: sobre los individuos activos debe reproducir res$ind$coord
chk <- coord_atipicos_sup_acm(res.MCA.limpio2,
                             car_data[!resultado2$es_atipico, , drop = FALSE])
all.equal(unname(chk$coord), unname(res.MCA.limpio2$ind$coord[, 1:5]))

# gráfico (activos + atípicos proyectados)
library(factoextra); library(ggrepel)
df_sup <- data.frame(sup$coord, prop = sup$diagnostico$prop_no_soportada,
                    id = rownames(sup$coord))

fviz_mca_ind(res.MCA.limpio2, col.ind = "blue", label = "ind",
            labelsize = 3, title = "ACM sin atípicos + atípicos proyectados") +
 geom_point(data = df_sup, aes(x = Dim.1, y = Dim.2),
              colour = "red", shape = 19, size = 2.5, inherit.aes = FALSE) +
 geom_text_repel(data = df_sup, aes(x = Dim.1, y = Dim.2, label = id),
                 colour = "red", size = 3, inherit.aes = FALSE) +
 theme_minimal()


#mejorar vista
fviz_mca_ind(
  res.MCA.limpio2,
  title = "ACM | Car",
  col.ind = "blue",
  mean.point = FALSE,
  alpha.ind = 0.8,
  geom.ind = c("point", "text"),
  pointsize = 1.4,
  labelsize = 3,
  ggtheme = theme_bw(base_size = 11)
) +
  geom_point(
    data = df_sup, 
    aes(x = Dim.1, y = Dim.2),
    colour = "red", 
    shape = 19, 
    size = 2.5, 
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = df_sup, 
    aes(x = Dim.1, y = Dim.2, label = id),
    colour = "red", 
    size = 3, 
    inherit.aes = FALSE
  ) +
  labs(
    subtitle = "Individuos activos (azul) y atípicos proyectados (rojo)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

















