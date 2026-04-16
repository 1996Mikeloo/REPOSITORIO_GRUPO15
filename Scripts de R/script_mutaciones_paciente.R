########## 1. Instalación y carga de librerías
# openxlsx: lectura de Excel
install.packages("openxlsx")
library(openxlsx)

# dplyr y janitor: limpieza de datos
install.packages("dplyr")
install.packages("janitor")
library(dplyr)
library(janitor)

# ggplot2: visualización
auto_pkgs <- c("ggplot2")
install.packages(setdiff(auto_pkgs, installed.packages()[,"Package"]))
library(ggplot2)

########## 2. Importar datos desde Excel
datoscol <- openxlsx::read.xlsx(
  xlsxFile = "/Users/sara/Desktop/ExcelsR/ddPCR_database_090525.xlsx",
  sheet    = 1,
  colNames = TRUE
) %>%
  janitor::clean_names()

########## 3. Cálculo de correlación de Spearman
# Extraer columnas de interés ALL
vaf1 <- datoscol$vafngsbulk   # VAF NGS bulk
vaf2 <- datoscol$vafbulk   # VAF PCR bulk
vaf3 <- datoscol$vafnk
vaf4 <- datoscol$vaflt
vaf5 <- datoscol$vaflb
vaf6 <- datoscol$vafmo
vaf7 <- datoscol$vafgran

# Test de Spearman 
shapiro.test(vaf1)
shapiro.test(vaf2)
median(datoscol$patient_age, na.rm = TRUE)
range(datoscol$patient_age, na.rm= TRUE)

# Rango (mínimo y máximo)
rango  <- range(x, na.rm = TRUE)
spearman_test <- cor.test(
  x      = vaf1,
  y      = vaf2,
  method = "spearman",
  exact  = FALSE
)
print(spearman_test)

rango  <- range(x, na.rm = TRUE)
spearman_test <- cor.test(
  x      = vaf2,
  y      = vaf7,
  method = "spearman",
  exact  = FALSE
)
print(spearman_test)

## Gráfico de correlación Spearman 
ggplot(datoscol, aes(x = vaf1, y = vaf2)) +
  geom_point(color = "blue", alpha = 0.6) +       # puntos en azul
  geom_smooth(method = "loess", se = TRUE,
              color = "red", fill = "purple",    # línea en rojo y banda en rosa
              size = 1) +
  # Anotar coeficiente de Spearman (rho)
  annotate(
    "text",
    x = min(vaf2, na.rm = TRUE),
    y = max(vaf1, na.rm = TRUE),
    label = paste0("rho = ", round(spearman_test$estimate, )),
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 6
  ) +
  labs(
    title = NULL,
    x = "VAF NGS",
    y = "VAF ddPCR"
  ) +
  theme_minimal() +
  theme(
    axis.title.x  = element_text(size = 16),  # ↑ tamaño del título X
    axis.title.y  = element_text(size = 16),  # ↑ tamaño del título Y
    axis.text.x   = element_text(size = 14),  # ↑ tamaño de los valores X
    axis.text.y   = element_text(size = 14)   # ↑ tamaño de los valores Y
  )



########## 4. Diagrama circular de pacientes con LNH
# Datos de ejemplo: CHIP vs No CHIP
lnh_counts <- data.frame(
  group = c("CHIP", "No CHIP"),
  count = c(10, 15)
) %>%
  mutate(
    pct   = count / sum(count),
    label = paste0(group, "
", count, " (", round(pct * 100, 1), "%)")
  )

# Gráfico circular con colores personalizados
# Define un vector de colores (ajusta nombres según tus preferencias)
custom_colors1 <- c("coral3","darkolivegreen3")  # CHIP y No CHIP respectivamente

ggplot(lnh_counts, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 7,
    fontface = "bold"
  ) +
  scale_fill_manual(values = custom_colors1) +  # aplicar colores personalizados
  labs(
    title = "Pacientes con LNH",
    fill  = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title     = element_text(hjust = 0.5)
  )



########## 4. Diagrama circular de genes mutados
genetype_counts <- data.frame(
  group = c("ASXL1", "CBL", "DNMT3A", "GNB1", "JAK2", "PPM1D", "TET2", "TP53"),
  count = c(2, 1, 3, 1, 1, 5, 4, 4)
) %>%
  mutate(
    pct   = count / sum(count),
    label = paste0(group, "
", count, " (", round(pct * 100, 1), "%)")
  )

ggplot(genetype_counts, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(x = 1.85, label = label),
    position = position_stack(vjust = 0.5),
    size = 6.5,
    fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2") +      # paleta Set2
  labs(
    title = NULL,
    fill  = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title     = element_text(hjust = 0.5)
  )


########## 5. Distribución de mutaciones por paciente
# 5.1 Calcular mutaciones por paciente
a_mut_per_patient <- datoscol %>%
  count(datoscol$tumor_sample_barcode, name = "n_mutations")

# 5.2 Calcular distribución: cuántos pacientes tienen X mutaciones
dist_mut_patients <- a_mut_per_patient %>%
  count(n_mutations, name = "n_patients") %>%
  arrange(n_mutations)

# 5.3 Barplot de distribución
ggplot(dist_mut_patients, aes(
  x = factor(n_mutations),  # número de mutaciones
  y = n_patients            # número de pacientes
)) +
  geom_col(fill = "darkgreen") +
  geom_text(
    aes(label = n_patients),
    vjust = -0.5,
    size = 7
  ) +
  labs(
    title = NULL,
    x     = "Número de mutaciones",
    y     = "Número de pacientes",
    size = 7
  ) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),  # tamaño del título del eje X
    axis.title.y = element_text(size = 20)   # tamaño del título del eje Y
  )
