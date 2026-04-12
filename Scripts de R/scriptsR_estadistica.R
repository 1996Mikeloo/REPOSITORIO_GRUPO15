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


########## 6. Co-mutation plot con maftools
# 6.1 Instalar y cargar maftools (vía BiocManager si es necesario)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("maftools", quietly = TRUE)) {
  BiocManager::install("maftools", ask = FALSE, update = FALSE)
}
library(maftools)

# 6.2 Preparar datos para el formato MAF
# Renombrar columnas y añadir la columna t_vaf requerida por maftools
# Asegúrate de que la columna 'vafbulk' existe en 'datoscol' después de clean_names()
# Si se llama diferente (p.ej., 'vaf_bulk'), ajústalo aquí.
datosmaf <- datoscol %>%
  rename(
    Hugo_Symbol            = hugo_symbol, 
    Tumor_Sample_Barcode   = tumor_sample_barcode,
    Chromosome             = chromosome,
    Start_Position         = start_position,
    End_Position           = end_position,
    Reference_Allele       = reference_allele,
    Tumor_Seq_Allele2      = tumor_seq_allele2,
    Variant_Classification = variant_classification
  )
  
# 6.3 Crear objeto MAF
maf_obj <- read.maf(maf = datosmaf, verbose = FALSE)

# 6.5 Obtener resumen de muestras (incluye VAF promedio si t_vaf está presente y es numérico)
ss <- getSampleSummary(maf_obj)

library(RColorBrewer)
# Extrae las categorías únicas de tipo de mutación
muts <- unique(datosmaf$Variant_Classification)
# Pide tantos colores como categorías tengas
brewer_cols <- brewer.pal(n = length(muts), name = "Set2")
# Ponles nombres iguales a las categorías
names(brewer_cols) <- muts
# 6.6 Dibujar OncoPlot con VAF
    oncoplot(
      maf                     = maf_obj,
      removeNonMutated        = TRUE,
      showTumorSampleBarcodes = TRUE,
      drawColBar              = TRUE,
      colors                  = brewer_cols,  
      showPct                 = TRUE,
      showTitle               = FALSE,
      legendFontSize          = 2,
      fontSize = 0.9,
      SampleNamefontSize = 1.4
          )

    
    
    library(dplyr)
    library(ggplot2)
    
############# 7. Grafico de barras con VAF de cada mutación ordenadas por VAF    
# 1) Ordenar las mutaciones de mayor a menor FractionalAbundance
    #    y convertir 'mutation' en factor con ese orden
    plot_data <- datoscol %>%
      arrange(desc(vafngsbulk)) %>%
      mutate(
        mutation = factor(mutation, levels = unique(mutation))
      )
    
    # 2) Elegir una paleta de colores para genes (una por gen)
    #    Aquí usamos Set2 de RColorBrewer, ajusta 'n' al número de genes únicos
    library(RColorBrewer)
    genes <- unique(plot_data$hugo_symbol)
    paleta_genes <- brewer.pal(n = length(genes), name = "Set2")
    names(paleta_genes) <- genes
    
    # 3) Crear el barplot
    ggplot(plot_data, aes(x = mutation, y = vafngsbulk, fill = hugo_symbol)) +
      geom_col() +
      scale_fill_manual(values = paleta_genes, name = "Gen") +
      labs(
        title = NULL,
        x     = "Mutación",
        y     = "Fractional Abundance (%)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 80, hjust = 1, size=10),
        plot.title  = element_text(hjust = 0.5)
      )

################# 8. Graficos Timepoints 
    library(dplyr)
    library(tidyr)
    library(ggplot2)
 
    datoscol5p <- openxlsx::read.xlsx(
      xlsxFile = "/Users/sara/Desktop/ExcelsR/ddPCR_database_090525.xlsx",
      sheet    = 2,
      colNames = TRUE
    ) %>%
      janitor::clean_names()
        
       
    # ——— 1) Define tus columnas de población y la columna de timepoint ———
    time_col  <- "sample_time" 
    id_cols   <- c("pid", "hugo_symbol", "mutation")  
    pop_cols  <- c("vafnk", "vaflt", "vaflb", "vafmo", "vafgran")
    
    
    # ——— 2) Pivotar a formato largo ———
    df_long2 <- datoscol5p %>%
      select(all_of(c(id_cols, time_col, pop_cols))) %>%
      pivot_longer(
        cols      = all_of(pop_cols),
        names_to  = "Population",
        values_to = "VAF"
      ) %>%
      mutate(
        TimeF = factor(.data[[time_col]], levels = sort(unique(.data[[time_col]])))
      )
    
    # ——— 3) Hacer el mini-panel con facet (o color) ———
    ggplot(df_long2, aes(x = TimeF, y = VAF, group = Population, color = Population)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      facet_wrap(~ pid + hugo_symbol + mutation, ncol = 3, scales = "free_y", strip.position = "top") +
      scale_color_brewer(palette = "Set2") +
      labs(
        x     = NULL,
        y     = "Fractional Allele Frequency (%)",
        color = "Población"
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position   = "right",
        strip.text        = element_text(face = "bold", size = 12),
        panel.spacing     = unit(0.5, "lines"),
        axis.text.x       = element_text(size = 10),
        axis.title.y      = element_text(size = 12)
      )
    
    geom_point(size = 2) +
      facet_wrap(~ patient_id + mutation_gene + mutation,
                 ncol       = 2,          # dos columnas de paneles
                 scales     = "free_y",   # eje Y libre en cada panel
                 strip.position = "top") +

      
      
      library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(openxlsx)
    library(janitor)
    
    # 1. Leer datos y limpiar nombres
    datoscol5p <- read.xlsx(
      xlsxFile = "/Users/sara/Desktop/ExcelsR/ddPCR_database_090525.xlsx",
      sheet    = 4,
      colNames = TRUE
    ) %>%
      clean_names()
    
    # 2. Columnas de interés
    time_col  <- "sample_time"
    id_cols   <- c("pid", "hugo_symbol", "mutation")
    pop_cols  <- c("vafnk", "vaflt", "vaflb", "vafmo", "vafgran")
    
    # 3. Pivotar a formato largo
    df_long2 <- datoscol5p %>%
      select(all_of(c(id_cols, time_col, pop_cols))) %>%
      pivot_longer(
        cols      = all_of(pop_cols),
        names_to  = "Population",
        values_to = "VAF"
      ) %>%
      mutate(
        TimeF = factor(.data[[time_col]], levels = sort(unique(.data[[time_col]]))),
        # Renombrar categorías de población para la leyenda
        Population = recode(Population,
                            vafnk   = "NK",
                            vaflt   = "LT",
                            vaflb   = "LB",
                            vafmo   = "Mono",
                            vafgran = "Gran")
      )
    
    # 4. Obtener todos los pacientes únicos
    pacientes <- unique(df_long2$pid)
    
    # 5. Crear gráfico por paciente
    for (paciente in pacientes) {
      p <- df_long2 %>%
        filter(pid == paciente) %>%
        ggplot(aes(x = TimeF, y = VAF, group = Population, color = Population)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        facet_wrap(~ hugo_symbol + mutation, ncol = 2, scales = "free_y") +
        scale_color_brewer(palette = "Set2", name = "Población") +
        labs(
          title = paste("PID:", paciente),
          x = "Timepoint",
          y = "Fractional Allele Frequency (%)"
        ) +
        theme_minimal(base_size = 15) +
        theme(
          legend.position = "right",
          strip.text = element_text(face = "bold", size = 12),
          panel.spacing = unit(0.5, "lines"),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12)
        )
      
      print(p)  # para visualizar en RStudio
      
      # Guardar automáticamente si quieres:
      # ggsave(paste0("PID_", paciente, ".pdf"), plot = p, width = 10, height = 6)
    }
    