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
    
