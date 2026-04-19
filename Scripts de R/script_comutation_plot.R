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
    

    
