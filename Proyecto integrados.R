#Luis Andres Casillas Casillas A01645008
#Nubia Selene Garcidueñas Barajas A01352303

library(seqinr)

setwd("D:/TEC DE MONTERREY/SEMESTRE2/ANÁLISIS DE BIOLOGÍA COMPUTACIONAL/Evidencia 2, 2/20 Paises/20 Paises")


# Crear una lista vacía para almacenar las variables de las variantes
variantes <- list()

carpeta <- "20 paises"
archivos <- list.files(path = carpeta, pattern = "\\.fasta$", full.names = TRUE)
for (archivo in archivos) {
  # Excluir los archivos "secuencia_alineada.fasta" y "secuencias.fasta"
  if (!basename(archivo) %in% c("secuencia_alineada.fasta", "secuencias.fasta")) {
    # Extraer el nombre de la variante del archivo
    nombre <- tools::file_path_sans_ext(basename(archivo))
    # Leer el archivo y almacenarlo en la lista de variantes
    variantes[[nombre]] <- read.fasta(archivo)
  }
}

# Crear una lista vacía para almacenar las longitudes de las secuencias
longitudes <- list()

# Iterar sobre cada variante en la lista
for (nombre_variante in names(variantes)) {
  # Obtener la secuencia de la variante actual
  secuencia <- variantes[[nombre_variante]]
  # Calcular la longitud de la secuencia y almacenarla en la lista de longitudes
  longitudes[[nombre_variante]] <- length(secuencia[[1]])
}
# Ahora la lista "longitudes" contiene la longitud de cada secuencia en cada variante

print(longitudes)

# Función para comparar secuencias

  # Crear una lista vacía para almacenar los nucleótidos de las variantes
  nucleotidos <- list()
  
  # Iterar sobre los archivos de secuencias
  for (archivo in archivos) {
    # Excluir los archivos "secuencia_alineada.fasta" y "secuencias.fasta"
    if (!basename(archivo) %in% c("secuencia_alineada.fasta", "secuencias.fasta")) {
      # Extraer el nombre de la variante del archivo
      nombre <- tools::file_path_sans_ext(basename(archivo))
      # Leer el archivo de secuencia
      secuencia <- read.fasta(archivo)
      # Contar los nucleótidos y almacenarlos en la lista de nucleótidos
      nucleotidos[[nombre]] <- table(unlist(secuencia))
    }
  }

# Combinar los datos de las variantes en un solo data frame
nucleotido_data <- do.call(rbind, lapply(names(nucleotidos), function(variante) {
  data.frame(Variante = rep(variante, length(nucleotidos[[variante]])), 
             Nucleotido = names(nucleotidos[[variante]]), 
             Cantidad = as.numeric(nucleotidos[[variante]]),
             stringsAsFactors = FALSE)
}))

# Cargar la librería ggplot2
library(ggplot2)

# Crear el gráfico de barras
ggplot(nucleotido_data, aes(x = Variante, y = Cantidad, fill = Nucleotido)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Secuencias de SARS-CoV-2",
       x = "Variantes/Países",
       y = "Cantidad",
       fill = "Nucleótidos") +
  theme_minimal() +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

library(seqinr)
library(ape)
library(gridExtra)

setwd("D:/TEC DE MONTERREY/SEMESTRE2/ANÁLISIS DE BIOLOGÍA COMPUTACIONAL/Evidencia 2, 2/20 Paises/20 Paises")

Argentina <- read.fasta("Argentina.fasta")
Australia <- read.fasta("Australia.fasta")
Belgium <- read.fasta("Belgium.fasta")
Canada <- read.fasta("Canada.fasta")
Chile <- read.fasta("Chile.fasta")
China <- read.fasta("China.fasta")
Colombia <- read.fasta("Colombia.fasta")
Ecuador <- read.fasta("Ecuador.fasta")
Georgia <- read.fasta("Georgia.fasta")
Greece <- read.fasta("Greece.fasta")
India <- read.fasta("India.fasta")
Ireland <- read.fasta("Ireland.fasta")
Italy <- read.fasta("Italy.fasta")
Malaysia <- read.fasta("Malaysia.fasta")
NewZeland <- read.fasta("NewZeland.fasta")
Poland <- read.fasta("Poland.fasta")
Romania <- read.fasta("Romania.fasta")
Russia <- read.fasta("Russia.fasta")
Thailand <- read.fasta("Thailand.fasta")
UK <- read.fasta("UK.fasta")

virus <- c("China", "Argentina", "Canada", "Chile", "Colombia",
           "Ecuador", "Georgia", "Greece", "India", "Ireland",
           "Italy", "Malaysia", "NewZeland", "Poland", "Romania",
           "Russia", "Thailand", "UK", "Australia", "Belgium")

secuencias <- c(Argentina, Australia, Belgium, Canada, Chile, China, Colombia, Ecuador,
                Georgia, Greece, India, Ireland, Italy, Malaysia, NewZeland, Poland,
                Romania, Russia, Thailand, UK)

write.dna(secuencias,  file ="secuencias.fasta", format = "fasta", append = FALSE,
          nbcol = 6, colsep = " ", colw = 10)

secuencia_no_alineada <- readDNAStringSet("secuencias.fasta", format = "fasta")
secuencia_no_alineada <- OrientNucleotides(secuencia_no_alineada)

# Alineamiento
secuencia_alineada <- AlignSeqs(secuencia_no_alineada)

# Resultado
BrowseSeqs(secuencia_alineada)
ConsensusSequence(secuencia_no_alineada)

writeXStringSet(secuencia_alineada, file="secuencia_alineada.fasta")

alineada <- read.alignment("secuencia_alineada.fasta", format = "fasta")
matriz_distancia <- dist.alignment(alineada, matrix = "similarity")

temp <- as.data.frame(as.matrix(matriz_distancia))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) + scale_color_viridis()

arbol <- nj(matriz_distancia)
arbol <- ladderize(arbol)

# Arbol
plot(arbol, cex = 0.6)
title("Arbol Filogenetico")

library(ggtree)
p1 <- ggtree(arbol, branch.length='none', layout='circular') + geom_tiplab()

library(ggmsa)
p2 <- ggmsa(secuencia_no_alineada, 1, 20, color = "Chemistry_AA")

combined_plot <- grid.arrange(p1, p2, nrow = 1)
print(combined_plot)

#¿Son muy diferentes las variantes entre cada país? ¿Es diferente el SARS-CoV-2 entre las diferentes poblaciones: Asiática, Hispana, Europea o Africana?
#Las variantes del SARS-COv-2 identificadas en los países más afectados son muy similares entre sí a través de los continentes.
#Lo que indica que las variaciones genéticas entre las poblaciones de Asia, Hispanoamérica, Europa y Africa son escasas.
#Esto sugiere que no hay diferencias importantes en la estructura genética del virus en diverasa regiones globales.

#Interpretación y conclusión
#El análisis de ADN del virus revela un alto grado de similitud entre las variantes globales, con diferencias regionales más específicas.
#India y Reino Unido parecen ser centros iniciales de propagación, y las conexiones entre países como Colombia y Ecuador, o Bélgica y Polonia, reflejan patrones de movilidad y comercio.
#En general, las similitudes indican una rápida dispersión internacional del virus.
