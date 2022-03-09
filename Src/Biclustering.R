# Title: Biclustering 
# Autor: Victor Jesus Enriquez Castro
# Mail: victorec@lcg.unam.mx

# Agregamos las librerias que utilizaremos
library(QUBIC)
library(RColorBrewer)

# Una vez que descargamos el dataset que utilizaremos le damos lectura a los datos
data  <- read.table("../Data/Datos_Genes.tsv", header = FALSE, sep = "\t", quote="", row.names = 1)
matriz_1 <- data.matrix(data,rownames.force = NA)

# Hacemos el biclustering de nuestros datos en este caso utilizando el paquete QUBIC
bclusters <- biclust::biclust(matriz_1, method = BCQU())

# Visualizamos el numero de biclusters y sus dimensiones
summary(bclusters)

# Generamos Heatmaps que nos muestren de una manera mas legible los biclusters
par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(matriz_1, bclsuters, number = c(17), showlabel = TRUE, col = paleta)

par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(matriz_1, bclsuters, number = c(19), showlabel = TRUE, col = paleta)

par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(matriz_1, bclsuters, number = c(30), showlabel = TRUE, col = paleta)

par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(matriz_1, bclsuters, number = c(31), showlabel = TRUE, col = paleta)

# Generamos redes que ilustres los biclusters individuales
red <- qunetwork(matriz_1, bclusters, number = 17, group = 17, method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], layout = "spring", minimum = 0.6,
                 color = cbind(rainbow(length(red[[2]]) - 1), "blue"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = 19, group = 19, method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], layout = "spring", minimum = 0.6,
                 color = cbind(rainbow(length(red[[2]]) - 1), "blue"), edge.label = FALSE)

red <- quredwork(matriz_1, bclusters, number = 30, group = 30, method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], layout = "spring", minimum = 0.6,
                 color = cbind(rainbow(length(red[[2]]) - 1), "blue"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = 31, group = 31, method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], layout = "spring", minimum = 0.6,
                 color = cbind(rainbow(length(red[[2]]) - 1), "blue"), edge.label = FALSE)

# Generamos redes que ilustren los pares de biclusters 
red <- qunetwork(matriz_1, bclusters, number = c(17,19), group = c(17, 19), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = c(17,30), group = c(17, 30), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = c(17,31), group = c(17, 31), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = c(19,30), group = c(19, 30), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = c(19,31), group = c(19, 31), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

red <- qunetwork(matriz_1, bclusters, number = c(30,31), group = c(30, 31), method = "spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(red[[1]], groups = red[[2]], legend.cex = 0.5, layout = "spring",
                 minimum = 0.6, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

