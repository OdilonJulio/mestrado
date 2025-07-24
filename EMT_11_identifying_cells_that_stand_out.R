

###### 1 - Análise de Variabilidade Global

# Calcular a variabilidade (desvio padrão) por célula
cell_variability <- apply(t_matrix_R30@transcriptogramS2[, -c(1,2)], 2, sd)  # Exclui colunas 'Protein' e 'Position'

# Ordenar células por variabilidade
top_cells <- sort(cell_variability, decreasing = TRUE)

# Visualizar as 10 células mais variáveis
head(top_cells, 10)

# Gráfico de variabilidade
ggplot(data.frame(Cell=names(cell_variability), Variability=cell_variability), 
       aes(x=reorder(Cell, Variability), y=Variability)) +
  geom_point() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Células", y = "Desvio Padrão", title = "Variabilidade de Expressão por Célula")

###### 2 - Análise de Expressão Média

# Calcular média de expressão por célula
cell_means <- colMeans(t_matrix_R30@transcriptogramS2[, -c(1,2)])

# Identificar outliers
mean_outliers <- cell_means[cell_means > mean(cell_means) + 2*sd(cell_means) | 
                              cell_means < mean(cell_means) - 2*sd(cell_means)]

# Visualização
boxplot(cell_means, main="Distribuição de Expressão Média por Célula")
points(rep(1, length(mean_outliers)), mean_outliers, col="red", pch=19)

###### 3 - Análise de PCA (Componentes Principais)

# Executar PCA (se ainda não tiver sido feito)
pca_data <- t_matrix_R30@transcriptogramS2[, -c(1,2)]
pca_result <- prcomp(t(pca_data), scale. = TRUE)

# Visualizar células no espaço PCA
library(ggfortify)
autoplot(pca_result, label = TRUE, label.size = 3,
         main = "Células no Espaço PCA") +
  theme_minimal()

# Identificar células extremas nos componentes principais
pca_scores <- pca_result$x
outlier_cells <- rownames(pca_scores)[abs(pca_scores[,1]) > 3*sd(pca_scores[,1]) | 
                                        abs(pca_scores[,2]) > 3*sd(pca_scores[,2])]

###### 4 - Análise de Clusterização

# Clusterização hierárquica
dist_matrix <- dist(t(pca_data))
hclust_result <- hclust(dist_matrix)

# Dendrograma
plot(hclust_result, cex = 0.6, main = "Clusterização de Células")
rect.hclust(hclust_result, k = 5, border = 2:6)  # Agrupa em 5 clusters

# Identificar clusters pequenos ou isolados
cluster_groups <- cutree(hclust_result, k = 5)
small_clusters <- which(table(cluster_groups) < 0.1*length(cluster_groups))
outlier_cluster_cells <- names(cluster_groups)[cluster_groups %in% small_clusters]

###### 5 - Análise de Genes Marcadores

# Supondo que você tenha uma lista de genes marcadores
marker_genes <- c("GENE1", "GENE2", "GENE3")  # Substitua pelos seus genes

# Encontrar células com alta expressão de múltiplos marcadores
marker_expression <- colSums(t_matrix_R30@transcriptogramS2[rownames(t_matrix_R30@transcriptogramS2) %in% marker_genes, -c(1,2)])
top_marker_cells <- names(sort(marker_expression, decreasing = TRUE))[1:10]

###### 6 - Análise Combinada (Score de Outlier)

# Criar um score combinado de outlier
outlier_score <- scale(cell_variability) + scale(cell_means) + 
  scale(abs(pca_scores[,1])) + scale(abs(pca_scores[,2]))

# Identificar as células mais extremas
top_outliers <- names(sort(outlier_score, decreasing = TRUE))[1:10]

###### 7 - Visualização de Células Selecionadas

# Supondo que 'outlier_cells' contém suas células de interesse
plot_data <- t_matrix_R30@transcriptogramS2[, c("Position", outlier_cells[1:3])]

ggplot(melt(plot_data, id.vars = "Position"), 
       aes(x = Position, y = value, color = variable)) +
  geom_line() +
  labs(title = "Padrão de Expressão de Células que se Destacam",
       y = "Expressão", color = "Célula") +
  theme_minimal()


###### 8 - Soma total da expressão por célula


mat <- t_matrix_R30@transcriptogramS2
expr_only <- mat[, !colnames(mat) %in% c("Protein", "Position")]

total_expr <- colSums(expr_only)
highlight_cell <- names(total_expr)[which.max(total_expr)]

barplot(sort(total_expr, decreasing = TRUE), las = 2, main = "Total expression per cell")
abline(v = which(names(sort(total_expr, decreasing = TRUE)) == highlight_cell), col = "red", lwd = 2)


###### 9 - Variância por célula

cell_var <- apply(expr_only, 2, var)
highlight_cell_var <- names(cell_var)[which.max(cell_var)]

plot(cell_var, main = "Variance per cell", ylab = "Variance", xlab = "Cell index")
points(which.max(cell_var), max(cell_var), col = "red", pch = 19)

###### 10 - Distância da média global

mean_profile <- rowMeans(expr_only)
dists <- apply(expr_only, 2, function(cell) sqrt(sum((cell - mean_profile)^2)))
highlight_cell_dist <- names(dists)[which.max(dists)]

hist(dists, breaks = 50, main = "Euclidean distance from mean cell", xlab = "Distance")
abline(v = dists[highlight_cell_dist], col = "red", lwd = 2)

###### 11 - PCA para detectar outliers

pca_res <- pca_result_R30$pca_result
pc1_scores <- pca_res$x[,1]

highlight_cell_pca <- names(pc1_scores)[which.max(abs(pc1_scores))]

plot(pc1_scores, main = "PCA1 scores by cell", ylab = "PC1", xlab = "Cell")
abline(v = which.max(abs(pc1_scores)), col = "blue")

###### 12 - Z-score por célula

z_matrix <- scale(t(expr_only))  # z-score por gene (linhas: células)
cell_z_var <- apply(z_matrix, 1, var)

highlight_cell_z <- rownames(z_matrix)[which.max(cell_z_var)]


###### 13 - Visualizar a célula destacada

cell_to_plot <- highlight_cell_var  # ou qualquer outro dos métodos acima

plot(
  t_matrix_R30@transcriptogramS2$Position,
  expr_only[[cell_to_plot]],
  type = "l", col = "blue", lwd = 1,
  main = paste("Expression profile of standout cell:", cell_to_plot),
  xlab = "Gene position", ylab = "Expression"
)

