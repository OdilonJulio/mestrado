# normalization.R
# Normaliza os dados e ajusta os batches.

library(Seurat)

# Carregar objeto combinado
load("filtered_combined_seurat.RData")

# Extrair a lista de camadas
layers <- Layers(filtered_combined_seurat)

# Criar uma lista para armazenar as matrizes de expressão de cada camada
count_matrices <- lapply(layers, function(layer_name) {
  # Acessar os dados de contagem para cada camada
  GetAssayData(filtered_combined_seurat, assay = "RNA", slot = layer_name)
})

# Combinar as matrizes de contagem
combined_expression_matrix <- do.call(cbind, count_matrices)

# Verificar as dimensões da matriz combinada
dim(combined_expression_matrix)


combined_expression_matrix_full <- as.matrix(combined_expression_matrix)

# Normalização Total Count

# 1. Calcular a soma de cada coluna (Total Count por célula)
column_sums <- colSums(combined_expression_matrix_full)

# 2. Dividir cada elemento pela soma da respectiva coluna
normalized_matrix <- sweep(combined_expression_matrix_full, 2, column_sums, FUN = "/")

# 3. Verificar se a soma de cada coluna é 1
column_sums_normalized <- colSums(normalized_matrix)

# Salvar matriz normalizada
save(normalized_matrix, file = "normalized_matrix.RData")
