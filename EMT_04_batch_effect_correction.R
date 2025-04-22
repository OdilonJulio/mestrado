# batch_effect_correction.R
# Realiza a análise de transcriptograma.

# Carregar matriz normalizada
load("normalized_matrix.RData")

# Filtrar colunas para "notreated_batch1" e "notreated_batch2"
cols_batch1 <- grep("^notreated-batch1", colnames(normalized_matrix))
cols_batch2 <- grep("^notreated-batch2", colnames(normalized_matrix))

# Submatrizes para cada batch
sub_matrix_batch1 <- normalized_matrix[, cols_batch1]
sub_matrix_batch2 <- normalized_matrix[, cols_batch2]

# Calcular as médias por linha
mean_batch1 <- rowMeans(sub_matrix_batch1)
mean_batch2 <- rowMeans(sub_matrix_batch2)

# Criar um data frame com os resultados
result <- data.frame(
  Gene = rownames(normalized_matrix),
  Mean_Batch1 = mean_batch1,
  Mean_Batch2 = mean_batch2,
  Ratio_Batch2_Batch1 = ifelse(mean_batch1 == 0, 0, mean_batch2 / mean_batch1) # Divisão das médias
)

## Multiplicando todas as colunas de normalized_matrix que possuem "batch1" no colname.

colnames_batch1 <- grep("batch1", colnames(normalized_matrix))

# criando cópia de normalized_matrix

copy_matrix <- normalized_matrix

# Multiplicando a razão pelas colunas selecionadas.

copy_matrix[, colnames_batch1] <- sweep(normalized_matrix[, colnames_batch1], 1, result$Ratio_Batch2_Batch1, "*")


# Salvar matriz com efeito de lote corrigido.
save(copy_matrix, file = "copy_matrix.RData")
