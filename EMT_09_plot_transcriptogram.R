# plot_transcriptogram.R
# Gera gráficos de diferenças entre transcriptogramas originais e reconstruídos.

library(ggplot2)

# Carregar dados reconstruídos
load("reconstructed_r0.RData")
load("reconstructed_r30.RData")

# Função para plotar diferenças
plot_transcriptogram <- function(df_original, df_reconstruido, max_pcs = 10) {
  df_original <- as.data.frame(df_original)
  df_reconstruido <- as.data.frame(df_reconstruido)
  max_pcs <- min(max_pcs, ncol(df_original))
  
  diff_data <- data.frame(PCs = integer(), Difference = numeric())
  
  for (num_pcs in 1:max_pcs) {
    diff <- abs(df_original[, num_pcs] - df_reconstruido[, num_pcs])
    total_diff <- sum(diff, na.rm = TRUE)
    diff_data <- rbind(diff_data, data.frame(PCs = num_pcs, Difference = total_diff))
  }
  
  ggplot(diff_data, aes(x = factor(PCs), y = Difference)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Diferença por PC entre Original e Reconstruído",
         x = "Número da PC",
         y = "Soma das Diferenças Absolutas")
}

# Plotar diferenças
plot_transcriptogram(df_t_matrix_batch_effect_correction_R0_for_PCA, reconstructed_r30, max_pcs = 10000)