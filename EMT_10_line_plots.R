# line_plots.R
# Gera gráficos de linhas para os transcriptogramas.

library(ggplot2)

# Carregar dados
load("reconstructed_r0.RData")
load("reconstructed_r30.RData")

# Gráfico de linhas para R0
ggplot(df_t_matrix_batch_effect_correction_R0_for_PCA, aes(x = Day, y = PC1, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(x = "Dia", y = "Componente Principal 1", title = "Gráfico de Linha de PCA (R0)")

# Gráfico de linhas para R30
plot(t(df_t_matrix_batch_effect_correction_R30_for_PCA)[, 2], type = "p", pch = 16, 
     col = "blue", cex = 0.1, xlab = "Índice", ylab = "Valor", main = "Gráfico de Pontos com Linhas (R30)")
lines(t(reconstructed_r30)[, 2], type = "o", cex = 0.1, col = "red")