# reconstruction_metrics_R2.R
# Implementação exata da fórmula LaTeX para cálculo do R²

# Função que calcula o R² conforme definição matemática
calculate_transcriptogram_R2 <- function(T_matrix, T_reconstructed) {
  # Verificação de dimensões
  if (!identical(dim(T_matrix), dim(T_reconstructed))) {
    stop("As matrizes devem ter dimensões idênticas")
  }
  
  # Número de células e genes (notação igual ao LaTeX)
  N_células <- nrow(T_matrix)
  N_genes <- ncol(T_matrix)
  
  # Cálculo de T_med,i (média por gene - Eq. 2 do LaTeX)
  T_med_i <- colMeans(T_matrix)
  
  # Proteção contra divisão por zero (limiar 1e-10)
  T_med_i_sq <- ifelse(T_med_i^2 == 0, 1e-10, T_med_i^2)
  
  # Cálculo do erro quadrático normalizado (Eq. 1 do LaTeX)
  squared_error <- (T_matrix - T_reconstructed)^2
  normalized_error <- sweep(squared_error, 2, T_med_i_sq, "/")
  
  # Cálculo final do R²
  R2 <- sum(normalized_error) / (N_células * N_genes)
  
  return(R2)
}

# -------------------------------------------------------------------------
# Adaptação para seu pipeline específico (com seus nomes de objetos)
# -------------------------------------------------------------------------

# 1. Preparação dos dados (igual ao seu script anterior)
load("t_matrix_R0.RData")
load("t_matrix_R30.RData")
load("reconstructed_transcriptograms_final.RData")

# 2. Preparar matrizes (funções do seu script)
prepare_original_matrix <- function(transcriptogram_obj) {
  df <- as.data.frame(transcriptogram_obj@transcriptogramS2)
  mat <- as.matrix(df[, -c(1, 2)])  # Remove colunas não-numéricas
  t_mat <- t(mat)
  colnames(t_mat) <- df$Protein
  return(t_mat)
}

prepare_reconstructed_matrix <- function(reconstructed_matrix, original_cell_names) {
  t_reconstructed <- t(reconstructed_matrix)
  rownames(t_reconstructed) <- original_cell_names
  return(t_reconstructed)
}

# 3. Aplicação às suas matrizes específicas
T_matrix_R0 <- prepare_original_matrix(t_matrix_R0)
T_reconstructed_R0 <- prepare_reconstructed_matrix(reconstructed_r0, rownames(T_matrix_R0))

T_matrix_R30 <- prepare_original_matrix(t_matrix_R30)
T_reconstructed_R30 <- prepare_reconstructed_matrix(reconstructed_r30, rownames(T_matrix_R30))

# 4. Cálculo do R² (fiel ao LaTeX)
R2_R0 <- calculate_transcriptogram_R2(T_matrix_R0, T_reconstructed_R0)
R2_R30 <- calculate_transcriptogram_R2(T_matrix_R30, T_reconstructed_R30)

# 5. Saída formatada
cat("-------------------------------------\n")
cat("Resultados do R² (definição LaTeX):\n")
cat("-------------------------------------\n")
cat(sprintf("Raio 0:  R² = %.6f\n", R2_R0))
cat(sprintf("Raio 30: R² = %.6f\n", R2_R30))
cat("-------------------------------------\n")
cat("Nota: Valores mais próximos de 0 indicam melhor reconstrução.\n")

# 6. Salvar resultados (opcional)
metrics <- list(
  R0 = list(R2 = R2_R0),
  R30 = list(R2 = R2_R30)
)
save(metrics, file = "R2_metrics_LaTeX_definition.RData")