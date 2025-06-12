# transcriptogram_reconstruction_validation_optimized.R
# Versão otimizada para alto desempenho

library(transcriptogramer)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(Matrix)

# 1. Configurar paralelização
registerDoParallel(cores = detectCores() - 2)

# 2. Carregar dados necessários
load("reconstructed_matrix_R0.RData")
load("reconstructed_matrix_R30.RData")
load("t_matrix_R0.RData")
load("t_matrix_R30.RData")
load("pca_result_R0.RData")
load("pca_result_R30.RData")

# 3. Funções otimizadas
prepare_original_matrix <- function(transcriptogram_obj) {
  df <- transcriptogram_obj@transcriptogramS2
  rownames(df) <- df$Protein
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  t(as.matrix(df))
}

calculate_r_squared <- function(original, reconstructed) {
  gene_means <- rowMeans(original, na.rm = TRUE)
  mode_value <- median(original) / 1000
  zero_means <- gene_means == 0
  gene_means[zero_means] <- mode_value
  
  squared_diff <- (original - reconstructed)^2
  normalized_diff <- squared_diff / (gene_means^2)
  
  list(
    global = mean(normalized_diff, na.rm = TRUE),
    by_gene = rowMeans(normalized_diff, na.rm = TRUE),
    by_cell = colMeans(normalized_diff, na.rm = TRUE),
    adjusted_zero_genes = sum(zero_means),
    adjustment_value = mode_value
  )
}

reconstruct_with_n_pcs <- function(pca_result, original_matrix, n_pcs) {
  rotation <- pca_result$pca_result$rotation[, 1:n_pcs, drop = FALSE]
  components <- pca_result$pca_result$x[, 1:n_pcs, drop = FALSE]
  center <- pca_result$pca_result$center
  
  reconstructed <- components %*% t(rotation)
  
  if (!is.null(center)) {
    reconstructed <- reconstructed + rep(center, each = nrow(components))
  }
  
  t(reconstructed)
}

analyze_pc_components <- function(pca_result, original_matrix, radius_label) {
  max_pcs <- ncol(pca_result$pca_result$x)
  
  results <- foreach(n = 1:max_pcs, .combine = rbind) %dopar% {
    reconstructed <- reconstruct_with_n_pcs(pca_result, original_matrix, n)
    r2 <- calculate_r_squared(original_matrix, reconstructed)
    data.frame(n_pcs = n, r_squared = r2$global)
  }
  
  results <- results[order(results$n_pcs), ]
  
  plot <- ggplot(results, aes(x = n_pcs, y = r_squared)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue") +
    labs(title = paste("Erro de reconstrução por número de PCs -", radius_label),
         x = "Número de Componentes Principais",
         y = "Erro R²",
         caption = paste("Análise para", radius_label)) +
    theme_minimal()
  
  list(results = results, plot = plot)
}

# 4. Preparar matrizes originais (convertidas para formato esparso)
original_matrix_R0 <- as(prepare_original_matrix(t_matrix_R0), "dgCMatrix")
original_matrix_R30 <- as(prepare_original_matrix(t_matrix_R30), "dgCMatrix")

# 5. Executar análises em paralelo
analysis_results <- foreach(r = list(list(pca_result_R0, original_matrix_R0, "R0"),
                           list(pca_result_R30, original_matrix_R30, "R30")),
                          .combine = list) %dopar% {
  analyze_pc_components(r[[1]], r[[2]], r[[3]])
}

analysis_R0 <- analysis_results[[1]]
analysis_R30 <- analysis_results[[2]]

# 6. Cálculo final otimizado
calculate_full_recon <- function(original, reconstructed) {
  chunk_size <- 1000
  n_genes <- nrow(original)
  chunks <- split(1:n_genes, ceiling(seq_along(1:n_genes)/chunk_size))
  
  results <- foreach(chunk = chunks, .combine = c) %dopar% {
    calculate_r_squared(original[chunk, ], reconstructed[chunk, ])$global
  }
  
  list(global = mean(results))
}

full_recon_r2_R0 <- calculate_full_recon(original_matrix_R0, reconstructed_matrix_R0)
full_recon_r2_R30 <- calculate_full_recon(original_matrix_R30, reconstructed_matrix_R30)

# 7. Resultados e gráficos (mesmo código original)
cat("=== Resultados para R0 ===\n")
cat("Erro R² global (todos PCs):", full_recon_r2_R0$global, "\n")
print(head(analysis_R0$results))

cat("\n=== Resultados para R30 ===\n")
cat("Erro R² global (todos PCs):", full_recon_r2_R30$global, "\n")
print(head(analysis_R30$results))

print(analysis_R0$plot)
print(analysis_R30$plot)

# 8. Salvar resultados
save(analysis_R0, analysis_R30, full_recon_r2_R0, full_recon_r2_R30,
     file = "transcriptogram_reconstruction_validation.RData")

ggsave("reconstruction_error_R0.png", analysis_R0$plot, width = 8, height = 6)
ggsave("reconstruction_error_R30.png", analysis_R30$plot, width = 8, height = 6)

# 9. Limpar memória
rm(list = ls()[!ls() %in% c("analysis_R0", "analysis_R30", "full_recon_r2_R0", "full_recon_r2_R30")])
gc(full = TRUE)
