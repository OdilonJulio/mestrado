library(transcriptogramer)
library(dplyr)
library(parallel)
library(progress)

# Função para carregar os dados necessários
load_required_data <- function() {
  # Lista de arquivos necessários
  required_files <- c(
    "reconstructed_matrix_R0.RData",
    "reconstructed_matrix_R30.RData",
    "t_matrix_R0.RData",
    "t_matrix_R30.RData",
    "pca_result_R0.RData",
    "pca_result_R30.RData"
  )
  
  # Verificar se todos os arquivos existem
  missing_files <- setdiff(required_files, list.files())
  if (length(missing_files) > 0) {
    stop("Arquivos necessários não encontrados: ", paste(missing_files, collapse = ", "))
  }
  
  # Carregar cada arquivo
  loaded_objects <- list()
  for (file in required_files) {
    load(file, envir = .GlobalEnv)
    loaded_objects <- c(loaded_objects, file)
    gc()
  }
  
  message("Arquivos carregados com sucesso: ", paste(loaded_objects, collapse = ", "))
}

# Função para preparar a matriz original
prepare_original_matrix <- function(t_matrix) {
  # Extrair a matriz de expressão (excluindo colunas não numéricas)
  expr_mat <- as.matrix(t_matrix@transcriptogramS2[, !names(t_matrix@transcriptogramS2) %in% c("Protein", "Position")])
  
  # Definir nomes das linhas como os IDs das proteínas
  rownames(expr_mat) <- t_matrix@transcriptogramS2$Protein
  
  # Transpor a matriz para ter genes nas linhas e células nas colunas
  t(expr_mat)
}

# Função para reconstruir usando n componentes principais
reconstruct_with_pcs <- function(pca_result, n_pcs) {
  # Extrair os componentes necessários
  rotation <- pca_result$rotation[, 1:n_pcs, drop = FALSE]
  components <- pca_result$x[, 1:n_pcs, drop = FALSE]
  center <- pca_result$center
  
  # Reconstruir a matriz
  reconstructed <- tcrossprod(components, rotation)
  
  # Adicionar o centro de volta, se existir
  if (!is.null(center)) {
    reconstructed <- sweep(reconstructed, 2, center, "+")
  }
  
  reconstructed
}

# Função para calcular R² entre matriz original e reconstruída
calculate_r_squared <- function(original, reconstructed) {
  # Calcular médias por gene (linha)
  gene_means <- rowMeans(original, na.rm = TRUE)
  
  # Substituir zeros pelo valor modal (para evitar divisão por zero)
  zero_genes <- gene_means == 0
  if (any(zero_genes)) {
    mode_val <- as.numeric(names(which.max(table(original)))) / 1000
    gene_means[zero_genes] <- mode_val
  }
  
  # Calcular erro quadrático médio normalizado
  squared_error <- (original - reconstructed)^2
  normalized_error <- squared_error / (gene_means^2)
  
  # Retornar métricas
  list(
    global_r2 = mean(normalized_error),
    gene_r2 = rowMeans(normalized_error),
    cell_r2 = colMeans(normalized_error),
    adjusted_zero_genes = sum(zero_genes),
    adjustment_value = ifelse(any(zero_genes), mode_val, NA)
  )
}

# Função principal para reconstruir com diferentes números de PCs
reconstruct_transcriptograms <- function(pca_result, original_matrix, label, max_pcs = 10) {
  results <- list()
  
  # Barra de progresso
  pb <- progress_bar$new(
    format = paste(label, "[:bar] :percent | Tempo restante: :eta"),
    total = max_pcs,
    clear = FALSE
  )
  
  for (n in 1:max_pcs) {
    # Reconstruir com 1 até n PCs
    recon_matrix <- reconstruct_with_pcs(pca_result, 1:n)
    
    # Calcular métricas de qualidade
    metrics <- calculate_r_squared(original_matrix, recon_matrix)
    
    # Armazenar resultados
    results[[paste0("PCs_1_to_", n)]] <- list(
      reconstructed_matrix = recon_matrix,
      metrics = metrics
    )
    
    pb$tick()
  }
  
  results
}

# Função principal
main <- function() {
  # 1. Carregar os dados
  message("\nCarregando dados...")
  load_required_data()
  
  # 2. Preparar matrizes originais
  message("\nPreparando matrizes originais...")
  original_R0 <- prepare_original_matrix(t_matrix_R0)
  original_R30 <- prepare_original_matrix(t_matrix_R30)
  
  # 3. Limpar objetos não necessários
  rm(t_matrix_R0, t_matrix_R30)
  gc()
  
  # 4. Reconstruir transcriptogramas para R0
  message("\nReconstruindo transcriptogramas para R0...")
  reconstructions_R0 <- reconstruct_transcriptograms(
    pca_result_R0$pca_result,
    original_R0,
    "R0",
    max_pcs = 10
  )
  
  # 5. Reconstruir transcriptogramas para R30
  message("\nReconstruindo transcriptogramas para R30...")
  reconstructions_R30 <- reconstruct_transcriptograms(
    pca_result_R30$pca_result,
    original_R30,
    "R30",
    max_pcs = 10
  )
  
  # 6. Salvar resultados
  message("\nSalvando resultados...")
  save(reconstructions_R0, file = "transcriptogram_reconstructions_R0_1_to_10_PCs.RData")
  save(reconstructions_R30, file = "transcriptogram_reconstructions_R30_1_to_10_PCs.RData")
  
  message("\nProcesso concluído com sucesso!")
  message("Arquivos salvos:")
  message("- transcriptogram_reconstructions_R0_1_to_10_PCs.RData")
  message("- transcriptogram_reconstructions_R30_1_to_10_PCs.RData")
}

# Executar a análise
main()
