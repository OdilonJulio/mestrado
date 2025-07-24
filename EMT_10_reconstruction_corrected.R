library(transcriptogramer)
library(ggplot2)
library(gridExtra)

# 1. Função para carregar dados
load_required_data <- function() {
  required_files <- c("t_matrix_R30.RData", "pca_result_R30.RData")
  missing_files <- setdiff(required_files, list.files())
  
  if (length(missing_files) > 0) {
    stop("Arquivos faltando: ", paste(missing_files, collapse = ", "))
  }
  
  load("t_matrix_R30.RData", envir = .GlobalEnv)
  load("pca_result_R30.RData", envir = .GlobalEnv)
  message("Dados carregados com sucesso")
}

# 2. Função para obter dados originais da célula
get_cell_data <- function(t_matrix, cell_name) {
  cell_columns <- setdiff(colnames(t_matrix@transcriptogramS2), c("Protein", "Position"))
  
  if (!cell_name %in% cell_columns) {
    stop("Célula não encontrada. Exemplos: ", paste(head(cell_columns), collapse = ", "))
  }
  
  data.frame(
    Position = t_matrix@transcriptogramS2$Position,
    Expression = t_matrix@transcriptogramS2[[cell_name]],
    Gene = rownames(t_matrix@transcriptogramS2)
  )
}

# 3. Função de reconstrução corrigida
reconstruct_expression <- function(pca_result, cell_name, n_pcs) {
  # Verificar se a célula existe
  if (!cell_name %in% rownames(pca_result$x)) {
    stop("Célula não encontrada no PCA")
  }
  
  # Obter componentes
  rotation <- pca_result$rotation[, 1:n_pcs, drop = FALSE]
  scores <- pca_result$x[cell_name, 1:n_pcs, drop = FALSE]
  
  # Reconstrução correta
  reconstructed <- scores %*% t(rotation)
  
  # Adicionar média se existir
  if (!is.null(pca_result$center)) {
    reconstructed <- reconstructed + pca_result$center
  }
  
  # Retornar como vetor nomeado
  setNames(as.numeric(reconstructed), colnames(reconstructed))
}

# 4. Função para plotagem
create_comparison_plot <- function(original, reconstructed, title) {
  df <- data.frame(
    Position = original$Position,
    Original = original$Expression,
    Reconstructed = reconstructed
  )
  
  ggplot(df, aes(x = Position)) +
    geom_line(aes(y = Original, color = "Original"), linewidth = 0.2, alpha = 0.8) +
    geom_line(aes(y = Reconstructed, color = "Reconstructed"), linewidth = 0.15, alpha = 0.6) +
    scale_color_manual(values = c("Original" = "black", "Reconstructed" = "red")) +
    labs(title = title, x = "Gene", y = "Expression", color = "") +
    theme_minimal(base_size = 10) +  # Fonte menor para mais espaço
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 0),
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      panel.grid.major = element_line(linewidth = 0.05),  # Grades muito sutis
      panel.grid.minor = element_blank()
    )
    
}

# 5. Função principal
analyze_cell <- function(cell_name, max_pcs = 5) {
  # Carregar dados
  load_required_data()
  
  # Obter dados originais
  original_data <- get_cell_data(t_matrix_R30, cell_name)
  
  # Preparar lista para resultados
  results <- list()
  plots <- list()
  
  # Processar cada número de PCs
  for (n in 1:max_pcs) {
    message("Processando ", n, " PC(s)...")
    
    # Reconstruir expressão
    reconstructed <- reconstruct_expression(pca_result_R30$pca_result, cell_name, n)
    
    # Verificar se a reconstrução funcionou
    if (length(reconstructed) != nrow(original_data)) {
      stop("Erro na reconstrução: tamanho incompatível")
    }
    
    # Criar plot
    plots[[n]] <- create_comparison_plot(
      original_data,
      reconstructed,
      paste("Using", n, ifelse(n == 1, "PC", "PCs"))
    )
    
    # Armazenar resultados
    results[[paste0("PC", n)]] <- reconstructed
  }
  
  # Criar gráfico combinado
  combined_plot <- grid.arrange(
    grobs = plots, 
    ncol = 2, 
    top = paste("Analysis for cell:", cell_name)
  )
  
  # Salvar o gráfico combinado
  output_dir <- "images"
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  filename <- paste0(output_dir, "/TranscriptogramReconstructed_", 
                     gsub("[^[:alnum:]]", "_", cell_name), "_", max_pcs, "PCs.png")
  
  ggsave(
    filename = filename,
    plot = combined_plot,
    width = 12, 
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  # Mostrar o gráfico na sessão R
  print(combined_plot)
  
  # Retornar resultados
  list(
    original = original_data,
    reconstructions = results,
    plots = plots,
    combined_plot = combined_plot
  )
}

# Verificar células disponíveis
load_required_data()
message("\nCélulas na matriz de expressão (10 primeiras):")
print(head(setdiff(colnames(t_matrix_R30@transcriptogramS2), c("Protein", "Position")), 10))

message("\nCélulas no PCA (10 primeiras):")
print(head(rownames(pca_result_R30$pca_result$x), 10))

# Executar análise (substitua pelo nome correto)

results <- analyze_cell("TGFbeta1-8day-batch1_CCATGTCAGTCATGCT-1", max_pcs = 6)
results <- analyze_cell("notreated-batch1_CTACACCCATGAACCT-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-8day-batch1_GCTTGAAGTCGCATCG-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-8day-batch1_CACTCCAAGCTGCAAG-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-1day-batch2_ACTTCCGTCGAGAACG-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-8day-batch1_AACACGTCAGATCGGA-1", max_pcs = 6)
results <- analyze_cell("notreated-batch2_GATGGAGGTTATGTGC-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-4day-batch1_CACACAACAGGAACGT-1", max_pcs = 6)
results <- analyze_cell("notreated-batch1_CAGCATACACCTCGTT-1", max_pcs = 6)
results <- analyze_cell("TGFbeta1-4day-batch1_GAAATGAGTAAGTGGC-1", max_pcs = 6)
