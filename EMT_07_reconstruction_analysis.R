# reconstruction_analysis_final.R
library(ggplot2)
library(patchwork)

## 1. Funções Auxiliares -------------------------------------------------

# Função para extrair condições dos nomes das amostras
extract_condition <- function(sample_names, conditions) {
  sapply(sample_names, function(name) {
    matched <- grep(paste(conditions, collapse = "|"), name, value = TRUE)
    if (length(matched) > 0) sub("_.*", "", matched[1]) else "Unknown"
  })
}

# Função para calcular centróides
calculate_centroid <- function(pca_data) {
  if (!"Day" %in% colnames(pca_data)) {
    stop("A coluna 'Day' não existe no dataframe fornecido")
  }
  centroid <- aggregate(. ~ Day, data = pca_data, FUN = mean)
  return(centroid)
}

# Função para expandir centros de massa
expand_center_of_mass <- function(pca_data, centroid_data) {
  counts <- table(pca_data$Day)
  valid_days <- names(counts)[names(counts) %in% centroid_data$Day]
  row_indices <- rep(match(valid_days, centroid_data$Day), counts[valid_days])
  expanded_centroid <- centroid_data[row_indices, -1, drop = FALSE]
  return(t(expanded_centroid))
}

# Função para reconstruir transcriptogramas com verificação dimensional
reconstruct_transcriptogram <- function(pca_result, center_of_mass) {
  rotation <- pca_result$pca_result$rotation
  components <- pca_result$pca_result$x
  reconstructed <- components %*% t(rotation)
  return(reconstructed + center_of_mass)
}

## 2. Carregamento de Dados ----------------------------------------------

cat("=== CARREGAMENTO DE DADOS ===\n")

# Verificar arquivos necessários
required_files <- c("pca_result_R0.RData", "pca_result_R30.RData",
                    "t_matrix_R0.RData", "t_matrix_R30.RData")
missing_files <- setdiff(required_files, list.files())
if (length(missing_files)) {
  stop(paste("Arquivos faltantes:", paste(missing_files, collapse = ", ")))
}

# Carregar dados em ambientes separados
load_data <- function(file) {
  env <- new.env()
  load(file, envir = env)
  return(env)
}

env_pca_r0 <- load_data("pca_result_R0.RData")
env_pca_r30 <- load_data("pca_result_R30.RData")
env_t0 <- load_data("t_matrix_R0.RData")
env_t30 <- load_data("t_matrix_R30.RData")

# Extrair objetos
pca_result_R0 <- env_pca_r0$pca_result_R0
pca_result_R30 <- env_pca_r30$pca_result_R30
T_matrix_R0 <- env_t0$t_matrix_R0
T_matrix_R30 <- env_t30$t_matrix_R30

# Verificar dimensões
cat("\nDimensões das matrizes carregadas:\n")
cat("T_matrix_R0:", dim(T_matrix_R0), "\n")
cat("T_matrix_R30:", dim(T_matrix_R30), "\n")

## 3. Preparação dos Dados -----------------------------------------------

# Preparar dados PCA
prepare_pca_data <- function(pca_result, conditions) {
  pca_df <- as.data.frame(pca_result$pca_result$x)
  pca_df$Day <- extract_condition(rownames(pca_df), conditions)
  return(pca_df)
}

pca_r0 <- prepare_pca_data(pca_result_R0, c("notreated-batch1", "notreated-batch2"))
pca_r30 <- prepare_pca_data(pca_result_R30, c("notreated-batch1", "notreated-batch2"))

## 4. Processamento Principal --------------------------------------------

tryCatch({
  cat("\n=== PROCESSAMENTO INICIADO ===\n")
  
  # Calcular centróides
  centroid_r0 <- calculate_centroid(pca_r0)
  centroid_r30 <- calculate_centroid(pca_r30)
  
  # Expandir centros de massa
  center_of_mass_r0 <- expand_center_of_mass(pca_r0, centroid_r0)
  center_of_mass_r30 <- expand_center_of_mass(pca_r30, centroid_r30)
  
  cat("\nDimensões após expansão:\n")
  cat("center_of_mass_r0:", dim(center_of_mass_r0), "\n")
  cat("center_of_mass_r30:", dim(center_of_mass_r30), "\n")
  
  # Função segura para reconstrução
  safe_reconstruct <- function(pca_result, com_matrix, ref_matrix) {
    # Ajustar orientação da matriz de centros de massa
    if (ncol(com_matrix) != ncol(pca_result$pca_result$rotation)) {
      com_matrix <- t(com_matrix)
    }
    
    # Reconstruir
    reconstructed <- reconstruct_transcriptogram(pca_result, com_matrix)
    
    # Verificar e ajustar dimensões
    if (!identical(dim(reconstructed), dim(ref_matrix))) {
      if (identical(dim(reconstructed), rev(dim(ref_matrix)))) {
        reconstructed <- t(reconstructed)
      } else {
        stop("Dimensões incompatíveis após reconstrução")
      }
    }
    return(reconstructed)
  }
  
  # Reconstruir transcriptogramas
  cat("\nReconstruindo R0...\n")
  reconstructed_r0 <- safe_reconstruct(pca_result_R0, center_of_mass_r0, T_matrix_R0)
  
  cat("\nReconstruindo R30...\n")
  reconstructed_r30 <- safe_reconstruct(pca_result_R30, center_of_mass_r30, T_matrix_R30)
  
  # Verificação final
  cat("\n=== VERIFICAÇÃO FINAL ===\n")
  cat("reconstructed_r0:", dim(reconstructed_r0), "\n")
  cat("reconstructed_r30:", dim(reconstructed_r30), "\n")
  
  # Salvar resultados
  save(reconstructed_r0, reconstructed_r30, 
       file = "reconstructed_transcriptograms_final.RData")
  
  cat("\n=== PROCESSAMENTO CONCLUÍDO COM SUCESSO ===\n")
  
}, error = function(e) {
  cat("\n=== ERRO ===\n")
  cat(conditionMessage(e), "\n")
  
  # Diagnóstico
  cat("\nEstado atual:\n")
  if (exists("center_of_mass_r0")) cat("center_of_mass_r0:", dim(center_of_mass_r0), "\n")
  if (exists("center_of_mass_r30")) cat("center_of_mass_r30:", dim(center_of_mass_r30), "\n")
  if (exists("reconstructed_r0")) cat("reconstructed_r0:", dim(reconstructed_r0), "\n")
  if (exists("reconstructed_r30")) cat("reconstructed_r30:", dim(reconstructed_r30), "\n")
})