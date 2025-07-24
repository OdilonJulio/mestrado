# --- Requisitos ---
library(dplyr)
library(readr)

# --- Carregar objetos necessários ---
load("~/mestrado/pca_result_R30.RData")  # pca_result_R30
load("~/mestrado/t_matrix_R30.RData")    # contém os nomes dos genes

# --- Carregar dados de médias por intervalo (gerado no script "curva de médias") ---
mean_df <- read_csv("PC_means_by_interval_PC2_to_PC9.csv")  # Deve conter colunas: interval, PC1_mean, PCn_index, PCn_mean

# --- Configuração ---
rotation_matrix <- pca_result_R30$pca_result$rotation   # genes x PCs
gene_names <- t_matrix_R30@transcriptogramS2[, 1]       # nomes dos genes (posição original dos genes no transcriptograma)

# --- Identificar PCs presentes na média ---
pcs_presentes <- unique(mean_df$PCn_index)
num_pcs <- length(pcs_presentes)

# --- Construir matriz onde:
# - Primeira coluna: coeficientes de PC1 (PC1_mean)
# - Outras colunas: médias de coeficientes de PC2, PC3, ..., PCn em cada intervalo

# Pivot para formato largo
wide_matrix <- mean_df %>%
  dplyr::select(interval, PC1_mean, PCn_index, PCn_mean) %>%
  tidyr::pivot_wider(names_from = PCn_index, values_from = PCn_mean) %>%
  arrange(PC1_mean)

# Substituir NAs por zero
wide_matrix[is.na(wide_matrix)] <- 0


original_means <- pca_result_R30$pca_result$center  # Médias originais por gene

# --- Reconstrução dos transcriptogramas por intervalo ---
reconstructed_list <- lapply(1:nrow(wide_matrix), function(i) {
  row <- wide_matrix[i, ]
  pc1_coeff <- row$PC1_mean
  
  # Inicializa vetor com a contribuição da PC1
  reconstructed_vector <- pc1_coeff * rotation_matrix[, 1] + original_means  # Adiciona a média # PC1 é sempre a 1ª coluna
  
  # Para as demais PCs
  for (pc_index in pcs_presentes) {
    # Posição correta da PC no rotation matrix
    pc_number <- as.integer(gsub("PC", "", pc_index))
    pc_coeff <- row[[pc_index]]
    
    # Adiciona a contribuição da PCn
    reconstructed_vector <- reconstructed_vector + (pc_coeff * rotation_matrix[, pc_number])
  }
  
  data.frame(
    Gene = gene_names,
    Expression = as.numeric(reconstructed_vector),
    PC1_value = pc1_coeff,
    IntervalID = row$interval
  )
})

# --- Unir tudo em um único data.frame ---
reconstructed_from_means_df <- bind_rows(reconstructed_list)

# --- Salvar resultado ---
write.csv(reconstructed_from_means_df, "transcriptograms_reconstructed_by_PC1_weighted.csv", row.names = FALSE)

### MOVIE ###

library(scales)
library(ggplot2)

# --- Diretório para salvar figuras ---
output_dir <- "frames_pc1"
dir.create(output_dir, showWarnings = FALSE)

# --- Normalizar valores para melhor visualização (opcional) ---
# Você pode desativar essa parte se quiser valores reais
normalize_expression <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# --- Loop: gerar imagem para cada intervalo ---
unique_vals <- sort(unique(reconstructed_from_means_df$PC1_value))

for (val in unique_vals) {
  plot_df <- reconstructed_from_means_df %>%
    filter(PC1_value == val) %>%
    mutate(NormExpression = normalize_expression(Expression))
  
  p <- ggplot(plot_df, aes(x = seq_along(NormExpression), y = NormExpression)) +
    geom_line(color = "steelblue", size = 0.8) +
    labs(
      title = paste0("Transcriptogram (PC1 = ", signif(val, 4), ")"),
      x = "Gene Position",
      y = "Normalized Expression"
    ) +
    theme_minimal(base_size = 14) +
    ylim(0, 1)
  
  # Formatar nome do arquivo com padding (ordenação lexicográfica correta)
  filename <- sprintf("%s/pc1_%06d.png", output_dir, round(val * 100000))
  ggsave(filename, plot = p, width = 10, height = 4)
}

library(gifski)

# Lista todos os arquivos PNG no diretório, em ordem
png_files <- list.files("frames_pc1", pattern = ".*\\.png$", full.names = TRUE)

# Certifique-se de que estão na ordem correta (pelo nome)
png_files <- sort(png_files)

# Crie o GIF
gifski(png_files, 
       gif_file = "images/transcriptogram_animation.gif",
       width = 1000, 
       height = 400,
       delay = 0.1,  # tempo entre frames (em segundos)
       loop = TRUE)  # se deve repetir a animação




#### BIOLOGIC #####

library(magick)
library(ggplot2)
library(gganimate)
library(cowplot)

# 1. Carregar os GIFs existentes
gif_transcriptogram <- image_read("images/transcriptogram_animation.gif")
gif_biologic <- image_read("images/Biologic.jpg")  # Converta para GIF se necessário

# 2. Criar o gráfico de enriquecimento biológico como imagem
# Primeiro, gerar o gráfico ggplot
biological_data <- data.frame(
  pathway = c("TGF-beta", "Wnt", "Metabolism", "Cell Cycle", "ECM"),
  activation = c(0.8, 0.6, 0.4, 0.2, 0.5)
)

bio_plot <- ggplot(biological_data, aes(x = pathway, y = activation, fill = pathway)) +
  geom_col() +
  scale_fill_viridis_d() +
  labs(title = "Pathway Activation During EMT",
       x = "Biological Pathways",
       y = "Activation Level") +
  theme_minimal()

# Salvar como imagem temporária
temp_file <- tempfile(fileext = ".png")
ggsave(temp_file, bio_plot, width = 8, height = 4, dpi = 600)
bio_image <- image_read(temp_file)
file.remove(temp_file)

# 3. Redimensionar e combinar as imagens
combined <- image_composite(
  image_scale(gif_transcriptogram, "800x600"),
  image_scale(bio_image, "800x400"),
  offset = "+0+600"  # Posiciona abaixo do transcriptograma
)

# 4. Criar animação final
final_gif <- image_animate(combined, fps = 10) %>% 
  image_annotate("EMT Progression Analysis", 
                 size = 30, 
                 color = "black",
                 location = "+10+10")

# 5. Salvar o resultado
image_write(final_gif, "images/combined_emt_analysis.gif")
