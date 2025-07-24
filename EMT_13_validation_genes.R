############################################
### 1. CARREGAMENTO DE DADOS E PACOTES ####
############################################

# Carrega dados prévios
load("~/mestrado/pca_result_R30.RData")
load("~/mestrado/t_matrix_R30.RData")
load("~/mestrado/t_matrix_R0.RData") # Para cálculo do SD

# Instala e carrega pacotes necessários
required_packages <- c("ggplot2", "ggrepel", "dplyr", "RColorBrewer", "patchwork", "biomaRt", "matrixStats")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

############################################
### 2. PREPARAÇÃO DOS DADOS DE PCA #########
############################################

# Extrai informações básicas
gene_names <- t_matrix_R30@transcriptogramS2[, 1]
gene_positions <- t_matrix_R30@transcriptogramS2[, 2]
pca_rotation <- pca_result_R30[["pca_result"]][["rotation"]]
rownames(pca_rotation) <- gene_names

# Filtra para as 6 primeiras PCs
n_pcs <- 6
top_percent <- 0.10

# Carrega genes extremos pré-identificados
top_extreme_table <- read.csv("genes_destaque_extremos_por_PC.csv", stringsAsFactors = FALSE) %>%
  filter(PC %in% paste0("PC", 1:n_pcs)) %>%
  mutate(Type = ifelse(Direction == "Up", "Peak", "Valley"))

############################################
### 3. CÁLCULO DO DESVIO PADRÃO (SD) ######
############################################

# Função para calcular SD de forma eficiente
calculate_gene_sd <- function() {
  # Extrai dados de expressão (ignora as 2 primeiras colunas)
  exp_data <- as.matrix(t_matrix_R0@transcriptogramS2[, -c(1:2)])
  
  data.frame(
    Gene = t_matrix_R0@transcriptogramS2[, 1],
    Position = t_matrix_R0@transcriptogramS2[, 2],
    SD = rowSds(exp_data, na.rm = TRUE) # Usando matrixStats para velocidade
  )
}

gene_sd <- calculate_gene_sd()

# Seleciona top 10% genes com maior SD
high_sd_genes <- gene_sd %>%
  arrange(desc(SD)) %>%
  slice_head(prop = top_percent) %>%
  mutate(Type = "High_SD")

############################################
### 4. ANOTAÇÃO DOS GENES ##################
############################################

# Converte IDs para nomes de genes
annotate_genes <- function(gene_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  getBM(
    attributes = c("ensembl_peptide_id", "external_gene_name"),
    filters = "ensembl_peptide_id",
    values = gene_ids,
    mart = mart
  )
}

gene_ids <- unique(c(top_extreme_table$GeneID, high_sd_genes$Gene))
gene_annotations <- annotate_genes(gene_ids)

# Adiciona nomes aos dados
top_extreme_table <- left_join(top_extreme_table, gene_annotations, 
                               by = c("GeneID" = "ensembl_peptide_id"))
high_sd_genes <- left_join(high_sd_genes, gene_annotations,
                           by = c("Gene" = "ensembl_peptide_id"))

# Salva genes com alto SD
write.csv(high_sd_genes %>% dplyr::select(GeneID = Gene, GeneName = external_gene_name, Position, SD),
          "top_high_SD_genes.csv", row.names = FALSE)

############################################
### 5. IDENTIFICAÇÃO DE GENES COMUNS #######
############################################

find_common_genes <- function() {
  inner_join(
    top_extreme_table %>% dplyr::select(GeneID, GeneName = external_gene_name, PC, Type, Loading),
    high_sd_genes %>% dplyr::select(Gene, GeneName = external_gene_name),
    by = c("GeneID" = "Gene", "GeneName")
  ) %>%
    left_join(gene_sd, by = c("GeneID" = "Gene")) %>%
    mutate(Common = TRUE)
}

common_genes <- find_common_genes()

# Salva genes comuns
write.csv(common_genes %>% dplyr::select(GeneID, GeneName, PC, Loading, SD),
          "genes_in_both_groups.csv", row.names = FALSE)

############################################
### 6. PREPARAÇÃO DO GRÁFICO PRINCIPAL ####
############################################

# Prepara dados de PCA para plotagem
prepare_pc_data <- function() {
  pc_means <- colMeans(pca_rotation[, 1:n_pcs])
  
  data.frame(
    Gene = rep(gene_names, n_pcs),
    Position = rep(gene_positions, n_pcs),
    PC = factor(rep(paste0("PC", 1:n_pcs), each = length(gene_names))),
    Loading = as.vector(pca_rotation[, 1:n_pcs]) - rep(pc_means, each = length(gene_names))
  )
}

pc_df <- prepare_pc_data()

# Configurações visuais
pc_colors <- brewer.pal(n_pcs, "Set1")
names(pc_colors) <- paste0("PC", 1:n_pcs)

############################################
### 7. CONSTRUÇÃO DO GRÁFICO ##############
############################################

# create_main_plot <- function() {
  # 1. Configurações iniciais ----
  y_range <- range(pc_df$Loading)
  y_limits <- y_range + c(-1, 1) * diff(y_range) * 0.05
  
  # 2. Preparação dos dados para labels ----
  label_genes <- bind_rows(
    top_extreme_table %>% 
      group_by(PC) %>% 
      slice_max(abs(Loading), n = 3),
    high_sd_genes %>% 
      arrange(desc(SD)) %>% 
      slice_head(n = 20),
    common_genes
  ) %>%
    distinct(GeneID, .keep_all = TRUE)
  
  # 3. Construção do gráfico base ----
  p <- ggplot() +
    # Linhas das PCs
    geom_line(
      data = pc_df, 
      aes(x = Position, y = Loading, color = PC), 
      alpha = 0.7, 
      linewidth = 0.5
    ) +
    
    # Linha de referência zero
    geom_hline(
      yintercept = 0, 
      linetype = "dashed", 
      color = "gray50"
    ) +
    
    # Pontos de SD (eixo secundário)
    geom_point(
      data = high_sd_genes,
      aes(x = Position, y = max(y_limits) * 1.05, size = SD),
      color = "darkgreen", 
      alpha = 0.3, 
      shape = 16
    ) +
    
    # Extremos de PCA
    geom_point(
      data = top_extreme_table,
      aes(x = Position, y = Loading, fill = PC),
      shape = ifelse(top_extreme_table$Type == "Peak", 24, 25),
      color = "black", 
      size = 2
    ) +
    
    # Genes comuns
    geom_point(
      data = common_genes,
      aes(x = Position, y = Loading),
      shape = 23, 
      fill = "gold", 
      color = "black", 
      size = 3
    ) +
    
    # Escalas e cores
    scale_color_manual(values = pc_colors) +
    scale_fill_manual(values = pc_colors) +
    scale_y_continuous(
      limits = y_limits,
      sec.axis = sec_axis(
        trans = ~./max(y_limits) * max(high_sd_genes$SD),
        name = "Standard Deviation"
      )
    ) +
    
    # Labels e tema
    labs(
      title = "Principal Components 1-6 with Highlighted Genes",
      subtitle = "Triangles: Top 5% PCA loadings (▲ = positive, ▼ = negative)\nGreen dots: Top 5% highly variable genes | Stars: Genes in both groups",
      x = "Gene Position",
      y = "PCA Rotation (centered)",
      color = "Principal Component",
      fill = "Principal Component",
      size = "Standard Deviation"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 20),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  # 4. Adicionar labels ----
  p_with_labels <- p +
    geom_text_repel(
      data = label_genes,
      aes(x = Position, y = Loading, label = GeneName),
      size = 3, 
      max.overlaps = Inf, 
      box.padding = 0.5,
      segment.color = "gray50"
    ) +
    
    geom_text_repel(
      data = common_genes,
      aes(x = Position, y = Loading, label = GeneName),
      size = 3.5, 
      fontface = "bold", 
      color = "darkorange",
      max.overlaps = Inf, 
      box.padding = 0.7
    )
  
  return(p_with_labels)
} # sem genes table
create_main_plot <- function() {
  # 1. Configurações iniciais
  y_range <- range(pc_df$Loading)
  y_limits <- y_range + c(-1, 1) * diff(y_range) * 0.05
  
  # 2. Preparação dos dados para labels
  label_genes <- bind_rows(
    top_extreme_table %>% 
      group_by(PC) %>% 
      slice_max(abs(Loading), n = 3),
    high_sd_genes %>% 
      arrange(desc(SD)) %>% 
      slice_head(n = 20),
    common_genes
  ) %>%
    distinct(GeneID, .keep_all = TRUE)
  
  # 3. Preparar tabela de genes comuns
  table_data <- common_genes %>%
    mutate(Direction = ifelse(Loading > 0, "Up", "Down")) %>%
    arrange(PC, desc(abs(Loading))) %>%
    dplyr::select(GeneName, Direction) %>%
    distinct(GeneName, .keep_all = TRUE)
  
  table_grob <- gridExtra::tableGrob(
    table_data,
    rows = NULL,
    cols = c("Gene", "Direction"),
    theme = gridExtra::ttheme_minimal(
      base_size = 8,
      padding = unit(c(2, 2), "mm"),
      core = list(
        fg_params = list(hjust = 0, x = 0.05),
        bg_params = list(fill = c("white", "#f7f7f7"))
      )
    )
  )
  
  # 4. Construção do gráfico base
  p <- ggplot() +
    # Linhas das PCs
    geom_line(
      data = pc_df, 
      aes(x = Position, y = Loading, color = PC), 
      alpha = 0.7, 
      linewidth = 0.5
    ) +
    
    # Linha de referência zero
    geom_hline(
      yintercept = 0, 
      linetype = "dashed", 
      color = "gray50"
    ) +
    
    # Pontos de SD
    geom_point(
      data = high_sd_genes,
      aes(x = Position, y = max(y_limits) * 1.05, size = SD),
      color = "darkgreen", 
      alpha = 0.3, 
      shape = 16
    ) +
    
    # Extremos de PCA
    geom_point(
      data = top_extreme_table,
      aes(x = Position, y = Loading, fill = PC),
      shape = ifelse(top_extreme_table$Type == "Peak", 24, 25),
      color = "black", 
      size = 2
    ) +
    
    # Genes comuns
    geom_point(
      data = common_genes,
      aes(x = Position, y = Loading),
      shape = 23, 
      fill = "gold", 
      color = "black", 
      size = 3
    ) +
    
    # Escalas e cores
    scale_color_manual(values = pc_colors) +
    scale_fill_manual(values = pc_colors) +
    scale_y_continuous(
      limits = y_limits,
      sec.axis = sec_axis(
        trans = ~./max(y_limits) * max(high_sd_genes$SD),
        name = "Standard Deviation"
      )
    ) +
    
    # Labels e tema
    labs(
      title = "Principal Components 1-6 with Highlighted Genes",
      subtitle = "Triangles: Top 10% PCA loadings (▲ = positive, ▼ = negative)\nGreen dots: Top 10% highly variable genes \nDiamond: Genes in both groups",
      x = "Gene Position",
      y = "PCA Rotation (centered)",
      color = "Principal Component",
      fill = "Principal Component",
      size = "Standard Deviation"
    ) +
    
    # Adiciona a tabela
    
    # Adiciona a tabela no canto inferior direito
    annotation_custom(
      grob = table_grob,
      xmin = max(pc_df$Position) - 0.05*diff(range(pc_df$Position)),  # Direita
      xmax = max(pc_df$Position),
      ymin = min(y_limits),
      ymax = min(y_limits) + 0.05*diff(range(y_limits))  # Ajuste de altura
    ) +
    
    # Título da tabela posicionado acima
    annotate(
      "text",
      x = max(pc_df$Position) - 0.125*diff(range(pc_df$Position)),  # Centralizado acima
      y = min(y_limits) + 0.4*diff(range(y_limits)),  # Imediatamente acima
      label = "Common Genes \n(High SD + PCA extremes)",
      hjust = 0.1,  # Centralizado
      vjust = 1,
      size = 4,
      fontface = "bold"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 20),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  # 5. Adicionar labels
  p_with_labels <- p +
    geom_text_repel(
      data = label_genes,
      aes(x = Position, y = Loading, label = GeneName),
      size = 3, 
      max.overlaps = Inf, 
      box.padding = 0.5,
      segment.color = "gray50"
    ) +
    
    geom_text_repel(
      data = common_genes,
      aes(x = Position, y = Loading, label = GeneName),
      size = 3.5, 
      fontface = "bold", 
      color = "darkorange",
      max.overlaps = Inf, 
      box.padding = 0.7
    )
  
  return(p_with_labels)
}


# Gera e salva o gráfico
final_plot <- create_main_plot()
ggsave("images/PCs_linhas_com_destaques.png", final_plot, width = 16, height = 10, dpi = 600, units = "in")

