# quality_control.R
# Realiza o controle de qualidade (QC) e filtra os dados.

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

# Carregar objetos Seurat
load("seurat_objects.RData")

# Identificar doublets e adicionar porcentagem de genes mitocondriais
seurat_objects <- lapply(seurat_objects, function(seurat) {
  sce <- as.SingleCellExperiment(seurat)
  sce <- scDblFinder(sce)
  seurat$doublet <- sce$scDblFinder.class
  seurat <- subset(seurat, subset = doublet == "singlet")
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  return(seurat)
})

# Aplicar controle de qualidade para cada objeto Seurat individual
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  subset(seurat, subset = nFeature_RNA > 500 & percent.mt < 20)
})

# Combinar os objetos filtrados
filtered_combined_seurat <- merge(
  filtered_seurat_objects[[1]],
  y = filtered_seurat_objects[-1],
  add.cell.ids = sample_names,
  project = "FilteredCombined"
)

# Salvar objeto combinado
save(filtered_combined_seurat, file = "filtered_combined_seurat.RData")