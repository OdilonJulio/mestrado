# plot_PCs_reconstructed.R

library(dplyr)
library(ggplot2)
library(patchwork)

load_data <- function() {
  files <- c(
    "reconstructed_matrix_R0.RData", 
    "reconstructed_matrix_R30.RData",
    "t_matrix_R0.RData", 
    "t_matrix_R30.RData",
    "pca_result_R0.RData", 
    "pca_result_R30.RData",
    "transcriptogram_validation_results_parallel.RData"
  )
  
  loaded_files <- character(0)
  
  for (file in files) {
    if (!file.exists(file)) {
      warning("File not found: ", file)
      next
    }
    if (!grepl("\\.RData$", file, ignore.case = TRUE)) {
      warning("Skipping non-RData file: ", file)
      next
    }
    
    message("Loading ", file)
    load(file, envir = .GlobalEnv)
    loaded_files <- c(loaded_files, file)
    gc()
  }
  
  message("Successfully loaded ", length(loaded_files), " of ", length(files), " files")
  return(invisible(loaded_files))
}

load_data()

# Função para preparar dados com correção da cumulativa
prepare_percent_non_cumulative_r2 <- function(pca_res, r2_res, condition_label, max_pcs = NULL) {
  all_sdev <- pca_res$pca_result$sdev
  total_var <- sum(all_sdev^2)
  
  if (is.null(max_pcs)) {
    max_pcs <- length(all_sdev)
  }
  
  sdev <- all_sdev[1:max_pcs]
  
  # Corrigir início da curva cumulativa para alinhar com R²
  cumulative_var <- (cumsum(sdev^2) / total_var) * 100
  total_r2 <- sum(r2_res$r_squared[1:max_pcs])
  r2_percent <- (r2_res$r_squared[1:max_pcs] / total_r2) * 100
  
  # Ajuste para começar no ponto de R² da PC1
  cumulative_var <- cumulative_var - cumulative_var[1] + r2_percent[1]
  
  df_elbow <- data.frame(PC = 1:max_pcs, CumVar = cumulative_var)
  fit1 <- lm(CumVar ~ PC, df_elbow[1:2, ])
  fit2 <- lm(CumVar ~ PC, tail(df_elbow, 2))
  elbow_x <- (coef(fit2)[1] - coef(fit1)[1]) / (coef(fit1)[2] - coef(fit2)[2])
  elbow_y <- coef(fit1)[2] * elbow_x + coef(fit1)[1]
  
  list(
    data = data.frame(
      PC = 1:max_pcs,
      R2_percent = r2_percent,
      CumVar_percent = cumulative_var,
      Condition = condition_label
    ),
    elbow_x = elbow_x,
    elbow_y = elbow_y
  )
}


res_R0 <- prepare_percent_non_cumulative_r2(pca_result_R0, analysis_R0$results, "R0", max_pcs = 50)
res_R30 <- prepare_percent_non_cumulative_r2(pca_result_R30, analysis_R30$results, "R30", max_pcs = 50)
res_R0_all <- prepare_percent_non_cumulative_r2(pca_result_R0, analysis_R0$results, "R0")
res_R30_all <- prepare_percent_non_cumulative_r2(pca_result_R30, analysis_R30$results, "R30")

df_plot <- bind_rows(
  res_R0$data %>% mutate(Metric = "R² per PC"),
  res_R30$data %>% mutate(Metric = "R² per PC"),
  res_R0$data %>% mutate(Metric = "Cumulative Variance", R2_percent = CumVar_percent),
  res_R30$data %>% mutate(Metric = "Cumulative Variance", R2_percent = CumVar_percent)
)

df_plot_all <- bind_rows(
  res_R0_all$data %>% mutate(Metric = "R² per PC"),
  res_R30_all$data %>% mutate(Metric = "R² per PC"),
  res_R0_all$data %>% mutate(Metric = "Cumulative Variance", R2_percent = CumVar_percent),
  res_R30_all$data %>% mutate(Metric = "Cumulative Variance", R2_percent = CumVar_percent)
)

elbows <- data.frame(
  x = c(res_R0$elbow_x, res_R30$elbow_x),
  y = c(res_R0$elbow_y, res_R30$elbow_y),
  label = c(
    paste0("Elbow R0 (PC ≈ ", round(res_R0$elbow_x, 1), ")"),
    paste0("Elbow R30 (PC ≈ ", round(res_R30$elbow_x, 1), ")")
  )
)

elbows_all <- data.frame(
  x = c(res_R0_all$elbow_x, res_R30_all$elbow_x),
  y = c(res_R0_all$elbow_y, res_R30_all$elbow_y)
)

# ----- Plot 1: Até 50 PCs -----
gg_50 <- ggplot() +
  geom_line(data = df_plot %>% filter(Metric == "Cumulative Variance"),
            aes(x = PC, y = R2_percent, color = Condition), size = 1) +
  geom_point(data = df_plot %>% filter(Metric == "R² per PC"),
             aes(x = PC, y = R2_percent, color = Condition, shape = Condition), size = 2) +
  geom_line(data = df_plot %>% filter(Metric == "R² per PC"),
            aes(x = PC, y = R2_percent, color = Condition), linetype = "dotted") +
  geom_point(data = elbows, aes(x = x, y = y), shape = 8, color = "black", size = 3) +
  geom_text(data = elbows, aes(x = x, y = y, label = label),
            vjust = -1, size = 3, color = "black") +
  labs(
    title = "Cumulative Variance (%) + R² per PC (%) — First 50 PCs",
    x = "Number of Principal Components",
    y = "Cumulative Variance (%) / R² per PC (%)",
    color = "Condition",
    shape = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(limits = c(0, 100))

# ----- Plot 2: Todas as PCs -----
gg_all <- ggplot() +
  geom_line(data = df_plot_all %>% filter(Metric == "Cumulative Variance"),
            aes(x = PC, y = R2_percent, color = Condition), size = 1) +
  geom_point(data = df_plot_all %>% filter(Metric == "R² per PC"),
             aes(x = PC, y = R2_percent, color = Condition, shape = Condition), size = 1.5) +
  geom_line(data = df_plot_all %>% filter(Metric == "R² per PC"),
            aes(x = PC, y = R2_percent, color = Condition), linetype = "dotted") +
  geom_point(data = elbows_all, aes(x = x, y = y), shape = 8, color = "black", size = 3) +
  labs(
    title = "Cumulative Variance (%) + R² per PC (%) — All PCs",
    x = "Number of Principal Components",
    y = "Cumulative Variance (%) / R² per PC (%)",
    caption = paste0(
      "Elbow R0 ≈ PC ", round(res_R0_all$elbow_x, 1),
      " | Elbow R30 ≈ PC ", round(res_R30_all$elbow_x, 1)
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 10, face = "italic", hjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 100))

# ----- Combinar lado a lado e salvar -----
final_plot <- gg_50 + gg_all + plot_layout(ncol = 2)

ggsave("images/r2_vs_cumulative_percent_elbow_combined.png", final_plot, width = 14, height = 6, dpi = 600)
print(final_plot)
