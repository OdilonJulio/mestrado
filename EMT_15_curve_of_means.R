library(tidyr)
library(dplyr)
library(ggplot2)

# Load PCA result
load("~/mestrado/pca_result_R30.RData")

# --- Configuration ---
num_pcs <- 8         # Number of PCs beyond PC1 (PC2 to PC(num_pcs + 1))
num_intervals <- 100  # Number of intervals along PC1

# --- Extract PCA matrix ---
pc_matrix <- pca_result_R30$pca_result$x

# --- Define intervals based on PC1 range ---
x_breaks <- seq(min(pc_matrix[, 1]), max(pc_matrix[, 1]), length.out = num_intervals + 1)
interval_levels <- levels(cut(pc_matrix[, 1], breaks = x_breaks, include.lowest = TRUE))

# --- Function to generate PC1 vs PCn summary data ---
generate_pc1_vs_pcn_data <- function(n) {
  data <- data.frame(PC1 = pc_matrix[, 1], PCn = pc_matrix[, n])
  data$interval <- cut(data$PC1, breaks = x_breaks, include.lowest = TRUE)
  
  grouped <- data %>%
    group_by(interval) %>%
    summarise(
      PC1_mean = mean(PC1, na.rm = TRUE),
      PCn_mean = mean(PCn, na.rm = TRUE),
      .groups = "drop"
    )
  
  grouped$PCn_index <- paste0("PC", n)
  return(grouped)
}

# --- Generate data for all PCs ---
all_paths_data <- lapply(2:(num_pcs + 1), generate_pc1_vs_pcn_data)
# Reorder factor by numeric PC index
all_paths_df <- bind_rows(all_paths_data) %>%
  filter(!is.na(PCn_index))  # Remove facetas indesejadas

all_paths_df$PCn_index <- factor(
  all_paths_df$PCn_index,
  levels = paste0("PC", sort(2:(num_pcs + 1)))
)


# --- Determine global Y-axis limits based on PC with largest variation ---
amplitude_by_pc <- all_paths_df %>%
  group_by(PCn_index) %>%
  summarise(
    ymin = min(PCn_mean, na.rm = TRUE),
    ymax = max(PCn_mean, na.rm = TRUE),
    amplitude = ymax - ymin
  ) %>%
  arrange(desc(amplitude))

ymin <- amplitude_by_pc$ymin[1]
ymax <- amplitude_by_pc$ymax[1]
max_pc <- amplitude_by_pc$PCn_index[1]

# --- Plot with smooth interpolation ---
plot_paths <- ggplot(all_paths_df, aes(x = PC1_mean, y = PCn_mean)) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue", span = 0.2, na.rm = TRUE) +
  facet_wrap(~PCn_index, scales = "fixed", ncol = 4) +
  labs(
    title = paste0("PC1 vs PC2 to PC", num_pcs + 1,
                   " — Intervals: ", num_intervals,
                   "\nY-axis scaled by largest variation: ", max_pc),
    x = "PC1 (mean per interval)",
    y = "PCn (mean per interval)"
  ) +
  ylim(ymin, ymax) +
  theme_minimal(base_size = 13)

print(plot_paths)

# --- Save high-resolution plot ---
ggsave(
  filename = paste0("images/PC1_vs_PC2_to_PC", num_pcs + 1, "_smoothed_Yscaled_", max_pc, ".png"),
  plot = plot_paths,
  width = 20,
  height = ceiling(num_pcs / 4) * 3.5,
  dpi = 300,
  units = "in",
  limitsize = FALSE
)

write.csv(
  all_paths_df,
  file = paste0("PC_means_by_interval_PC2_to_PC", num_pcs + 1, ".csv"),
  row.names = FALSE
)


# Instale se necessário
# install.packages("gganimate")
library(gganimate)

# Criar a animação
anim <- ggplot(all_paths_df, aes(x = PC1_mean, y = PCn_mean)) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue", span = 0.2, na.rm = TRUE) +
  labs(
    title = paste0("PC1 vs ", "{closest_state} — Intervals: ", num_intervals),
    x = "PC1 (mean per interval)",
    y = "PCn (mean per interval)"
  ) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_minimal(base_size = 14) +
  transition_states(PCn_index, transition_length = 2, state_length = 1) +
  ease_aes('cubic-in-out')

# Exportar como vídeo MP4 (requer ffmpeg) ou como GIF
animate(anim, renderer = gifski_renderer(), width = 1000, height = 700, duration = 10, fps = 2)
anim_save("images/PC1_vs_PCn_animation.gif", animation = last_animation())
