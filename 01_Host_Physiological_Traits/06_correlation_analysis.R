library(tidyverse)

## =========================================================
## SETTINGS
## =========================================================
file_L   <- "input_files/plasma_llactate.tsv"
file_D   <- "input_files/plasma_dlactate.tsv"
file_HSI <- "input_files/hepatosomatic_index.tsv"
file_GLU <- "input_files/plasma_glucose.tsv"

OUTPUT_DIR <- "results/correlation_analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

trait_cols <- c(
  "Plasma_D_Lactate",
  "Plasma_L_Lactate",
  "HSI",
  "Plasma_Glucose"
)

diet_cols <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)

## =========================================================
## READ AND HARMONIZE DATA
## =========================================================

L_dat <- read_tsv(file_L, show_col_types = FALSE) %>%
  transmute(
    sampleID_STPN2309,
    diet = factor(diet),
    day = factor(day),
    rearing_tank,
    Plasma_L_Lactate
  )

D_dat <- read_tsv(file_D, show_col_types = FALSE) %>%
  transmute(
    sampleID_STPN2309,
    diet = factor(diet),
    day = factor(day),
    rearing_tank,
    Plasma_D_Lactate
  )

HSI_dat <- read_tsv(file_HSI, show_col_types = FALSE) %>%
  rename(sampleID_STPN2309 = SampleID_STPN2309) %>%
  transmute(
    sampleID_STPN2309,
    diet = factor(diet),
    day = factor(day),
    rearing_tank,
    HSI
  )

GLU_dat <- read_tsv(file_GLU, show_col_types = FALSE) %>%
  rename(sampleID_STPN2309 = SampleID_STPN2309) %>%
  mutate(
    Plasma_Glucose = rowMeans(
      across(starts_with("replicate_")),
      na.rm = TRUE
    )
  ) %>%
  transmute(
    sampleID_STPN2309,
    diet = factor(diet),
    day = factor(day),
    rearing_tank,
    Plasma_Glucose
  )

## =========================================================
## MERGE DATASETS AND SELECT COMPLETE CASES
## =========================================================

id_cols <- c("sampleID_STPN2309", "diet", "day", "rearing_tank")

merged_dat <- L_dat %>%
  inner_join(D_dat, by = id_cols) %>%
  inner_join(HSI_dat, by = id_cols) %>%
  inner_join(GLU_dat, by = id_cols)

host_traits <- merged_dat %>%
  filter(diet != "Fasted") %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC"))
  ) %>%
  select(
    sampleID_STPN2309,
    diet,
    day,
    rearing_tank,
    all_of(trait_cols)
  ) %>%
  drop_na(all_of(trait_cols))

cat("Number of complete non-fasted samples:", nrow(host_traits), "\n")

write_tsv(
  host_traits,
  file.path(OUTPUT_DIR, "host_traits_complete_cases.tsv")
)

## =========================================================
## PAIRWISE SPEARMAN CORRELATIONS
## =========================================================

pair_list <- combn(trait_cols, 2, simplify = FALSE)

cor_results <- map_dfr(pair_list, function(pair) {
  
  pair_dat <- host_traits %>%
    select(all_of(pair)) %>%
    drop_na()
  
  test <- suppressWarnings(
    cor.test(
      pair_dat[[1]],
      pair_dat[[2]],
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    trait_1 = pair[1],
    trait_2 = pair[2],
    n = nrow(pair_dat),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance_raw = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    significance_BH = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_value)

cat("\nPairwise Spearman correlations:\n")
print(cor_results, n = Inf)

write_tsv(
  cor_results,
  file.path(OUTPUT_DIR, "host_trait_pairwise_spearman_correlations.tsv")
)

## =========================================================
## SAVE RHO AND P-VALUE MATRICES
## =========================================================

host_trait_matrix <- host_traits %>%
  select(all_of(trait_cols))

rho_mat <- cor(
  host_trait_matrix,
  method = "spearman",
  use = "complete.obs"
)

p_mat <- matrix(
  NA_real_,
  nrow = length(trait_cols),
  ncol = length(trait_cols),
  dimnames = list(trait_cols, trait_cols)
)

p_adj_BH_mat <- p_mat

diag(p_mat) <- 0
diag(p_adj_BH_mat) <- 0

for (i in seq_len(nrow(cor_results))) {
  
  trait_1 <- cor_results$trait_1[i]
  trait_2 <- cor_results$trait_2[i]
  
  p_mat[trait_1, trait_2] <- cor_results$p_value[i]
  p_mat[trait_2, trait_1] <- cor_results$p_value[i]
  
  p_adj_BH_mat[trait_1, trait_2] <- cor_results$p_adj_BH[i]
  p_adj_BH_mat[trait_2, trait_1] <- cor_results$p_adj_BH[i]
}

write_tsv(
  as.data.frame(rho_mat) %>% rownames_to_column("trait"),
  file.path(OUTPUT_DIR, "host_trait_spearman_rho_matrix.tsv")
)

write_tsv(
  as.data.frame(p_mat) %>% rownames_to_column("trait"),
  file.path(OUTPUT_DIR, "host_trait_spearman_pvalue_matrix.tsv")
)

write_tsv(
  as.data.frame(p_adj_BH_mat) %>% rownames_to_column("trait"),
  file.path(OUTPUT_DIR, "host_trait_spearman_BH_pvalue_matrix.tsv")
)

## =========================================================
## EXTRACT AND SAVE GLUCOSE VS HSI STATISTICS
## =========================================================

glu_hsi_stats <- cor_results %>%
  filter(
    trait_1 == "HSI",
    trait_2 == "Plasma_Glucose"
  )

cat("\nPlasma Glucose vs HSI:\n")
print(glu_hsi_stats)

write_tsv(
  glu_hsi_stats,
  file.path(OUTPUT_DIR, "plasma_glucose_vs_HSI_statistics.tsv")
)

## =========================================================
## PLOT: PLASMA GLUCOSE VS HSI ONLY
## =========================================================

plot_dat_glu_hsi <- host_traits %>%
  select(
    sampleID_STPN2309,
    diet,
    day,
    rearing_tank,
    Plasma_Glucose,
    HSI
  ) %>%
  drop_na()

write_tsv(
  plot_dat_glu_hsi,
  file.path(OUTPUT_DIR, "data_used_plasma_glucose_vs_HSI.tsv")
)

ct_glu_hsi <- cor.test(
  plot_dat_glu_hsi$Plasma_Glucose,
  plot_dat_glu_hsi$HSI,
  method = "spearman",
  exact = FALSE
)

rho_glu_hsi <- unname(ct_glu_hsi$estimate)
p_glu_hsi <- ct_glu_hsi$p.value

p_text <- if (p_glu_hsi < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", signif(p_glu_hsi, 2))
}

annotation_text <- paste0(
  "Spearman rho = ", sprintf("%.2f", rho_glu_hsi),
  "\n", p_text,
  "\nn = ", nrow(plot_dat_glu_hsi)
)

p_glu_hsi_plot <- ggplot(
  plot_dat_glu_hsi,
  aes(x = Plasma_Glucose, y = HSI)
) +
  geom_point(
    aes(color = diet),
    size = 3.5,
    alpha = 0.65
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "firebrick",
    fill = "grey75",
    linewidth = 1.1
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = annotation_text,
    hjust = -0.05,
    vjust = 1.2,
    size = 5
  ) +
  scale_color_manual(
    values = diet_cols,
    name = "Diet",
    breaks = c("NC", "HC")
  ) +
  labs(
    x = "Plasma Glucose (g/L)",
    y = "Hepato-Somatic Index (%)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major = element_line(color = "grey82", linewidth = 0.45),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 15)
  )

print(p_glu_hsi_plot)

ggsave(
  file.path(OUTPUT_DIR, "plasma_glucose_vs_HSI.pdf"),
  p_glu_hsi_plot,
  width = 5,
  height = 4.5
)

ggsave(
  file.path(OUTPUT_DIR, "plasma_glucose_vs_HSI.tiff"),
  p_glu_hsi_plot,
  width = 5,
  height = 4.5,
  dpi = 1200,
  compression = "lzw"
)

cat("\nSaved outputs to:", OUTPUT_DIR, "\n")
