############################
## LOAD LIBRARIES
############################
library(microeco)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(data.table)
library(eulerr)

set.seed(123)
############################
## OUTPUT DIRECTORIES
############################
dir.create("03_16S_diversity_analysis", showWarnings = FALSE, recursive = TRUE)

############################
## GLOBAL OPTIONS
############################
FIG_WIDTH  <- 6
FIG_HEIGHT <- 4
FIG_DPI    <- 1200

diet_colors <- c(
  NC = "#44a7c4",
  HC = "#f99943"
)

############################
## HELPER FUNCTIONS
############################

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("ns")
}

style_axes <- function(p) {
  p +
    theme(
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.key.size = unit(1.0, "cm"),
      legend.text = ggtext::element_markdown(size = 12),
      legend.title = element_text(size = 13)
    )
}

style_pcoa <- function(p) {
  p +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      fill  = guide_legend(override.aes = list(size = 4))
    )
}

save_plot_dual <- function(plot_obj, filename_base,
                           width = FIG_WIDTH,
                           height = FIG_HEIGHT,
                           dpi = FIG_DPI,
                           outdir = "03_16S_diversity_analysis") {
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".tiff")),
    plot = plot_obj,
    device = "tiff",
    dpi = dpi,
    width = width,
    height = height,
    units = "in",
    compression = "lzw"
  )
  
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    device = "pdf",
    dpi = dpi,
    width = width,
    height = height,
    units = "in"
  )
}

make_permanova_label <- function(df) {
  p_val <- df$`Pr(>F)`[1]
  r2    <- df$R2[1]
  
  paste0(
    "PERMANOVA\n",
    "p = ", format.pval(p_val, digits = 3, eps = 0.001),
    "\nR² = ", round(r2, 3)
  )
}

standardize_alpha_data_stat <- function(df, group_col, dataset_subset, dataset_name) {
  df %>%
    mutate(
      Group = .data[[group_col]],
      Group_Type = group_col,
      dataset_subset = dataset_subset,
      dataset_name = dataset_name
    ) %>%
    select(Group, Group_Type, dataset_subset, everything(), -all_of(group_col))
}

standardize_alpha_diff <- function(df, dataset_subset, dataset_name, group_var) {
  df %>%
    mutate(
      dataset_subset = dataset_subset,
      dataset_name   = dataset_name,
      group_var      = group_var
    )
}

standardize_permanova <- function(df, dataset_subset, distance, group_var) {
  as.data.frame(df) %>%
    mutate(
      dataset_subset = dataset_subset,
      distance       = distance,
      group_var      = group_var
    )
}

create_beta_plot <- function(
    beta_obj,
    permanova_df,
    title,
    plot_color,
    colors = NULL,
    loading_arrow = TRUE,
    point_size = 3,
    plot_shape = NULL,
    plot_type = "point"
) {
  label_txt <- make_permanova_label(permanova_df)
  
  p <- beta_obj$plot_ordination(
    plot_color = plot_color,
    plot_shape = plot_shape,
    loading_arrow = loading_arrow,
    choices = c(1, 2),
    point_size = point_size,
    plot_type = plot_type
  ) +
    labs(title = title) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label_txt,
      hjust = 1.05,
      vjust = 1.15,
      size = 3.5,
      fontface = "bold"
    )
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  style_pcoa(p)
}


#########################################################
## GLOBAL MID-GUT ALPHA DIVERSITY
#########################################################

load("RData/16s_dataset_rarefied.RData")
mt_rarefied$sample_sums() %>% range()

t_alpha_global <- trans_alpha$new(
  dataset = mt_rarefied,
  group = "sample_type"
)

global_alpha_data_stat <- t_alpha_global$data_stat %>%
  standardize_alpha_data_stat(
    group_col = "sample_type",
    dataset_subset = "global",
    dataset_name = "global_midgut"
  )

write.csv(
  global_alpha_data_stat,
  file = file.path("03_16S_diversity_analysis", "global_alpha_div_data_stat.csv"),
  row.names = FALSE
)

t_alpha_global$cal_diff(
  measure = c("Observed", "Pielou", "Shannon", "PD"),
  method = "wilcox",
  p_adjust_method = "BH",
  return_model = TRUE
)

global_alpha_diff_results <- t_alpha_global$res_diff %>%
  standardize_alpha_diff(
    dataset_subset = "global",
    dataset_name = "global_midgut",
    group_var = "sample_type"
  )

write.csv(
  global_alpha_diff_results,
  file = file.path("03_16S_diversity_analysis", "global_alpha_div_diff_results.csv"),
  row.names = FALSE
)

global_observed_plot <- t_alpha_global$plot_alpha(
  measure = "Observed",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
)

global_observed_plot <- style_axes(global_observed_plot)

save_plot_dual(
  plot_obj = global_observed_plot,
  filename_base = "global_alpha_div_observed"
)

global_shannon_plot <- t_alpha_global$plot_alpha(
  measure = "Shannon",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
)

global_shannon_plot <- style_axes(global_shannon_plot)

save_plot_dual(
  plot_obj = global_shannon_plot,
  filename_base = "global_alpha_div_shannon"
)

############################
## GLOBAL - JACCARD PCoA
############################

t_beta_global_jaccard <- trans_beta$new(
  dataset = mt_rarefied,
  measure = "jaccard",
  group = "sample_type"
)

t_beta_global_jaccard$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_global_jaccard$cal_manova(manova_all = TRUE)

global_jaccard_permanova <- standardize_permanova(
  df = t_beta_global_jaccard$res_manova,
  dataset_subset = "global",
  distance = "jaccard",
  group_var = "sample_type"
)

write.csv(
  global_jaccard_permanova,
  file = file.path("03_16S_diversity_analysis", "global_jaccard_permanova_results.csv"),
  row.names = FALSE
)

global_jaccard_plot <- create_beta_plot(
  beta_obj = t_beta_global_jaccard,
  permanova_df = global_jaccard_permanova,
  title = "Jaccard distance (PCoA)",
  plot_color = "sample_type",
  colors = NULL,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = global_jaccard_plot,
  filename_base = "global_jaccard_PCoA"
)

############################
## GLOBAL - BRAY PCoA
############################

t_beta_global_bray <- trans_beta$new(
  dataset = mt_rarefied,
  measure = "bray",
  group = "sample_type"
)

t_beta_global_bray$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_global_bray$cal_manova(manova_all = TRUE)

global_bray_permanova <- standardize_permanova(
  df = t_beta_global_bray$res_manova,
  dataset_subset = "global",
  distance = "bray",
  group_var = "sample_type"
)

write.csv(
  global_bray_permanova,
  file = file.path("03_16S_diversity_analysis", "global_bray_permanova_results.csv"),
  row.names = FALSE
)

global_bray_plot <- create_beta_plot(
  beta_obj = t_beta_global_bray,
  permanova_df = global_bray_permanova,
  title = "Bray-Curtis distance (PCoA)",
  plot_color = "sample_type",
  colors = NULL,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = global_bray_plot,
  filename_base = "global_bray_PCoA"
)



#########################################################
## DIGESTA ALPHA DIVERSITY
#########################################################

load("RData/16s_dataset_dig_rarefied.RData")
mt_dig_rarefied$sample_sums() %>% range()

t_alpha_digesta <- trans_alpha$new(
  dataset = mt_dig_rarefied,
  group = "Diet"
)

digesta_alpha_data_stat <- t_alpha_digesta$data_stat %>%
  standardize_alpha_data_stat(
    group_col = "Diet",
    dataset_subset = "digesta",
    dataset_name = "digesta"
  )

write.csv(
  digesta_alpha_data_stat,
  file = file.path("03_16S_diversity_analysis", "digesta_alpha_div_data_stat.csv"),
  row.names = FALSE
)

t_alpha_digesta$cal_diff(
  measure = c("Observed", "Pielou", "Shannon", "PD"),
  method = "wilcox",
  p_adjust_method = "BH",
  return_model = TRUE
)

digesta_alpha_diff_results <- t_alpha_digesta$res_diff %>%
  standardize_alpha_diff(
    dataset_subset = "digesta",
    dataset_name = "digesta",
    group_var = "Diet"
  )

write.csv(
  digesta_alpha_diff_results,
  file = file.path("03_16S_diversity_analysis", "digesta_alpha_div_diff_results.csv"),
  row.names = FALSE
)

digesta_observed_plot <- t_alpha_digesta$plot_alpha(
  measure = "Observed",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
) +
  scale_color_manual(values = diet_colors)

digesta_observed_plot <- style_axes(digesta_observed_plot)

save_plot_dual(
  plot_obj = digesta_observed_plot,
  filename_base = "digesta_alpha_div_observed"
)

digesta_shannon_plot <- t_alpha_digesta$plot_alpha(
  measure = "Shannon",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
) +
  scale_color_manual(values = diet_colors)

digesta_shannon_plot <- style_axes(digesta_shannon_plot)

save_plot_dual(
  plot_obj = digesta_shannon_plot,
  filename_base = "digesta_alpha_div_shannon"
)

############################
## DIGESTA - JACCARD PCoA
############################

t_beta_digesta_jaccard <- trans_beta$new(
  dataset = mt_dig_rarefied,
  measure = "jaccard",
  group = "Diet"
)

t_beta_digesta_jaccard$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_digesta_jaccard$cal_manova(manova_all = TRUE)

digesta_jaccard_permanova <- standardize_permanova(
  df = t_beta_digesta_jaccard$res_manova,
  dataset_subset = "digesta",
  distance = "jaccard",
  group_var = "Diet"
)

write.csv(
  digesta_jaccard_permanova,
  file = file.path("03_16S_diversity_analysis", "digesta_jaccard_permanova_results.csv"),
  row.names = FALSE
)

digesta_jaccard_plot <- create_beta_plot(
  beta_obj = t_beta_digesta_jaccard,
  permanova_df = digesta_jaccard_permanova,
  title = "Jaccard distance (PCoA)",
  plot_color = "Diet",
  colors = diet_colors,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = digesta_jaccard_plot,
  filename_base = "digesta_jaccard_PCoA"
)

############################
## DIGESTA - BRAY PCoA
############################

t_beta_digesta_bray <- trans_beta$new(
  dataset = mt_dig_rarefied,
  measure = "bray",
  group = "Diet"
)

t_beta_digesta_bray$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_digesta_bray$cal_manova(manova_all = TRUE)

digesta_bray_permanova <- standardize_permanova(
  df = t_beta_digesta_bray$res_manova,
  dataset_subset = "digesta",
  distance = "bray",
  group_var = "Diet"
)

write.csv(
  digesta_bray_permanova,
  file = file.path("03_16S_diversity_analysis", "digesta_bray_permanova_results.csv"),
  row.names = FALSE
)

digesta_bray_plot <- create_beta_plot(
  beta_obj = t_beta_digesta_bray,
  permanova_df = digesta_bray_permanova,
  title = "Bray-Curtis distance (PCoA)",
  plot_color = "Diet",
  colors = diet_colors,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = digesta_bray_plot,
  filename_base = "digesta_bray_PCoA"
)

#########################################################
## MUCUS ALPHA DIVERSITY
#########################################################

load("RData/16s_dataset_mucus_rarefied.RData")
mt_mucus_rarefied$sample_sums() %>% range()

t_alpha_mucus <- trans_alpha$new(
  dataset = mt_mucus_rarefied,
  group = "Diet"
)

mucus_alpha_data_stat <- t_alpha_mucus$data_stat %>%
  standardize_alpha_data_stat(
    group_col = "Diet",
    dataset_subset = "mucus",
    dataset_name = "mucus"
  )

write.csv(
  mucus_alpha_data_stat,
  file = file.path("03_16S_diversity_analysis", "mucus_alpha_div_data_stat.csv"),
  row.names = FALSE
)

t_alpha_mucus$cal_diff(
  measure = c("Observed", "Pielou", "Shannon", "PD"),
  method = "wilcox",
  p_adjust_method = "BH",
  return_model = TRUE
)

mucus_alpha_diff_results <- t_alpha_mucus$res_diff %>%
  standardize_alpha_diff(
    dataset_subset = "mucus",
    dataset_name = "mucus",
    group_var = "Diet"
  )

write.csv(
  mucus_alpha_diff_results,
  file = file.path("03_16S_diversity_analysis", "mucus_alpha_div_diff_results.csv"),
  row.names = FALSE
)

mucus_observed_plot <- t_alpha_mucus$plot_alpha(
  measure = "Observed",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
) +
  scale_color_manual(values = diet_colors)

mucus_observed_plot <- style_axes(mucus_observed_plot)

save_plot_dual(
  plot_obj = mucus_observed_plot,
  filename_base = "mucus_alpha_div_observed"
)

mucus_shannon_plot <- t_alpha_mucus$plot_alpha(
  measure = "Shannon",
  xtext_angle = 0,
  xtext_size = 15,
  add = "jitter",
  add_sig_text_size = 15,
  y_start = -0.1,
  y_increase = -0.1,
  order_x_mean = FALSE
) +
  scale_color_manual(values = diet_colors)

mucus_shannon_plot <- style_axes(mucus_shannon_plot)

save_plot_dual(
  plot_obj = mucus_shannon_plot,
  filename_base = "mucus_alpha_div_shannon"
)

############################
## MUCUS - JACCARD PCoA
############################

t_beta_mucus_jaccard <- trans_beta$new(
  dataset = mt_mucus_rarefied,
  measure = "jaccard",
  group = "Diet"
)

t_beta_mucus_jaccard$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_mucus_jaccard$cal_manova(manova_all = TRUE)

mucus_jaccard_permanova <- standardize_permanova(
  df = t_beta_mucus_jaccard$res_manova,
  dataset_subset = "mucus",
  distance = "jaccard",
  group_var = "Diet"
)

write.csv(
  mucus_jaccard_permanova,
  file = file.path("03_16S_diversity_analysis", "mucus_jaccard_permanova_results.csv"),
  row.names = FALSE
)

mucus_jaccard_plot <- create_beta_plot(
  beta_obj = t_beta_mucus_jaccard,
  permanova_df = mucus_jaccard_permanova,
  title = "Jaccard distance (PCoA)",
  plot_color = "Diet",
  colors = diet_colors,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = mucus_jaccard_plot,
  filename_base = "mucus_jaccard_PCoA"
)

############################
## MUCUS - BRAY PCoA
############################

t_beta_mucus_bray <- trans_beta$new(
  dataset = mt_mucus_rarefied,
  measure = "bray",
  group = "Diet"
)

t_beta_mucus_bray$cal_ordination(method = "PCoA", ncomp = 2)
t_beta_mucus_bray$cal_manova(manova_all = TRUE)

mucus_bray_permanova <- standardize_permanova(
  df = t_beta_mucus_bray$res_manova,
  dataset_subset = "mucus",
  distance = "bray",
  group_var = "Diet"
)

write.csv(
  mucus_bray_permanova,
  file = file.path("03_16S_diversity_analysis", "mucus_bray_permanova_results.csv"),
  row.names = FALSE
)

mucus_bray_plot <- create_beta_plot(
  beta_obj = t_beta_mucus_bray,
  permanova_df = mucus_bray_permanova,
  title = "Bray-Curtis distance (PCoA)",
  plot_color = "Diet",
  colors = diet_colors,
  loading_arrow = TRUE,
  point_size = 3
)

save_plot_dual(
  plot_obj = mucus_bray_plot,
  filename_base = "mucus_bray_PCoA"
)

#########################################################
## MERGED CSV TABLES
#########################################################

all_alpha_data_stat <- bind_rows(
  global_alpha_data_stat,
  digesta_alpha_data_stat,
  mucus_alpha_data_stat
)

write.csv(
  all_alpha_data_stat,
  file = file.path("03_16S_diversity_analysis", "all_alpha_div_data_stat.csv"),
  row.names = FALSE
)

all_alpha_diff_results <- bind_rows(
  global_alpha_diff_results,
  digesta_alpha_diff_results,
  mucus_alpha_diff_results
)

write.csv(
  all_alpha_diff_results,
  file = file.path("03_16S_diversity_analysis", "all_alpha_div_diff_results.csv"),
  row.names = FALSE
)

all_beta_permanova_results <- bind_rows(
  global_jaccard_permanova,
  global_bray_permanova,
  digesta_jaccard_permanova,
  digesta_bray_permanova,
  mucus_jaccard_permanova,
  mucus_bray_permanova
)

write.csv(
  all_beta_permanova_results,
  file = file.path("03_16S_diversity_analysis", "all_beta_permanova_results.csv"),
  row.names = FALSE
)

############################
## SESSION INFO
############################

writeLines(
  capture.output(sessioninfo::session_info()),
  con = file.path("03_16S_diversity_analysis", "session_info.txt")
)

#########################################################
## COMBINED 3-COLUMN FIGURE
## Column 1 = Sample type : Shannon / Jaccard / Bray-Curtis
## Column 2 = Digesta     : Shannon / Jaccard / Bray-Curtis
## Column 3 = Mucus       : Shannon / Jaccard / Bray-Curtis
#########################################################


# =========================================================
# USER-CONTROLLED PARAMETERS
# =========================================================

# ---- General export ----
COMBO_WIDTH  <- 17
COMBO_HEIGHT <- 16
COMBO_DPI    <- 1200

# ---- Column titles ----
COLUMN_TITLE_SIZE <- 40
COLUMN_TITLE_FACE <- "bold"

# ---- Alpha diversity plot controls ----
ALPHA_AXIS_TEXT_SIZE    <- 30
ALPHA_AXIS_TITLE_SIZE   <- 30
ALPHA_LEGEND_TEXT_SIZE  <- 25
ALPHA_LEGEND_TITLE_SIZE <- 20
ALPHA_SIG_TEXT_SIZE     <- 20
ALPHA_POINT_SIZE        <- 3

# ---- PCoA plot controls ----
PCOA_TITLE_SIZE         <- 24
PCOA_AXIS_TEXT_SIZE     <- 25
PCOA_AXIS_TITLE_SIZE    <- 25
PCOA_LEGEND_TEXT_SIZE   <- 25
PCOA_LEGEND_TITLE_SIZE  <- 25
PCOA_POINT_SIZE         <- 2

# ---- PERMANOVA annotation inside PCoA ----
PCOA_PERMANOVA_TEXT_SIZE <- 8
PCOA_PERMANOVA_HJUST     <- 0.5
PCOA_PERMANOVA_VJUST     <- 1.5
PCOA_PERMANOVA_X         <- Inf
PCOA_PERMANOVA_Y         <- Inf

# ---- Relative heights inside each column ----
COLUMN_HEIGHTS <- c(0.12, 1, 1, 1)  # title, shannon, jaccard, bray

# =========================================================
# HELPER: column title as a standalone patchwork element
# =========================================================
make_column_title <- function(title_text,
                              size = COLUMN_TITLE_SIZE,
                              face = COLUMN_TITLE_FACE) {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label = title_text,
      gp = grid::gpar(fontsize = size, fontface = face)
    )
  )
}

# =========================================================
# HELPER: custom alpha styling
# =========================================================
style_alpha_custom <- function(p,
                               axis_text_size = ALPHA_AXIS_TEXT_SIZE,
                               axis_title_size = ALPHA_AXIS_TITLE_SIZE,
                               legend_text_size = ALPHA_LEGEND_TEXT_SIZE,
                               legend_title_size = ALPHA_LEGEND_TITLE_SIZE) {
  p +
    theme(
      axis.text.x  = element_text(size = axis_text_size, colour = "black"),
      axis.text.y  = element_text(size = axis_text_size, colour = "black"),
      axis.title.x = element_text(size = axis_title_size, colour = "black"),
      axis.title.y = element_text(size = axis_title_size, colour = "black"),
      legend.text  = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      plot.title   = element_text(hjust = 0.5, face = "bold")
    )
}

# =========================================================
# HELPER: custom PCoA styling
# =========================================================
style_pcoa_custom <- function(p,
                              title_size = PCOA_TITLE_SIZE,
                              axis_text_size = PCOA_AXIS_TEXT_SIZE,
                              axis_title_size = PCOA_AXIS_TITLE_SIZE,
                              legend_text_size = PCOA_LEGEND_TEXT_SIZE,
                              legend_title_size = PCOA_LEGEND_TITLE_SIZE) {
  p +
    theme(
      plot.title   = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.text.x  = element_text(size = axis_text_size),
      axis.text.y  = element_text(size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size),
      legend.text  = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      fill  = guide_legend(override.aes = list(size = 4))
    )
}

# =========================================================
# HELPER: make alpha plot with full control
# =========================================================
make_alpha_panel <- function(alpha_obj,
                             measure = "Shannon",
                             add_color_scale = NULL,
                             axis_text_size = ALPHA_AXIS_TEXT_SIZE,
                             axis_title_size = ALPHA_AXIS_TITLE_SIZE,
                             legend_text_size = ALPHA_LEGEND_TEXT_SIZE,
                             legend_title_size = ALPHA_LEGEND_TITLE_SIZE,
                             sig_text_size = ALPHA_SIG_TEXT_SIZE,
                             y_start = -0.1,
                             y_increase = -0.1) {
  
  p <- alpha_obj$plot_alpha(
    measure = measure,
    xtext_angle = 0,
    xtext_size = axis_text_size,
    add = "jitter",
    add_sig_text_size = sig_text_size,
    y_start = y_start,
    y_increase = y_increase,
    order_x_mean = FALSE
  )
  
  if (!is.null(add_color_scale)) {
    p <- p + add_color_scale
  }
  
  p <- style_alpha_custom(
    p,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    legend_text_size = legend_text_size,
    legend_title_size = legend_title_size
  )
  
  p
}

# =========================================================
# HELPER: make beta plot with PERMANOVA text + title control
# =========================================================
make_beta_panel <- function(beta_obj,
                            permanova_df,
                            plot_color,
                            title = NULL,
                            colors = NULL,
                            plot_shape = NULL,
                            loading_arrow = TRUE,
                            point_size = PCOA_POINT_SIZE,
                            plot_type = "point",
                            permanova_text_size = PCOA_PERMANOVA_TEXT_SIZE,
                            axis_text_size = PCOA_AXIS_TEXT_SIZE,
                            axis_title_size = PCOA_AXIS_TITLE_SIZE,
                            legend_text_size = PCOA_LEGEND_TEXT_SIZE,
                            legend_title_size = PCOA_LEGEND_TITLE_SIZE,
                            title_size = PCOA_TITLE_SIZE,
                            permanova_x_npc = 0.03,
                            permanova_y_npc = 0.97) {
  
  label_txt <- make_permanova_label(permanova_df)
  
  p <- beta_obj$plot_ordination(
    plot_color = plot_color,
    plot_shape = plot_shape,
    loading_arrow = loading_arrow,
    choices = c(1, 2),
    point_size = point_size,
    plot_type = plot_type
  )
  
  # ---- Extract plotted point coordinates from ggplot object ----
  pb <- ggplot_build(p)
  
  # Usually the first layer contains the sample points
  point_data <- pb$data[[1]]
  
  x_rng <- range(point_data$x, na.rm = TRUE)
  y_rng <- range(point_data$y, na.rm = TRUE)
  
  # ---- Convert npc-like relative positions to real data coordinates ----
  # permanova_x_npc: 0 = far left, 1 = far right
  # permanova_y_npc: 0 = bottom,   1 = top
  permanova_x <- x_rng[1] + permanova_x_npc * diff(x_rng)
  permanova_y <- y_rng[1] + permanova_y_npc * diff(y_rng)
  
  p <- p +
    annotate(
      "text",
      x = permanova_x,
      y = permanova_y,
      label = label_txt,
      hjust = 0,
      vjust = 1,
      size = permanova_text_size,
      fontface = "bold"
    )
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  p <- p +
    theme(
      plot.margin = ggplot2::margin(t = 12, r = 8, b = 8, l = 8)
    )
  
  p <- style_pcoa_custom(
    p,
    title_size = title_size,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    legend_text_size = legend_text_size,
    legend_title_size = legend_title_size
  )
  
  p
}

# =========================================================
# REBUILD THE 9 PLOTS
# =========================================================

# ---- Global column ----
p_global_shannon <- make_alpha_panel(
  alpha_obj = t_alpha_global,
  measure = "Shannon",
  add_color_scale = NULL
)

p_global_jaccard <- make_beta_panel(
  beta_obj = t_beta_global_jaccard,
  permanova_df = global_jaccard_permanova,
  plot_color = "sample_type",
  title = "Jaccard",
  colors = NULL
)

p_global_bray <- make_beta_panel(
  beta_obj = t_beta_global_bray,
  permanova_df = global_bray_permanova,
  plot_color = "sample_type",
  title = "Bray-Curtis",
  colors = NULL
)

# ---- Digesta column ----
p_digesta_shannon <- make_alpha_panel(
  alpha_obj = t_alpha_digesta,
  measure = "Shannon",
  add_color_scale = scale_color_manual(values = diet_colors)
)

p_digesta_jaccard <- make_beta_panel(
  beta_obj = t_beta_digesta_jaccard,
  permanova_df = digesta_jaccard_permanova,
  plot_color = "Diet",
  title = "Jaccard",
  colors = diet_colors
)

p_digesta_bray <- make_beta_panel(
  beta_obj = t_beta_digesta_bray,
  permanova_df = digesta_bray_permanova,
  plot_color = "Diet",
  title = "Bray-Curtis",
  colors = diet_colors
)

# ---- Mucus column ----
p_mucus_shannon <- make_alpha_panel(
  alpha_obj = t_alpha_mucus,
  measure = "Shannon",
  add_color_scale = scale_color_manual(values = diet_colors)
)

p_mucus_jaccard <- make_beta_panel(
  beta_obj = t_beta_mucus_jaccard,
  permanova_df = mucus_jaccard_permanova,
  plot_color = "Diet",
  title = "Jaccard",
  colors = diet_colors
)

p_mucus_bray <- make_beta_panel(
  beta_obj = t_beta_mucus_bray,
  permanova_df = mucus_bray_permanova,
  plot_color = "Diet",
  title = "Bray-Curtis",
  colors = diet_colors
)

# =========================================================
# BUILD EACH COLUMN
# =========================================================

col_global <- (
  make_column_title("Sample type") /
    p_global_shannon /
    p_global_jaccard /
    p_global_bray
) + plot_layout(heights = COLUMN_HEIGHTS)

col_digesta <- (
  make_column_title("Digesta") /
    p_digesta_shannon /
    p_digesta_jaccard /
    p_digesta_bray
) + plot_layout(heights = COLUMN_HEIGHTS)

col_mucus <- (
  make_column_title("Mucus") /
    p_mucus_shannon /
    p_mucus_jaccard /
    p_mucus_bray
) + plot_layout(heights = COLUMN_HEIGHTS)

# =========================================================
# FINAL COMBINED FIGURE
# =========================================================

combined_3col_figure <- col_global | col_digesta | col_mucus

combined_3col_figure <- combined_3col_figure +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

# Show in R
combined_3col_figure

# =========================================================
# SAVE
# =========================================================

ggsave(
  filename = file.path("03_16S_diversity_analysis", "combined_shannon_jaccard_bray_3columns.tiff"),
  plot = combined_3col_figure,
  device = "tiff",
  dpi = COMBO_DPI,
  width = COMBO_WIDTH,
  height = COMBO_HEIGHT,
  units = "in",
  compression = "lzw"
)

ggsave(
  filename = file.path("03_16S_diversity_analysis", "combined_shannon_jaccard_bray_3columns.pdf"),
  plot = combined_3col_figure,
  device = "pdf",
  dpi = COMBO_DPI,
  width = COMBO_WIDTH,
  height = COMBO_HEIGHT,
  units = "in"
)





## dbRDA: DIGESTA + MUCUS

library(microeco)
library(ggplot2)

#
# Settings
# 
RDS_FILE <- "rds/metadata_for_omics_data.rds"

DBRDA_OUTDIR <- file.path(
  "03_16S_diversity_analysis",
  "dbRDA"
)

dir.create(DBRDA_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Choose: "jaccard" for presence/absence or "bray" for abundance
DISTANCE <- "jaccard"

# missing environmental values are imputed by microeco/mice.
COMPLETE_NA <- TRUE

# Recommended because traits are on very different scales
STANDARDIZE_ENV <- TRUE

ENV_VARS <- c(
  "plasma_glucose",
  "HSI",
  "body_weight",
  "plasma_DLactate",
  "plasma_LLactate"
)

DIET_COLORS <- c(
  "NC" = "#56B4E9",
  "HC" = "#E69F00"
)

set.seed(123)


# -----------------------------
# Load new metadata RDS
# -----------------------------

omics_data <- readRDS(RDS_FILE)

stopifnot(
  all(c("metadata", "host_traits") %in% names(omics_data))
)

metadata <- as.data.frame(omics_data$metadata)
host_traits <- as.data.frame(omics_data$host_traits)

stopifnot(
  !is.null(rownames(metadata)),
  !is.null(rownames(host_traits))
)


# -----------------------------
# Helper functions
# -----------------------------

to_numeric_safe <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  suppressWarnings(as.numeric(x))
}


make_dbrda_sample_table <- function(
    microtable_object,
    metadata,
    host_traits,
    microbiome_id_column,
    env_vars
) {
  
  # Match metadata and host traits using S_001, S_002, ...
  shared_fish_ids <- intersect(
    rownames(metadata),
    rownames(host_traits)
  )
  
  if (length(shared_fish_ids) < 3) {
    stop("Too few matched S_### identifiers between metadata and host_traits.")
  }
  
  meta_matched <- metadata[
    shared_fish_ids,
    ,
    drop = FALSE
  ]
  
  traits_matched <- host_traits[
    shared_fish_ids,
    ,
    drop = FALSE
  ]
  
  # Keep metadata values when present.
  # Fill missing values using host_traits.
  for (trait in env_vars) {
    
    if (!trait %in% colnames(meta_matched) &&
        !trait %in% colnames(traits_matched)) {
      stop("Environmental variable not found: ", trait)
    }
    
    meta_values <- if (trait %in% colnames(meta_matched)) {
      to_numeric_safe(meta_matched[[trait]])
    } else {
      rep(NA_real_, nrow(meta_matched))
    }
    
    trait_values <- if (trait %in% colnames(traits_matched)) {
      to_numeric_safe(traits_matched[[trait]])
    } else {
      rep(NA_real_, nrow(meta_matched))
    }
    
    meta_matched[[trait]] <- ifelse(
      is.na(meta_values),
      trait_values,
      meta_values
    )
  }
  
  # Convert grouping variables to factors
  factor_vars <- intersect(
    c(
      "diet",
      "day_categ",
      "rearing_tank",
      "group",
      "phase",
      "diet_phase"
    ),
    colnames(meta_matched)
  )
  
  for (var in factor_vars) {
    meta_matched[[var]] <- droplevels(
      as.factor(meta_matched[[var]])
    )
  }
  
  if (!microbiome_id_column %in% colnames(meta_matched)) {
    stop(
      "Column not found in metadata: ",
      microbiome_id_column
    )
  }
  
  # Remove fish without a microbiome sample ID
  microbiome_ids <- trimws(
    as.character(meta_matched[[microbiome_id_column]])
  )
  
  keep <- !is.na(microbiome_ids) & nzchar(microbiome_ids)
  
  meta_matched <- meta_matched[
    keep,
    ,
    drop = FALSE
  ]
  
  microbiome_ids <- microbiome_ids[keep]
  
  if (anyDuplicated(microbiome_ids)) {
    duplicated_ids <- unique(
      microbiome_ids[duplicated(microbiome_ids)]
    )
    
    stop(
      "Duplicated microbiome sample IDs found: ",
      paste(head(duplicated_ids, 10), collapse = ", ")
    )
  }
  
  # Change S_### row names to the actual microbiome sample IDs
  rownames(meta_matched) <- microbiome_ids
  
  # Obtain the real sample order from the OTU table
  otu_sample_ids <- colnames(microtable_object$otu_table)
  
  if (is.null(otu_sample_ids)) {
    stop("No sample IDs found in colnames(microtable_object$otu_table).")
  }
  
  missing_metadata <- setdiff(
    otu_sample_ids,
    rownames(meta_matched)
  )
  
  if (length(missing_metadata) > 0) {
    stop(
      "These microbiome samples have no matched metadata in ",
      microbiome_id_column,
      ":\n",
      paste(head(missing_metadata, 20), collapse = "\n")
    )
  }
  
  # Reorder metadata exactly to OTU-table sample order
  sample_table <- meta_matched[
    otu_sample_ids,
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(sample_table),
      otu_sample_ids
    )
  )
  
  sample_table
}


run_dbrda <- function(
    microtable_object,
    location_name,
    microbiome_id_column,
    metadata,
    host_traits,
    env_vars,
    distance,
    output_dir,
    complete_na = TRUE,
    standardize_env = TRUE
) {
  
  # Do not overwrite the original rarefied object
  mt_dbrda <- microtable_object$clone(deep = TRUE)
  
  # Build correctly matched sample metadata
  sample_table_dbrda <- make_dbrda_sample_table(
    microtable_object = mt_dbrda,
    metadata = metadata,
    host_traits = host_traits,
    microbiome_id_column = microbiome_id_column,
    env_vars = env_vars
  )
  
  mt_dbrda$sample_table <- sample_table_dbrda
  mt_dbrda$tidy_dataset()
  
  stopifnot(
    identical(
      rownames(mt_dbrda$sample_table),
      colnames(mt_dbrda$otu_table)
    )
  )
  
  message(
    location_name,
    ": ",
    nrow(mt_dbrda$sample_table),
    " samples matched successfully."
  )
  
  # Calculate the distance only if it is not already present
  if (is.null(mt_dbrda$beta_diversity) ||
      !distance %in% names(mt_dbrda$beta_diversity)) {
    
    mt_dbrda$cal_betadiv(measure = distance)
  }
  
  # dbRDA
  t_env <- trans_env$new(
    dataset = mt_dbrda,
    env_cols = env_vars,
    standardize = standardize_env,
    complete_na = complete_na
  )
  
  t_env$cal_ordination(
    method = "dbRDA",
    use_measure = distance
  )
  
  # Permutation tests
  t_env$cal_ordination_anova(
    permutations = 999
  )
  
  t_env$cal_ordination_envfit(
    permutations = 999
  )
  
  # Adjust arrow lengths for display only
  t_env$trans_ordination(
    adjust_arrow_length = TRUE,
    max_perc_env = 0.8
  )
  
  distance_title <- switch(
    distance,
    "jaccard" = "Jaccard",
    "bray" = "Bray-Curtis",
    distance
  )
  
  p <- t_env$plot_ordination(
    plot_color = "diet",
    env_text_color = "black",
    point_size = 5,
    env_text_size = 5
  ) +
    scale_color_manual(values = DIET_COLORS) +
    labs(
      title = paste0(
        "dbRDA of ",
        distance_title,
        " dissimilarities: ",
        location_name
      ),
      color = "Diet"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 6)
      )
    )
  
  print(p)
  
  file_prefix <- file.path(
    output_dir,
    paste0(
      "dbRDA_",
      tolower(location_name),
      "_",
      distance
    )
  )
  
  # Save plot
  ggsave(
    filename = paste0(file_prefix, ".tiff"),
    plot = p,
    width = 7,
    height = 5,
    units = "in",
    dpi = 1200,
    compression = "lzw"
  )
  
  ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = p,
    width = 7,
    height = 5,
    units = "in"
  )
  
  # Save exact metadata and statistical outputs
  write.csv(
    data.frame(
      sample_id = rownames(mt_dbrda$sample_table),
      mt_dbrda$sample_table,
      check.names = FALSE
    ),
    file = paste0(file_prefix, "_sample_metadata.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      sample_id = rownames(t_env$data_env),
      t_env$data_env,
      check.names = FALSE
    ),
    file = paste0(file_prefix, "_environmental_data_used.csv"),
    row.names = FALSE
  )
  
  capture.output(
    t_env$res_ordination_R2,
    file = paste0(file_prefix, "_R2.txt")
  )
  
  capture.output(
    t_env$res_ordination_terms,
    file = paste0(file_prefix, "_anova_terms.txt")
  )
  
  capture.output(
    t_env$res_ordination_axis,
    file = paste0(file_prefix, "_anova_axes.txt")
  )
  
  capture.output(
    t_env$res_ordination_envfit,
    file = paste0(file_prefix, "_envfit.txt")
  )
  
  saveRDS(
    list(
      microtable = mt_dbrda,
      trans_env = t_env,
      sample_table = mt_dbrda$sample_table,
      environmental_data = t_env$data_env
    ),
    file = paste0(file_prefix, "_results.rds")
  )
  
  invisible(
    list(
      microtable = mt_dbrda,
      trans_env = t_env,
      plot = p
    )
  )
}


# Run dbRDA: Digesta

dbrda_digesta <- run_dbrda(
  microtable_object = mt_dig_rarefied,
  location_name = "Digesta",
  microbiome_id_column = "SampleID_digesta_microbiota_PGTB",
  metadata = metadata,
  host_traits = host_traits,
  env_vars = ENV_VARS,
  distance = DISTANCE,
  output_dir = DBRDA_OUTDIR,
  complete_na = COMPLETE_NA,
  standardize_env = STANDARDIZE_ENV
)


# Run dbRDA: Mucus

dbrda_mucus <- run_dbrda(
  microtable_object = mt_mucus_rarefied,
  location_name = "Mucus",
  microbiome_id_column = "SampleID_mucus_microbiota_PGTB",
  metadata = metadata,
  host_traits = host_traits,
  env_vars = ENV_VARS,
  distance = DISTANCE,
  output_dir = DBRDA_OUTDIR,
  complete_na = COMPLETE_NA,
  standardize_env = STANDARDIZE_ENV
)





# Euler diagram
ANALYSIS_OUTDIR <- "03_16S_diversity_analysis"
EULER_OUTDIR <- file.path(ANALYSIS_OUTDIR, "Euler_diagram")

dir.create(ANALYSIS_OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(EULER_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Helper: define group order consistently
get_groups_2groups <- function(meta, group_var) {
  observed_groups <- unique(as.character(meta[[group_var]]))
  
  if (length(observed_groups) != 2) {
    stop(paste0("'", group_var, "' must contain exactly 2 groups. Found: ",
                paste(observed_groups, collapse = ", ")))
  }
  
  if (group_var == "Diet") {
    groups <- c("NC", "HC")
  } else if (group_var == "sample_type") {
    groups <- c("Digesta", "Mucus")
  } else {
    groups <- sort(observed_groups)
  }
  
  if (!all(groups %in% observed_groups)) {
    stop(paste0(
      "Expected groups not found in '", group_var, "'. Found: ",
      paste(observed_groups, collapse = ", ")
    ))
  }
  
  groups
}

# =========================================================
# Helper: build presence/absence sets for 2 groups
# =========================================================
make_sets_2groups <- function(otu, meta, group_var) {
  groups <- get_groups_2groups(meta, group_var)
  
  sets <- lapply(groups, function(g) {
    sam <- rownames(meta)[meta[[group_var]] == g]
    rownames(otu)[rowSums(otu[, sam, drop = FALSE] > 0) > 0]
  })
  
  names(sets) <- groups
  sets
}

# =========================================================
# Helper: compute count regions for 2 groups
# =========================================================
calc_count_regions <- function(sets) {
  g1 <- names(sets)[1]
  g2 <- names(sets)[2]
  
  only1  <- setdiff(sets[[g1]], sets[[g2]])
  only2  <- setdiff(sets[[g2]], sets[[g1]])
  shared <- intersect(sets[[g1]], sets[[g2]])
  
  out <- c(length(only1), length(shared), length(only2))
  names(out) <- c(paste0(g1, "_only"), "Shared", paste0(g2, "_only"))
  out
}

# Helper: compute numratio for 2 groups
# numratio = percentage of feature number
calc_numratio_2groups <- function(sets) {
  g1 <- names(sets)[1]
  g2 <- names(sets)[2]
  
  only1  <- setdiff(sets[[g1]], sets[[g2]])
  only2  <- setdiff(sets[[g2]], sets[[g1]])
  shared <- intersect(sets[[g1]], sets[[g2]])
  
  total_features <- length(unique(c(sets[[g1]], sets[[g2]])))
  
  res <- c(
    length(only1) / total_features,
    length(shared) / total_features,
    length(only2) / total_features
  )
  
  names(res) <- c(paste0(g1, "_only"), "Shared", paste0(g2, "_only"))
  res
}

# Helper: compute seqratio for 2 groups
# seqratio = percentage of feature abundance
# IMPORTANT: use the same group order as the Euler sets
calc_seqratio_2groups <- function(otu, meta, groups, group_var) {
  if (length(groups) != 2) {
    stop("groups must contain exactly 2 values.")
  }
  
  otu_group <- lapply(groups, function(g) {
    sam <- rownames(meta)[meta[[group_var]] == g]
    rowSums(otu[, sam, drop = FALSE])
  })
  
  otu_group <- do.call(cbind, otu_group)
  colnames(otu_group) <- groups
  
  g1 <- otu_group[, 1]
  g2 <- otu_group[, 2]
  
  shared_taxa  <- rownames(otu_group)[g1 > 0 & g2 > 0]
  g1_only_taxa <- rownames(otu_group)[g1 > 0 & g2 == 0]
  g2_only_taxa <- rownames(otu_group)[g2 > 0 & g1 == 0]
  
  total_abund <- sum(otu_group)
  
  res <- c(
    sum(g1[g1_only_taxa]) / total_abund,
    sum(g1[shared_taxa] + g2[shared_taxa]) / total_abund,
    sum(g2[g2_only_taxa]) / total_abund
  )
  
  names(res) <- c(paste0(groups[1], "_only"), "Shared", paste0(groups[2], "_only"))
  res
}

# Helper: create explicit labels for regions
make_region_labels <- function(groups) {
  g1 <- groups[1]
  g2 <- groups[2]
  c(g1, paste0(g1, " \u2229 ", g2), g2)
}

# Helper: save Euler plot with counts + numratio + seqratio
save_euler_plot <- function(
    fit,
    numratio,
    seqratio,
    filename_base,
    fills,
    main_title,
    groups,
    outdir = EULER_OUTDIR,
    width = 7,
    height = 5,
    dpi = 1200
) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Full output path, without extension
  output_base <- file.path(outdir, filename_base)
  
  region_labels <- make_region_labels(groups)
  
  numratio_txt <- paste(
    paste0(region_labels, ": ", sprintf("%.2f%%", 100 * numratio)),
    collapse = "   |   "
  )
  
  seqratio_txt <- paste(
    paste0(region_labels, ": ", sprintf("%.2f%%", 100 * seqratio)),
    collapse = "   |   "
  )
  
  draw_plot <- function() {
    
    p <- plot(
      fit,
      fills = fills,
      alpha = 0.6,
      quantities = list(type = "counts", cex = 3.2),
      labels = list(cex = 3),
      shape = "ellipse",
      edges = TRUE
    )
    
    print(p)
    
    grid::pushViewport(
      grid::viewport(
        x = 0.5, y = 1,
        width = 1, height = 0.1,
        just = c("center", "top")
      )
    )
    
    grid::grid.text(
      main_title,
      gp = grid::gpar(fontsize = 30, fontface = "bold")
    )
    
    grid::popViewport()
    
    grid::grid.text(
      paste0("Numratio (% of feature number): ", numratio_txt),
      x = 0.5, y = 0.085,
      gp = grid::gpar(fontsize = 11)
    )
    
    grid::grid.text(
      paste0("Seqratio (% of feature abundance): ", seqratio_txt),
      x = 0.5, y = 0.035,
      gp = grid::gpar(fontsize = 11)
    )
  }
  
  tiff(
    filename = paste0(output_base, ".tiff"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  grid::grid.newpage()
  draw_plot()
  dev.off()
  
  pdf(
    file = paste0(output_base, ".pdf"),
    width = width,
    height = height
  )
  grid::grid.newpage()
  draw_plot()
  dev.off()
  
  message("Saved Euler plots to: ", normalizePath(outdir))
}
# Helper: run one complete analysis
run_euler_pipeline <- function(otu, meta, group_var, fills,
                               main_title, filename_base,
                               width = 7, height = 5) {
  
  sets <- make_sets_2groups(otu, meta, group_var)
  groups <- names(sets)
  
  message("Groups in ", main_title, ": ", paste(groups, collapse = " vs "))
  
  fit  <- euler(sets)
  
  count_regions <- calc_count_regions(sets)
  numratio      <- calc_numratio_2groups(sets)
  seqratio      <- calc_seqratio_2groups(otu, meta, groups, group_var)
  
  region_labels <- make_region_labels(groups)
  
  cat("\n=============================\n")
  cat(main_title, "\n")
  cat("=============================\n")
  print(count_regions)
  
  cat("\nNumratio (% of feature number):\n")
  print(setNames(round(100 * as.numeric(numratio), 2), region_labels))
  
  cat("\nSeqratio (% of feature abundance):\n")
  print(setNames(round(100 * as.numeric(seqratio), 2), region_labels))
  cat("\n")
  
  save_euler_plot(
    fit = fit,
    numratio = numratio,
    seqratio = seqratio,
    filename_base = filename_base,
    fills = fills,
    main_title = main_title,
    groups = groups,
    width = width,
    height = height
  )
  
  invisible(list(
    sets = sets,
    fit = fit,
    count_regions = count_regions,
    numratio = numratio,
    seqratio = seqratio,
    region_labels = region_labels
  ))
}

# 1) SAMPLE TYPE from mt_rarefied
load("RData/16s_dataset_rarefied.RData")

otu_sample  <- mt_rarefied$otu_table
meta_sample <- mt_rarefied$sample_table

res_sample_type <- run_euler_pipeline(
  otu = otu_sample,
  meta = meta_sample,
  group_var = "sample_type",
  fills = c("#8DA0CB", "#FC8D62"),   # Digesta, Mucus
  main_title = "Sample type",
  filename_base = "Euler_sample_type",
  width = 7,
  height = 5
)

# =========================================================
# 2) DIGESTA: NC vs HC
# =========================================================
load("RData/16s_dataset_dig_rarefied.RData")

otu_dig  <- mt_dig_rarefied$otu_table
meta_dig <- mt_dig_rarefied$sample_table

res_digesta <- run_euler_pipeline(
  otu = otu_dig,
  meta = meta_dig,
  group_var = "Diet",
  fills = c("#56B4E9", "#E69F00"),   # NC, HC
  main_title = "Digesta",
  filename_base = "Euler_digesta_diet",
  width = 7,
  height = 5
)

# =========================================================
# 3) MUCUS: NC vs HC
# =========================================================
load("RData/16s_dataset_mucus_rarefied.RData")

otu_muc  <- mt_mucus_rarefied$otu_table
meta_muc <- mt_mucus_rarefied$sample_table

res_mucus <- run_euler_pipeline(
  otu = otu_muc,
  meta = meta_muc,
  group_var = "Diet",
  fills = c("#56B4E9", "#E69F00"),   # NC, HC
  main_title = "Mucus",
  filename_base = "Euler_mucus_diet",
  width = 7,
  height = 5
)

# =========================================================
# Optional: collect all results in one summary table
# =========================================================
summary_table <- rbind(
  data.frame(
    Comparison = "Sample type",
    Region = res_sample_type$region_labels,
    Count = as.numeric(res_sample_type$count_regions),
    Numratio_percent = round(100 * as.numeric(res_sample_type$numratio), 2),
    Seqratio_percent = round(100 * as.numeric(res_sample_type$seqratio), 2)
  ),
  data.frame(
    Comparison = "Digesta",
    Region = res_digesta$region_labels,
    Count = as.numeric(res_digesta$count_regions),
    Numratio_percent = round(100 * as.numeric(res_digesta$numratio), 2),
    Seqratio_percent = round(100 * as.numeric(res_digesta$seqratio), 2)
  ),
  data.frame(
    Comparison = "Mucus",
    Region = res_mucus$region_labels,
    Count = as.numeric(res_mucus$count_regions),
    Numratio_percent = round(100 * as.numeric(res_mucus$numratio), 2),
    Seqratio_percent = round(100 * as.numeric(res_mucus$seqratio), 2)
  )
)

write.csv(summary_table,
          file = file.path("03_16S_diversity_analysis", 
                           "Euler_numratio_seqratio_summary.csv"), row.names = FALSE)

summary_table




