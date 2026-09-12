library(phyloseq)
library(microViz)
library(dplyr)
library(forcats)
library(scales)
library(ggplot2)
library(microeco)
library(file2meco)
library(RColorBrewer)
library(cowplot)
library(ggside)


# set output directory
rdata_path <- file.path("RData", "16s_dataset_microViz.RData")
output_dir <- "02_microbiota_composition_microViz"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

# Confirm that the input RData file exists
if (!file.exists(rdata_path)) {
  stop("Input file not found: ", rdata_path)
}

# load dataset
load(rdata_path)

# Convert microeco object to phyloseq
ps_raw <- meco2phyloseq(mt_abund)

# Fix unknown taxonomy BEFORE prepending ranks
ps_fixed <- ps_raw %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")
ps_fixed

group_levels <- c("Digesta_NC", "Digesta_HC", "Mucus_NC", "Mucus_HC")

group_labels <- c(
  "Digesta_NC" = "Digesta NC",
  "Digesta_HC" = "Digesta HC",
  "Mucus_NC"   = "Mucus NC",
  "Mucus_HC"   = "Mucus HC"
)

group_cols <- setNames(
  RColorBrewer::brewer.pal(n = 4, name = "Dark2"),
  group_levels
)

clr_pca <- ps_fixed %>%
  tax_transform(
    trans = "clr",
    rank = "Genus",
    keep_counts = TRUE
  ) %>%
  ord_calc(method = "PCA")

pca <- clr_pca %>%
  ord_plot(
    colour = "sample_type_diet",
    shape = "sample_type",
    alpha = 0.8,
    size = 3.5,
    stroke = 1.3
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = group_levels,
    labels = group_labels,
    name = "Compartment + Diet"
  ) +
  scale_fill_manual(
    values = group_cols,
    breaks = group_levels,
    labels = group_labels,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      Digesta = 22,
      Mucus   = 24
    ),
    name = "Compartment"
  ) +
  ggside::geom_xsidedensity(
    aes(fill = sample_type_diet),
    alpha = 0.45,
    show.legend = FALSE
  ) +
  ggside::geom_ysidedensity(
    aes(fill = sample_type_diet),
    alpha = 0.45,
    show.legend = FALSE
  ) +
  ggside::theme_ggside_void() +
  guides(
    colour = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(size = 5, stroke = 1.4, alpha = 1)
    ),
    shape = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(size = 5, stroke = 1.4, alpha = 1)
    )
  ) +
  theme_classic(base_size = 14) +
  ggside::theme_ggside_void() +
  theme(
    # Larger density panels: top x-density and right y-density
    ggside.panel.scale.x = 0.11,
    ggside.panel.scale.y = 0.11,
    ggside.panel.spacing.x = unit(1, "pt"),
    ggside.panel.spacing.y = unit(1, "pt"),
    
    # Stack colour legend, then shape legend
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.justification = "left",
    legend.spacing.y = unit(1, "mm"),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.55, "cm"),
    
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold")
  )


print(clr_pca)
pca

iris_plot <- clr_pca %>%
  ord_plot_iris(
    tax_level = "Genus",
    n_taxa = 10,
    anno_colour = "sample_type_diet",
    anno_colour_style = list(size = 3),
    anno_binary = "sample_type_diet"
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = group_levels,
    labels = group_labels,
    name = "Compartment + Diet"
  ) +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.7, "cm")
  )
iris_plot

combined_main <- cowplot::plot_grid(
  pca,
  iris_plot,
  nrow = 1,
  rel_widths = c(3, 4)
)

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Midgut-Associated Microbiota Composition",
    fontfamily = "mono",
    fontface = "bold",
    size = 22,
    x = 0.5,
    hjust = 0.5
  )

combined <- cowplot::plot_grid(
  title,
  combined_main,
  ncol = 1,
  rel_heights = c(0.06, 1)
)
combined

ps_fixed %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 5, facet_by = "sample_type_diet",
    label = NULL, bar_outline_colour = NA
  ) +
  coord_flip() +
  theme(axis.ticks.y = element_blank())


# Labels for nicer facet names
group_levels <- c("Digesta_NC", "Digesta_HC", "Mucus_NC", "Mucus_HC")

group_labels <- c(
  "Digesta_NC" = "Digesta NC",
  "Digesta_HC" = "Digesta HC",
  "Mucus_NC"   = "Mucus NC",
  "Mucus_HC"   = "Mucus HC"
)

# Make a copy only for plotting
ps_clean <- ps_raw %>%
  tax_fix() %>%
  phyloseq_validate(remove_undetected = TRUE)

ps_bar <- ps_clean

sd <- data.frame(sample_data(ps_bar))
sd$sample_type_diet_label <- factor(
  sd$sample_type_diet,
  levels = group_levels,
  labels = group_labels[group_levels]
)

sample_data(ps_bar) <- sample_data(sd)

bar_phylum <- ps_bar %>%
  comp_barplot(
    tax_level = "Phylum",
    n_taxa = 5,
    facet_by = "sample_type_diet_label",
    label = NULL,
    bar_outline_colour = NA
  ) +
  coord_flip() +
  labs(
    title = "Midgut-Associated Microbiota Relative Composition",
    x = NULL,
    y = "Relative abundance"
  ) +
  theme(
    plot.title = element_text(
      family = "mono",
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      #family = "mono",
      size = 16,
      face = "bold"
    ),
    legend.title = element_text(
      #family = "mono",
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      #family = "mono",
      size = 12
    ),
    axis.title = element_text(
      #family = "mono",
      size = 14,
      face = "bold"
    ),
    axis.text.x = element_text(
      #family = "mono",
      size = 12
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

bar_phylum


bar_genus <- ps_bar %>%
  comp_barplot(
    tax_level = "Genus",
    n_taxa = 10,
    facet_by = "sample_type_diet_label",
    label = NULL,
    bar_outline_colour = NA
  ) +
  coord_flip() +
  labs(
    title = "Midgut-Associated Microbiota Relative Composition",
    x = NULL,
    y = "Relative abundance"
  ) +
  theme(
    plot.title = element_text(
      family = "mono",
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      #family = "mono",
      size = 16,
      face = "bold"
    ),
    legend.title = element_text(
      #family = "mono",
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      #family = "mono",
      size = 12
    ),
    axis.title = element_text(
      #family = "mono",
      size = 14,
      face = "bold"
    ),
    axis.text.x = element_text(
      #family = "mono",
      size = 12
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

bar_genus


ggsave(
  filename = file.path(output_dir, "combined_plot.pdf"),
  plot = combined,
  width = 18,
  height = 7,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "combined_plot.tiff"),
  plot = combined,
  width = 18,
  height = 7,
  units = "in",
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "phylum_relative_composition.tiff"),
  plot = bar_phylum,
  width = 13,
  height = 5,
  units = "in",
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "phylum_relative_composition.pdf"),
  plot = bar_phylum,
  width = 13,
  height = 5,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "genus_relative_composition.tiff"),
  plot = bar_genus,
  width = 15,
  height = 6,
  units = "in",
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "genus_relative_composition.pdf"),
  plot = bar_genus,
  width = 15,
  height = 6,
  units = "in"
)

