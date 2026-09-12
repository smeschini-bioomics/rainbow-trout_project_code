# libraries
library(dplyr)
library(stringr)
library(ggplot2)

# load dataset
load("01_MSDAP_DEA/2026-04-29_16-53-21/dataset.RData")

# Output directory
OUTPUT_DIR <- "02_DEA_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Prepare data using qvalue directly
de_dot <- dataset$de_proteins %>%
  mutate(
    contrast_clean = str_extract(
      contrast,
      regex(
        "NC_day_(\\d+)\\s+vs\\s+HC_day_\\1",
        ignore_case = TRUE
      )
    ),
    day = as.integer(
      str_match(contrast_clean, "day_(\\d+)")[, 2]
    ),
    direction = case_when(
      qvalue < 0.05 & foldchange.log2 > 0 ~ "HC > NC",
      qvalue < 0.05 & foldchange.log2 < 0 ~ "HC < NC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(contrast_clean),
    !is.na(day),
    is.finite(foldchange.log2),
    is.finite(qvalue),
    !is.na(direction)
  )

# Count significant proteins per algorithm/day/direction
dot_data <- de_dot %>%
  group_by(day, contrast_clean, dea_algorithm, direction) %>%
  summarise(
    n_proteins = n_distinct(protein_id),
    .groups = "drop"
  ) %>%
  mutate(
    day = factor(day, levels = c(1, 2, 3, 4, 10, 15, 22)),
    direction = factor(direction, levels = c("HC > NC", "HC < NC")),
    dea_algorithm = factor(
      dea_algorithm,
      levels = c("deqms", "ebayes", "msempire", "msqrob", "msqrobsum")
    )
  )

cat("Significant DEA rows:", nrow(de_dot), "\n")
print(table(de_dot$day, de_dot$direction))
#------------------------------------------------------------
# Plot
#------------------------------------------------------------
p <- ggplot(
  dot_data,
  aes(x = direction, y = dea_algorithm, size = n_proteins, colour = direction)
) +
  geom_point(alpha = 0.80) +
  geom_text(
    aes(
      label = n_proteins,
      vjust = ifelse(direction == "HC > NC", -0.8, 1.5)
    ),
    colour = "black",
    size = 5,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ day,
    nrow = 1,
    labeller = labeller(day = function(d) paste("Day", d))
  ) +
  scale_size_continuous(
    name = "Number of significant proteins\n(qvalue < 0.05)",
    range = c(5, 20),
    breaks = c(1, 5, 20, 100, 400, 900)
  ) +
  scale_colour_manual(
    name = "Direction\n(HC vs NC)",
    values = c("HC < NC" = "#44a7c4", "HC > NC" = "#f99943"),
    labels = c("HC < NC" = "HC < NC", "HC > NC" = "HC > NC")
  ) +
  coord_cartesian(clip = "off") +
  guides(
    colour = guide_legend(override.aes = list(size = 6))
  ) +
  labs(
    title = "Midgut Proteome DEA",
    x = NULL,
    y = "DEA Algorithm"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(size = 22,
                              face = "bold",
                              family = "mono",
                              hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 18),
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.margin = margin(5.5, 25, 5.5, 5.5)
  )

print(p)


# Save outputs
saveRDS(
  p,
  file = file.path(
    OUTPUT_DIR,
    "dea_dotplot_qvalue.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "dea_dotplot_qvalue.tiff"
  ),
  plot = p,
  width = 13,
  height = 5,
  units = "in",
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "dea_dotplot_qvalue.pdf"
  ),
  plot = p,
  width = 13,
  height = 5,
  units = "in"
)
