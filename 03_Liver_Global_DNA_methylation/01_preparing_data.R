library(readr)
library(dplyr)

methylome_rel_con <- read_tsv("input_files/methylome.tsv", show_col_types = FALSE) %>%
  mutate(
    total_cytosine = dC + mdC + hmdC,
    dC_rel   = dC   / total_cytosine,
    mdC_rel  = mdC  / total_cytosine,
    hmdC_rel = hmdC / total_cytosine
  )

write_tsv(methylome_rel_con, "input_files/methylome_for_dirichlet.tsv")

cat("Done. Common samples:", nrow(methylome_rel_con),
    "\nSaved: methylome_for_dirichlet.tsv\n")
