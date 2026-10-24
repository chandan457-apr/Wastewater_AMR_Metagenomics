
# ================================================================
# Figure3 MAG-level ARG/MGE analysis
# for 33 HQ MAGs.
#
# Structure:
#   A: proportion of ARG carrying contigs for each ARG
#   B: proportion of ARG carrying contigs for each drug class
#   C: relationship between ARG proportion and MGE co occurrence
# ================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)

# -------------------------------
# 1. Read inputs
# -------------------------------
args <- read_tsv(
  file.path(data_dir, "ARG_hits_identity90.tsv"),
  show_col_types = FALSE
)

transposons <- read_tsv(
  file.path(data_dir, "Transposon_hits_identity90.tsv"),
  show_col_types = FALSE
)

plasmids <- read_tsv(
  file.path(data_dir, "Plasmid_hits_score90.tsv"),
  show_col_types = FALSE
)

mag_summary <- read_tsv(
  file.path(data_dir, "MAG_summary_identity90.tsv"),
  show_col_types = FALSE
)

args <- args %>%
  mutate(
    MAG_ID = as.character(MAG_ID),
    Contig_ID = as.character(Contig_ID)
  )

transposons <- transposons %>%
  mutate(
    MAG_ID = as.character(MAG_ID),
    Contig_ID = as.character(Contig_ID)
  )

plasmids <- plasmids %>%
  mutate(
    MAG_ID = as.character(MAG_ID),
    Contig_ID = as.character(Contig_ID)
  )

# -------------------------------
# 2. Same contig ARG-MGE table
# -------------------------------
arg_keys <- args %>%
  distinct(MAG_ID, Contig_ID)

plasmid_keys <- plasmids %>%
  distinct(MAG_ID, Contig_ID) %>%
  mutate(Plasmid = TRUE)

transposon_keys <- transposons %>%
  distinct(MAG_ID, Contig_ID) %>%
  mutate(Transposon = TRUE)

arg_mge <- arg_keys %>%
  left_join(plasmid_keys, by = c("MAG_ID", "Contig_ID")) %>%
  left_join(transposon_keys, by = c("MAG_ID", "Contig_ID")) %>%
  mutate(
    Plasmid = replace_na(Plasmid, FALSE),
    Transposon = replace_na(Transposon, FALSE),
    MGE_Associated = Plasmid | Transposon
  ) %>%
  left_join(args, by = c("MAG_ID", "Contig_ID")) %>%
  left_join(
    plasmids %>%
      select(MAG_ID, Contig_ID, Plasmid_Score),
    by = c("MAG_ID", "Contig_ID")
  )

# -------------------------------
# 3. Figure 3A
#    Proportion of ARG carrying
#    contigs for each ARG
# -------------------------------
n_arg_contigs <- n_distinct(args$Contig_ID)

arg_panel <- args %>%
  distinct(Contig_ID, ARG) %>%
  count(ARG, name = "ARG_contigs") %>%
  mutate(
    proportion = ARG_contigs / n_arg_contigs,
    label = paste0(
      ARG_contigs, "/",
      n_arg_contigs,
      " (",
      percent(proportion, accuracy = 0.1),
      ")"
    )
  ) %>%
  arrange(desc(proportion))

fig3A <- arg_panel %>%
  mutate(
    ARG = fct_reorder(ARG, proportion)
  ) %>%
  ggplot(aes(x = proportion, y = ARG)) +
  geom_col(
    width = 0.68,
    fill = "grey35"
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.08,
    size = 3.0
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(arg_panel$proportion) * 1.25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "a",
    subtitle = "Proportion of ARG carrying contigs",
    x = "Proportion of ARG-positive contigs",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(color = "black")
  )

# -------------------------------
# 4. Figure 3B
#    Proportion of ARG carrying
#    contigs for each drug class
# -------------------------------
class_panel <- args %>%
  distinct(Contig_ID, Drug_Class) %>%
  count(Drug_Class, name = "ARG_contigs") %>%
  mutate(
    proportion = ARG_contigs / n_arg_contigs,
    label = paste0(
      ARG_contigs, "/",
      n_arg_contigs,
      " (",
      percent(proportion, accuracy = 0.1),
      ")"
    )
  ) %>%
  arrange(desc(proportion))

fig3B <- class_panel %>%
  mutate(
    Drug_Class = fct_reorder(Drug_Class, proportion)
  ) %>%
  ggplot(aes(x = proportion, y = Drug_Class)) +
  geom_col(
    width = 0.68,
    fill = "grey35"
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.08,
    size = 3.0
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(class_panel$proportion) * 1.25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "b",
    subtitle = "Proportion of ARG carrying contigs",
    x = "Proportion of ARG-positive contigs",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(color = "black")
  )

# -------------------------------
# 5. Figure 3C
#    dot plot:
#    x = proportion ARG carrying
#        contigs for class
#    y = proportion co occurring
#        with MGE
#    size = unique ARG count
# -------------------------------
dot_panel <- arg_mge %>%
  distinct(
    Contig_ID,
    Drug_Class,
    ARG,
    MGE_Associated
  ) %>%
  group_by(Drug_Class) %>%
  summarise(
    ARG_contigs = n_distinct(Contig_ID),
    MGE_contigs = n_distinct(
      Contig_ID[MGE_Associated]
    ),
    unique_ARGs = n_distinct(ARG),
    .groups = "drop"
  ) %>%
  mutate(
    ARG_proportion = ARG_contigs / n_arg_contigs,
    MGE_proportion = if_else(
      ARG_contigs > 0,
      MGE_contigs / ARG_contigs,
      0
    )
  )

fig3C <- ggplot(
  dot_panel,
  aes(
    x = ARG_proportion,
    y = MGE_proportion,
    size = unique_ARGs
  )
) +
  geom_point(
    shape = 21,
    fill = "grey35",
    colour = "black",
    stroke = 0.4
  ) +
  geom_text_repel(
    aes(label = Drug_Class),
    size = 3.0,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(dot_panel$ARG_proportion) * 1.20),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(c(dot_panel$MGE_proportion, 0.05)) * 1.25),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_size_continuous(
    range = c(3.5, 8),
    breaks = pretty_breaks(n = 4)
  ) +
  labs(
    title = "c",
    x = "Proportion of ARG carrying contigs",
    y = "Proportion co-occurring with MGE",
    size = "Unique ARGs"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

# -------------------------------
# 6. Combine as a 3 panel
#    Figure 3
# -------------------------------
figure3 <- (fig3A | fig3B) / fig3C +
  plot_annotation(
    title = "Figure 3 | ARGs on HQ MAGs and their association with MGEs",
    subtitle = paste0(
      "33 HQ MAGs; BacAnt ARG/transposon identity ≥90%; ",
      "geNomad plasmid score ≥0.9; same-contig ARG–MGE association"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        size = 10.5,
        colour = "grey25"
      )
    )
  )

# -------------------------------
# 7. Save formats
# -------------------------------
ggsave(
  file.path(fig_dir, "Figure3.png"),
  figure3,
  width = 11,
  height = 8.7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure3.pdf"),
  figure3,
  width = 11,
  height = 8.7,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# -------------------------------
# 8. Save the plotting tables
# -------------------------------
write_tsv(
  arg_panel,
  file.path(fig_dir, "Figure3A_ARG_proportions.tsv")
)

write_tsv(
  class_panel,
  file.path(fig_dir, "Figure3B_drug_class_proportions.tsv")
)

write_tsv(
  dot_panel,
  file.path(fig_dir, "Figure3C_ARG_MGE_dotplot.tsv")
)

write_tsv(
  arg_mge,
  file.path(fig_dir, "Figure3_ARG_MGE_contig_table.tsv")
)

