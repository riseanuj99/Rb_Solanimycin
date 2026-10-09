#!/usr/bin/env Rscript

# ================================================================
# MANUSCRIPT FIGURES AND TABLES
# Conservation of the canonical Dickeya Vfm system in Rouxiella
#
# INPUT:
# /scratch/al98750/Roux/06_QS_Vfm/05_results/
#
# OUTPUT:
# /scratch/al98750/Roux/06_QS_Vfm/07_manuscript/
# ================================================================


# ================================================================
# 1. PACKAGES
# ================================================================

suppressPackageStartupMessages({

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(patchwork)

})


# ================================================================
# 2. PATHS
# ================================================================

BASE_DIR <- "/scratch/al98750/Roux/06_QS_Vfm"

RESULT_DIR <- file.path(
  BASE_DIR,
  "05_results"
)

OUT_DIR <- file.path(
  BASE_DIR,
  "07_manuscript"
)

FIG_DIR <- file.path(
  OUT_DIR,
  "figures"
)

TABLE_DIR <- file.path(
  OUT_DIR,
  "tables"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ================================================================
# 3. INPUT FILES
# ================================================================

REFERENCE_FILE <- file.path(
  RESULT_DIR,
  "Ddadantii_exact_VfmA_Z_reference.tsv"
)

SUMMARY_FILE <- file.path(
  RESULT_DIR,
  "EXACT_VfmA_Z_cluster_summary.tsv"
)

MATRIX_FILE <- file.path(
  RESULT_DIR,
  "EXACT_VfmA_Z_presence_matrix.tsv"
)

LONG_FILE <- file.path(
  RESULT_DIR,
  "EXACT_VfmA_Z_best_hits_all_genomes.tsv"
)

RB_PRIMARY_FILE <- file.path(
  RESULT_DIR,
  "EXACT_20GA0316_VfmA_Z.tsv"
)

DSOLANI_FILE <- file.path(
  RESULT_DIR,
  "EXACT_Dsolani_MK10_VfmA_Z.tsv"
)


required_files <- c(

  REFERENCE_FILE,
  SUMMARY_FILE,
  MATRIX_FILE,
  LONG_FILE,
  RB_PRIMARY_FILE,
  DSOLANI_FILE

)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {

  stop(
    paste(
      "Missing required files:",
      paste(missing_files, collapse = "\n")
    )
  )

}


# ================================================================
# 4. READ DATA
# ================================================================

reference <- read_tsv(
  REFERENCE_FILE,
  show_col_types = FALSE
)

summary_df <- read_tsv(
  SUMMARY_FILE,
  show_col_types = FALSE
)

presence <- read_tsv(
  MATRIX_FILE,
  show_col_types = FALSE
)

long_hits <- read_tsv(
  LONG_FILE,
  show_col_types = FALSE
)

primary_hits <- read_tsv(
  RB_PRIMARY_FILE,
  show_col_types = FALSE
)

dsolani_hits <- read_tsv(
  DSOLANI_FILE,
  show_col_types = FALSE
)


# ================================================================
# 5. DEFINE GENE ORDER
# ================================================================

gene_order <- reference$gene


# Canonical order should be:
#
# vfmY K L M N O P Q R S T U V W X J I H G F E D C B Z A
#

print(gene_order)


# ================================================================
# 6. GENOME ORDER
# ================================================================

POSITIVE_CONTROL <- "Dsolani_MK10"

PRIMARY_RB <- "GCF_020740305.1"


other_rb <- setdiff(
  summary_df$genome,
  c(
    POSITIVE_CONTROL,
    PRIMARY_RB
  )
)


genome_order <- c(

  POSITIVE_CONTROL,
  PRIMARY_RB,
  other_rb

)


# ================================================================
# 7. DISPLAY LABELS
# ================================================================

make_genome_label <- function(x) {

  case_when(

    x == POSITIVE_CONTROL ~
      "D. solani MK10",

    x == PRIMARY_RB ~
      "R. badensis 20GA0316",

    TRUE ~ x

  )

}


genome_labels <- sapply(
  genome_order,
  make_genome_label
)


# ================================================================
# 8. CLASSIFICATION COLORS
# ================================================================

status_colors <- c(

  "HIGH"      = "#0072B2",

  "CANDIDATE" = "#56B4E9",

  "WEAK"      = "#E69F00",

  "ABSENT"    = "#E5E5E5"

)


# ================================================================
# 9. PREPARE PRESENCE MATRIX FOR HEATMAP
# ================================================================

presence_long <- presence %>%

  pivot_longer(

    cols = -genome,

    names_to = "gene",

    values_to = "classification"

  ) %>%

  mutate(

    gene = factor(
      gene,
      levels = gene_order
    ),

    classification = factor(

      classification,

      levels = c(
        "HIGH",
        "CANDIDATE",
        "WEAK",
        "ABSENT"
      )

    ),

    genome_label = make_genome_label(genome),

    genome_label = factor(

      genome_label,

      levels = rev(genome_labels)

    )

  )


# ================================================================
# 10. FIGURE 1A
# HEATMAP OF ALL 26 VFM GENES
# ================================================================

p_heatmap <- ggplot(

  presence_long,

  aes(
    x = gene,
    y = genome_label,
    fill = classification
  )

) +

  geom_tile(

    color = "white",

    linewidth = 0.35

  ) +

  scale_fill_manual(

    values = status_colors,

    drop = FALSE

  ) +

  scale_x_discrete(

    labels = function(x) {

      parse(
        text = paste0(
          "italic(",
          x,
          ")"
        )
      )

    }

  ) +

  labs(

    x = NULL,

    y = NULL,

    fill = "Homology\nclassification"

  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.text.x = element_text(

      angle = 45,

      hjust = 1,

      vjust = 1,

      size = 9

    ),

    axis.text.y = element_text(
      size = 9
    ),

    axis.ticks = element_blank(),

    axis.line = element_blank(),

    legend.title = element_text(
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      5,
      5,
      5,
      5
    )

  )


# ================================================================
# 11. PREPARE ACCEPTED-HOMOLOG BAR PLOT
#
# Accepted =
# HIGH + CANDIDATE
# ================================================================

bar_long <- summary_df %>%

  mutate(

    genome_label = make_genome_label(genome),

    genome_label = factor(

      genome_label,

      levels = rev(genome_labels)

    )

  ) %>%

  select(

    genome,

    genome_label,

    HIGH,

    CANDIDATE

  ) %>%

  pivot_longer(

    cols = c(
      HIGH,
      CANDIDATE
    ),

    names_to = "classification",

    values_to = "count"

  ) %>%

  mutate(

    classification = factor(

      classification,

      levels = c(
        "CANDIDATE",
        "HIGH"
      )

    )

  )


summary_plot <- summary_df %>%

  mutate(

    genome_label = make_genome_label(genome),

    genome_label = factor(

      genome_label,

      levels = rev(genome_labels)

    )

  )


# ================================================================
# 12. FIGURE 1B
# NUMBER OF ACCEPTED VFM HOMOLOGS
# ================================================================

p_bar <- ggplot(

  bar_long,

  aes(

    x = count,

    y = genome_label,

    fill = classification

  )

) +

  geom_col(

    width = 0.78

  ) +

  geom_text(

    data = summary_plot,

    aes(

      x = accepted_total + 0.45,

      y = genome_label,

      label = accepted_total

    ),

    inherit.aes = FALSE,

    size = 3

  ) +

  scale_fill_manual(

    values = status_colors[
      c(
        "CANDIDATE",
        "HIGH"
      )
    ],

    drop = FALSE

  ) +

  scale_x_continuous(

    limits = c(
      0,
      28
    ),

    breaks = c(
      0,
      5,
      10,
      15,
      20,
      25
    ),

    expand = expansion(
      mult = c(
        0,
        0
      )
    )

  ) +

  labs(

    x = "Accepted Vfm homologs",

    y = NULL,

    fill = NULL

  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.text.y = element_blank(),

    axis.ticks.y = element_blank(),

    axis.line.y = element_blank(),

    legend.position = "none",

    plot.margin = margin(
      5,
      5,
      5,
      0
    )

  )


# ================================================================
# 13. COMBINE FIGURE 1
# ================================================================

figure1 <- (

  p_heatmap |

  p_bar

) +

  plot_layout(

    widths = c(
      3.4,
      1
    )

  ) +

  plot_annotation(

    tag_levels = "A"

  )


# ================================================================
# 14. SAVE FIGURE 1
# ================================================================

ggsave(

  filename = file.path(
    FIG_DIR,
    "Figure_Vfm_conservation.pdf"
  ),

  plot = figure1,

  width = 14,

  height = 7.5,

  units = "in"

)


ggsave(

  filename = file.path(
    FIG_DIR,
    "Figure_Vfm_conservation_600dpi.png"
  ),

  plot = figure1,

  width = 14,

  height = 7.5,

  units = "in",

  dpi = 600

)


# ================================================================
# 15. SUPPLEMENTAL FIGURE
#
# Frequency of accepted homologs among the 17 Rouxiella genomes
# ================================================================

rb_gene_frequency <- presence_long %>%

  filter(
    genome != POSITIVE_CONTROL
  ) %>%

  filter(
    classification %in%
      c(
        "HIGH",
        "CANDIDATE"
      )
  ) %>%

  count(

    gene,

    classification,

    name = "genome_count"

  ) %>%

  complete(

    gene = factor(
      gene,
      levels = gene_order
    ),

    classification = factor(

      c(
        "HIGH",
        "CANDIDATE"
      ),

      levels = c(
        "CANDIDATE",
        "HIGH"
      )

    ),

    fill = list(
      genome_count = 0
    )

  )


p_frequency <- ggplot(

  rb_gene_frequency,

  aes(

    x = gene,

    y = genome_count,

    fill = classification

  )

) +

  geom_col(

    width = 0.75

  ) +

  scale_fill_manual(

    values = status_colors[
      c(
        "CANDIDATE",
        "HIGH"
      )
    ]

  ) +

  scale_x_discrete(

    labels = function(x) {

      parse(
        text = paste0(
          "italic(",
          x,
          ")"
        )
      )

    }

  ) +

  scale_y_continuous(

    limits = c(
      0,
      17
    ),

    breaks = c(
      0,
      5,
      10,
      15,
      17
    ),

    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )

  ) +

  labs(

    x = NULL,

    y = "R. badensis genomes with accepted homolog",

    fill = "Classification"

  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.text.x = element_text(

      angle = 45,

      hjust = 1,

      vjust = 1

    ),

    legend.position = "top"

  )


ggsave(

  filename = file.path(
    FIG_DIR,
    "FigureS_Vfm_gene_frequency.pdf"
  ),

  plot = p_frequency,

  width = 10,

  height = 5,

  units = "in"

)


ggsave(

  filename = file.path(
    FIG_DIR,
    "FigureS_Vfm_gene_frequency_600dpi.png"
  ),

  plot = p_frequency,

  width = 10,

  height = 5,

  units = "in",

  dpi = 600

)


# ================================================================
# 16. TABLE 1
# FULL GENOME-LEVEL SUMMARY
# ================================================================

table1 <- summary_df %>%

  transmute(

    Genome = make_genome_label(genome),

    `High-confidence homologs` = HIGH,

    `Candidate homologs` = CANDIDATE,

    `Weak hits` = WEAK,

    `No detectable hit` = ABSENT,

    `Accepted homologs` = accepted_total,

    `Accepted homologs on top contig` =
      accepted_on_top_contig,

    `Top-contig span (bp)` =
      top_contig_span_bp,

    `Gene-order concordance` =
      gene_order_concordance,

    `VfmE/H/I same locus` =
      VfmEHI_same_locus,

    `Coherent Vfm locus` =
      coherent_Vfm_locus

  )


write_csv(

  table1,

  file.path(
    TABLE_DIR,
    "Table1_Vfm_genome_summary.csv"
  )

)


write_tsv(

  table1,

  file.path(
    TABLE_DIR,
    "Table1_Vfm_genome_summary.tsv"
  )

)


# ================================================================
# 17. COMPACT MAIN-TEXT TABLE
# ================================================================

rb_summary <- summary_df %>%

  filter(
    genome != POSITIVE_CONTROL
  )


compact_table <- tibble(

  Comparison = c(

    "D. solani MK10",

    "R. badensis 20GA0316",

    "R. badensis, all 17 genomes"

  ),

  HIGH = c(

    summary_df$HIGH[
      summary_df$genome ==
        POSITIVE_CONTROL
    ],

    summary_df$HIGH[
      summary_df$genome ==
        PRIMARY_RB
    ],

    paste0(
      min(rb_summary$HIGH),
      "-",
      max(rb_summary$HIGH)
    )

  ),

  CANDIDATE = c(

    summary_df$CANDIDATE[
      summary_df$genome ==
        POSITIVE_CONTROL
    ],

    summary_df$CANDIDATE[
      summary_df$genome ==
        PRIMARY_RB
    ],

    paste0(
      min(rb_summary$CANDIDATE),
      "-",
      max(rb_summary$CANDIDATE)
    )

  ),

  `Accepted homologs` = c(

    summary_df$accepted_total[
      summary_df$genome ==
        POSITIVE_CONTROL
    ],

    summary_df$accepted_total[
      summary_df$genome ==
        PRIMARY_RB
    ],

    paste0(
      min(rb_summary$accepted_total),
      "-",
      max(rb_summary$accepted_total)
    )

  ),

  `VfmE/H/I same locus` = c(

    "Yes",

    "No",

    paste0(
      sum(
        rb_summary$VfmEHI_same_locus
      ),
      "/17"
    )

  ),

  `Coherent Vfm locus` = c(

    "Yes",

    "No",

    paste0(
      sum(
        rb_summary$coherent_Vfm_locus
      ),
      "/17"
    )

  )

)


write_csv(

  compact_table,

  file.path(
    TABLE_DIR,
    "Table1_compact_Vfm_comparison.csv"
  )

)


write_tsv(

  compact_table,

  file.path(
    TABLE_DIR,
    "Table1_compact_Vfm_comparison.tsv"
  )

)


# ================================================================
# 18. SUPPLEMENTAL TABLE:
# EXACT VFM GENE REFERENCE
# ================================================================

write_csv(

  reference,

  file.path(
    TABLE_DIR,
    "TableS1_Ddadantii_Vfm26_reference.csv"
  )

)


# ================================================================
# 19. SUPPLEMENTAL TABLE:
# PRESENCE / CLASSIFICATION MATRIX
# ================================================================

write_csv(

  presence,

  file.path(
    TABLE_DIR,
    "TableS2_Vfm26_classification_matrix.csv"
  )

)


# ================================================================
# 20. SUPPLEMENTAL TABLE:
# ALL BEST HITS
# ================================================================

write_csv(

  long_hits,

  file.path(
    TABLE_DIR,
    "TableS3_Vfm_best_hits_all_genomes.csv"
  )

)


# ================================================================
# 21. SUPPLEMENTAL TABLE:
# 20GA0316 DETAILED RESULTS
# ================================================================

write_csv(

  primary_hits,

  file.path(
    TABLE_DIR,
    "TableS4_20GA0316_Vfm_best_hits.csv"
  )

)


# ================================================================
# 22. SUPPLEMENTAL TABLE:
# D. SOLANI POSITIVE CONTROL
# ================================================================

write_csv(

  dsolani_hits,

  file.path(
    TABLE_DIR,
    "TableS5_Dsolani_Vfm_best_hits.csv"
  )

)


# ================================================================
# 23. SUPPLEMENTAL TABLE:
# GENE FREQUENCY ACROSS ROUXIELLA
# ================================================================

gene_frequency_table <- rb_gene_frequency %>%

  mutate(
    gene = as.character(gene)
  ) %>%

  pivot_wider(

    names_from = classification,

    values_from = genome_count,

    values_fill = 0

  ) %>%

  mutate(

    accepted_total = HIGH + CANDIDATE

  ) %>%

  arrange(

    factor(
      gene,
      levels = gene_order
    )

  )


write_csv(

  gene_frequency_table,

  file.path(
    TABLE_DIR,
    "TableS6_Vfm_gene_frequency_Rouxiella.csv"
  )

)


# ================================================================
# 24. PRINT IMPORTANT RESULTS
# ================================================================

cat("\n")
cat("============================================================\n")
cat(" MANUSCRIPT OUTPUT COMPLETE\n")
cat("============================================================\n")
cat("\n")

cat(
  "Figures:\n",
  FIG_DIR,
  "\n\n",
  sep = ""
)

cat(
  "Tables:\n",
  TABLE_DIR,
  "\n\n",
  sep = ""
)


cat("Compact comparison:\n")
print(compact_table)


cat("\n20GA0316 summary:\n")

print(

  summary_df %>%

    filter(
      genome == PRIMARY_RB
    )

)


cat("\nD. solani MK10 summary:\n")

print(

  summary_df %>%

    filter(
      genome == POSITIVE_CONTROL
    )

)


cat("\nDone.\n")