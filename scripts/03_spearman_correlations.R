# Purpose: Run stratified Spearman correlations for Figure 5a.
# Inputs: results/prepared_data.rds and results/growth_summary.csv.
# Outputs: Figure 5a audit/results CSV files and the heatmap SVG.

source("scripts/00_helpers.R")
suppressPackageStartupMessages(library(ggplot2))

x <- prepare_pairs()
dir.create("results/figure_5a", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  x |>
    distinct(Treatment, sampling_time, animal_id, sample_id),
  "results/figure_5a/pairing_audit.csv"
)

res <- x |>
  group_by(Treatment, sampling_time) |>
  group_modify(~ {
    d <- .x
    keep <- d |>
      group_by(genus) |>
      summarise(
        prev = mean(rel_t0 > 0 | rel_t72 > 0, na.rm = TRUE),
        mean_rel = mean((rel_t0 + rel_t72) / 2, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(prev >= 0.20, mean_rel >= 0.0001) |>
      pull(genus)

    purrr::map_dfr(keep, function(g) {
      z <- d |>
        filter(genus == g, !is.na(delta_CLR), !is.na(auc_logCFU))

      if (
        nrow(z) < 5 ||
          !is.finite(sd(z$delta_CLR)) ||
          !is.finite(sd(z$auc_logCFU)) ||
          sd(z$delta_CLR) == 0 ||
          sd(z$auc_logCFU) == 0
      ) {
        return(tibble(
          genus = g,
          rho = NA_real_,
          p = NA_real_,
          n = nrow(z),
          status = "insufficient_or_constant"
        ))
      }

      test <- suppressWarnings(
        cor.test(
          z$delta_CLR,
          z$auc_logCFU,
          method = "spearman",
          exact = FALSE
        )
      )

      tibble(
        genus = g,
        rho = unname(test$estimate),
        p = test$p.value,
        n = nrow(z),
        status = "ok"
      )
    })
  }) |>
  ungroup() |>
  group_by(Treatment, sampling_time) |>
  mutate(padj = p.adjust(p, method = "BH")) |>
  ungroup() |>
  mutate(
    genus_clean = clean_genus(genus),
    sig_label = signif_mark(p, padj)
  )

readr::write_csv(
  res,
  "results/figure_5a/spearman_all_results.csv"
)

selected <- res |>
  filter(status == "ok", plot_taxon(genus_clean)) |>
  group_by(genus) |>
  filter(any(p < 0.05, na.rm = TRUE)) |>
  ungroup() |>
  mutate(comparison = paste(Treatment, sampling_time, sep = " / "))

if (nrow(selected)) {
  fig <- ggplot(
    selected,
    aes(
      x = comparison,
      y = forcats::fct_reorder(genus_clean, p, .fun = min, .desc = TRUE),
      fill = rho
    )
  ) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = sig_label), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    labs(x = NULL, y = NULL, fill = "Spearman rho") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    "results/figure_5a/spearman_heatmap.svg",
    fig,
    width = 13,
    height = max(6, 0.20 * n_distinct(selected$genus) + 2)
  )
} else {
  warning("No nominally associated genera for Figure 5a")
}
