# Purpose: Fit per-treatment linear mixed-effects models for Figure 5b.
# Inputs: results/prepared_data.rds and results/growth_summary.csv.
# Outputs: Figure 5b model-results CSV and heatmap SVG.

source("scripts/00_helpers.R")
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(ggplot2)
})

x <- prepare_pairs()
dir.create("results/figure_5b", recursive = TRUE, showWarnings = FALSE)

fit_genus <- function(d, g) {
  z <- d |>
    filter(
      genus == g,
      !is.na(delta_CLR),
      !is.na(auc_logCFU),
      !is.na(animal_id)
    )

  empty <- function(status) {
    tibble(
      genus = g,
      beta_std = NA_real_,
      se = NA_real_,
      statistic = NA_real_,
      p = NA_real_,
      n = nrow(z),
      n_animals = n_distinct(z$animal_id),
      model_status = status,
      singular = NA
    )
  }

  if (nrow(z) < 8) {
    return(empty("too_few_observations"))
  }
  if (
    !is.finite(sd(z$delta_CLR)) ||
      !is.finite(sd(z$auc_logCFU)) ||
      sd(z$delta_CLR) == 0 ||
      sd(z$auc_logCFU) == 0
  ) {
    return(empty("zero_variance"))
  }

  z <- z |>
    mutate(
      delta_CLR_z = as.numeric(scale(delta_CLR)),
      auc_z = as.numeric(scale(auc_logCFU)),
      sampling_time = factor(
        sampling_time,
        levels = c("pre_treatment", "2d_treatment", "6w_treatment")
      ),
      animal_id = factor(animal_id)
    )

  fit <- tryCatch(
    lmer(
      auc_z ~ delta_CLR_z + sampling_time + (1 | animal_id),
      data = z,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(empty("model_failed"))
  }

  t <- broom.mixed::tidy(fit, effects = "fixed") |>
    filter(term == "delta_CLR_z")

  if (!nrow(t)) {
    return(empty("term_not_found"))
  }

  tibble(
    genus = g,
    beta_std = t$estimate,
    se = t$std.error,
    statistic = t$statistic,
    p = t$p.value,
    n = nrow(z),
    n_animals = n_distinct(z$animal_id),
    model_status = "ok",
    singular = isSingular(fit)
  )
}

res <- x |>
  group_by(Treatment) |>
  group_modify(~ {
    d <- .x
    keep <- d |>
      group_by(genus) |>
      summarise(
        prev_pair = mean(!is.na(delta_CLR)),
        prev_nonzero = mean(rel_t0 > 0 | rel_t72 > 0, na.rm = TRUE),
        mean_rel = mean((rel_t0 + rel_t72) / 2, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(
        prev_pair >= 0.20,
        prev_nonzero >= 0.20,
        mean_rel >= 0.0001
      ) |>
      pull(genus)

    purrr::map_dfr(keep, ~ fit_genus(d, .x))
  }) |>
  ungroup() |>
  group_by(Treatment) |>
  mutate(padj = p.adjust(p, method = "BH")) |>
  ungroup() |>
  mutate(
    genus_clean = clean_genus(genus),
    sig_label = signif_mark(p, padj)
  )

readr::write_csv(
  res,
  "results/figure_5b/mixed_models_all_results.csv"
)

selected <- res |>
  filter(model_status == "ok", plot_taxon(genus_clean)) |>
  group_by(genus) |>
  filter(any(p < 0.05, na.rm = TRUE)) |>
  ungroup()

if (nrow(selected)) {
  fig <- ggplot(
    selected,
    aes(
      x = Treatment,
      y = forcats::fct_reorder(genus_clean, p, .fun = min, .desc = TRUE),
      fill = beta_std
    )
  ) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = sig_label), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      na.value = "grey90"
    ) +
    labs(x = NULL, y = NULL, fill = "Standardized beta") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    "results/figure_5b/mixed_models_heatmap.svg",
    fig,
    width = 9,
    height = max(6, 0.20 * n_distinct(selected$genus) + 2)
  )
} else {
  warning("No nominally associated genera for Figure 5b")
}
