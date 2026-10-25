# Purpose: Define shared functions used by the Figure 5 analyses.
# Inputs: Abundance tables and results/prepared_data.rds.
# Outputs: Reusable functions; this script does not write files directly.

suppressPackageStartupMessages(library(tidyverse))


# Data preparation helpers -------------------------------------------------

read_taxa <- function(path) {
  stopifnot(file.exists(path))
  x <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(x) <- "numeric"

  if (anyNA(x) || any(x < 0)) {
    stop("Invalid abundance table: ", path)
  }

  x
}

extract_genus <- function(x) {
  g <- stringr::str_extract(x, "g__[^;]+")
  g <- stringr::str_remove(g, "^g__")
  g[is.na(g) | g == "" | g == "__"] <- NA_character_
  out <- ifelse(is.na(g), paste0("Unresolved_", make.names(x)), g)
  stringr::str_replace_all(out, "[^A-Za-z0-9_\\-\\.]", "_")
}

collapse_genus <- function(x) {
  tibble(
    genus = extract_genus(rownames(x)),
    as.data.frame(x, check.names = FALSE)
  ) |>
    group_by(genus) |>
    summarise(
      across(everything(), sum),
      .groups = "drop"
    ) |>
    column_to_rownames("genus") |>
    as.matrix()
}

union_taxa <- function(a, b) {
  taxa <- sort(union(rownames(a), rownames(b)))

  expand <- function(x) {
    out <- matrix(
      0,
      length(taxa),
      ncol(x),
      dimnames = list(taxa, colnames(x))
    )
    out[rownames(x), ] <- x
    out
  }

  cbind(expand(a), expand(b))
}


# Metadata helpers ---------------------------------------------------------

recode_treatment <- function(x) {
  x <- tolower(as.character(x))

  dplyr::case_when(
    x %in% c("control", "ctrl") ~ "Control",
    x %in% c(
      "oxytetraciclina",
      "oxytetracycline",
      "oxitetraciclina"
    ) ~ "Oxytetracycline",
    x %in% c("penicilina", "penicillin", "penicilin") ~ "Penicilin",
    x %in% c("tilosina", "tylosine") ~ "Tylosine",
    TRUE ~ x
  )
}

recode_time <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    x %in% c("R1", "r1", "pre", "pre_treatment", "pretreatment") ~
      "pre_treatment",
    x %in% c("R2", "r2", "2d", "2d_treatment", "2d_post_treatment") ~
      "2d_treatment",
    x %in% c("R4", "r4", "6w", "6w_treatment", "6w_post_treatment") ~
      "6w_treatment",
    TRUE ~ x
  )
}

animal_number <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    stringr::str_detect(x, "\\.") ~
      stringr::str_extract(x, "(?<=\\.)[0-9]+$"),
    stringr::str_detect(x, "_") ~
      stringr::str_extract(x, "(?<=_)[0-9]+$"),
    TRUE ~ stringr::str_extract(x, "[0-9]+$")
  )
}


# Analysis and plotting helpers -------------------------------------------

clr <- function(counts, pseudocount = 0.5) {
  z <- log(counts + pseudocount)
  sweep(z, 2, colMeans(z), "-")
}

clean_genus <- function(x) {
  x <- stringr::str_remove(as.character(x), "^g__")
  x <- gsub("uncultured_bacterium", "uncultured", x)
  x <- gsub("_group$", "", x)
  stringr::str_squish(gsub("_", " ", x))
}

plot_taxon <- function(x) {
  !grepl(
    "^Unresolved|^uncultured|Ambiguous|Unknown|Incertae",
    x,
    ignore.case = TRUE
  )
}

signif_mark <- function(p, q) {
  dplyr::case_when(
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    p < 0.05 ~ "·",
    TRUE ~ ""
  )
}

prepare_pairs <- function() {
  dat <- readRDS("results/prepared_data.rds")
  growth <- readr::read_csv(
    "results/growth_summary.csv",
    show_col_types = FALSE
  ) |>
    transmute(
      Treatment = recode_treatment(treatment),
      sampling_time = recode_time(sampling_time),
      animal_id = animal_number(sample_id),
      sample_id,
      auc_logCFU
    )

  m <- dat$metadata |>
    transmute(
      SampleID,
      Treatment = as.character(Treatment),
      sampling_time = as.character(Time),
      animal_id = animal_number(Sheep),
      Time_Listeria = as.character(Time_Listeria)
    ) |>
    filter(Time_Listeria %in% c("t0", "t72"))

  if (anyDuplicated(m$SampleID)) {
    stop("Duplicated sample IDs in metadata")
  }
  if (!all(m$SampleID %in% colnames(dat$counts))) {
    stop("Missing sample counts")
  }

  cmat <- clr(dat$counts[, m$SampleID, drop = FALSE], 0.5)
  rmat <- dat$rel_abundance[, m$SampleID, drop = FALSE]

  tolong <- function(mat, name) {
    as.data.frame(mat, check.names = FALSE) |>
      tibble::rownames_to_column("genus") |>
      tidyr::pivot_longer(
        -genus,
        names_to = "SampleID",
        values_to = name
      )
  }

  a <- tolong(cmat, "CLR") |>
    left_join(m, by = "SampleID") |>
    select(-SampleID) |>
    pivot_wider(names_from = Time_Listeria, values_from = CLR)

  b <- tolong(rmat, "rel") |>
    left_join(m, by = "SampleID") |>
    select(-SampleID) |>
    pivot_wider(
      names_from = Time_Listeria,
      values_from = rel,
      names_prefix = "rel_"
    )

  joined <- inner_join(
    a,
    b,
    by = c("genus", "Treatment", "sampling_time", "animal_id")
  ) |>
    mutate(delta_CLR = t72 - t0) |>
    filter(!is.na(delta_CLR)) |>
    inner_join(
      growth,
      by = c("Treatment", "sampling_time", "animal_id")
    )

  if (anyDuplicated(
    joined[c("genus", "Treatment", "sampling_time", "animal_id")]
  )) {
    stop(
      "Non-unique paired genus/animal/time combinations; ",
      "check sample matching"
    )
  }
  if (!nrow(joined)) {
    stop("No microbiome/growth matches")
  }

  joined
}
