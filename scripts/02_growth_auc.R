# Purpose: Calculate Listeria growth AUC using trapezoidal integration.
# Inputs: data/growth.csv with logCFU measurements at 0, 24, 48, and 72 h.
# Outputs: results/growth_summary.csv.

source("scripts/00_helpers.R")

x <- readr::read_csv2(
  "data/growth.csv",
  show_col_types = FALSE
)

times <- c(0, 24, 48, 72)
cols <- c(
  "sample_id",
  "treatment",
  "sampling_time",
  paste0("logCFU_t", times)
)

if (!all(cols %in% names(x))) {
  stop("Missing growth columns")
}

growth_key <- c("sample_id", "treatment", "sampling_time")

if (anyNA(x[growth_key])) {
  stop("Missing values in growth sample identifiers")
}
if (anyDuplicated(x[growth_key])) {
  stop("Duplicate growth records for sample_id, treatment and sampling_time")
}

auc <- function(z) {
  if (anyNA(z)) {
    return(NA_real_)
  }

  sum(
    diff(times) *
      (head(z, -1) + tail(z, -1)) / 2
  )
}

x <- x |>
  mutate(
    across(
      all_of(paste0("logCFU_t", times)),
      as.numeric
    )
  ) |>
  rowwise() |>
  mutate(
    auc_logCFU = auc(
      c_across(all_of(paste0("logCFU_t", times)))
    ),
    final_logCFU = logCFU_t72,
    delta_logCFU_72 = logCFU_t72 - logCFU_t0
  ) |>
  ungroup()

dir.create("results", showWarnings = FALSE)
readr::write_csv(
  x,
  "results/growth_summary.csv"
)

message("AUC computed for ", nrow(x), " growth records")
