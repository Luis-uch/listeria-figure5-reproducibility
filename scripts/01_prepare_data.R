# Purpose: Aggregate abundance tables at genus level and align sample metadata.
# Inputs: Metadata and t0/t72 count and relative-abundance tables in data/.
# Outputs: results/prepared_data.rds and results/sample_counts.csv.

source("scripts/00_helpers.R")

input <- "data"
metadata <- readr::read_tsv(
  file.path(input, "mapping_file.tsv"),
  show_col_types = FALSE
) |>
  mutate(
    SampleID = as.character(SampleID),
    Sheep = as.character(Sheep)
  )

required <- c("SampleID", "Treatment", "Time", "Time_Listeria", "Sheep")

if (!all(required %in% names(metadata))) {
  stop("Missing metadata columns")
}

counts <- union_taxa(
  collapse_genus(
    read_taxa(
      file.path(input, "Microbial Metagenomics - t0_abundance.csv")
    )
  ),
  collapse_genus(
    read_taxa(
      file.path(input, "Microbial Metagenomics - t72_abundance.csv")
    )
  )
)

rel <- union_taxa(
  collapse_genus(
    read_taxa(
      file.path(input, "Microbial Metagenomics - t0_rel_abundance.csv")
    )
  ),
  collapse_genus(
    read_taxa(
      file.path(input, "Microbial Metagenomics - t72_rel_abundance.csv")
    )
  )
)

if (anyDuplicated(colnames(counts)) || anyDuplicated(metadata$SampleID)) {
  stop("Duplicate SampleID")
}
if (!setequal(colnames(counts), colnames(rel))) {
  stop("Count/relative-abundance sample mismatch")
}

metadata <- metadata |>
  filter(SampleID %in% colnames(counts)) |>
  arrange(Time_Listeria, Treatment, Sheep, Time)

if (!setequal(metadata$SampleID, colnames(counts))) {
  stop("Samples missing from metadata")
}

counts <- counts[, metadata$SampleID, drop = FALSE]
rel <- rel[, metadata$SampleID, drop = FALSE]

if (!identical(rownames(counts), rownames(rel))) {
  stop("Taxa mismatch")
}

dir.create("results", showWarnings = FALSE)
saveRDS(
  list(
    metadata = metadata,
    counts = counts,
    rel_abundance = rel
  ),
  "results/prepared_data.rds"
)

readr::write_csv(
  metadata |>
    count(Treatment, Time, Time_Listeria, name = "n_samples"),
  "results/sample_counts.csv"
)

message(
  "Prepared ",
  ncol(counts),
  " samples and ",
  nrow(counts),
  " genera"
)
