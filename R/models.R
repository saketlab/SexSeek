#' Feature names, in the order the coefficient table stores them.
#'
#' Library-normalised expression and depth covariates still depend on assay
#' and tissue. Transfer to a different setting requires external validation.
#' @noRd
.model_features <- c(
  "(Intercept)",
  "log_spec_cpm",
  "log_inact_cpm",
  "frac_spec_detected",
  "inact_detected",
  "log10_depth"
)

#' Fitted logistic coefficients, one row per species per feature
#'
#' Empty until `data-raw/fit_models.R` has been run against the benchmark
#' cohort. A species with no rows falls back to the ratio model.
#'
#' @return A data frame with columns `scientific_name`, `feature`, `estimate`.
#' @export
SexModels <- function() {
  if (is.null(.cache$models)) {
    df <- .read_tsv(.extdata("models.tsv"))
    if (nrow(df)) df$estimate <- as.numeric(df$estimate)
    .cache$models <- df
  }
  .cache$models
}

#' @noRd
.coefs_for <- function(species) {
  m <- SexModels()
  rows <- m[m$scientific_name == species, , drop = FALSE]
  if (nrow(rows) == 0) {
    return(NULL)
  }
  b <- rows$estimate[match(.model_features, rows$feature)]
  if (anyNA(b)) {
    return(NULL)
  }
  stats::setNames(b, .model_features)
}

#' Turn the raw per-unit measurements into model features.
#' @noRd
.features <- function(f) {
  c(
    "(Intercept)" = 1,
    log_spec_cpm = log1p(f$spec_cpm),
    log_inact_cpm = log1p(f$inact_cpm),
    frac_spec_detected = if (f$n_spec_panel > 0) f$n_spec / f$n_spec_panel else 0,
    inact_detected = as.numeric(f$inact_cpm > 0),
    log10_depth = log10(max(f$total, 1))
  )
}

#' Probability that the unit is the heterogametic sex.
#'
#' Expressed as P(heterogametic) rather than P(male) so the same fitted model
#' works for a ZW species, where the heterogametic sex is female.
#' @noRd
.p_het <- function(f, coefs) {
  stats::plogis(sum(coefs * .features(f)[names(coefs)]))
}
