#' Provisional decision constants.
#'
#' All cutoffs are provisional. Requiring multiple markers is a safeguard,
#' not a calibrated false-positive rate: ambient RNA can affect several genes.
#' Reference distributions must be fitted by species, chemistry and tissue.
#'
#' Note what is NOT here: the preflightx read-level numbers (R_MALE = 0.50,
#' R_FEMALE = 0.20, an 8-read floor). Those were fitted on reads, on two species,
#' and do not transfer to transcript-level data.
#' @noRd
.thresholds <- list(
  min_marker_genes = 2L,
  # library size cancels in Y/(Y+X); gene length, capture and tissue effects do not
  min_gametolog_frac = 0.05,
  # shares the X reference with the gametolog fraction, so not independent
  min_inact_frac = 0.05,
  # floor against stray counts
  min_cpm = 1,
  min_total_counts = 1000,
  # the gap between p_low and p_high abstains as uncertain
  p_high = 0.90,
  p_low = 0.10
)

#' Estimate genetic sex from expression counts
#'
#' Scores two complementary lines of evidence, expression of sex-specific
#' chromosome markers and expression of the X-inactivation marker, and reports
#' the state they jointly imply. It does not force a male/female answer.
#'
#' Absence of Y signal is **not** evidence of female: loss of Y in aged male
#' blood, shallow sequencing and incomplete annotation all look identical to a
#' female sample. Most species in the shipped panel have no annotated
#' X-inactivation marker at all, so for those a female verdict is never
#' available and the honest answer is `uncertain`.
#'
#' For single-cell input the default is to aggregate, and the biological unit is
#' the **donor**, not the library. Pass `group` for anything that may hold more
#' than one donor. A male cell can carry zero Y UMIs through dropout alone, and
#' ambient RNA in a mixed-sex run carries Y transcripts into female droplets,
#' which is a direct false-male mechanism. `per_cell = TRUE` is therefore
#' **unvalidated**: it has not been tested against cell-level donor labels or
#' doublet calls, and its output should be read as evidence, not as a call.
#'
#' @param x A counts matrix, sparse matrix, data frame, file path, or a
#'   supported container (Seurat, SingleCellExperiment).
#' @param species Species name, scientific or common. Detected from the gene
#'   identifiers when `NULL`.
#' @param group Optional grouping of columns into units, as a vector the length
#'   of `ncol(x)` or the name of a column in the object's metadata. Defaults to
#'   treating the whole matrix as one unit.
#' @param per_cell Report one row per column instead of aggregating.
#' @param model Which model settles the male/female question. `"ratio"` uses
#'   marker detection and normalised expression directly and is always
#'   available. `"logistic"` uses coefficients fitted per species (see
#'   [SexModels()]), and falls back to `"ratio"` for any species that has
#'   none. Structural decisions (insufficient depth, both signals present,
#'   no inactivation marker in the panel) are made before either model runs.
#' @param ... Passed to methods.
#' @param input_scale `"auto"` checks every stored value for fractional counts;
#'   `"counts"` requires integer-like values; `"normalised"` explicitly declares
#'   normalised or estimated expression and caps confidence even if rounded.
#'   Integer-like values alone cannot establish raw-count provenance.
#' @param annotation Optional sex-chromosome table from [SexChromosomeGenes()],
#'   or `"genevintage"` to fetch the newest annotation for the resolved species.
#'   `NULL` uses shipped markers without downloading. A supplied table must
#'   belong to the input species; use the full chromosome inventory, including
#'   zero-expression genes. Curated markers absent from it are excluded.
#' @return A data frame, one row per unit, with verdict, heuristic confidence,
#'   marker CPM, gametolog and inactivation fractions, depth and diagnostic flags.
#'   Confidence categories are not calibrated probabilities. `inact_measured`
#'   is `FALSE` when the species panel has no inactivation marker or the input
#'   lacks its row; `inact_score` is then 0 without being measured. `per_cell = TRUE`
#'   adds `per_cell_unvalidated` and caps confidence at medium.
#' @examples
#' path <- system.file("extdata", "example_counts.tsv.gz", package = "SexSeek")
#'
#' # A file path works directly; species is detected from the gene IDs.
#' # Without `group` the whole matrix is one unit: two donors mixed.
#' EstimateSex(path)[, c("unit", "verdict", "flags")]
#'
#' # One unit per sample.
#' counts <- as.matrix(read.delim(path, row.names = 1))
#' res <- EstimateSex(counts, group = colnames(counts))
#' res[, c("unit", "verdict", "confidence", "gametolog_frac", "inact_frac")]
#' @export
EstimateSex <- function(x, ...) UseMethod("EstimateSex")

#' @rdname EstimateSex
#' @export
EstimateSex.default <- function(x, species = NULL, group = NULL,
                                per_cell = FALSE, model = c("ratio", "logistic"),
                                annotation = NULL,
                                input_scale = c("auto", "counts", "normalised"), ...) {
  m <- .as_counts(x, ...)
  .estimate_on_matrix(m,
    species = species, group = group, per_cell = per_cell,
    model = match.arg(model), annotation = annotation,
    input_scale = match.arg(input_scale)
  )
}

#' @rdname EstimateSex
#' @export
EstimateSex.character <- function(x, ...) {
  if (length(x) != 1L || is.na(x)) stop("Give a single file path.", call. = FALSE)
  if (grepl("\\.rds(\\.gz)?$", x, ignore.case = TRUE) && file.exists(x)) {
    object <- readRDS(x)
    if (is.character(object)) {
      stop("RDS must contain expression data or a supported container, not a file path.",
        call. = FALSE
      )
    }
    return(EstimateSex(object, ...))
  }
  EstimateSex.default(x, ...)
}

#' The one place the call is actually made. Every method funnels here.
#' @noRd
.estimate_on_matrix <- function(m, species, group, per_cell,
                                model = "ratio", annotation = NULL,
                                input_scale = "auto") {
  fractional <- .looks_normalised(m)
  if (input_scale == "counts" && fractional) {
    stop("`input_scale = 'counts'` requires integer-like counts; fractional values found.",
      call. = FALSE
    )
  }
  s <- .scoring_setup(rownames(m), species, annotation)
  units <- .unit_index(m, group, per_cell)
  if (!is.null(s$reason)) {
    return(.unknown_frame(units, s$species, s$system, s$reason, s$note))
  }
  reg <- s$reg
  panel <- s$panel
  annotation <- s$annotation
  excluded <- sum(is.na(panel$gene_id))
  result <- .score_units(m, units, panel, reg, model, annotation,
    looks_normalised = input_scale == "normalised" || fractional
  )
  if (excluded > 0L) {
    result <- .cap_flag(result, "annotation_markers_missing")
    result$notes <- trimws(paste(
      result$notes, excluded,
      "curated marker(s) absent or ambiguous in the supplied annotation."
    ))
  }
  if (per_cell) result <- .cap_flag(result, "per_cell_unvalidated")
  result
}

#' Registry row and (annotated) scoring panel, or why scoring cannot run
#' @noRd
.scoring_setup <- function(features, species, annotation = NULL) {
  sp <- .resolve_or_detect(species, features)
  reg <- if (is.null(sp$species)) NULL else .species_row(sp$species)
  give_up <- function(reason, note) {
    list(
      species = reg$scientific_name %||% sp$species,
      system = reg$system %||% NA_character_, reason = reason, note = note
    )
  }
  if (is.null(reg)) {
    return(give_up("unrecognised", "Species could not be determined from the gene identifiers."))
  }
  if (identical(reg$status, "unsupported")) {
    return(give_up("unsupported", if (nzchar(reg$notes %||% "")) reg$notes else "Species is not supported."))
  }
  panel <- .scoring_panel(reg$scientific_name)
  if (is.null(panel)) {
    return(give_up("no_panel", "No usable marker tier for this species."))
  }
  annotation <- .prepare_annotation(annotation, reg$scientific_name)
  if (!is.null(annotation)) {
    chr <- ifelse(panel$role %in% c("inactivation", "male_specific"),
      "X", panel$role
    )
    panel <- .annotate_markers(panel, annotation, chr)
  }
  list(reg = reg, panel = panel, annotation = annotation)
}

#' Cap high confidence at medium and append `flag`
#' @noRd
.cap_flag <- function(result, flag) {
  result$confidence[result$confidence == "high"] <- "medium"
  result$flags <- vapply(result$flags, function(f) {
    paste(c(f[nzchar(f)], flag), collapse = ";")
  }, character(1))
  result
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

#' @noRd
.resolve_or_detect <- function(species, features) {
  if (!is.null(species)) {
    hit <- ResolveSpecies(species)
    if (is.null(hit)) {
      stop("Unknown species '", species,
        "'. See SexSpecies() for what is supported.",
        call. = FALSE
      )
    }
    return(list(species = hit, confidence = "high", method = "given"))
  }
  DetectSpecies(features)
}

#' Marker rows for scoring, preferring core over secondary.
#'
#' A species is scored at one tier only so the selected marker inventory
#' and its gametolog references describe the same evidence set.
#' @noRd
.scoring_panel <- function(species) {
  # registry species without a panel return unknown, not an error
  all_p <- SexPanels()
  p <- all_p[all_p$scientific_name == species, , drop = FALSE]
  if (nrow(p) == 0) {
    return(NULL)
  }
  for (tier in c("core", "secondary")) {
    keep <- p[p$tier == tier, , drop = FALSE]
    # an X-only tier cannot call anything
    if (any(keep$role %in% c("Y", "W", "male_specific", "inactivation"))) {
      return(keep)
    }
  }
  NULL
}

#' @noRd
.specific_panel <- function(panel, system) {
  roles <- if (identical(system, "ZW")) "W" else c("Y", "male_specific")
  panel[panel$role %in% roles, , drop = FALSE]
}

#' Map panel genes onto matrix rows.
#'
#' Identifier first, symbol second: Ensembl renames symbols between releases
#' (chicken HINTW is HINTLA1 in release 116) and some orthologs have no symbol
#' at all, so a symbol match is a fallback, never the primary key.
#' @noRd
.match_rows <- function(features, panel) {
  hit <- .match_marker_rows(features, panel)
  unique(hit[!is.na(hit)])
}

#' Match one row per marker, retaining missing entries for paired matching.
#' @noRd
.match_marker_rows <- function(features, panel) {
  ids <- .strip_version(features)
  hit <- match(.strip_version(panel$gene_id), ids)
  hit[is.na(panel$gene_id)] <- NA_integer_

  miss <- is.na(hit) & !is.na(panel$gene_name) & nzchar(panel$gene_name)
  if (any(miss)) {
    by_sym <- match(
      .norm_symbol(panel$gene_name[miss]),
      .norm_symbol(features)
    )
    hit[miss] <- by_sym
  }
  hit
}

#' Columns to aggregate over.
#' @noRd
.unit_index <- function(m, group, per_cell) {
  n <- ncol(m)
  if (per_cell) {
    labs <- colnames(m)
    if (is.null(labs)) labs <- as.character(seq_len(n))
    if (anyNA(labs) || any(!nzchar(labs)) || anyDuplicated(labs)) {
      stop("Per-column calls require unique, non-missing column names.", call. = FALSE)
    }
    return(split(seq_len(n), factor(labs, levels = labs)))
  }
  if (is.null(group)) {
    return(list(all = seq_len(n)))
  }
  if (anyNA(group) || any(!nzchar(as.character(group)))) {
    stop("`group` cannot contain missing or empty labels.", call. = FALSE)
  }
  if (length(group) != n) {
    stop("`group` must have one entry per column of the matrix (", n, ").",
      call. = FALSE
    )
  }
  split(seq_len(n), factor(group))
}

#' Rows of `features` that scoring reads
#'
#' `spec` and `inact` are deduplicated row indices; `iy` and `ix` hold one
#' entry per gametolog pair, `NA` when unmatched.
#' @noRd
.scoring_rows <- function(features, panel, reg, annotation = NULL) {
  spec <- .specific_panel(panel, reg$system)
  inact <- panel[panel$role == "inactivation", , drop = FALSE]
  # only XY species have pairs; others fall back to detection plus the CPM floor
  gp <- GametologPairs()
  gp <- gp[gp$scientific_name == reg$scientific_name & gp$y_tier %in% panel$tier, , drop = FALSE]
  yp <- .annotate_markers(data.frame(
    gene_id = gp$y_gene_id, gene_name = gp$y_gene_name,
    stringsAsFactors = FALSE
  ), annotation, rep("Y", nrow(gp)))
  xp <- .annotate_markers(data.frame(
    gene_id = gp$x_gene_id, gene_name = gp$x_gene_name,
    stringsAsFactors = FALSE
  ), annotation, rep("X", nrow(gp)))
  list(
    spec = .match_rows(features, spec), inact = .match_rows(features, inact),
    iy = .match_marker_rows(features, yp), ix = .match_marker_rows(features, xp)
  )
}

#' @noRd
.score_units <- function(m, units, panel, reg, model, annotation = NULL,
                         looks_normalised = .looks_normalised(m)) {
  r <- .scoring_rows(rownames(m), panel, reg, annotation)
  i_spec <- r$spec
  i_inact <- r$inact
  iy <- r$iy
  ix <- r$ix
  complete <- !is.na(iy) & !is.na(ix) & iy != ix
  n_pairs <- sum(complete)
  i_gy <- unique(iy[complete])
  i_gx <- unique(ix[complete])
  # X reference survives a filtered-out Y row
  i_xref <- unique(ix[!is.na(ix)])
  n_matched <- length(i_spec) + length(i_inact)

  # XIST filtered from the matrix is fixable; XIST absent from the species is not
  species_has_inact <- any(panel$role == "inactivation")
  matrix_has_inact <- length(i_inact) > 0

  coefs <- if (identical(model, "logistic")) {
    .coefs_for(reg$scientific_name)
  } else {
    NULL
  }

  het <- if (identical(reg$heterogametic, "female")) "female" else "male"
  homo <- if (het == "male") "female" else "male"

  out <- lapply(names(units), function(u) {
    cols <- units[[u]]
    # missing diagnostic values abstain below
    total <- sum(Matrix::colSums(m[, cols, drop = FALSE], na.rm = TRUE),
      na.rm = TRUE
    )

    spec_cpm <- .cpm(m, i_spec, cols, total)
    inact_cpm <- .cpm(m, i_inact, cols, total)
    n_spec <- .n_detected(m, i_spec, cols)

    gy <- if (length(i_gy)) sum(m[i_gy, cols, drop = FALSE], na.rm = TRUE) else 0
    gx <- if (length(i_gx)) sum(m[i_gx, cols, drop = FALSE], na.rm = TRUE) else 0
    gfrac <- if (length(i_gy) && (gy + gx) > 0) gy / (gy + gx) else NA_real_

    inact_ct <- if (length(i_inact)) {
      sum(m[i_inact, cols, drop = FALSE], na.rm = TRUE)
    } else {
      0
    }
    xref <- if (length(i_xref)) sum(m[i_xref, cols, drop = FALSE], na.rm = TRUE) else 0
    ifrac <- if (length(i_xref) && (inact_ct + xref) > 0) {
      inact_ct / (inact_ct + xref)
    } else {
      NA_real_
    }

    f <- list(
      gfrac = gfrac, ifrac = ifrac, n_pairs = n_pairs,
      spec_cpm = spec_cpm, inact_cpm = inact_cpm, n_spec = n_spec,
      n_spec_panel = length(i_spec), total = total,
      # homogametic calls need an inactivation marker
      can_call_homo = matrix_has_inact, het = het, homo = homo,
      species_has_inact = species_has_inact,
      missing_marker_values = anyNA(m[unique(c(i_spec, i_inact, i_gy, i_xref)),
        cols,
        drop = FALSE
      ])
    )
    call <- .decide(f, model, coefs)
    call <- .cap(call, reg, looks_normalised, f)

    data.frame(
      unit = u, species = reg$scientific_name, system = reg$system,
      panel_status = reg$status, verdict = call$verdict,
      confidence = call$confidence,
      model = call$model,
      p_heterogametic = call$p,
      gametolog_frac = if (is.null(f$gfrac)) NA_real_ else f$gfrac,
      inact_frac = if (is.null(f$ifrac)) NA_real_ else f$ifrac,
      n_gametolog_pairs = n_pairs,
      y_score = spec_cpm, inact_score = inact_cpm,
      # inact_score is 0 here without being measured
      inact_measured = matrix_has_inact,
      qc_score = total, n_y_core_detected = n_spec,
      n_panel_genes_matched = n_matched,
      flags = call$flags, notes = call$notes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' @noRd
.cpm <- function(m, rows, cols, total) {
  if (length(rows) == 0 || total == 0) {
    return(0)
  }
  sum(m[rows, cols, drop = FALSE], na.rm = TRUE) / total * 1e6
}

#' @noRd
.n_detected <- function(m, rows, cols) {
  if (length(rows) == 0) {
    return(0L)
  }
  sum(Matrix::rowSums(m[rows, cols, drop = FALSE], na.rm = TRUE) > 0,
    na.rm = TRUE
  )
}

#' The decision table.
#'
#' Both signals present is an unresolved biological state, including possible
#' mixtures. Neither present is not evidence for the homogametic sex.
#' @noRd
.decide <- function(f, model, coefs) {
  out <- function(verdict, confidence, flags, notes, m = "ratio", p = NA_real_) {
    list(
      verdict = verdict, confidence = confidence, flags = flags,
      notes = notes, model = m, p = p
    )
  }

  # structural gates run first; a binary model cannot express these states
  if (isTRUE(f$missing_marker_values)) {
    return(out(
      "unknown", "low", "missing_marker_values",
      "Missing marker measurements cannot be treated as zero expression."
    ))
  }
  if (f$total < .thresholds$min_total_counts) {
    return(out("unknown", "low", "low_depth", "Too few counts to call."))
  }

  # no matched pairs (ZW, or X partners missing): detection alone
  spec_pos <- if (!is.na(f$gfrac)) {
    f$n_spec >= .thresholds$min_marker_genes &&
      f$spec_cpm >= .thresholds$min_cpm &&
      f$gfrac >= .thresholds$min_gametolog_frac
  } else {
    f$n_spec >= .thresholds$min_marker_genes &&
      f$spec_cpm >= .thresholds$min_cpm
  }
  inact_pos <- f$can_call_homo && if (!is.na(f$ifrac)) {
    f$inact_cpm >= .thresholds$min_cpm &&
      f$ifrac >= .thresholds$min_inact_frac
  } else {
    f$inact_cpm >= .thresholds$min_cpm
  }

  if (spec_pos && inact_pos) {
    return(out(
      "possible_mixed", "low", "both_signals",
      "Both signals present; possible mixture, ambient RNA or altered chromosome regulation."
    ))
  }

  # model settles het vs homo within what the gates left
  if (identical(model, "logistic")) {
    if (is.null(coefs)) {
      r <- .decide_ratio(f, spec_pos, inact_pos)
      r$flags <- paste(
        c(r$flags[nzchar(r$flags)], "no_fitted_model"),
        collapse = ";"
      )
      r$notes <- trimws(paste(
        r$notes,
        "No fitted model for this species; fell back to the ratio model."
      ))
      return(r)
    }
    result <- .decide_logistic(f, coefs)
    # A fitted probability cannot manufacture positive biological evidence.
    unsupported_call <- (result$verdict == f$het && !spec_pos) ||
      (result$verdict == f$homo && !inact_pos)
    if (unsupported_call) {
      result$verdict <- "uncertain"
      result$confidence <- "low"
      result$flags <- "model_without_positive_evidence"
      result$notes <- "Model prediction lacks the required positive marker evidence."
    }
    return(result)
  }
  .decide_ratio(f, spec_pos, inact_pos)
}

#' Detection-and-expression rule. Always available, needs no training.
#' @noRd
.decide_ratio <- function(f, spec_pos, inact_pos) {
  out <- function(verdict, confidence, flags, notes) {
    list(
      verdict = verdict, confidence = confidence, flags = flags,
      notes = notes, model = "ratio", p = NA_real_
    )
  }
  if (spec_pos) {
    return(out(f$het, "high", "", ""))
  }
  if (inact_pos) {
    return(out(f$homo, "high", "", ""))
  }
  if (!f$can_call_homo) {
    if (isTRUE(f$species_has_inact)) {
      # usually gene filtering dropped the lncRNA
      return(out(
        "uncertain", "low", "inactivation_marker_missing_from_input",
        paste0(
          "No sex-specific signal, and the inactivation marker row is absent ",
          "from this matrix, so ", f$homo, " cannot be called. Re-run with an ",
          "unfiltered gene set to recover it."
        )
      ))
    }
    return(out(
      "uncertain", "low", "no_inactivation_marker",
      paste0(
        "No sex-specific signal, and this species has no annotated ",
        "inactivation marker, so ", f$homo, " cannot be called positively."
      )
    ))
  }
  out(
    "uncertain", "low", "no_signal",
    "Neither signal detected; could be loss of Y, low depth or erosion."
  )
}

#' Fitted logistic rule, with an abstain band rather than a 0.5 cut.
#' @noRd
.decide_logistic <- function(f, coefs) {
  p <- .p_het(f, coefs)
  out <- function(verdict, confidence, flags, notes) {
    list(
      verdict = verdict, confidence = confidence, flags = flags,
      notes = notes, model = "logistic", p = p
    )
  }
  if (p >= .thresholds$p_high) {
    return(out(f$het, "high", "", ""))
  }
  if (p <= .thresholds$p_low) {
    # without an inactivation marker, low P only means no Y seen
    if (!f$can_call_homo) {
      flag <- if (isTRUE(f$species_has_inact)) {
        "inactivation_marker_missing_from_input"
      } else {
        "no_inactivation_marker"
      }
      return(out(
        "uncertain", "low", flag,
        paste0(
          "Low P(", f$het, "), but no inactivation marker is available, so ",
          f$homo, " cannot be called positively."
        )
      ))
    }
    return(out(f$homo, "high", "", ""))
  }
  out(
    "uncertain", "medium", "abstain_band",
    paste0("P(", f$het, ") = ", signif(p, 3), " falls in the abstain band.")
  )
}

#' @noRd
.unknown_frame <- function(units, species, system, flag, note) {
  data.frame(
    unit = names(units),
    species = if (is.null(species)) NA_character_ else species,
    system = system, panel_status = flag, verdict = "unknown",
    confidence = "low", model = NA_character_, p_heterogametic = NA_real_,
    gametolog_frac = NA_real_, inact_frac = NA_real_, n_gametolog_pairs = 0L,
    y_score = NA_real_, inact_score = NA_real_, inact_measured = FALSE,
    qc_score = NA_real_, n_y_core_detected = NA_integer_,
    n_panel_genes_matched = 0L, flags = flag, notes = note,
    stringsAsFactors = FALSE
  )
}

#' Detect fractional values throughout the matrix without densifying sparse input.
#' This detects non-integer values, not their experimental provenance.
#' @noRd
.looks_normalised <- function(m, chunk_size = 1000000L) {
  x <- if (methods::is(m, "sparseMatrix")) m@x else as.vector(m)
  if (!length(x)) {
    return(FALSE)
  }
  for (start in seq.int(1, length(x), by = chunk_size)) {
    block <- x[seq.int(start, min(length(x), start + chunk_size - 1))]
    if (any(abs(block - round(block)) > 1e-8, na.rm = TRUE)) {
      return(TRUE)
    }
  }
  FALSE
}

#' Cap the reported confidence by what the evidence can actually support.
#'
#' These are heuristic categories, not calibrated probabilities. Registry
#' status describes the marker panel's history, not validation of the current
#' ratio cutoffs. Experimental panels and fractional input lower the category.
#' @noRd
.cap <- function(call, reg, looks_normalised, f) {
  order_of <- c(low = 1L, medium = 2L, high = 3L)
  cap <- 3L
  extra <- character()

  # experimental panels are uncalibrated orthology guesses
  if (!identical(reg$status, "validated")) {
    cap <- min(cap, 2L)
    extra <- c(extra, "experimental_panel")
  }
  if (looks_normalised) {
    cap <- min(cap, 2L)
    extra <- c(extra, "normalised_input")
  }
  # one marker row is one mismap from a wrong call
  if (isTRUE(f$n_spec_panel <= 1)) {
    cap <- min(cap, 2L)
    extra <- c(extra, "thin_panel_match")
  }

  if (order_of[[call$confidence]] > cap) {
    call$confidence <- names(order_of)[cap]
  }
  if (length(extra)) {
    call$flags <- paste(c(call$flags[nzchar(call$flags)], extra), collapse = ";")
  }
  call
}
