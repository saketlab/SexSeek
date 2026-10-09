.verdict_col <- c(
  female = "#c0392b", male = "#2c6fbb", uncertain = "grey45",
  possible_mixed = "#8e44ad", unknown = "grey70"
)

#' circle female, triangle male, cross when not recorded
#' @importFrom ggplot2 .data
#' @noRd
.shape_scale <- function() {
  ggplot2::scale_shape_manual(
    name = "recorded", values = c(female = 16, male = 17), na.value = 4
  )
}

#' Plot the evidence behind each call
#'
#' One point per unit: sex-specific marker expression against inactivation
#' marker expression, both as log10(CPM + 1). Colour is the call; the shape is
#' the recorded sex when it is given (circle female, triangle male, cross when
#' not recorded). When no unit has a measured inactivation marker (none in the
#' species panel, or its row absent from the input), the plot has one axis,
#' marker expression per unit.
#'
#' The result is a ggplot, so scales and theme can be replaced, e.g.
#' `+ ggplot2::scale_shape_manual(values = c(female = 1, male = 2))`.
#'
#' @param calls Output of [EstimateSex()] or [CompareSex()].
#' @param recorded Optional recorded sex named by unit; see [CompareSex()].
#' @param main Plot title.
#' @return A ggplot. Its `data` is `calls`, with the recorded columns when
#'   `recorded` was given.
#' @examples
#' ex <- ExampleData("mouse_bulk")
#' calls <- EstimateSex(ex$counts, species = ex$species, group = colnames(ex$counts))
#' PlotSexEvidence(calls, ex$recorded, main = ex$study)
#' @export
PlotSexEvidence <- function(calls, recorded = NULL, main = NULL) {
  if (!is.null(recorded)) calls <- CompareSex(calls, recorded)
  need <- c("unit", "verdict", "y_score", "inact_score", "inact_measured")
  if (!all(need %in% names(calls))) {
    stop("'calls' must be the output of EstimateSex().", call. = FALSE)
  }
  if (all(is.na(calls$y_score))) {
    stop("No marker scores to plot: the species has no usable panel.", call. = FALSE)
  }
  x <- ggplot2::aes(log10(.data$y_score + 1))
  xlab <- "Sex-specific markers, log10(CPM + 1)"
  if (any(calls$inact_measured)) {
    p <- ggplot2::ggplot(calls, x) +
      ggplot2::aes(y = log10(.data$inact_score + 1)) +
      ggplot2::labs(x = xlab, y = "Inactivation marker, log10(CPM + 1)")
  } else {
    o <- order(match(calls$verdict, names(.verdict_col)), calls$y_score)
    p <- ggplot2::ggplot(calls, x) +
      ggplot2::aes(y = .data$unit) +
      ggplot2::scale_y_discrete(limits = calls$unit[o]) +
      ggplot2::labs(x = xlab, y = NULL)
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$verdict), size = 3) +
    ggplot2::scale_colour_manual(
      name = "call", values = .verdict_col,
      breaks = intersect(names(.verdict_col), calls$verdict)
    ) +
    ggplot2::labs(title = main) +
    ggplot2::theme_bw()
  if (!"recorded_sex" %in% names(calls)) {
    return(p)
  }
  p + ggplot2::aes(shape = .data$recorded_sex) + .shape_scale()
}

#' Plot the fraction of cells carrying a marker count
#'
#' Bars of [MarkerCellFraction()] per unit. When `recorded` is given, each bar
#' is topped by the recorded sex (circle female, triangle male). The result is
#' a ggplot, so scales and theme can be replaced.
#'
#' @inheritParams MarkerCellFraction
#' @param recorded Optional recorded sex named by unit; see [CompareSex()].
#' @param main Plot title.
#' @return A ggplot whose `data` is the [MarkerCellFraction()] table, with
#'   `recorded_sex` when `recorded` was given.
#' @examples
#' ex <- ExampleData("mouse_sc")
#' PlotCellFraction(ex$counts, ex$species, ex$recorded, main = ex$study)
#' @export
PlotCellFraction <- function(x, species, recorded = NULL, group = NULL, main = NULL) {
  frac <- MarkerCellFraction(x, species, group = group)
  if (!is.null(recorded)) frac$recorded_sex <- .recorded_sex(.recorded_for(frac$unit, recorded))
  o <- if (is.null(recorded)) order(frac$fraction) else order(frac$recorded_sex, frac$fraction)
  p <- ggplot2::ggplot(frac, ggplot2::aes(.data$unit, .data$fraction)) +
    ggplot2::geom_col(fill = "grey80") +
    ggplot2::scale_x_discrete(limits = frac$unit[o]) +
    ggplot2::scale_y_continuous(limits = c(0, 1.1)) +
    ggplot2::labs(x = NULL, y = "cells with a sex-specific marker count", title = main) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (is.null(recorded)) {
    return(p)
  }
  p + ggplot2::geom_point(ggplot2::aes(y = .data$fraction + 0.07, shape = .data$recorded_sex), size = 3) +
    .shape_scale()
}
