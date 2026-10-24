.verdict_col <- c(
  female = "#c0392b", male = "#2c6fbb", uncertain = "grey45",
  possible_mixed = "#8e44ad", unknown = "grey70"
)

#' Shape scale: circle is female, triangle is male, cross is not recorded.
#' @importFrom ggplot2 .data
#' @noRd
.shape_scale <- function() {
  ggplot2::scale_shape_manual(
    name = "recorded", values = c(female = 16, male = 17), na.value = 4
  )
}

#' Plot evidence behind each call
#'
#' Plot sex-specific markers against the inactivation marker,
#' as log10(CPM + 1).
#'
#' @param calls Output of [EstimateSex()] or [CompareSex()].
#' @param recorded Optional recorded sex, named by unit. See [CompareSex()].
#' @param main Plot title.
#' @return A ggplot. Colour shows the call. Shape shows the recorded sex (circle
#'   is female, triangle is male, cross is not recorded). Its `data` is `calls`,
#'   with the recorded columns when you give `recorded`.
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

#' Plot the fraction of cells with a marker count
#'
#' One bar of [MarkerCellFraction()] for each unit. If you give `recorded`, a
#' point above each bar shows the recorded sex (circle is female, triangle is male).
#'
#' @inheritParams MarkerCellFraction
#' @param recorded Optional recorded sex, named by unit. See [CompareSex()].
#' @param main Plot title.
#' @return A ggplot. Its `data` is the [MarkerCellFraction()] table, with
#'   `recorded_sex` when you give `recorded`.
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
