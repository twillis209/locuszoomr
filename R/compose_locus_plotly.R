
#' Stack one or more scatter panels above a shared gene track
#'
#' locus_plotly() composes exactly one scatter with one gene track. Rather
#' than change its interface - it is exported and its current form is under
#' review upstream - this composes N panels for zoom()'s use only. If that
#' review settles, locus_plotly() could be refactored to call this and the
#' duplication collapsed.
#'
#' The gene track is built from `loci[[1]]`: every trait shares one window,
#' so they share one set of genes.
#'
#' @param loci List of 'locus' objects, all covering the same window.
#' @param ylabs Character vector of y axis labels, same length as `loci`.
#' @param heights Relative panel heights, defaulting to equal scatter panels
#'   over a gene track taking 0.3.
#' @param maxrows,width,filter_gene_biotype,pcutoff Passed through.
#' @param dynamic,scrollZoom As in [locus_plotly()].
#' @param ... Passed to [scatter_plotly()].
#' @return A 'plotly' object.
#' @noRd
compose_locus_plotly <- function(loci, ylabs, heights = NULL,
                                 maxrows = 12, width = 600,
                                 filter_gene_biotype = NULL,
                                 pcutoff = 5e-8,
                                 dynamic = TRUE, scrollZoom = FALSE, ...) {
  n <- length(loci)
  if (is.null(heights)) {
    heights <- c(rep(0.7 / n, n), 0.3)
  }
  g <- genetrack_ly(loci[[1]], filter_gene_biotype = filter_gene_biotype,
                    maxrows = maxrows, width = width, blanks = "show",
                    dynamic = FALSE)
  gt <- attr(g, "genetrack_data")

  panels <- lapply(seq_len(n), function(i) {
    scatter_plotly(loci[[i]], pcutoff = pcutoff, ylab = ylabs[i],
                   showlegend = (i == 1L), ...)
  })

  sp <- plotly::subplot(c(panels, list(g)), shareX = TRUE, nrows = n + 1L,
                        heights = heights, titleY = TRUE, margin = 0)
  sp <- remap_overlaying_yaxes(sp)

  out <- if (!dynamic || is.null(gt)) {
    sp
  } else {
    # Same fail-open contract as locus_plotly(): the R side of
    # add_genetrack_relayout() can stop() if it cannot resolve the gene
    # track traces after subplot(), and that must degrade to a static plot
    # rather than taking the whole app down.
    tryCatch(
      add_genetrack_relayout(sp, gt$TX, gt$EX, gt$cfg),
      error = function(e) {
        warning("dynamic gene track disabled: ", conditionMessage(e),
                call. = FALSE)
        sp
      })
  }
  if (scrollZoom) out <- plotly::config(out, scrollZoom = TRUE)
  out
}


#' Repoint secondary y axes at their own panel after subplot()
#'
#' subplot() renumbers the y axes it merges but not their `overlaying`
#' references, so every recombination axis after panel 1's still points at
#' "y" and its line is drawn inside panel 1. Rematch each overlay to the base
#' axis sharing its domain; ambiguous matches are left alone. Single-panel
#' plots are unaffected, which is why locus_plotly() never needed this.
#'
#' @param p A 'plotly' object returned by [plotly::subplot()].
#' @return `p` with each `overlaying` reference repointed.
#' @noRd
remap_overlaying_yaxes <- function(p) {
  lay <- p$x$layout
  nms <- grep("^yaxis[0-9]*$", names(lay), value = TRUE)
  base <- nms[vapply(nms, function(n) is.null(lay[[n]]$overlaying), logical(1))]
  for (n in setdiff(nms, base)) {
    dom <- lay[[n]]$domain
    if (is.null(dom)) next
    hit <- base[vapply(base,
                       function(b) isTRUE(all.equal(lay[[b]]$domain, dom)),
                       logical(1))]
    if (length(hit) == 1L) {
      p$x$layout[[n]]$overlaying <- sub("^yaxis", "y", hit)
    }
  }
  p
}
