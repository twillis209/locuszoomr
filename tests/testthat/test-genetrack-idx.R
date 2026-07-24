skip_if_no_ensdb <- function() {
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  # locus() requires the ens_db package to be *attached* (it checks
  # `.packages()`), not merely installed, so attach it here.
  library(EnsDb.Hsapiens.v75)
}

make_locus <- function() {
  data(SLE_gwas_sub, package = "locuszoomr", envir = environment())
  locus(SLE_gwas_sub, gene = "IRF5", flank = c(7e4, 2e5), LD = "r2",
        ens_db = "EnsDb.Hsapiens.v75")
}

test_that("resolve_genetrack_idx finds the gene track in a standalone plot", {
  skip_if_no_ensdb()
  loc <- make_locus()
  built <- plotly::plotly_build(genetrack_ly(loc))
  idx <- resolve_genetrack_idx(built)

  expect_equal(idx$lineTrace, 0L)
  expect_equal(idx$labelTrace, 1L)
  expect_equal(idx$xaxis, "x")
  expect_equal(idx$yaxis, "y")
  # every exon shape belongs to the gene panel; there are no others
  expect_equal(length(idx$shapeIdx), length(built$x$layout$shapes))
})

test_that("resolve_genetrack_idx survives subplot axis renaming", {
  skip_if_no_ensdb()
  loc <- make_locus()
  built <- plotly::plotly_build(locus_plotly(loc))
  idx <- resolve_genetrack_idx(built)

  # gene traces sit after the scatter panel's LD-group traces
  expect_gt(idx$lineTrace, 0L)
  expect_equal(idx$labelTrace, idx$lineTrace + 1L)

  # x is shared between panels, y is not
  expect_equal(idx$xaxis, "x")
  expect_equal(idx$yaxis, "y2")

  # the p-value threshold line belongs to the scatter panel and must be excluded
  expect_lt(length(idx$shapeIdx), length(built$x$layout$shapes))
  excluded <- setdiff(seq_along(built$x$layout$shapes) - 1L, idx$shapeIdx)
  for (i in excluded) {
    expect_false(identical(built$x$layout$shapes[[i + 1L]]$yref, "y2"))
  }
})

test_that("resolve_genetrack_idx errors when the gene track is absent", {
  p <- plotly::plot_ly(x = 1:3, y = 1:3, type = "scatter", mode = "markers")
  expect_error(resolve_genetrack_idx(plotly::plotly_build(p)),
               "gene track traces")
})
