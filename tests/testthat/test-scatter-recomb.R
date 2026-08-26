skip_if_no_ensdb <- function() {
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  library(EnsDb.Hsapiens.v75)
}

a_locus <- function() {
  data(SLE_gwas_sub, envir = environment())
  locus(SLE_gwas_sub, gene = "IRF5", flank = 1e5,
        ens_db = "EnsDb.Hsapiens.v75")
}

# Same shape link_recomb() attaches: see parse_ucsc_track().
with_recomb <- function(loc, n = 20L) {
  brk <- seq(loc$xrange[1], loc$xrange[2], length.out = n + 1L)
  loc$recomb <- data.frame(chrom = paste0("chr", loc$seqname),
                           start = brk[-length(brk)], end = brk[-1],
                           value = seq(1, 60, length.out = n))
  loc
}

traces <- function(p) suppressWarnings(plotly::plotly_build(p))$x$data

recomb_trace <- function(p) {
  d <- traces(p)
  hit <- vapply(d, function(t) identical(t$name, "recombination"), logical(1))
  list(idx = which(hit), n = length(d), trace = d[hit][[1]])
}

test_that("the recombination line is drawn on top of the points", {
  skip_if_no_ensdb()
  # With webGL every trace shares one canvas, so paint order is trace order.
  # Added first - as it was - the line vanishes under a dense window of
  # points, which is the whole reason it moved.
  r <- recomb_trace(scatter_plotly(with_recomb(a_locus())))
  expect_equal(r$idx, r$n)
  expect_equal(r$trace$yaxis, "y2")
})

test_that("the recombination line is drawn on top with beta symbols too", {
  skip_if_no_ensdb()
  loc <- with_recomb(a_locus())
  set.seed(1)
  loc$data$beta <- rnorm(nrow(loc$data))
  r <- recomb_trace(scatter_plotly(loc, beta = "beta"))
  expect_equal(r$idx, r$n)
})

test_that("recomb_col survives the move to the end of the trace list", {
  skip_if_no_ensdb()
  # The line carries no `color` of its own, so plotly's `colors` mapping hands
  # it one from the default palette - by trace position, which the reorder
  # changed. That lands on marker$color, which a mode = "lines" trace never
  # draws, so it is line$color that has to keep reading recomb_col.
  loc <- with_recomb(a_locus())
  expect_equal(recomb_trace(scatter_plotly(loc))$trace$line$color, "blue")
  expect_equal(
    recomb_trace(scatter_plotly(loc, recomb_col = "green"))$trace$line$color,
    "green")
})

test_that("the scatter colours are unaffected by the presence of recombination", {
  skip_if_no_ensdb()
  # `colors = scheme` moved from the recombination trace to the scatter one.
  # The colours the points end up with must not have moved with it.
  loc <- a_locus()
  cols <- function(p) {
    d <- traces(p)
    d <- d[!vapply(d, function(t) identical(t$name, "recombination"),
                   logical(1))]
    vapply(d, function(t) unlist(t$marker$color)[1], character(1))
  }
  expect_equal(cols(scatter_plotly(with_recomb(loc))), cols(scatter_plotly(loc)))

  # ...including on the LD path, which builds a six-level scheme.
  set.seed(1)
  loc$data$ld <- runif(nrow(loc$data))
  expect_equal(cols(scatter_plotly(with_recomb(loc))), cols(scatter_plotly(loc)))
})
