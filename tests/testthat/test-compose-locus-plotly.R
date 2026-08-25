skip_if_no_ensdb <- function() {
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  # locus() requires the ens_db package to be *attached* (it checks
  # `.packages()`), not merely installed, so attach it here. Matches the
  # pattern in test-genetrack-dynamic.R etc.
  library(EnsDb.Hsapiens.v75)
}

two_loci <- function() {
  data(SLE_gwas_sub, envir = environment())
  a <- locus(SLE_gwas_sub, gene = "IRF5", flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v75")
  b <- a
  set.seed(42)
  # `logP`, not `p`: scatter_plotly() plots data[, loc$yvar], and locus()
  # fixed yvar to "logP" at construction. Perturbing `p` here would be
  # inert - both panels would plot identical values, and the tests could
  # not tell composing two distinct panels from composing one twice.
  b$data$logP <- runif(nrow(b$data), 0, 8)   # a second, unrelated trait
  list(a, b)
}

# A minimal stand-in for what link_recomb() attaches: a data frame of
# chrom/start/end/value (see parse_ucsc_track()). scatter_plotly() reads only
# $start and $value, but its presence is what switches that panel onto a
# secondary y axis - which is precisely the thing that historically moved the
# gene track from y2 to y3.
add_fake_recomb <- function(loc, n = 20L) {
  xr <- loc$xrange
  brk <- seq(xr[1], xr[2], length.out = n + 1L)
  loc$recomb <- data.frame(chrom = paste0("chr", loc$seqname),
                           start = brk[-length(brk)],
                           end = brk[-1],
                           value = seq(1, 60, length.out = n))
  loc
}

test_that("compose_locus_plotly returns a plotly object", {
  skip_if_no_ensdb()
  p <- compose_locus_plotly(two_loci(), c("A", "B"), dynamic = FALSE)
  expect_s3_class(p, "plotly")
})

test_that("the gene track re-layout resolves across three panels", {
  skip_if_no_ensdb()
  p <- compose_locus_plotly(two_loci(), c("A", "B"), dynamic = TRUE)
  h <- p$jsHooks$render
  expect_gt(length(h), 0)
  idx <- h[[1]]$data$idx

  # Deliberately NOT expect_equal(idx$lineTrace, 5L) and friends: absolute
  # trace numbers depend on how many colour groups scatter_plotly() emitted,
  # which changes with the LD/significance levels present in the window, so
  # pinning them would break on data changes that are not regressions. The
  # invariants below are the ones the design actually rests on.

  # subplot() numbers y axes by row, so two scatter panels put the gene
  # track - always the last row - on y3. This is the whole technical risk of
  # the feature: the JS re-packer addresses the gene panel by axis name.
  expect_equal(idx$yaxis, "y3")
  # ...while x stays shared across all three rows (shareX = TRUE).
  expect_equal(idx$xaxis, "x")
  # add_segments() then add_text() are adjacent, as in the standalone and
  # single-panel cases (see test-genetrack-idx.R).
  expect_equal(idx$labelTrace, idx$lineTrace + 1L)

  # The shape partition: everything selected must belong to the gene panel,
  # and nothing left out may. Both scatter panels contribute a p-value
  # threshold line, so the excluded set is non-empty.
  shapes <- p$x$layout$shapes
  expect_lt(length(idx$shapeIdx), length(shapes))
  excluded <- setdiff(seq_along(shapes) - 1L, idx$shapeIdx)
  expect_gt(length(excluded), 0)
  for (i in excluded) {
    expect_false(identical(shapes[[i + 1L]]$yref, idx$yaxis))
  }
})

test_that("the gene track still resolves with two traits and recombination", {
  skip_if_no_ensdb()
  # Two traits + recombination is a DEFAULT configuration - zoom()'s recomb
  # checkbox is ticked whenever a recomb object is supplied - yet it adds a
  # secondary y axis to panel 1, pushing every axis after it along. Never
  # rendered before this test.
  loci <- two_loci()
  loci[[1]] <- add_fake_recomb(loci[[1]])
  p <- compose_locus_plotly(loci, c("A", "B"), dynamic = TRUE)
  idx <- p$jsHooks$render[[1]]$data$idx

  # panel 1 takes y (scatter) and y2 (recombination), panel 2 takes y3,
  # gene track y4.
  expect_equal(idx$yaxis, "y4")
  expect_equal(idx$xaxis, "x")
  expect_equal(idx$labelTrace, idx$lineTrace + 1L)

  shapes <- p$x$layout$shapes
  excluded <- setdiff(seq_along(shapes) - 1L, idx$shapeIdx)
  for (i in excluded) {
    expect_false(identical(shapes[[i + 1L]]$yref, idx$yaxis))
  }
})

test_that("a single locus composes as two panels", {
  skip_if_no_ensdb()
  p <- compose_locus_plotly(two_loci()[1], "A", dynamic = FALSE)
  expect_s3_class(p, "plotly")
})
