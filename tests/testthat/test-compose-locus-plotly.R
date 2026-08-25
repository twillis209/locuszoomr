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
  b$data$p <- runif(nrow(b$data))     # a second, unrelated trait
  list(a, b)
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
  expect_true(is.numeric(idx$lineTrace))
  expect_true(is.numeric(idx$labelTrace))
  expect_false(identical(idx$lineTrace, idx$labelTrace))
  expect_true(nzchar(idx$yaxis))
})

test_that("a single locus composes as two panels", {
  skip_if_no_ensdb()
  p <- compose_locus_plotly(two_loci()[1], "A", dynamic = FALSE)
  expect_s3_class(p, "plotly")
})
