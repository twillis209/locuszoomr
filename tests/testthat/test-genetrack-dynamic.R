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

make_dense_locus <- function() {
  data(SLE_gwas_sub, package = "locuszoomr", envir = environment())
  locus(SLE_gwas_sub, gene = "UBE2L3", flank = 5e5,
        ens_db = "EnsDb.Hsapiens.v75")
}

hex_re <- "^#[0-9A-Fa-f]{6}$"

test_that("cfg on the genetrack_data attribute carries hex colours, not R colour names", {
  skip_if_no_ensdb()
  loc <- make_locus()

  p <- genetrack_ly(loc)
  cfg <- attr(p, "genetrack_data")$cfg
  expect_match(cfg$geneCol, hex_re)
  expect_match(cfg$exonCol, hex_re)
  expect_match(cfg$exonBorder, hex_re)
  # the raw defaults are R colour names, not valid CSS - if cfg ever carried
  # them unconverted this would catch it
  expect_false(cfg$geneCol %in% c("blue4", "skyblue"))
  expect_false(cfg$exonCol == "blue4")
  expect_false(cfg$exonBorder == "blue4")

  # a non-default, caller-supplied colour name must also arrive hexed
  p2 <- genetrack_ly(loc, gene_col = "red", exon_col = "red", exon_border = "red")
  cfg2 <- attr(p2, "genetrack_data")$cfg
  expect_match(cfg2$geneCol, hex_re)
  expect_match(cfg2$exonCol, hex_re)
  expect_match(cfg2$exonBorder, hex_re)
  expect_equal(cfg2$geneCol, col2hex("red"))
  expect_equal(cfg2$exonCol, col2hex("red"))
  expect_equal(cfg2$exonBorder, col2hex("red"))
})

test_that("the payload carries the unfiltered gene tables, not the maxrows-truncated panel", {
  skip_if_no_ensdb()
  loc <- make_dense_locus()

  filtered <- genetrack_ly(loc, maxrows = 8, plot = FALSE)
  p <- genetrack_ly(loc, maxrows = 8)
  gt <- attr(p, "genetrack_data")

  expect_gt(nrow(gt$TX), nrow(filtered$TX))
  expect_gt(nrow(gt$EX), nrow(filtered$EX))
  expect_gt(max(gt$TX$row), gt$cfg$maxrows)
})

test_that("dynamic = FALSE injects no JS hook; dynamic = TRUE injects LZR.attach", {
  skip_if_no_ensdb()
  loc <- make_locus()

  g_off <- genetrack_ly(loc, dynamic = FALSE)
  expect_true(is.null(g_off$jsHooks$render) || length(g_off$jsHooks$render) == 0)

  g_on <- genetrack_ly(loc, dynamic = TRUE)
  expect_equal(length(g_on$jsHooks$render), 1)
  expect_true(grepl("LZR.attach", g_on$jsHooks$render[[1]]$code, fixed = TRUE))

  p_off <- locus_plotly(loc, dynamic = FALSE)
  expect_true(is.null(p_off$jsHooks$render) || length(p_off$jsHooks$render) == 0)

  p_on <- locus_plotly(loc, dynamic = TRUE)
  expect_equal(length(p_on$jsHooks$render), 1)
  expect_true(grepl("LZR.attach", p_on$jsHooks$render[[1]]$code, fixed = TRUE))
})

test_that("a locus with no matching genes returns without error and without a JS hook", {
  skip_if_no_ensdb()
  loc <- make_locus()

  g <- genetrack_ly(loc, filter_gene_name = "NO_SUCH_GENE", dynamic = TRUE)
  expect_s3_class(g, "plotly")
  expect_true(is.null(g$jsHooks$render) || length(g$jsHooks$render) == 0)

  # plotly::subplot() warns "Can only have one: config" here because the
  # blank no-genes gene panel and the scatter panel each carry their own
  # plotly::config() call; that collision is pre-existing subplot()
  # behaviour, unrelated to the dynamic wiring under test here.
  p <- suppressWarnings(
    locus_plotly(loc, filter_gene_name = "NO_SUCH_GENE", dynamic = TRUE)
  )
  expect_s3_class(p, "plotly")
  expect_true(is.null(p$jsHooks$render) || length(p$jsHooks$render) == 0)
})
