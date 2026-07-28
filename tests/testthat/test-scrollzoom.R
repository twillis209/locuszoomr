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

test_that("scrollZoom is off by default on both entry points", {
  skip_if_no_ensdb()
  loc <- make_locus()

  expect_false(isTRUE(genetrack_ly(loc)$x$config$scrollZoom))
  expect_false(isTRUE(locus_plotly(loc)$x$config$scrollZoom))
})

test_that("scrollZoom = TRUE reaches the widget config on both entry points", {
  skip_if_no_ensdb()
  loc <- make_locus()

  expect_true(genetrack_ly(loc, scrollZoom = TRUE)$x$config$scrollZoom)
  expect_true(locus_plotly(loc, scrollZoom = TRUE)$x$config$scrollZoom)
})

test_that("enabling scrollZoom does not clobber the existing modebar config", {
  skip_if_no_ensdb()
  loc <- make_locus()

  # plotly::config() merges rather than replaces, but that is the whole risk
  # here: a replacing implementation would silently drop the modebar
  # customisation that both entry points rely on.
  for (p in list(genetrack_ly(loc, scrollZoom = TRUE),
                 locus_plotly(loc, scrollZoom = TRUE))) {
    expect_false(p$x$config$displaylogo)
    expect_true("select2d" %in% p$x$config$modeBarButtonsToRemove)
    expect_true("lasso2d" %in% p$x$config$modeBarButtonsToRemove)
  }
})

test_that("scrollZoom composes with the dynamic gene track handler", {
  skip_if_no_ensdb()
  loc <- make_locus()

  # Scroll-zoom fires plotly_relayout, so the re-layout handler must survive
  # the extra config() call for the gene tracks to re-pack while scrolling.
  p <- locus_plotly(loc, scrollZoom = TRUE, dynamic = TRUE)
  expect_length(p$jsHooks$render, 1)
  expect_match(p$jsHooks$render[[1]]$code, "LZR.attach", fixed = TRUE)
  expect_true(p$x$config$scrollZoom)

  g <- genetrack_ly(loc, scrollZoom = TRUE, dynamic = TRUE)
  expect_length(g$jsHooks$render, 1)
  expect_true(g$x$config$scrollZoom)
})

test_that("locus_plotly does not warn about mismatched configs", {
  skip_if_no_ensdb()
  loc <- make_locus()

  # plotly::subplot() calls ensure_one(plots, "config") and warns "Can only
  # have one: config" when its inputs disagree. Setting scrollZoom on
  # genetrack_ly()'s config even to FALSE makes it differ from
  # scatter_plotly()'s and triggers that warning on every call, so the flag
  # must only ever be added when it is TRUE.
  expect_no_warning(locus_plotly(loc))
  expect_no_warning(locus_plotly(loc, scrollZoom = TRUE))
  expect_no_warning(locus_plotly(loc, scrollZoom = TRUE, dynamic = FALSE))
})

test_that("scrollZoom is applied on the static path too", {
  skip_if_no_ensdb()
  loc <- make_locus()

  p <- locus_plotly(loc, scrollZoom = TRUE, dynamic = FALSE)
  expect_true(p$x$config$scrollZoom)
  expect_null(p$jsHooks$render)
})
