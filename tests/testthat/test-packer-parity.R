test_that("pack_rows assigns the lowest non-overlapping row", {
  # two overlapping intervals plus one that clears the first
  expect_equal(pack_rows(c(0, 0.5, 2), c(1, 1.5, 3)), c(1L, 2L, 1L))
})

test_that("pack_rows treats touching endpoints as non-overlapping", {
  expect_equal(pack_rows(c(0, 1), c(1, 2)), c(1L, 1L))
})

test_that("pack_rows handles an empty input", {
  expect_equal(pack_rows(numeric(0), numeric(0)), integer(0))
})

test_that("the JavaScript packer agrees with the R packer on every fixture", {
  skip_on_cran()
  node <- Sys.which("node")
  skip_if(node == "", "node not available")

  cases <- jsonlite::fromJSON(
    testthat::test_path("fixtures", "packer-cases.json"),
    simplifyDataFrame = FALSE
  )

  harness <- testthat::test_path("packer-parity.js")
  skip_if(!file.exists(harness), "parity harness not available")

  # system.file() resolves correctly both in the source tree (under
  # devtools::load_all()) and in an installed package, where inst/js is
  # flattened to js/ — R CMD check runs tests against the installed copy, so
  # a path built from inst/js literally (as the JS harness used to do)
  # cannot find the file there.
  relayout_js <- system.file("js", "genetrack-relayout.js", package = "locuszoomr")
  skip_if(relayout_js == "", "genetrack-relayout.js not found via system.file()")

  js <- jsonlite::fromJSON(system2(node, c(harness, relayout_js), stdout = TRUE),
                           simplifyDataFrame = FALSE)

  for (i in seq_along(cases)) {
    r_rows <- pack_rows(as.numeric(cases[[i]]$min), as.numeric(cases[[i]]$max))
    js_rows <- as.integer(unlist(js[[i]]$rows))
    expect_equal(js_rows, r_rows, info = cases[[i]]$name)
  }
})
