# The flicker at the end of a scroll gesture was zoom() re-rendering the whole
# locus widget for a window whose data it had already loaded, redrawing the
# same points in place. needs_refetch() is the guard that stops that.

test_that("zooming in needs no refetch", {
  expect_false(needs_refetch(c(1200, 1800), c(1000, 2000)))
})

test_that("panning inside the loaded window needs no refetch", {
  expect_false(needs_refetch(c(1100, 1900), c(1000, 2000)))
})

test_that("a view identical to the loaded window needs no refetch", {
  expect_false(needs_refetch(c(1000, 2000), c(1000, 2000)))
})

test_that("touching either edge exactly needs no refetch", {
  expect_false(needs_refetch(c(1000, 1500), c(1000, 2000)))
  expect_false(needs_refetch(c(1500, 2000), c(1000, 2000)))
})

test_that("panning off the left edge needs a refetch", {
  expect_true(needs_refetch(c(900, 1900), c(1000, 2000)))
})

test_that("panning off the right edge needs a refetch", {
  expect_true(needs_refetch(c(1100, 2100), c(1000, 2000)))
})

test_that("zooming out past both edges needs a refetch", {
  expect_true(needs_refetch(c(500, 2500), c(1000, 2000)))
})

test_that("a single base past an edge is enough to need a refetch", {
  expect_true(needs_refetch(c(999L, 2000L), c(1000L, 2000L)))
  expect_true(needs_refetch(c(1000L, 2001L), c(1000L, 2000L)))
})
