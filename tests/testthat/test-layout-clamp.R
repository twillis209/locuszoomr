# The JS label anchor must mirror mapRow()'s edge clamping (R/genetracks.R):
# when a gene straddles the plot edge, its label is pulled inside the boundary
# so it stays visible, but only if the gene is wide enough to hold it.
#
# Unlike R, the JS clamps to the visible window itself rather than to xlim
# widened by 4%: plotly clips at the axis range, so R's deliberate overshoot
# would leave the label off-screen whenever halfw is smaller than it.
#
# Geometry of the fixtures, all with x0 = 10, x1 = 20, pxPerData = 60 and the
# harness's 6px-per-character stub:
#   label     = "--" + 4-char name = 6 chars = 36px = 0.6 data units
#   halfw     = 0.3

run_clamp_cases <- function() {
  skip_on_cran()
  node <- Sys.which("node")
  skip_if(node == "", "node not available")

  harness <- testthat::test_path("layout-clamp.js")
  skip_if(!file.exists(harness), "layout-clamp harness not available")

  relayout_js <- system.file("js", "genetrack-relayout.js",
                             package = "locuszoomr")
  skip_if(relayout_js == "", "genetrack-relayout.js not found")

  out <- system2(node, c(harness, relayout_js), stdout = TRUE)
  res <- jsonlite::fromJSON(out, simplifyDataFrame = FALSE)
  stats::setNames(res, vapply(res, function(x) x$name, character(1)))
}

anchor <- function(res, case, gene) {
  a <- Filter(function(x) x$gene == gene, res[[case]]$anchors)
  if (length(a) != 1L) stop("expected exactly one anchor for ", gene)
  a[[1]]$mid
}

test_that("a gene overhanging the left edge has its label pulled inside", {
  res <- run_clamp_cases()
  # unclamped midpoint would be (2 + 17) / 2 = 9.5, whose label spans
  # 9.2..9.8 and so starts left of the visible window. Clamped anchor sits
  # half a label width inside it: 10 + 0.3.
  expect_equal(anchor(res, "left edge, gene wide enough to hold its label",
                      "AAAA"), 10.3)
})

test_that("a gene overhanging the right edge has its label pulled inside", {
  res <- run_clamp_cases()
  # unclamped midpoint (13 + 28) / 2 = 20.5, label spans 20.2..20.8, past the
  # right edge of the visible window. Clamped anchor: 20 - 0.3.
  expect_equal(anchor(res, "right edge, gene wide enough to hold its label",
                      "BBBB"), 19.7)
})

test_that("a gene too narrow to hold its label is left alone", {
  res <- run_clamp_cases()
  # Clamping would need the label to fit between the boundary and the gene
  # end: 10 + 0.6 = 10.6, but the gene ends at 9.8. R declines to move it,
  # and so must the JS -- otherwise the label would float off its own gene.
  expect_equal(anchor(res, "left edge, gene too narrow to hold its label",
                      "CCCC"), 5.9)
})

test_that("a fully visible gene keeps its plain midpoint", {
  res <- run_clamp_cases()
  expect_equal(anchor(res, "fully visible gene is untouched", "DDDD"), 14)
})

test_that("a gene with no label is never clamped", {
  res <- run_clamp_cases()
  # halfw is 0 for an unlabelled gene, which would otherwise satisfy the
  # clamp conditions trivially and shift the footprint for no visible reason.
  expect_equal(anchor(res, "gene with no label is never clamped", ""), 9.5)
})
