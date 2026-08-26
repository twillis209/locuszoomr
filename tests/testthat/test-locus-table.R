# DT's format*() family works by appending a columnDefs entry carrying a JS
# `render` function, targeted at 0-based column indices. So "is this column
# formatted?" is answerable by collecting the targets of every entry that has
# a render.
formatted_cols <- function(w, nms) {
  defs <- w$x$options$columnDefs
  tgt <- unlist(lapply(defs, function(d) if (!is.null(d$render)) d$targets))
  sort(nms[unlist(tgt) + 1L])
}

a_frame <- function() {
  data.frame(rsID = c("rs1", "rs2"),
             CHR = c(19L, 19L),
             BP = c(44905791, 44908822),
             P = c(1e-9, 2.5e-4),
             beta = c(0.1234, -0.5678),
             stringsAsFactors = FALSE)
}

test_that("coordinates are left unformatted and measurements are not", {
  d <- a_frame()
  w <- locus_table(d, coord_cols = c("CHR", "BP"))
  # 3 significant figures would render 44905791 as 4.49e+07 and 19 as 19.0.
  expect_equal(formatted_cols(w, colnames(d)), c("P", "beta"))
})

test_that("a character chromosome column changes nothing", {
  # detect_cols() can resolve `chrom` to a character column ("X", "chr7"),
  # which is not numeric and so was never a formatSignif target anyway.
  d <- a_frame()
  d$CHR <- as.character(d$CHR)
  w <- locus_table(d, coord_cols = c("CHR", "BP"))
  expect_equal(formatted_cols(w, colnames(d)), c("P", "beta"))
})

test_that("a frame with no measurements does not error", {
  # formatSignif() rejects an empty column selection, and rsID/CHR/BP is a
  # legitimate input - locus() requires nothing else.
  d <- a_frame()[, c("rsID", "CHR", "BP")]
  w <- expect_silent(locus_table(d, coord_cols = c("CHR", "BP")))
  expect_equal(length(formatted_cols(w, colnames(d))), 0L)
})

test_that("the underlying values are unchanged", {
  # The formatting is a render hook: the data shipped to the browser must
  # still be the real coordinates, so sorting and filtering stay numeric.
  d <- a_frame()
  w <- locus_table(d, coord_cols = c("CHR", "BP"))
  expect_true(any(grepl("44905791", w$x$data$BP, fixed = TRUE)))
})
