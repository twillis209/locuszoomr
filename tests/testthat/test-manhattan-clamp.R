# A p-value that underflowed to zero in the source is clamped to the smallest
# representable double before plotting, and -log10(5e-324) is 323.3. That is an
# artefact of the clamp, not a measurement, and one such SNP sets the y axis for
# the whole manhattan plot: in a 21M row Alzheimer GWAS a single clamped point
# stretched the axis to 323 while the largest real value was 114.5.

mh_data <- function(n = 400, clamp_at = NULL, top = 1e-115) {
  d <- data.frame(
    chrom = rep(c("1", "2"), each = n / 2),
    pos   = c(seq_len(n / 2), seq_len(n / 2)) * 1000,
    p     = rep(0.5, n),
    rsid  = paste0("rs", seq_len(n)),
    stringsAsFactors = FALSE
  )
  d$p[10] <- top                       # strongest genuine hit
  if (!is.null(clamp_at)) d$p[clamp_at] <- 0
  d$p[d$p < 5e-324] <- 5e-324          # what zoom() does before plotting
  d
}

mh <- function(d) {
  suppressMessages(
    manhattan(d, "chrom", "pos", "p", "rsid", pcutoff = 5e-8, npoints = NA))
}

test_that("without clamped points nothing changes", {
  m <- mh(mh_data())
  expect_setequal(unique(m$data$col), c(1, 2, 3))
  expect_equal(max(m$data$logP), 115, tolerance = 1)
})

test_that("a clamped point does not stretch the y axis", {
  m <- mh(mh_data(clamp_at = 200))
  # would be 323.3 if the clamped value drove the axis
  expect_equal(max(m$data$logP), 115, tolerance = 1)
  expect_lt(max(m$data$logP), 200)
})

test_that("clamped points are pegged to the largest real value", {
  m <- mh(mh_data(clamp_at = 200))
  real_top <- m$data$logP[m$data$rsid == "rs10"]
  expect_equal(m$data$logP[m$data$rsid == "rs200"], real_top)
})

test_that("clamped points get their own colour level, overriding significance", {
  m <- mh(mh_data(clamp_at = 200))
  expect_equal(m$data$col[m$data$rsid == "rs200"], 4)
  # the genuine top hit stays on the significant level
  expect_equal(m$data$col[m$data$rsid == "rs10"], 3)
})

test_that("the colour scheme resolves with and without clamped points", {
  full <- c(c("royalblue", "skyblue", "red"), "black")
  for (d in list(mh_data(), mh_data(clamp_at = 200))) {
    lv <- as.numeric(levels(as.factor(mh(d)$data$col)))
    expect_false(anyNA(full[lv]))
  }
})

test_that("several clamped points are all pegged and recoloured", {
  m <- mh(mh_data(clamp_at = c(200, 201, 202)))
  expect_equal(sum(m$data$col == 4), 3)
  expect_equal(length(unique(m$data$logP[m$data$col == 4])), 1)
})

test_that("clamping is reported rather than silent", {
  expect_message(
    manhattan(mh_data(clamp_at = 200), "chrom", "pos", "p", "rsid",
              pcutoff = 5e-8, npoints = NA),
    "below floating point precision")
})
