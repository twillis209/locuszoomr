# The property that matters: the indexed path must return exactly the rows
# the full scan returned. Everything else is an optimisation detail.
naive_rows <- function(d, chrom, pos, seqname, xrange) {
  which(as.character(d[, chrom]) == as.character(seqname) &
          d[, pos] >= xrange[1] & d[, pos] <= xrange[2])
}

a_frame <- function(n = 500L, seed = 1L) {
  set.seed(seed)
  data.frame(CHR = sample(c(1:3, "X"), n, replace = TRUE),
             BP = sample.int(1e6, n),
             P = runif(n),
             rsID = paste0("rs", seq_len(n)),
             stringsAsFactors = FALSE)
}

test_that("the index partitions the rows into contiguous per-chromosome runs", {
  d <- a_frame()
  idx <- build_locus_index(d, "CHR", "BP")
  ix <- idx$index
  expect_equal(nrow(idx$data), nrow(d))
  expect_equal(ix$first[1], 1L)
  expect_equal(ix$last[nrow(ix)], nrow(d))
  # no gaps, no overlaps
  expect_equal(ix$first[-1], ix$last[-nrow(ix)] + 1L)
  # every row's chromosome agrees with the range it falls in
  for (k in seq_len(nrow(ix))) {
    expect_true(all(as.character(idx$data[ix$first[k]:ix$last[k], "CHR"]) ==
                      ix$chrom[k]))
  }
  # ...and positions ascend within each range, which is what findInterval needs
  for (k in seq_len(nrow(ix))) {
    expect_false(is.unsorted(idx$data[ix$first[k]:ix$last[k], "BP"]))
  }
})

test_that("indexed lookup returns the same rows as a full scan", {
  d <- a_frame(2000L)
  idx <- build_locus_index(d, "CHR", "BP")
  set.seed(99)
  for (i in 1:50) {
    chr <- sample(c(1:3, "X", "22"), 1)          # "22" is absent on purpose
    a <- sample.int(1e6, 1); b <- a + sample.int(3e5, 1)
    got <- locus_subset(idx, "BP", chr, c(a, b))
    want <- d[naive_rows(d, "CHR", "BP", chr, c(a, b)), ]
    expect_equal(sort(got$rsID), sort(want$rsID),
                 info = paste("chr", chr, a, b))
  }
})

test_that("boundary positions are included, and left to locus() to trim", {
  # locus() filters with strict > and <, so the index must not pre-drop a
  # variant sitting exactly on an edge - locus() would then never see it and
  # the two paths could disagree in the other direction.
  d <- data.frame(CHR = "1", BP = c(100, 200, 300), P = 0.5,
                  rsID = c("a", "b", "c"), stringsAsFactors = FALSE)
  idx <- build_locus_index(d, "CHR", "BP")
  expect_equal(locus_subset(idx, "BP", "1", c(100, 300))$rsID,
               c("a", "b", "c"))
  expect_equal(locus_subset(idx, "BP", "1", c(150, 250))$rsID, "b")
})

test_that("empty results are empty, not errors", {
  d <- a_frame()
  idx <- build_locus_index(d, "CHR", "BP")
  expect_equal(nrow(locus_subset(idx, "BP", "Y", c(1, 1e6))), 0L)   # no such chr
  expect_equal(nrow(locus_subset(idx, "BP", "1", c(2e6, 3e6))), 0L) # past the end
  expect_equal(nrow(locus_subset(idx, "BP", "1", c(-100, 0))), 0L)  # before the start
  expect_equal(locus_rows(idx, "BP", "Y", c(1, 1e6)), integer(0))
})

test_that("a numeric chromosome column indexes the same as a character one", {
  # detect_cols() resolves `chrom` to whatever the file had; both are common,
  # and coords$chr arrives as character either way.
  d <- a_frame()
  d_num <- d[d$CHR != "X", ]
  d_chr <- d_num
  d_num$CHR <- as.integer(d_num$CHR)
  a <- locus_subset(build_locus_index(d_num, "CHR", "BP"), "BP", "2", c(1, 5e5))
  b <- locus_subset(build_locus_index(d_chr, "CHR", "BP"), "BP", "2", c(1, 5e5))
  expect_equal(sort(a$rsID), sort(b$rsID))
  expect_gt(nrow(a), 0L)
})

test_that("duplicate positions are not split across the boundary", {
  # Multi-allelic sites give repeated positions; findInterval must take all
  # of them or none.
  d <- data.frame(CHR = "1", BP = c(rep(100, 5), rep(200, 3)), P = 0.5,
                  rsID = paste0("r", 1:8), stringsAsFactors = FALSE)
  idx <- build_locus_index(d, "CHR", "BP")
  expect_equal(nrow(locus_subset(idx, "BP", "1", c(100, 100))), 5L)
  expect_equal(nrow(locus_subset(idx, "BP", "1", c(99, 201))), 8L)
  expect_equal(nrow(locus_subset(idx, "BP", "1", c(101, 199))), 0L)
})
