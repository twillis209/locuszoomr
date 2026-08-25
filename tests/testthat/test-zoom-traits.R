test_that("check_trait_cols accepts a frame with every required column", {
  d <- data.frame(CHR = "1", BP = 1L, P = 0.5, rsID = "rs1")
  expect_true(check_trait_cols(d, c("CHR", "BP", "P", "rsID")))
})

test_that("check_trait_cols names every missing column", {
  d <- data.frame(CHR = "1", BP = 1L)
  expect_error(check_trait_cols(d, c("CHR", "BP", "P", "rsID")),
               "P")
  expect_error(check_trait_cols(d, c("CHR", "BP", "P", "rsID")),
               "rsID")
})

test_that("check_trait_cols names the offending argument", {
  d <- data.frame(CHR = "1")
  expect_error(check_trait_cols(d, c("CHR", "BP"), name = "data2"),
               "data2")
})

build_df <- function(ids, pos) {
  data.frame(rsID = ids, BP = pos, stringsAsFactors = FALSE)
}

test_that("check_same_build passes when shared ids agree on position", {
  ids <- paste0("rs", 1:2000)
  a <- build_df(ids, seq_along(ids) * 100L)
  b <- build_df(ids, seq_along(ids) * 100L)
  expect_silent(check_same_build(a, b, "BP", "rsID"))
})

test_that("check_same_build warns when shared ids disagree on position", {
  ids <- paste0("rs", 1:2000)
  a <- build_df(ids, seq_along(ids) * 100L)
  b <- build_df(ids, seq_along(ids) * 100L + 5000L)
  expect_warning(check_same_build(a, b, "BP", "rsID"), "genome build")
})

test_that("check_same_build tolerates a small disagreement", {
  ids <- paste0("rs", 1:2000)
  p <- seq_along(ids) * 100L
  a <- build_df(ids, p)
  p2 <- p; p2[1:100] <- p2[1:100] + 7L      # 95% agree
  b <- build_df(ids, p2)
  expect_silent(check_same_build(a, b, "BP", "rsID"))
})

test_that("check_same_build skips when too few ids are shared", {
  a <- build_df(paste0("rs", 1:50), 1:50 * 100L)
  b <- build_df(paste0("rs", 1:50), 1:50 * 100L + 999L)
  expect_message(check_same_build(a, b, "BP", "rsID"), "too few")
})

test_that("trait_labels prefers explicit names", {
  expect_equal(trait_labels(c("AD", "Dizziness"), "x", "y"),
               c("AD", "Dizziness"))
})

test_that("trait_labels falls back to the deparsed argument names", {
  expect_equal(trait_labels(NULL, "ad", "dizzy"), c("ad", "dizzy"))
})

test_that("trait_labels rejects unusable deparsed expressions", {
  expect_equal(trait_labels(NULL, 'read_gwas("x.tsv")', "df2"),
               c("Trait 1", "df2"))
})

test_that("trait_labels rejects over-long deparsed expressions", {
  expect_equal(trait_labels(NULL, strrep("a", 30), "df2"),
               c("Trait 1", "df2"))
})

fake_locus <- function(ids) {
  structure(list(data = data.frame(rsID = ids, stringsAsFactors = FALSE),
                 labs = "rsID"),
            class = "locus")
}

test_that("join_ld attaches r2 by SNP id", {
  loc2 <- fake_locus(c("rs1", "rs2", "rs3"))
  ref <- data.frame(snp = c("rs1", "rs3"), ld = c(1, 0.4),
                    stringsAsFactors = FALSE)
  out <- join_ld(loc2, ref, "rsID")
  expect_equal(out$data$ld, c(1, NA, 0.4))
})

test_that("join_ld leaves unmatched SNPs as NA", {
  loc2 <- fake_locus(c("rsX", "rsY"))
  ref <- data.frame(snp = "rs1", ld = 1, stringsAsFactors = FALSE)
  out <- join_ld(loc2, ref, "rsID")
  expect_true(all(is.na(out$data$ld)))
})

test_that("join_ld is a no-op given an empty reference", {
  loc2 <- fake_locus(c("rs1", "rs2"))
  ref <- data.frame(snp = character(0), ld = numeric(0),
                    stringsAsFactors = FALSE)
  out <- join_ld(loc2, ref, "rsID")
  expect_false("ld" %in% colnames(out$data))
})
