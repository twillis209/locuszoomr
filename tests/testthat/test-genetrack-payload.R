test_that("genetrack_payload builds tx and ex in Mb with 0-based exon indices", {
  TX <- data.frame(
    gene_id      = c("ENSG1", "ENSG2"),
    gene_name    = c("AAA", "BBB"),
    gene_biotype = c("protein_coding", "lincRNA"),
    start        = c(1.0, 2.0),
    end          = c(1.5, 2.5),
    strand       = c("+", "-"),
    stringsAsFactors = FALSE
  )
  EX <- data.frame(
    gene_id = c("ENSG2", "ENSG1", "ENSG_absent"),
    start   = c(2.1, 1.1, 9.0),
    end     = c(2.2, 1.2, 9.1),
    stringsAsFactors = FALSE
  )
  p <- genetrack_payload(TX, EX, cfg = list(italics = FALSE))

  expect_equal(nrow(p$tx), 2L)
  expect_equal(p$tx$priority, 1:2)
  expect_equal(p$tx$start, c(1.0, 2.0))

  # strand arrows, trailing for +, leading for -
  expect_equal(p$tx$label, c("AAA&#8594;", "&#8592;BBB"))

  # exons index into tx 0-based, and exons of unknown genes are dropped
  expect_equal(nrow(p$ex), 2L)
  expect_equal(p$ex$gene_idx, c(1L, 0L))
})

test_that("genetrack_payload wraps names in italics when asked", {
  TX <- data.frame(
    gene_id = "ENSG1", gene_name = "AAA", gene_biotype = "protein_coding",
    start = 1, end = 2, strand = "+", stringsAsFactors = FALSE
  )
  EX <- data.frame(gene_id = character(), start = numeric(), end = numeric())
  p <- genetrack_payload(TX, EX, cfg = list(italics = TRUE))

  expect_equal(p$tx$label, "<i>AAA</i>&#8594;")
})

test_that("genetrack_payload emits an empty label for blank gene names", {
  TX <- data.frame(
    gene_id = "ENSG1", gene_name = "", gene_biotype = "protein_coding",
    start = 1, end = 2, strand = "+", stringsAsFactors = FALSE
  )
  EX <- data.frame(gene_id = character(), start = numeric(), end = numeric())
  p <- genetrack_payload(TX, EX, cfg = list(italics = FALSE))

  expect_equal(p$tx$label, "")
})
