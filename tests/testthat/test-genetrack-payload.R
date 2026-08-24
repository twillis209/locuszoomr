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

test_that("genetrack_payload's hover column matches genetrack_ly()'s hovertext format", {
  TX <- data.frame(
    gene_id      = c("ENSG1", "ENSG2"),
    gene_name    = c("AAA", "BBB"),
    gene_biotype = c("protein_coding", "lincRNA"),
    start        = c(1.0, 2.0),
    end          = c(1.5, 2.5),
    strand       = c("+", "-"),
    stringsAsFactors = FALSE
  )
  EX <- data.frame(gene_id = character(), start = numeric(), end = numeric())
  p <- genetrack_payload(TX, EX, cfg = list(italics = FALSE))

  # Same format as R/genetrack_ly.R:138-143: gene_name, fullname (absent here,
  # contributes nothing), Gene ID, Biotype, Start/End back-converted to bp.
  expect_equal(
    p$tx$hover,
    c(
      "AAA<br>Gene ID: ENSG1<br>Biotype: protein_coding<br>Start: 1e+06<br>End: 1500000",
      "BBB<br>Gene ID: ENSG2<br>Biotype: lincRNA<br>Start: 2e+06<br>End: 2500000"
    )
  )

  # hover is independent of label: italics/strand arrows must not leak in,
  # and a blank gene name (which blanks the label) must not blank the hover.
  expect_false(any(grepl("&#8594;|&#8592;", p$tx$hover)))
})

test_that("genetrack_payload's hover column tolerates a missing fullname column", {
  # genetrack_ly()'s TX has no `fullname` column in ordinary use (it comes
  # from ensembldb::genes(), which never sets one), so TX$fullname is NULL
  # there too; paste0() with a NULL argument contributes nothing rather than
  # producing character(0).
  TX <- data.frame(
    gene_id = "ENSG1", gene_name = "AAA", gene_biotype = "protein_coding",
    start = 1, end = 2, strand = "+", stringsAsFactors = FALSE
  )
  expect_null(TX$fullname)
  EX <- data.frame(gene_id = character(), start = numeric(), end = numeric())
  p <- genetrack_payload(TX, EX, cfg = list(italics = FALSE))

  expect_equal(p$tx$hover,
               "AAA<br>Gene ID: ENSG1<br>Biotype: protein_coding<br>Start: 1e+06<br>End: 2e+06")
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
