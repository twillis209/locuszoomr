# LDlinkR reports API-level failures by RETURNING a data frame that holds the
# message, rather than by raising a condition. try() therefore sees nothing
# wrong, and link_LD() went on to subset columns that were not there:
#
#   ldp[, "RS_Number"]  ->  Error: undefined columns selected
#
# This is the shape LDproxy() actually returned for a variant absent from the
# 1000G panel (captured from a live call, rs747519137, genome_build grch38).
ld_error_frame <- function(
    msg = "  error: rs747519137 Variant is not in 1000G reference panel.") {
  d <- data.frame(x = msg, stringsAsFactors = FALSE)
  colnames(d) <- "data_out[1, 1]"
  d
}

ld_good_frame <- function() {
  data.frame(RS_Number = c("rs1", "rs2"),
             Coord = c("chr1:100", "chr1:200"),
             R2 = c(1, 0.8),
             stringsAsFactors = FALSE)
}

test_that("a well formed LDproxy response is accepted", {
  expect_true(ld_response_ok(ld_good_frame(), c("RS_Number", "R2")))
  expect_true(ld_response_ok(ld_good_frame(), c("Coord", "R2")))
})

test_that("an LDlinkR error frame is rejected rather than subset", {
  expect_false(ld_response_ok(ld_error_frame(), c("RS_Number", "R2")))
  expect_false(ld_response_ok(ld_error_frame(), c("Coord", "R2")))
})

test_that("responses missing only the R2 column are rejected", {
  d <- ld_good_frame()
  d$R2 <- NULL
  expect_false(ld_response_ok(d, c("RS_Number", "R2")))
})

test_that("empty and non-data-frame responses are rejected", {
  expect_false(ld_response_ok(ld_good_frame()[0, ], c("RS_Number", "R2")))
  expect_false(ld_response_ok(NULL, c("RS_Number", "R2")))
  expect_false(ld_response_ok("error", c("RS_Number", "R2")))
})

test_that("the API message is recovered for reporting", {
  expect_match(ld_response_error(ld_error_frame()),
               "not in 1000G reference panel")
  # the "error:" prefix and surrounding whitespace are stripped
  expect_false(grepl("^\\s*error", ld_response_error(ld_error_frame())))
})

test_that("a good response yields no error message", {
  expect_true(is.na(ld_response_error(ld_good_frame())))
})
