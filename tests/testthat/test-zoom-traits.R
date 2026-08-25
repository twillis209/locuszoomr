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
