test_that("recomb_table() picks build-specific defaults", {
  expect_equal(recomb_table("hg38"), "recomb1000GAvg")
  expect_equal(recomb_table("hg19"), "hapMapRelease24CombinedRecombMap")
  expect_equal(recomb_table("hg38", "recomb1000GCEU"), "recomb1000GCEU")
  expect_error(recomb_table("mm10"), "specify `table`")
})


test_that("ucsc_track_url() builds a UCSC REST API query", {
  url <- ucsc_track_url("hg38", "recomb1000GAvg", "chr7", c(1234567, 1334567))
  expect_equal(
    url,
    paste0("https://api.genome.ucsc.edu/getData/track?genome=hg38",
           ";track=recomb1000GAvg;chrom=chr7;start=1234567;end=1334567"))
})


test_that("ucsc_track_url() does not use scientific notation", {
  url <- ucsc_track_url("hg38", "recomb1000GAvg", "chr1", c(1e6, 2.5e6))
  expect_match(url, "start=1000000;end=2500000", fixed = TRUE)
  expect_false(grepl("e+", url, fixed = TRUE))
})


test_that("ucsc_track_url() clamps negative start values", {
  url <- ucsc_track_url("hg38", "recomb1000GAvg", "chr1", c(-5000, 10000))
  expect_match(url, "start=0;end=10000", fixed = TRUE)
})


test_that("parse_ucsc_track() extracts the track data frame", {
  json <- '{"genome":"hg38","track":"recomb1000GAvg","chrom":"chr7",
    "recomb1000GAvg":[
      {"chrom":"chr7","start":1234567,"end":1239976,"value":1.22503},
      {"chrom":"chr7","start":1239976,"end":1245559,"value":1.22639}],
    "itemsReturned":2}'
  res <- parse_ucsc_track(json, "recomb1000GAvg")
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("chrom", "start", "end", "value"))
  expect_equal(nrow(res), 2L)
  expect_equal(res$start, c(1234567, 1239976))
  expect_equal(res$value, c(1.22503, 1.22639))
})


test_that("parse_ucsc_track() surfaces the UCSC error message", {
  json <- '{"error":"can not find track=nosuchtrack name for endpoint",
    "statusCode":400,"statusMessage":"Bad Request"}'
  expect_error(parse_ucsc_track(json, "nosuchtrack"),
               "can not find track")
})


test_that("parse_ucsc_track() returns 0 rows when the range holds no data", {
  json <- '{"genome":"hg38","track":"recomb1000GAvg","recomb1000GAvg":[],
    "itemsReturned":0}'
  res <- parse_ucsc_track(json, "recomb1000GAvg")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
  expect_equal(colnames(res), c("chrom", "start", "end", "value"))
})


test_that("parse_ucsc_track() errors if the track element is missing", {
  json <- '{"genome":"hg38","track":"recomb1000GAvg","itemsReturned":0}'
  expect_error(parse_ucsc_track(json, "recomb1000GAvg"),
               "no 'recomb1000GAvg' data")
})


test_that("parse_ucsc_track() warns when UCSC truncates the response", {
  json <- '{"genome":"hg38","track":"recomb1000GAvg",
    "recomb1000GAvg":[{"chrom":"chr7","start":1,"end":2,"value":1}],
    "maxItemsReached":true,"itemsReturned":1}'
  expect_warning(parse_ucsc_track(json, "recomb1000GAvg"), "truncated")
})


test_that("recomb_bw_url() maps builds to their gbdb folders", {
  expect_equal(
    recomb_bw_url("hg38", "recomb1000GAvg"),
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recomb1000GAvg.bw")
  expect_equal(
    recomb_bw_url("hg19", "hapMapRelease24CombinedRecombMap"),
    paste0("https://hgdownload.soe.ucsc.edu/gbdb/hg19/decode/",
           "hapMapRelease24CombinedRecombMap.bw"))
  expect_null(recomb_bw_url("mm10", "foo"))
})


test_that("link_recomb() attaches locally supplied GRanges data", {
  skip_if_not_installed("GenomicRanges")
  rec <- GenomicRanges::GRanges(
    seqnames = c("chr7", "chr7", "chr8"),
    ranges = IRanges::IRanges(start = c(1000, 2000, 1000),
                              end = c(1999, 2999, 1999)),
    score = c(1.5, 2.5, 9.9))
  loc <- structure(
    list(seqname = "7", xrange = c(1500, 2500), genome = "hg38"),
    class = "locus")
  out <- link_recomb(loc, recomb = rec)
  expect_equal(colnames(out$recomb), c("start", "end", "value"))
  expect_equal(out$recomb$value, c(1.5, 2.5))
})


test_that("query_recomb() falls back to the bigWig when the REST API fails", {
  fallback <- data.frame(chrom = "chr7", start = 1, end = 2, value = 1.5)
  local_mocked_bindings(
    query_recomb_api = function(...) stop("connection refused"),
    query_recomb_bw = function(...) fallback)
  expect_message(res <- query_recomb("hg38", c(1, 2), "7"),
                 "used bigWig fallback")
  expect_equal(res, fallback)
})


test_that("query_recomb() returns NULL when both routes fail", {
  local_mocked_bindings(
    query_recomb_api = function(...) stop("connection refused"),
    query_recomb_bw = function(...) stop("bigWig unreachable"))
  expect_message(res <- query_recomb("hg38", c(1, 2), "7"),
                 "connection refused")
  expect_null(res)
})


test_that("query_recomb() returns NULL for a build with no default table", {
  expect_message(res <- query_recomb("mm10", c(1, 2), "7"),
                 "no default recombination table")
  expect_null(res)
})


test_that("query_recomb() prefixes bare chromosome names", {
  local_mocked_bindings(
    query_recomb_api = function(gen, table, seqname, xrange) seqname,
    query_recomb_bw = function(...) stop("should not be called"))
  expect_equal(suppressMessages(query_recomb("hg38", c(1, 2), "7")), "chr7")
  expect_equal(suppressMessages(query_recomb("hg38", c(1, 2), "chr7")), "chr7")
})


test_that("link_recomb() queries UCSC over the network", {
  skip_on_cran()
  skip_if_offline("api.genome.ucsc.edu")
  loc <- structure(
    list(seqname = "7", xrange = c(1234567, 1334567), genome = "hg38"),
    class = "locus")
  out <- link_recomb(loc)
  expect_s3_class(out$recomb, "data.frame")
  expect_true(nrow(out$recomb) > 0)
  expect_true(all(c("start", "end", "value") %in% colnames(out$recomb)))
  expect_equal(out$recomb$start[1], 1234567)
  expect_true(all(out$recomb$value >= 0))
})


test_that("the bigWig fallback returns the same data as the REST API", {
  skip_on_cran()
  skip_if_offline("hgdownload.soe.ucsc.edu")
  skip_if_not_installed("rtracklayer")
  xrange <- c(1234567, 1334567)
  api <- query_recomb_api("hg38", "recomb1000GAvg", "chr7", xrange)
  bw <- query_recomb_bw("hg38", "recomb1000GAvg", "chr7", xrange)
  expect_equal(nrow(api), nrow(bw))
  expect_equal(api$start, bw$start)
  expect_equal(api$end, bw$end)
  # the REST API rounds values to 6 significant figures, the bigWig does not
  expect_equal(api$value, bw$value, tolerance = 1e-5)
})
