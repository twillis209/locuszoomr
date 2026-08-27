
# Position indexing for zoom(). At GWAS scale, filtering by chromosome and
# position scans tens of millions of rows per navigation; sorting once at
# startup makes every window a contiguous row range found by binary search.
# Confined to zoom(): locus() is unchanged and simply receives pre-windowed
# rows.


#' Sort a dataset by locus and record each chromosome's row range
#'
#' @param data Data frame of GWAS results.
#' @param chrom,pos Column names for chromosome and position.
#' @return A list with `data` (sorted) and `index`, a data frame of
#'   `chrom`/`first`/`last` giving each chromosome's contiguous row range.
#'   Chromosomes are compared as character, matching how `zoom()` carries
#'   `coords$chr`.
#' @noRd
build_locus_index <- function(data, chrom, pos) {
  # radix: orders of magnitude faster at this scale, and locale-independent
  o <- order(as.character(data[, chrom]), data[, pos], method = "radix")
  data <- data[o, , drop = FALSE]
  ch <- as.character(data[, chrom])
  first <- which(!duplicated(ch))
  list(data = data,
       index = data.frame(chrom = ch[first],
                          first = first,
                          last = c(first[-1] - 1L, length(ch)),
                          stringsAsFactors = FALSE))
}


#' Row numbers covering one window, by binary search
#'
#' Bounds are inclusive - laxer than `locus()`'s strict `>` / `<` - so
#' boundary variants are trimmed by `locus()` itself and the indexed path
#' matches the full scan exactly.
#'
#' @param idx The list returned by [build_locus_index()].
#' @param pos Column name for position.
#' @param seqname Chromosome, compared as character.
#' @param xrange Length-2 numeric window.
#' @return Integer vector of row numbers, empty if the window holds nothing.
#' @noRd
locus_rows <- function(idx, pos, seqname, xrange) {
  hit <- match(as.character(seqname), idx$index$chrom)
  if (is.na(hit)) return(integer(0))
  lo <- idx$index$first[hit]
  hi <- idx$index$last[hit]
  p <- idx$data[lo:hi, pos]
  # findInterval() returns the count of values <= x, hence +1 on the left
  i <- findInterval(xrange[1], p, left.open = TRUE) + 1L
  j <- findInterval(xrange[2], p)
  if (i > j) return(integer(0))
  seq.int(lo + i - 1L, lo + j - 1L)
}


#' The window's rows, as a data frame
#'
#' @inheritParams locus_rows
#' @return Data frame, possibly with zero rows.
#' @noRd
locus_subset <- function(idx, pos, seqname, xrange) {
  idx$data[locus_rows(idx, pos, seqname, xrange), , drop = FALSE]
}
