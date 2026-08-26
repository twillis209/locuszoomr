
# Position indexing for zoom()
#
# zoom() navigates by (chromosome, window) many times over the life of one
# session, and every navigation used to scan the whole dataset:
#
#   data[which(data[, chrom] == seqname), ]                    # locus():171
#   data[which(data[, pos] > x1 & data[, pos] < x2), ]         # locus():172
#
# On a 21M row GWAS each `==` allocates an 84 MB logical vector, `which()` an
# integer index, and the subset a fresh data frame - per trait, per
# navigation, plus the same scan again in output$table. Measured, the plot
# build takes 0.15s most of the time and ~1.6s whenever that allocation
# triggers a major collection, which is what a user perceives as lag.
#
# Sorting by (chromosome, position) once at startup makes every window a
# CONTIGUOUS range of rows, reachable by binary search. The scans go away and
# so does the garbage they generate.
#
# Deliberately confined to zoom(). locus() keeps its own filters - it has
# nowhere to hold an index across calls, and it is Myles's function.
# Handed a pre-windowed frame, its filters run over a few thousand rows and
# cost nothing.


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
  # method = "radix" rather than the default: on 21M rows the default shell
  # sort is minutes, radix is seconds. It also fixes the collation order for
  # character chromosomes, so the index does not depend on the locale.
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
#' Bounds are INCLUSIVE, which is deliberately laxer than `locus()`'s strict
#' `>` / `<`. The rows are handed to `locus()`, which re-applies its own
#' filter, so a variant sitting exactly on a boundary is dropped there rather
#' than here - and the indexed path returns exactly what the scan did.
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
  # findInterval() needs a sorted vector, which p is by construction.
  # +1 on the left because findInterval() returns the count of values <= x,
  # i.e. the index of the last value at or below the window's start.
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
