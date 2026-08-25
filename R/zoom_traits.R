#' Check a second trait's data frame carries the resolved columns
#'
#' Both traits share one set of column names resolved from the first, so a
#' mismatch has to be caught before the app launches; otherwise it surfaces
#' much later as an empty panel.
#'
#' @param data2 Data frame to check.
#' @param cols Character vector of required column names.
#' @param name Argument name to use in the error message.
#' @return `invisible(TRUE)`, or an error naming the missing columns.
#' @noRd
check_trait_cols <- function(data2, cols, name = "data2") {
  missing <- setdiff(cols, colnames(data2))
  if (length(missing) > 0) {
    stop("`", name, "` is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Warn when two traits look like different genome builds
#'
#' Matching on SNP id and comparing positions is what makes this a build
#' test rather than a coverage test: two GWAS may legitimately contain
#' different variants, but the same id at a different coordinate means
#' different builds. There is no liftover, and mismatched builds would
#' misalign the panels while looking entirely plausible.
#'
#' Warns rather than stops - a small disagreement can be legitimate, e.g.
#' multi-allelic normalisation, and the user may know better than this check.
#'
#' @param data,data2 The two datasets.
#' @param pos,labs Resolved position and SNP id column names.
#' @param min_shared Minimum shared ids needed to judge at all.
#' @param min_agree Proportion of shared ids that must agree on position.
#' @return `invisible(TRUE)`.
#' @noRd
check_same_build <- function(data, data2, pos, labs,
                             min_shared = 1000L, min_agree = 0.9) {
  i <- match(data2[, labs], data[, labs])
  ok <- !is.na(i)
  n <- sum(ok)
  if (n < min_shared) {
    message("Skipping genome build check: too few shared SNP ids (", n, ")")
    return(invisible(TRUE))
  }
  agree <- mean(data[i[ok], pos] == data2[ok, pos])
  if (agree < min_agree) {
    warning("`data` and `data2` may be on different genome builds: only ",
            round(agree * 100), "% of ", n,
            " shared SNP ids agree on position", call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve panel labels for the two traits
#'
#' `deparse(substitute(data))` gives a good label for `zoom(ad, dizzy)` but
#' an unusable one for `zoom(read_gwas("x.tsv"), df2)`, so anything that is
#' not a plain short name falls back to a positional label.
#'
#' @param trait_names User-supplied labels, or `NULL`.
#' @param expr1,expr2 Deparsed argument expressions.
#' @return Character vector of length 2.
#' @noRd
trait_labels <- function(trait_names, expr1, expr2) {
  if (!is.null(trait_names)) {
    return(rep_len(as.character(trait_names), 2L))
  }
  usable <- function(x) {
    length(x) == 1L && nchar(x) <= 20L && grepl("^[A-Za-z.][A-Za-z0-9._]*$", x)
  }
  c(if (usable(expr1)) expr1 else "Trait 1",
    if (usable(expr2)) expr2 else "Trait 2")
}
