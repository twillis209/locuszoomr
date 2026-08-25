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
