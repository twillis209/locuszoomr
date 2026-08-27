
#' Render a window's datapoints as a DT table
#'
#' Rounds measurements to 3 significant figures but leaves the coordinate
#' columns alone: they are labels, and `formatSignif()` would render position
#' 44905791 as 4.49e+07.
#'
#' @param d Data frame of datapoints, already subset and ordered.
#' @param coord_cols Column names to leave unformatted - the chromosome and
#'   position columns resolved by [zoom()].
#' @param digits Significant figures for the remaining numeric columns.
#' @return A 'datatables' htmlwidget.
#' @noRd
locus_table <- function(d, coord_cols, digits = 3) {
  cols <- colnames(d)[vapply(d, is.numeric, logical(1))]
  cols <- setdiff(cols, coord_cols)
  out <- datatable(d, rownames = FALSE)
  # formatSignif() errors on an empty selection (e.g. rsID/CHR/BP only)
  if (length(cols)) out <- formatSignif(out, cols, digits = digits)
  out
}
