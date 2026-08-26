
#' Render a window's datapoints as a DT table
#'
#' Rounds the measurements to 3 significant figures and leaves the coordinates
#' alone. Coordinates are numeric but they are labels, not measurements:
#' `formatSignif()` over every numeric column turns chromosome 19 into "19.0"
#' and position 44905791 into "4.49e+07", which is both unreadable and no
#' longer copy-pasteable back into the search box.
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
  # formatSignif() errors on an empty column selection, which is reachable:
  # a data frame of rsID/CHR/BP and nothing else leaves no measurements.
  if (length(cols)) out <- formatSignif(out, cols, digits = digits)
  out
}
