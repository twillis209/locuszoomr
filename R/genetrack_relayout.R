# Client-side gene track re-layout.
#
# The gene annotation panel is packed once by mapRow() at render time. These
# helpers ship the gene and exon tables into the widget so that JavaScript can
# re-pack them against the visible window whenever the user zooms or pans.


#' Build the gene track payload passed to JavaScript
#'
#' @param TX Data frame of transcripts, already ordered by [mapRow()] and with
#'   `start`/`end` in Mb.
#' @param EX Data frame of exons with `gene_id`, `start`, `end` in Mb.
#' @param cfg Named list of render settings. `italics` is required; the caller
#'   adds the remaining display settings.
#' @return A list with elements `tx`, `ex` and `cfg`.
#' @noRd
genetrack_payload <- function(TX, EX, cfg) {
  nm <- TX$gene_name
  if (isTRUE(cfg$italics)) nm <- paste0("<i>", nm, "</i>")
  pos <- as.character(TX$strand) == "+"
  label <- ifelse(pos, paste0(nm, "&#8594;"), paste0("&#8592;", nm))
  label[TX$gene_name == ""] <- ""

  tx <- data.frame(
    id       = as.character(TX$gene_id),
    name     = as.character(TX$gene_name),
    label    = label,
    start    = as.numeric(TX$start),
    end      = as.numeric(TX$end),
    strand   = as.character(TX$strand),
    priority = seq_len(nrow(TX)),
    stringsAsFactors = FALSE
  )

  gene_idx <- match(as.character(EX$gene_id), tx$id) - 1L
  keep <- !is.na(gene_idx)
  ex <- data.frame(
    gene_idx = as.integer(gene_idx[keep]),
    start    = as.numeric(EX$start[keep]),
    end      = as.numeric(EX$end[keep]),
    stringsAsFactors = FALSE
  )

  list(tx = tx, ex = ex, cfg = cfg)
}


#' Locate the gene track traces and shapes within a built plotly object
#'
#' `subplot()` renames axes and merges every panel's shapes into one array, so
#' positions cannot be assumed. The two gene track traces carry `meta` tags; the
#' exon shapes are then identified by sharing those traces' y axis reference.
#'
#' Note: `meta` is set as a plain (non-`I()`) value inside a data-bound
#' `add_segments()`/`add_text()` pipeline, so `plotly_build()` recycles it to
#' match each trace's row count rather than keeping it scalar; `add_segments()`
#' additionally interleaves `NA` at its line-break separators. Matching is
#' therefore done against any non-`NA` element rather than requiring a scalar.
#'
#' @param built A plotly object that has been through [plotly::plotly_build()].
#' @return A list of 0-based indices and axis identifiers.
#' @noRd
resolve_genetrack_idx <- function(built) {
  has_tag <- function(d, tag) {
    m <- d$meta
    if (is.null(m) || !is.character(m)) return(FALSE)
    any(m == tag, na.rm = TRUE)
  }

  li <- which(vapply(built$x$data, has_tag, logical(1),
                      tag = "locuszoomr_genetrack_lines"))
  ti <- which(vapply(built$x$data, has_tag, logical(1),
                      tag = "locuszoomr_genetrack_labels"))
  if (length(li) != 1L || length(ti) != 1L) {
    stop("could not locate gene track traces in the built plotly object",
         call. = FALSE)
  }

  xax <- built$x$data[[li]]$xaxis
  yax <- built$x$data[[li]]$yaxis
  if (is.null(xax)) xax <- "x"
  if (is.null(yax)) yax <- "y"

  shapes <- built$x$layout$shapes
  shp <- if (length(shapes) == 0L) integer(0) else {
    which(vapply(shapes, function(s) identical(s$yref, yax), logical(1)))
  }

  list(lineTrace  = as.integer(li - 1L),
       labelTrace = as.integer(ti - 1L),
       shapeIdx   = as.integer(shp - 1L),
       xaxis      = xax,
       yaxis      = yax)
}
