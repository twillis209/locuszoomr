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
