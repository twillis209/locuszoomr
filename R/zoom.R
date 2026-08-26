
#' Zoom browser to explore GWAS/eQTL results
#' 
#' Interactive genome browser to explore GWAS/eQTL results using a shiny
#' interface.
#' 
#' @details 
#' This launches a shiny app to explore the GWAS/eQTL results through visualising 
#' the Manhattan plot and exploring regional Manhattan plots of gene loci 
#' through selecting points or searching SNPs/genes. 
#' @param data Dataframe of GWAS results with columns for chromosome, position,
#'   p value and SNP rs IDs. Data.tables are coerced to dataframe.
#' @param ens_db Either a character string which specifies which Ensembl
#'   database package (version 86 and earlier for Homo sapiens) to query for
#'   gene and exon positions (see `ensembldb` Bioconductor package). Or an
#'   `ensembldb` object which can be obtained from the AnnotationHub database.
#'   See the vignette and the `AnnotationHub` Bioconductor package for how to
#'   create this object.
#' @param chrom Determines which column in `data` contains chromosome
#'   information If `NULL` tries to autodetect the column.
#' @param pos Determines which column in `data` contains position information.
#'   If `NULL` tries to autodetect the column.
#' @param p Determines which column in `data` contains SNP p-values. If `NULL`
#'   tries to autodetect the column.
#' @param labs Determines which column in `data` contains SNP rs IDs. If `NULL`
#'   tries to autodetect the column.
#' @param scheme Vector of 3 colours: 1st = normal points, 2nd = colour for
#'   significant points, 3rd = index SNP(s).
#' @param pcutoff Cut-off for p value significance. Defaults to p = 5e-08. Set
#'   to `NULL` to disable.
#' @param eqtl_gene Determines which column in `data` contains eQTL genes.
#' @param eqtl_beta Optional column name for beta coefficient to display upward
#'   triangles for positive beta and downward triangles for negative beta
#'   (significant SNPs only).
#' @param eqtl_scheme Colour scheme for eQTL genes.
#' @param add_hover Optional vector of column names in 'data' to add to the
#'   plotly hover text for scatter points.
#' @param mh_points Number of points to display in manhattan plot. Default is
#'   `1e5`.
#' @param recomb Optional `GRanges` class object of recombination data.
#' @param data2 Optional second dataframe of GWAS results, shown as an extra
#'   panel over the same window so two traits can be compared. Must use the
#'   same column names as `data`; there is no second column mapping. Both
#'   datasets must be on the same genome build, which is checked at startup.
#' @param trait_names Optional length-2 character vector labelling the two
#'   panels. Defaults to the deparsed argument names, so `zoom(ad, dizzy)`
#'   labels itself.
#' @param ld_token Personal access token for the LDlink API, available from
#'   <https://ldlink.nih.gov/?tab=apiaccess>. When empty the LD controls are
#'   hidden. LD is fetched on demand, not automatically: pressing "Get LD"
#'   pins the current index SNP as the reference variant and colours points by
#'   r^2 with it. The reference stays pinned while you pan and zoom, so the
#'   colouring keeps its meaning and repeat queries are served from the
#'   `memoise` cache rather than the API.
#' @param ld_pop 1000 Genomes population used for LD. Defaults to `"EUR"`. See
#'   `LDlinkR::LDproxy()` for the available codes.
#' @param AnnotationDb An `AnnotationDb` gene annotation database, specified
#'   either as a character string or as an `AnnotationDb` class object, used to
#'   obtain expanded gene names. The ensembl database specified in `ens_db` is
#'   queried first. Set to `NULL` to disable this feature.
#' @returns No return value. Opens an interactive shiny window.
#' @importFrom plotly plotlyOutput renderPlotly event_data config plotlyProxy
#' @importFrom plotly plotlyProxyInvoke layout
#' @importFrom shiny fluidPage tabsetPanel tabPanel fluidRow column actionButton 
#' @importFrom shiny icon uiOutput checkboxInput textOutput splitLayout req
#' @importFrom shiny textInput conditionalPanel h5 runApp debounce isolate
#' @importFrom shiny renderUI reactiveValues reactive observe observeEvent radioButtons
#' @importFrom shiny reactiveVal validate need renderText updateTextInput outputOptions
#' @importFrom shiny showNotification removeNotification
#' @importFrom shinyFeedback useShinyFeedback hideFeedback showFeedback
#' @importFrom shinyWidgets pickerInput pickerOptions dropdown
#' @importFrom shinycssloaders withSpinner
#' @importFrom htmltools tags br
#' @importFrom DT datatable formatSignif
#' @importFrom gtools mixedsort
#' @importFrom stats as.formula
#' @export

zoom <- function(data, ens_db,
                 chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                 scheme = c('royalblue', 'skyblue', 'red'),
                 pcutoff = 5e-8,
                 eqtl_gene = NULL,
                 eqtl_beta = NULL,
                 eqtl_scheme = c("#FF0000", "#00FFFF", "#FF9000", "#0080FF", "#FFFF00",
                                 "#0000FF", "#80DD00", "#8000FF", "#009900", "#FF00FF"),
                 add_hover = NULL,
                 mh_points = 1e5,
                 recomb = NULL,
                 ld_token = Sys.getenv("LDLINK_TOKEN"),
                 ld_pop = "EUR",
                 data2 = NULL,
                 trait_names = NULL,
                 AnnotationDb = "org.Hs.eg.db") {
  # Captured before `data` is reassigned below: once a formal's binding is
  # overwritten, substitute() returns the current value rather than the
  # caller's expression, which would deparse the entire dataset.
  #
  # Only a bare symbol is kept, and deparse() is never reached for anything
  # else. Two failure modes both follow from deparsing whatever arrives:
  #
  #  - do.call(zoom, list(data = big_df, ...)) without quote = TRUE splices
  #    the VALUE into the call, so substitute() hands back the data frame
  #    itself and deparse() serialises every row. A 50,000 row frame gives a
  #    12,505-element character vector; on a 21M row GWAS it is a
  #    multi-second stall on every call, single-trait ones included.
  #  - deparse() returns a VECTOR, so c() flattens the two results together.
  #    Once trait 1's expression needs more than one line (~60 chars),
  #    trait_expr[2] is a continuation line of trait 1 rather than trait 2's
  #    expression at all.
  #
  # Labels only ever want a plain short name, so bound the capture rather
  # than the string: anything else becomes "", which trait_labels() sends to
  # the positional fallback.
  sym <- function(e) if (is.name(e)) as.character(e) else ""
  trait_expr <- c(sym(substitute(data)), sym(substitute(data2)))
  data <- data.frame(data)
  # autodetect headings
  dc <- detect_cols(data, chrom, pos, p, labs)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  labs <- dc$labs
  # Resolved before any of the single-trait setup below, so a mismatched
  # second dataset fails now rather than as an empty panel later.
  trait_lab <- trait_labels(trait_names, trait_expr[1], trait_expr[2])
  if (!is.null(data2)) {
    data2 <- data.frame(data2)
    # add_hover is included because it is forwarded to BOTH panels:
    # scatter_plotly() does data[, i] for each name, so a column present in
    # `data` but not `data2` errors with "undefined columns selected" on the
    # first render - the late failure this startup check exists to prevent.
    check_trait_cols(data2, c(chrom, pos, p, labs, add_hover), "data2")
    check_same_build(data, data2, pos, labs)
    if (!is.null(eqtl_gene) || !is.null(eqtl_beta)) {
      # The two-trait panels share a single colour scheme sized for LD/
      # default colouring (3 levels), not for one-colour-per-eQTL-gene.
      # Rather than error, output$locus falls back to the default scheme
      # for both panels when data2 is set - so warn instead of silently
      # dropping the requested eQTL colouring. `eqtl_beta`'s up/down
      # triangles go the same way: compose_locus_plotly() is not passed
      # `beta` either.
      warning("eQTL colouring and beta-direction markers are disabled when ",
              "data2 is supplied: the two trait panels share one colour ",
              "scheme, so `eqtl_gene`/`eqtl_beta` are ignored", call. = FALSE)
    }
  }
  if (is.null(eqtl_gene)) {
    data[, labs] <- unique_snps(data, labs, chrom)
  } else {
    data[, labs] <- unique_snps(data, labs, eqtl_gene)
  }
  
  message("Generating Manhattan plot")
  # Union, not trait 1's alone: output$locus gates on `coords$chr %in%
  # chr_set`, so a chromosome present only in trait 2 would make clicking
  # its own manhattan do nothing at all. Trait 1 simply renders empty there,
  # which is honest and already handled.
  chr_set <- unique(data[, chrom])
  if (!is.null(data2)) chr_set <- union(chr_set, unique(data2[, chrom]))
  if (is.character(ens_db)) {
    if (!ens_db %in% (.packages())) {
      stop("Ensembl database not loaded. Try: library(", ens_db, ")",
           call. = FALSE)
    }
    edb <- get(ens_db)
  } else edb <- ens_db
  
  gene_db <- genes(edb, filter = AnnotationFilterList(
    SeqNameFilter(c(1:22, 'X', 'Y'))))
  gene_set <- unique(gene_db$gene_name)
  biotypes <- sort(unique(gene_db$gene_biotype))
  
  # lookup table for full length gene names using org.Hs.eg.db
  fullnames <- fullGeneNames(edb, AnnotationDb)
   
  if (!is.null(eqtl_gene)) {
    message("Setting eQTL colours")
    dat2 <- data[data[, p] < pcutoff, ]
    dat2 <- dat2[order(dat2[, chrom], dat2[, pos]), ]
    eqtl_set <- unique(dat2[, eqtl_gene])
    colours <- eqtl_scheme
    eqtl_colour <- rep_len(colours, length(eqtl_set))
    names(eqtl_colour) <- eqtl_set
    message(length(eqtl_set), " eQTL genes")
    rm(dat2, eqtl_set)
  }
  
  # apply min_p_snp to data for manhat?
  # smallest floating point
  data[which(data[, p] < 5e-324), p] <- 5e-324
  # Sort by (chromosome, position) once, so every window the user navigates
  # to is a contiguous run of rows found by binary search rather than a scan
  # of the whole dataset. See R/zoom_index.R. `data` is replaced by the
  # sorted copy: nothing downstream depends on the incoming row order
  # (manhattan() re-sorts, the index SNP is a which.max, the table orders by
  # p), and leaving both around would double the memory.
  data_idx <- build_locus_index(data, chrom, pos)
  data <- data_idx$data
  manhat <- manhattan(data, chrom, pos, p, labs, pcutoff = pcutoff,
                      npoints = mh_points)
  yrange <- range(manhat$data$logP, na.rm = TRUE)
  ymax <- yrange[2] + diff(yrange) * 0.05

  # Second genome-wide manhattan, so both traits can be scanned for
  # coinciding peaks before drilling in. Thinned independently: each trait's
  # own top mh_points SNPs are the interesting ones, and thinning trait 2 by
  # trait 1's selection would hide exactly the signals that differ.
  manhat2 <- NULL
  data2_idx <- NULL
  if (!is.null(data2)) {
    data2[which(data2[, p] < 5e-324), p] <- 5e-324
    data2_idx <- build_locus_index(data2, chrom, pos)
    data2 <- data2_idx$data
    manhat2 <- manhattan(data2, chrom, pos, p, labs, pcutoff = pcutoff,
                         npoints = mh_points)
    yrange2 <- range(manhat2$data$logP, na.rm = TRUE)
  }
  
  # LD needs a token, and its colouring would replace the eQTL gene colours
  # (scatter_plotly() switches to LD_scheme whenever an `ld` column exists),
  # so it is offered only when neither applies.
  show_ld <- nzchar(ld_token) && is.null(eqtl_gene)

  js <- '$(document).on("keydown", function(e) {
          if(e.keyCode == 13) {
            Shiny.onInputChange("enter", Math.random());
          }
        });'
  
  # https://shiny.posit.co/r/articles/build/packaging-javascript/
  
  ui <- fluidPage(
    tags$script(js),
    # 3 plotly scattergl figures gives error "too many active WebGL contexts"
    # see https://plotly.com/python/webgl-vs-svg/
    tags$script(src = "https://unpkg.com/virtual-webgl@1.0.6/src/virtual-webgl.js"),
    # includeScript(path = "/users/myles/dropbox/R scripts/locuszoomr/www/virtual-webgl.js"),
    useShinyFeedback(),
    tabsetPanel(
      tabPanel("Plot",
               fluidRow(
                 column(11,
                        # 300px for one trait; 220 each for two. Two full
                        # height strips would put 600px of manhattan above a
                        # locus panel that is itself 850px in two-trait mode,
                        # so the thing you navigated to would start below the
                        # fold. 220px still reads at genome scale.
                        withSpinner(
                          plotlyOutput("manhattan", width = "85vw",
                                       height = if (is.null(data2)) "300px" else "220px"),
                          type = 8, size = 0.7),
                        (if (!is.null(data2)) {
                          withSpinner(
                            plotlyOutput("manhattan2", width = "85vw",
                                         height = "220px"),
                            type = 8, size = 0.7)
                        } else NULL)
                 ),
                 column(1,
                        br(),
                        # One pair of buttons drives both strips - the point
                        # of stacking them is to read them together, so
                        # zooming one and not the other would defeat it.
                        actionButton("m_zoomin", NULL, icon = icon("magnifying-glass-plus")),
                        actionButton("m_zoomout", NULL, icon = icon("magnifying-glass-minus"))
                 )),
               fluidRow(
                 column(12,
                        uiOutput("ui_chrom")
                 )
               ),
               fluidRow(
                 column(3,
                        checkboxInput("show_chrom", "show chromosome")
                 )
               ),
               fluidRow(
                 column(4,
                        actionButton("left2", NULL, icon = icon("angles-left")),
                        actionButton("left", NULL, icon = icon("angle-left")),
                        actionButton("right", NULL, icon = icon("angle-right")),
                        actionButton("right2", NULL, icon = icon("angles-right")),
                        actionButton("zoomin", NULL, icon = icon("magnifying-glass-plus")),
                        actionButton("zoomout", NULL, icon = icon("magnifying-glass-minus"))
                        ),
                 column(3,
                        textOutput("pos"),
                        # LD lives behind the gear, so surface the pinned
                        # reference variant here: otherwise the colouring
                        # changes meaning with no visible reason why.
                        textOutput("ld_status"),
                        align = "centre", style='margin-top:7px;'),
                 column(4,
                        splitLayout(
                          textInput("tex", NULL, placeholder = "chr:start-end, rs or gene",
                                    width = "100%"),
                          actionButton("text_go", NULL, icon = icon("magnifying-glass"),
                                       class = "btn-success"),
                          cellWidths = c("75%", "25%")
                        )),
                 column(1,
                        dropdown(
                          # Always offered now. With `recomb` supplied the
                          # rates come from that object; without it
                          # link_recomb() queries the UCSC REST API per
                          # window, which costs about 0.2s and is memoised.
                          # Defaulted on only for supplied data, since that
                          # is free and was the previous behaviour, and off
                          # for the API so nothing goes over the network
                          # unasked.
                          checkboxInput("recomb", "show recombination rate",
                                        value = !is.null(recomb)),
                          # LD is an on-demand action, not a setting: each new
                          # reference variant costs an LDlink API call. Hidden
                          # without a token, and in eQTL mode, where the `ld`
                          # column would override the per-gene colours.
                          (if (show_ld) {
                            list(h5("Linkage disequilibrium"),
                                 actionButton("ld_get", "Get LD",
                                              icon = icon("circle-nodes"),
                                              class = "btn-primary btn-sm"),
                                 actionButton("ld_clear", "Clear",
                                              class = "btn-default btn-sm"),
                                 br(), br())
                          } else NULL),
                          pickerInput("biotype", h5("Select gene biotypes"),
                                      choices = biotypes, selected = biotypes,
                                      multiple = TRUE,
                                      options = pickerOptions(actionsBox = TRUE,
                                                              selectedTextFormat = 'count > 1')),
                          uiOutput("ui_genes"),
                          right = TRUE, icon = icon("gear")
                        ))
                 ),
                 fluidRow(
                   column(12,
                          # Taller with two traits: the gene track's share of
                          # the figure drops from 0.4 to 0.3 there, so at 600px
                          # its 12 rows would get ~15px each, below what the
                          # 9.8px labels need. 850px restores ~255px of gene
                          # track and still leaves ~300px per scatter panel.
                          # Single-trait geometry is unchanged.
                          plotlyOutput("locus", width = "95vw",
                                       height = if (is.null(data2)) 600 else 850),
                          br(), br()
                          # verbatimTextOutput("print")
                   )
                 )
      ),
      tabPanel("Table",
               fluidRow(
                 column(12, br(), DT::dataTableOutput("table"))))
    )
  )
  
  server <- function(input, output, session) {
    
    output$manhattan <- renderPlotly({
      p <- plotly_manhattan(manhat, labs, pcutline = NULL) %>%
        config(displayModeBar = FALSE)
      # Name the strip only when there are two, so the single-trait plot is
      # untouched.
      if (!is.null(data2)) {
        p <- p %>% layout(yaxis = list(title = paste0(trait_lab[1],
                                                      "  -log<sub>10</sub> P")))
      }
      p
    })

    output$manhattan2 <- renderPlotly({
      req(!is.null(manhat2))
      # Distinct source, so the click handler below can tell which trait was
      # clicked and look the SNP up in the right dataset.
      plotly_manhattan(manhat2, labs, pcutline = NULL,
                       source = "plotly_manh2") %>%
        layout(yaxis = list(title = paste0(trait_lab[2],
                                           "  -log<sub>10</sub> P"))) %>%
        config(displayModeBar = FALSE)
    })

    output$ui_chrom <- renderUI({
      req(input$show_chrom, coords$chr)
      fluidRow(
        column(11,
               withSpinner(
                 plotlyOutput("chrom", width = "85vw", height = "220px"),
                 type = 8, size = 0.7)
        ),
        column(1,
               br(),
               actionButton("chr_zoomin", NULL, icon = icon("magnifying-glass-plus")),
               actionButton("chr_zoomout", NULL, icon = icon("magnifying-glass-minus"))
        )
      )
    })
    
    output$chrom <- renderPlotly({
      req(coords$chr)
      # A whole chromosome is one contiguous run in the index, so this is a
      # slice rather than a scan of all 21M rows.
      chr_rows <- data_idx$index[match(as.character(coords$chr),
                                       data_idx$index$chrom), ]
      req(nrow(chr_rows) == 1L, !is.na(chr_rows$first))
      chr_manhat <- manhattan(data[chr_rows$first:chr_rows$last, ],
                              chrom, pos, p, labs, pcutoff = pcutoff,
                              npoints = 1e5)
      chr <- suppressWarnings(as.numeric(coords$chr))
      if ((!is.na(chr) && chr %% 2 == 0 || coords$chr == "Y")) {
        scheme[1] <- scheme[2]
      }
      yr <- range(chr_manhat$data$logP, na.rm = TRUE)
      isolate(chr_y$range <- yr)
      isolate(chr_y$max <- yr[2])
      isolate(xr <- view$xrange)
      
      plotly_manhattan(chr_manhat, labs, scheme = scheme,
                       source = "plotly_chrom") %>%
        layout(margin = list(t = 5),
               shapes = list(
                 list(type = "rect",
                      line = list(width = 1, color = "#00CD00"),
                      x0 = xr[1] / 1e6,
                      x1 = xr[2] / 1e6, y0 = 0, y1 = 1,
                      xref = "x", yref = "paper", layer = "below"))
               ) %>%
        config(displayModeBar = FALSE)
    })
    
    # `coords` is the window locus() was called for, i.e. the data currently
    # loaded, and is what output$locus depends on. `view` is the window
    # actually on screen, which drifts away from `coords` whenever the user
    # zooms or pans inside the loaded data. Keeping them apart is what stops a
    # scroll gesture from triggering a re-render that redraws the same points
    # in place; see the relayout observer below.
    coords <- reactiveValues(chr = NULL, xrange = NULL)
    view <- reactiveValues(chr = NULL, xrange = NULL)

    # Jump somewhere new: the plot has to move, so the data window and the
    # displayed window both change. Every navigation control goes through this.
    goto <- function(chr, xr) {
      xr <- as.integer(xr)
      coords$chr <- chr
      coords$xrange <- xr
      view$chr <- chr
      view$xrange <- xr
    }

    # hide picker at start
    output$coords_ok <- reactive({!is.null(coords$chr)})
    outputOptions(output, "coords_ok", suspendWhenHidden = FALSE)
    
    observe({
      s <- event_data("plotly_click", source = "plotly_manh")
      req(s)
      w <- which(data[, labs] == s$key)
      if (length(w) > 0) {
        goto(data[w[1], chrom], data[w[1], pos] + c(-5e5, 5e5))
      }
    })

    # Clicking the second trait's strip navigates the same way. The lookup
    # goes to data2, not data: a SNP can be present in one trait and absent
    # from the other, and resolving trait 2's key against trait 1 would
    # silently do nothing for exactly those SNPs.
    observe({
      s <- event_data("plotly_click", source = "plotly_manh2")
      req(s, !is.null(data2))
      w <- which(data2[, labs] == s$key)
      if (length(w) > 0) {
        goto(data2[w[1], chrom], data2[w[1], pos] + c(-5e5, 5e5))
      }
    })

    observe({
      s <- event_data("plotly_click", source = "plotly_chrom")
      req(s)
      w <- which(data[, labs] == s$key)
      if (length(w) > 0) {
        goto(data[w[1], chrom], data[w[1], pos] + c(-5e5, 5e5))
      }
    })
    
    # zoom manhattan y axis
    #
    # One pair of buttons drives both strips, but each keeps its OWN limit
    # and its own full range: two GWAS routinely differ by an order of
    # magnitude in power, so forcing a shared scale would flatten the weaker
    # trait to nothing. What is shared is the gesture, not the axis.
    m_ylim <- reactiveValues(max = yrange[2],
                             max2 = if (is.null(data2)) NULL else yrange2[2])

    # Push a y range to one strip. Factored out because there are now four
    # combinations of (zoom in, zoom out) x (trait 1, trait 2) and they
    # differ only in which limit and which output they touch.
    m_relayout <- function(id, lo, hi, title) {
      yr <- c(lo, hi)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy(id, session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = title,
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    }

    m_title <- function(i) {
      if (is.null(data2)) "-log<sub>10</sub> P" else
        paste0(trait_lab[i], "  -log<sub>10</sub> P")
    }

    observeEvent(input$m_zoomin, {
      m_ylim$max <- pmax(m_ylim$max * 0.88, 5)
      m_relayout("manhattan", yrange[1], m_ylim$max, m_title(1))
      if (!is.null(data2)) {
        m_ylim$max2 <- pmax(m_ylim$max2 * 0.88, 5)
        m_relayout("manhattan2", yrange2[1], m_ylim$max2, m_title(2))
      }
    })

    observeEvent(input$m_zoomout, {
      m_ylim$max <- pmin(m_ylim$max / 0.88, yrange[2])
      m_relayout("manhattan", yrange[1], m_ylim$max, m_title(1))
      if (!is.null(data2)) {
        m_ylim$max2 <- pmin(m_ylim$max2 / 0.88, yrange2[2])
        m_relayout("manhattan2", yrange2[1], m_ylim$max2, m_title(2))
      }
    })
    
    # zoom chrom y axis
    chr_y <- reactiveValues(max = 0, range = c(0, 0))
    
    observeEvent(input$chr_zoomin, {
      chr_y$max <- pmax(chr_y$max * 0.88, 5)
      yr <- c(chr_y$range[1], chr_y$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    observeEvent(input$chr_zoomout, {
      chr_y$max <- pmin(chr_y$max / 0.88, chr_y$range[2])
      yr <- c(chr_y$range[1], chr_y$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    input_biotype <- reactive({input$biotype}) %>% debounce(2000)
    
    genes <- reactiveValues(x = NULL)

    # The pinned LD reference variant, NULL when LD is off. Pinning rather than
    # following loc1$index_snp matters twice over: locus() recomputes the index
    # SNP for every window, so an unpinned reference would silently re-base the
    # colouring on a pan, and holding it fixed keeps link_LD()'s arguments
    # identical, so mem_LDproxy serves later windows from cache instead of
    # hitting the API again.
    ld_snp <- reactiveVal(NULL)
    # Published out of the render so the "Get LD" button knows what to pin.
    cur_index <- reactiveVal(NULL)

    output$locus <- renderPlotly({
      req(coords$chr %in% chr_set, coords$xrange)
      # Pre-windowed by binary search, so locus()'s own two filters run over a
      # few thousand rows instead of the whole dataset. Its strict > / <
      # bounds still apply, and locus_rows() is inclusive, so the result is
      # identical to passing the full frame.
      loc1 <- locus(data = locus_subset(data_idx, pos, coords$chr,
                                        coords$xrange),
                     xrange = coords$xrange,
                     seqname = coords$chr, ens_db = ens_db,
                     chrom = chrom, pos = pos, p = p, labs = labs)
      validate(need(loc1$data, "Locus contains no SNPs/datapoints"))
      validate(need(nrow(loc1$data) < 1.5e5, "Too many datapoints. Zoom in."))
      if (isTRUE(input$recomb)) {
        # recomb = NULL sends link_recomb() to the UCSC REST API for this
        # window. It returns the locus with $recomb left NULL on failure
        # rather than aborting, so the plot simply loses the track - which
        # would look like the checkbox does nothing, hence the notice.
        loc1 <- link_recomb(loc1, recomb = recomb)
        if (is.null(loc1$recomb)) {
          showNotification("No recombination data for this region",
                           type = "warning", duration = 5)
        }
      }
      isolate(cur_index(loc1$index_snp))

      # LD, when a reference variant has been pinned. link_LD() keys off
      # loc$index_snp, so override it rather than letting this window's own
      # index SNP take over. scatter_plotly() picks the colouring up on its
      # own once loc1$data has an `ld` column.
      pin <- ld_snp()
      if (!is.null(pin)) {
        loc1$index_snp <- pin
        # Capture rather than merely suppress link_LD's messages: when the API
        # declines, its reason ("Variant is not in 1000G reference panel") is
        # the actionable part, and would otherwise reach only the R console,
        # which nobody driving a browser is watching. The handler has to wrap
        # the call directly - a suppressMessages() inside would muffle each
        # message before this outer handler ever saw it.
        ld_msg <- NULL
        loc1 <- withCallingHandlers(
          link_LD(loc1, token = ld_token, pop = ld_pop),
          message = function(m) {
            txt <- conditionMessage(m)
            if (grepl("^LDproxy:", txt)) {
              ld_msg <<- trimws(sub("^LDproxy:", "", txt))
            }
            invokeRestart("muffleMessage")
          })
        # The blocking API call is done by here, so drop the "fetching"
        # notice whether it succeeded or not.
        removeNotification("ld_busy")
        if (!"ld" %in% colnames(loc1$data)) {
          # link_LD returns the locus untouched when the lookup fails, which
          # would otherwise look like the button did nothing at all.
          showNotification(
            paste0("LD failed for ", pin,
                   if (is.null(ld_msg)) "" else paste0(" - ", ld_msg)),
            type = "error", duration = 10)
        }
      }

      # Second trait, when supplied. Built from the same window, so it needs
      # no separate navigation state.
      loc2 <- NULL
      if (!is.null(data2)) {
        loc2 <- try(locus(data = locus_subset(data2_idx, pos, coords$chr,
                                              coords$xrange),
                          xrange = coords$xrange,
                          seqname = coords$chr, ens_db = ens_db,
                          chrom = chrom, pos = pos, p = p, labs = labs),
                    silent = TRUE)
        if (inherits(loc2, "try-error") || is.null(loc2$data)) {
          loc2 <- NULL
          # Falling back to the single-trait layout silently would make
          # "trait 2 has no data here" indistinguishable from "trait 2 has
          # no signal here" - the exact discrimination the second panel
          # exists to support. The single-panel fallback stays (an empty
          # plotly panel reads worse than none), but say why it happened.
          # try(silent = TRUE) also routes genuine locus() failures here, so
          # the wording covers both.
          showNotification(
            paste0(trait_lab[2], " has no datapoints in this window; ",
                   "showing ", trait_lab[1], " only"),
            type = "warning", duration = 5)
        } else {
          validate(need(nrow(loc2$data) < 1.5e5,
                        paste0("Too many datapoints in ", trait_lab[2],
                               ". Zoom in.")))
          # Recombination on the second panel too. The rate is a property of
          # the locus rather than of either trait, so this is the same line
          # drawn twice - but reading a peak against the rate is much easier
          # when the line sits in the panel you are looking at than when it
          # is one panel away. Costs nothing: link_recomb() memoises on
          # (genome, xrange, seqname, table), and both loci share all four,
          # so trait 2 is served from the cache trait 1 just populated.
          if (isTRUE(input$recomb)) {
            loc2 <- link_recomb(loc2, recomb = recomb)
          }
        }
      }
      # One pinned reference colours both panels: link_LD() already ran
      # against trait 1 above, so this is a match() rather than a second
      # API call.
      if (!is.null(loc2) && "ld" %in% colnames(loc1$data)) {
        ld_ref <- data.frame(snp = loc1$data[, labs],
                             ld = loc1$data$ld,
                             stringsAsFactors = FALSE)
        ld_ref <- ld_ref[!is.na(ld_ref$ld), ]
        loc2 <- join_ld(loc2, ld_ref, labs)
        # Panel 2 gets the same reference variant marked, not its own
        # lowest-p SNP: scatter_plotly() draws index_snp in a distinct
        # "index" style, and marking a variant that is not the LD reference
        # would contradict the colouring the panel is showing. Only on the
        # LD path - without LD, trait 2's own index SNP is the right mark.
        # `pin` can still be NULL here if the caller's own data happened to
        # carry an `ld` column, hence the guard.
        if (!is.null(pin)) loc2$index_snp <- pin
      }

      loc1$TX$fullname <- expandGenes(loc1$TX, fullnames)

      # req(nrow(loc1$data) > 0)
      # This block used to also count traces, to tell plotlyProxy() which ones
      # held the gene track. That proxy path is gone (see below), so only the
      # eQTL values consumed further down survive.
      if (!is.null(eqtl_gene) | !is.null(eqtl_beta)) {
        ind <- loc1$data[, p] < pcutoff
        eqtls <- loc1$data[ind, eqtl_gene]
        genes$x <- unique(eqtls)
      }

      if (!is.null(eqtl_gene)) {
        genes1 <- unique(eqtls)
        locscheme <- unname(c('grey', eqtl_colour[genes1]))
        if (!is.null(input$select_gene) && input$select_gene != "all") {
          # filter gene
          req(input$select_gene %in% unique(eqtls))  # stops double plot
          ok <- !ind | loc1$data[, eqtl_gene] == input$select_gene
          loc1$data <- loc1$data[ok, ]
          locscheme <- unname(c('grey', eqtl_colour[input$select_gene]))
        }
      } else locscheme <- c('grey', 'dodgerblue', 'red')
      
      isolate(width <- loc_width())
      # NOT isolated: the gene track is re-packed client-side by
      # inst/js/genetrack-relayout.js against a payload captured when the
      # widget is built, so changing the biotype filter has to rebuild the
      # widget to give the browser a fresh gene set. Pushing new genes in via
      # plotlyProxy() instead would leave that payload stale, and the next
      # zoom/pan would silently re-pack the pre-filter genes back in.
      biotype <- input_biotype()
      # maxrows: locus_plotly() defaults to 8, which drops a lot of genes on a
      # dense locus - a 1 Mb window round IRF5 needs 14 rows at this width, and
      # 21 in a narrow viewport. 12 rows needs roughly 19px each to stay clear
      # of the 9.8px labels, and both paths are sized to give it: single trait
      # is 0.4 * 600px = 240px, two traits 0.3 * 850px = 255px (the 850 comes
      # from the conditional height on plotlyOutput("locus") above, which
      # exists for exactly this reason - 0.3 * 600 would be ~15px a row).
      # Much above 12 and the labels start colliding with the row above.
      if (is.null(loc2)) {
        locus_plotly(loc1, filter_gene_biotype = biotype, pcutoff = pcutoff,
                     width = width, eqtl_gene = eqtl_gene, beta = eqtl_beta,
                     add_hover = add_hover, scheme = locscheme, maxrows = 12,
                     dynamic = TRUE, scrollZoom = TRUE)
      } else {
        # locscheme is sized for eQTL colouring (1 grey + one colour per
        # eQTL gene significant in this window, a per-window count unrelated
        # to 3). eqtl_gene is deliberately not forwarded to
        # compose_locus_plotly(), so its panels always use scatter_plotly()'s
        # default branch, which requires exactly the 3-tuple below -
        # anything else makes its factor(levels = scheme) call error.
        compose_locus_plotly(list(loc1, loc2), ylabs = trait_lab,
                             filter_gene_biotype = biotype, pcutoff = pcutoff,
                             width = width, maxrows = 12,
                             add_hover = add_hover,
                             scheme = if (is.null(eqtl_gene)) locscheme
                                      else c('grey', 'dodgerblue', 'red'),
                             dynamic = TRUE, scrollZoom = TRUE)
      }
    })
    
    output$ui_genes <- renderUI({
      # req(length(genes$x) > 1)
      g <- c("all", genes$x)
      isolate(ig <- input$select_gene)
      if (length(ig) == 0 || !ig %in% genes$x) ig <- "all"
      conditionalPanel("output.coords_ok",
                       radioButtons("select_gene", h5("eQTL genes"), 
                                    choices = g, selected = ig)
      )
    })
    
    outputOptions(output, "ui_genes", suspendWhenHidden = FALSE)
    
    # Nav buttons step relative to the window on screen, not the loaded one:
    # after a scroll-zoom those differ, and panning should move by what the
    # user can see.
    observeEvent(input$left2, {
      req(view$xrange)
      dif <- diff(view$xrange)
      goto(view$chr, pmax(view$xrange - dif, 0))
    })

    observeEvent(input$right2, {
      req(view$xrange)
      dif <- diff(view$xrange)
      goto(view$chr, view$xrange + dif)
    })

    observeEvent(input$left, {
      req(view$xrange)
      dif <- round(diff(view$xrange) / 2)
      goto(view$chr, pmax(view$xrange - dif, 0))
    })

    observeEvent(input$right, {
      req(view$xrange)
      dif <- round(diff(view$xrange) / 2)
      goto(view$chr, view$xrange + dif)
    })

    observeEvent(input$zoomin, {
      req(view$xrange)
      dif <- round(diff(view$xrange) / 4)
      goto(view$chr, view$xrange + c(dif, -dif))
    })

    observeEvent(input$zoomout, {
      req(view$xrange)
      dif <- round(diff(view$xrange) / 2)
      goto(view$chr, pmax(view$xrange + c(-dif, dif), 0))
    })
    
    output$pos <- renderText({
      req(view$chr %in% chr_set, view$xrange)
      paste0("chr ", view$chr, ": ", view$xrange[1], " - ",
             view$xrange[2])
    })

    # Pin the current window's index SNP and let output$locus do the fetch.
    # The first call blocks for several seconds on the LDlink API, so say so:
    # the notification is put up here, before the render is invalidated, and
    # torn down by the observer below once the new plot has been sent.
    observeEvent(input$ld_get, {
      snp <- cur_index()
      if (is.null(snp) || is.na(snp)) {
        showNotification("No index SNP in view", type = "warning")
        return()
      }
      showNotification(paste0("Fetching LD for ", snp, " (", ld_pop, ")"),
                       id = "ld_busy", duration = NULL)
      ld_snp(snp)
    })

    observeEvent(input$ld_clear, {
      ld_snp(NULL)
      removeNotification("ld_busy")
    })

    # Re-base LD onto a clicked point.
    #
    # Only while LD is already armed: unarmed, a stray click anywhere on the
    # plot would silently cost a multi-second API call. Arming happens even
    # when "Get LD" itself failed, which is what makes this usable at a locus
    # whose index SNP is absent from 1000G - press Get LD, get the rejection,
    # then click a common variant nearby.
    #
    # scatter_plotly() sets key = loc$labs on its point traces, so the clicked
    # SNP arrives directly. The gene track and recombination traces share
    # source = "plotly_locus" but set no key, so requiring one is what keeps a
    # click on a gene line from being taken as a reference variant.
    #
    # ld_snp() is read through isolate() so this observer depends only on the
    # click. Reading it reactively would re-enter on every re-base, see the
    # same stale event_data(), and only be stopped by the identical() guard
    # below - workable, but relying on a value comparison to terminate a loop
    # that never needs to start.
    observe({
      s <- event_data("plotly_click", source = "plotly_locus")
      req(s, !is.null(s$key))
      cur <- isolate(ld_snp())
      req(!is.null(cur))
      snp <- as.character(s$key)[1]
      req(!is.na(snp), nzchar(snp))
      if (identical(snp, cur)) return()
      showNotification(paste0("Fetching LD for ", snp, " (", ld_pop, ")"),
                       id = "ld_busy", duration = NULL)
      ld_snp(snp)
    })

    output$ld_status <- renderText({
      snp <- ld_snp()
      if (is.null(snp)) return("")
      # The hint earns its place: nothing else signals that the scatter is
      # clickable, and re-basing is the only way past an index SNP that the
      # reference panel does not contain.
      paste0("LD: ", snp, " (", ld_pop, ") - click a point to re-base")
    })
    
    # parse text box
    observeEvent(c(input$text_go, input$enter), {
      req(input$tex)
      chr <- NULL
      # Two cleaned forms, because they need different treatment. `tex` also
      # strips "chr" so "chr7:1-2" parses, but that would maim a gene whose
      # symbol contains it - CHRNA5 would become NA5 - so symbol and rsID
      # lookups use `q`, which only trims surrounding whitespace. Matching
      # those two branches against the raw input$tex, as they used to, meant a
      # stray leading space made " IRF5" fall through every branch and return
      # silently.
      q <- trimws(input$tex)
      tex <- gsub(" |chr", "", input$tex, ignore.case = TRUE)
      if (grepl(":", tex) && grepl("-", tex)) {
        # chr & range
        ss <- strsplit(tex, ":")[[1]]
        chr <- ss[1]
        xr <- as.integer(strsplit(ss[2], "-")[[1]])
      } else if (grepl(":", tex)) {
        # single position
        ss <- strsplit(tex, ":")[[1]]
        chr <- ss[1]
        xr <- as.integer(ss[2]) + c(-5e5, 5e5)
      } else if (any(w <- which(toupper(gene_set) == toupper(q)))) {
        gene <- gene_set[w]
        # compare the raw input, so a box holding " irf5" is still normalised
        if (!identical(input$tex, gene)) {
          updateTextInput(session, "tex", value = gene)
        }
        loc <- genes(edb, filter = AnnotationFilterList(
          GeneNameFilter(gene),
          SeqNameFilter(c(1:22, 'X', 'Y'))))
        if (length(loc) > 1) loc <- loc[1]
        chr <- names(seqlengths(loc))
        m <- mean(c(start(loc), end(loc)))
        xr <- as.integer(c(m - 5e5, m + 5e5))
      } else if (grepl("rs", q)) {
        w <- which(data[, labs] == q)
        if (length(w) > 0) {
          chr <- data[w[1], chrom]
          xr <- data[w[1], pos] + c(-5e5, 5e5)
        } else {
          showFeedback("tex", "not found")
          return()
        }
      } else return()
      xr <- as.integer(pmax(xr, 0))
      
      if (chr %in% chr_set) {
        goto(chr, xr)
        hideFeedback("tex")
      } else {
        showFeedback("tex", "not present")
      }
    })
    
    # Table tab
    #
    # Scoped to the window on screen rather than the whole dataset. On a
    # genome-wide GWAS `data` runs to tens of millions of rows: DT paginates
    # server-side so it would not ship all of that to the browser, but it
    # still sorts and filters the full frame on every interaction, and a table
    # of the locus being looked at is the more useful object anyway. Served
    # from the index, so opening the tab costs a binary search rather than a
    # scan of the whole dataset.
    #
    # `view`, not `coords`: the table should describe what is on screen, which
    # after a zoom-in is narrower than the window that was loaded.
    output$table <- DT::renderDataTable({
      validate(need(!is.null(view$chr) && !is.null(view$xrange),
                    "Select a locus to see its datapoints."))
      d <- locus_subset(data_idx, pos, view$chr, view$xrange)
      validate(need(nrow(d) > 0, "No datapoints in this window."))
      d <- d[order(d[, p]), ]
      locus_table(d, coord_cols = c(chrom, pos))
    })
    
    # detect change to x axis range
    #
    # Debounced because `scrollZoom = TRUE` turns one wheel gesture into a
    # burst of relayout events, and each one landing here would trigger a full
    # server round-trip: a fresh locus() call, a fresh ensembl query and a
    # complete re-render of the widget. The client-side re-pack in
    # inst/js/genetrack-relayout.js already redraws the gene track on every
    # one of those events, so the interaction stays responsive while the wheel
    # is turning; the server only needs to catch up once the user settles, to
    # pull in SNPs and genes outside the window originally fetched.
    locus_relayout <- reactive({
      event_data("plotly_relayout", source = "plotly_locus")
    }) %>% debounce(500)

    observeEvent(locus_relayout(), {
      req(coords$chr %in% chr_set, coords$xrange)
      s <- locus_relayout()
      req(c("xaxis.range[0]", "xaxis.range[1]") %in% names(s))
      xr <- as.integer(c(s$`xaxis.range[0]`, s$`xaxis.range[1]`) * 1e6)
      view$chr <- coords$chr
      view$xrange <- xr
      # Only go back to the server when the user has moved outside the data
      # locus() already returned. Zooming in, or panning within it, needs
      # nothing: plotly is already showing the right window client-side and the
      # gene track has been re-packed in the browser, so a re-render would
      # fetch the same rows again and redraw the same points in place - which
      # is exactly the flicker at the end of a scroll gesture.
      if (needs_refetch(xr, coords$xrange)) coords$xrange <- xr
    })
    
    loc_width <- reactiveVal(600)
    
    observe({
      loc_width(session$clientData$output_locus_width)
    })
    
    # The gene track used to be re-packed here, server-side, by pushing fresh
    # coordinates into the widget with plotlyProxy() whenever the plot width or
    # the biotype filter changed. That job now belongs entirely to
    # inst/js/genetrack-relayout.js, which re-packs in the browser on every
    # zoom, pan and resize. Two writers to the same traces and shapes would
    # race, and the browser would win with stale data: its gene payload is
    # captured when the widget is built, so anything proxied in afterwards is
    # discarded on the next re-pack. Width is handled client-side (the JS
    # measures the real rendered axis length, which is more accurate than the
    # `width` argument R packs against); biotype rebuilds the widget instead,
    # see output$locus above.

    # chrom highlight - tracks the window on screen, so it keeps up with
    # scroll-zooming even when no re-render happens
    observeEvent(view$xrange, {
      req(input$show_chrom, view$chr)
      plotlyProxy("chrom", session) %>%
        plotlyProxyInvoke("relayout",
                          list(shapes = list(
                            list(type = "rect",
                                 line = list(width = 1, color = "#00CD00"),
                                 x0 = view$xrange[1] / 1e6,
                                 x1 = view$xrange[2] / 1e6, y0 = 0, y1 = 1,
                                 xref = "x", yref = "paper", layer = "below")
                          )))
    })
    
    # output$print <- renderPrint({
    #   loc_width()
    # })
    
    # output$print <- renderPrint({
    #   req(coords$chr %in% chr_set, coords$xrange)
    #   s <- event_data("plotly_relayout", source = "plotly_locus")
    #   if (is.null(s)) return("Relayout events")
    #   str(s)
    # })
    
  }
  
  runApp(list(ui = ui, server = server))
}


manhattan <- function(data,
                      chrom = NULL, pos = NULL, p = NULL, labs = NULL,
                      pcutoff = 5e-08,
                      chromGap = NULL,
                      chromCols = c('royalblue', 'skyblue'),
                      sigCol = 'red',
                      npoints = 1e6) {
  # autodetect headings
  dc <- detect_cols(data, chrom, pos, p, labs)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  labs <- dc$labs
  
  if (!is.na(npoints) & nrow(data) > npoints) {
    index <- order(data[, p])
    if (npoints <= 1e5) {
      data <- data[index[seq_len(npoints)], ]
    } else {
      # thin points near x axis
      nplotly <- 1e5
      s1 <- seq_len(nplotly)
      s2len <- nrow(data) - nplotly
      s2 <- round(seq_len(npoints - nplotly) * s2len / (npoints - nplotly)) + nplotly
      data <- data[index[unique(c(s1, s2))], ]
    }
  }
  
  data$logP <- -log10(data[, p])
  # p-values that underflowed to zero in the source get clamped to the
  # smallest representable double before reaching here, and -log10 of that is
  # 323.3. That is an artefact of the clamp, not a measurement, and a single
  # such SNP sets the y axis for the entire plot: in a 21M row Alzheimer GWAS
  # one clamped point stretched the axis to 323 while the largest real value
  # was 114.5, so real signal occupied the bottom third of the panel.
  #
  # Peg them to the largest real value so the axis fits the data, and give
  # them their own colour level so they still read as off-scale rather than as
  # a genuine result equal to the strongest measured hit. Only `logP`, which
  # exists solely for plotting, is touched - the p-value column keeps whatever
  # the caller clamped it to.
  clamped <- !is.na(data[, p]) & data[, p] <= 5e-324
  real_max <- suppressWarnings(max(data$logP[!clamped], na.rm = TRUE))
  if (any(clamped) && is.finite(real_max)) {
    data$logP[clamped] <- real_max
    message(sum(clamped), " p-value(s) below floating point precision, ",
            "plotted at the largest measured value (", signif(real_max, 4), ")")
  }
  chrom_list <- mixedsort(unique(data[, chrom]), na.last = NA)
  chrom_list <- as.character(chrom_list)
  
  data[, chrom] <- factor(data[, chrom], levels = chrom_list)
  if (length(chrom_list) == 1) {
    data$genome_pos <- data[, pos]  # single chrom
  } else {
    maxpos <- tapply(data[, pos], data[, chrom], max, na.rm = TRUE)
    maxpos <- maxpos[chrom_list]  # reorder
    minpos <- tapply(data[, pos], data[, chrom], min, na.rm = TRUE)
    minpos <- minpos[chrom_list]  # reorder
    # calculate gap
    if (is.null(chromGap)) {
      chromGap <- sum(maxpos - minpos) / length(chrom_list) / 4.15
    }
    chrom_cumsum <- c(0, cumsum(maxpos - minpos + chromGap))
    chrom_cumsum2 <- chrom_cumsum - c(minpos, 0)
    chrom_cumsum <- chrom_cumsum[1:length(maxpos)]
    chrom_cumsum2 <- chrom_cumsum2[1:length(maxpos)]
    data$genome_pos <- data[, pos] + chrom_cumsum2[as.numeric(data[, chrom])]
  }
  data <- data[order(data$genome_pos), ]
  data$col <- ((as.numeric(data[, chrom]) - 1) %% length(chromCols)) + 1
  colScheme <- chromCols
  if (!is.na(sigCol)) {
    data$col[data[, p] < pcutoff] <- length(chromCols) + 1
    colScheme <- c(chromCols, sigCol)
  }
  # Applied after the significance level, which it deliberately overrides: a
  # clamped point is significant by definition, and the useful thing to convey
  # is that its value is a floor rather than a measurement. `data` has been
  # reordered above, so recompute rather than reusing the earlier vector.
  data$col[!is.na(data[, p]) & data[, p] <= 5e-324] <- length(chromCols) + 2L
  if (length(chrom_list) > 1) {
    xticks <- list(at = chrom_cumsum + 0.5 * (maxpos - minpos), 
                   labels = levels(data[, chrom]))
  } else xticks <- NULL
  
  ret <- list(data = data, xticks = xticks, pcutoff = pcutoff,
              chrom_list = chrom_list)
  class(ret) <- "manhattan"
  ret
}


plotly_manhattan <- function(obj,
                             labs,
                             scheme = c('royalblue', 'skyblue', 'red'),
                             clampCol = 'black',
                             xlab = "Chromosome",
                             pcutline = NULL,
                             source = "plotly_manh") {

  df <- obj$data
  df$col <- as.factor(df$col)
  # Appended rather than being part of `scheme` so that callers passing their
  # own three colours, as the chromosome panel does, keep working unchanged
  # and still get the extra level when clamped points are present. Levels are
  # only ever present in the factor when manhattan() assigned them, so the
  # subset below drops it again when there are none.
  scheme <- c(scheme, clampCol)
  scheme <- scheme[as.numeric(levels(df$col))]
  if (is.null(obj$xticks)) {
    # single chrom
    df$genome_pos <- df$genome_pos / 1e6
    if (xlab == "Chromosome") xlab <- paste(xlab, obj$chrom_list, "(Mb)")
  }
  xr <- range(df$genome_pos, na.rm = TRUE)
  xr <- xr + diff(xr) * c(-0.01, 0.01)
  yr <- range(df$logP, na.rm = TRUE)
  yr <- yr + diff(yr) * c(-0.05, 0.05)
  
  hline <- if (!is.null(pcutline)) {
    list(type = "line",
         line = list(width = 1, color = '#AAAAAA', dash = 'dash'),
         x0 = 0, x1 = 1, y0 = -log10(pcutline), y1 = -log10(pcutline),
         xref = "paper", layer = "below")
  } else NULL
  xlayout <- list(range = xr, title = xlab, ticks = "outside",
                  zeroline = FALSE, showline = TRUE, showgrid = FALSE)
  if (!is.null(obj$xticks)) {
    xlayout <- c(xlayout, list(tickvals = obj$xticks$at,
                               ticktext = obj$xticks$labels))
  }
  
  plot_ly(data = df, x = ~genome_pos, y = ~logP,
          color = ~col, colors = scheme,
          marker = list(size = 4, opacity = 0.8),
          text = as.formula(paste0('~', labs)),
          hoverinfo = 'text', key = as.formula(paste0('~', labs)),
          showlegend = FALSE,
          type = "scattergl", mode = "markers",
          source = source) %>%
    plotly::layout(xaxis = xlayout,
                   yaxis = list(range = yr,
                                title = "-log<sub>10</sub> P",
                                ticks = "outside",
                                zeroline = FALSE, showline = TRUE),
                   shapes = hline)
}


#' Does the window on screen need data the loaded window does not hold?
#'
#' `locus()` returns every datapoint inside the window it was called with, so
#' any view contained within that window can be drawn from what the browser
#' already has. Only a view reaching past either edge needs a fresh query.
#'
#' @param view_xr Numeric length-2, the window currently displayed.
#' @param loaded_xr Numeric length-2, the window `locus()` was called with.
#' @return `TRUE` if a re-query is needed.
#' @noRd
needs_refetch <- function(view_xr, loaded_xr) {
  view_xr[1] < loaded_xr[1] || view_xr[2] > loaded_xr[2]
}


unique_snps <- function(data, labs, append) {
  snps <- data[, labs]
  dups <- which(duplicated(snps))
  if (length(dups) > 0) {
    message("Duplicated SNPs found")
    snps[dups] <- make.unique(paste(snps[dups], data[dups, append], sep = "."))
  }
  snps
}


#' @importFrom ensembldb listColumns
fullGeneNames <- function(edb, AnnotationDb) {
  # check ens_db cols for 'description' first
  if (!"description" %in% listColumns(edb) | is.null(AnnotationDb)) return(NULL)
  
  if (!requireNamespace(AnnotationDb)) {
    stop("Gene annotation database '", AnnotationDb, "' is not installed")
  }
  if (is.character(AnnotationDb)) {
    AnnotationDb <- eval(str2lang(paste0(AnnotationDb, "::", AnnotationDb)))
  }
  alias <- AnnotationDbi::keys(AnnotationDb, "ALIAS")
  suppressMessages(
    AnnotationDbi::mapIds(AnnotationDb, alias,
                          "GENENAME", "ALIAS", multiVals = 'first')
  )
}


# Full genename lookup, returns hovertext
expandGenes <- function(TX, fullnames) {
  genelist <- TX$gene_name
  if ("description" %in% colnames(TX)) {
    # check ensembldb first
    out <- gsub(" \\[[^][]*]", "", TX$description)
  } else {
    if (is.null(fullnames)) return(genelist)
    out <- fullnames[genelist]
  }
  bad <- is.na(out) | out == "NULL" | out == ""
  out <- paste0("<br>", out)
  out[bad] <- ""
  out
}

