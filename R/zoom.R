
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
                 AnnotationDb = "org.Hs.eg.db") {
  data <- data.frame(data)
  # autodetect headings
  dc <- detect_cols(data, chrom, pos, p, labs)
  chrom <- dc$chrom
  pos <- dc$pos
  p <- dc$p
  labs <- dc$labs
  if (is.null(eqtl_gene)) {
    data[, labs] <- unique_snps(data, labs, chrom)
  } else {
    data[, labs] <- unique_snps(data, labs, eqtl_gene)
  }
  
  message("Generating Manhattan plot")
  chr_set <- unique(data[, chrom])
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
  manhat <- manhattan(data, chrom, pos, p, labs, pcutoff = pcutoff,
                      npoints = mh_points)
  yrange <- range(manhat$data$logP, na.rm = TRUE)
  ymax <- yrange[2] + diff(yrange) * 0.05
  
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
                        withSpinner(
                          plotlyOutput("manhattan", width = "85vw", height = "300px"),
                          type = 8, size = 0.7)
                 ),
                 column(1,
                        br(),
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
                          (if (!is.null(recomb)) {
                            checkboxInput("recomb", "show recombination rate", value = TRUE)
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
                          plotlyOutput("locus", width = "95vw", height = 600),
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
      plotly_manhattan(manhat, labs, pcutline = NULL) %>%
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
      chr_manhat <- manhattan(data[which(data[, chrom] == coords$chr), ],
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
    
    observe({
      s <- event_data("plotly_click", source = "plotly_chrom")
      req(s)
      w <- which(data[, labs] == s$key)
      if (length(w) > 0) {
        goto(data[w[1], chrom], data[w[1], pos] + c(-5e5, 5e5))
      }
    })
    
    # zoom manhattan y axis
    m_ylim <- reactiveValues(max = yrange[2])
    
    observeEvent(input$m_zoomin, {
      m_ylim$max <- pmax(m_ylim$max * 0.88, 5)
      yr <- c(yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
    })
    
    observeEvent(input$m_zoomout, {
      m_ylim$max <- pmin(m_ylim$max / 0.88, yrange[2])
      yr <- c(yrange[1], m_ylim$max)
      yr <- yr + diff(yr) * c(-0.05, 0.05)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout",
                          list(yaxis = list(range = yr,
                                            title = "-log<sub>10</sub> P",
                                            ticks = "outside",
                                            zeroline = FALSE, showline = TRUE)))
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
    
    output$locus <- renderPlotly({
      req(coords$chr %in% chr_set, coords$xrange)
      loc1 <- locus(data = data, xrange = coords$xrange,
                     seqname = coords$chr, ens_db = ens_db,
                     chrom = chrom, pos = pos, p = p, labs = labs)
      validate(need(loc1$data, "Locus contains no SNPs/datapoints"))
      validate(need(nrow(loc1$data) < 1.5e5, "Too many datapoints. Zoom in."))
      if (!is.null(recomb) && input$recomb) {
        loc1 <- link_recomb(loc1, recomb = recomb)
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
      # 21 in a narrow viewport. The gene panel here is 0.4 * 600px = 240px, so
      # 12 rows leaves ~19px each, comfortable for the 9.8px labels; much above
      # that and the labels start colliding with the row above.
      locus_plotly(loc1, filter_gene_biotype = biotype, pcutoff = pcutoff,
                   width = width, eqtl_gene = eqtl_gene, beta = eqtl_beta,
                   add_hover = add_hover, scheme = locscheme, maxrows = 12,
                   dynamic = TRUE, scrollZoom = TRUE)
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
    output$table <- DT::renderDataTable({
      cols <- colnames(data)[sapply(data, class) == "numeric"]
      datatable(data) %>% formatSignif(cols, digits = 3)
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
                             xlab = "Chromosome",
                             pcutline = NULL,
                             source = "plotly_manh") {
  
  df <- obj$data
  df$col <- as.factor(df$col)
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

