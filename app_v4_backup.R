library(shiny); library(bslib); library(DT)
library(ggplot2); library(dplyr); library(readr); library(tidyr)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
DATA_DIR <- if (dir.exists("data")) "data" else if (dir.exists(file.path("app","data"))) file.path("app","data") else "."
p  <- function(f) file.path(DATA_DIR, f)
rt <- function(f) suppressMessages(readr::read_tsv(f, show_col_types = FALSE))

datasets_all <- rt(p("datasets_master.tsv"))
datasets <- datasets_all %>% filter(assay %in% c("RNA-seq","ATAC-seq"))
proteome <- if (file.exists(p("proteome_master.tsv"))) rt(p("proteome_master.tsv")) else NULL
leads    <- if (file.exists(p("leads_to_follow.tsv"))) rt(p("leads_to_follow.tsv")) else NULL
tsne_df     <- rt(p("tsne_coords.tsv"))
de_df       <- rt(p("de_all_datasets.tsv"))
paths_df    <- rt(p("pathways.tsv"))
det_df      <- rt(p("det_all.tsv"))
tsne_det_df <- rt(p("tsne_det_coords.tsv"))
gene_model  <- rt(p("gene_model.tsv"))
gene_mut    <- rt(p("gene_mutations.tsv"))
gene_pal <- c(SCN1A = "#207394", SCN2A = "#CD5241")

keep_val <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && x != "All"

make_url <- function(pmid_doi, accession) {
  x <- as.character(pmid_doi)
  if (is.na(x) || x == "" || grepl("^PMID \\(", x)) {
    acc <- as.character(accession)
    if (grepl("^GSE", acc)) return(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", acc))
    sr <- sub("[;,-].*$", "", acc)
    if (grepl("^SRR|^SRX", sr)) return(paste0("https://www.ncbi.nlm.nih.gov/sra/?term=", sr))
    return(NA_character_)
  }
  if (grepl("^10\\.", x))   return(paste0("https://doi.org/", x))
  if (grepl("^PMID:", x))  return(paste0("https://pubmed.ncbi.nlm.nih.gov/", sub("PMID:","",x)))
  if (grepl("^PMC", x))    return(paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/", sub(" ","",x)))
  if (grepl("^eLife:", x)) return(paste0("https://elifesciences.org/articles/", sub("eLife:","",x)))
  NA_character_
}
datasets$url <- mapply(make_url, datasets$pmid_doi, datasets$accession)
datasets_all$url <- mapply(make_url, datasets_all$pmid_doi, datasets_all$accession)

org_ch  <- c("All", sort(unique(datasets$organism)))
chan_ch <- c("All","SCN1A","SCN2A")

ui <- page_navbar(
  title = "SCN1A / SCN2A Atlas",
  theme = bs_theme(version = 5, primary = "#207394"),

  nav_panel("Datasets", icon = icon("table"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("f_org", "Organism", org_ch),
        selectInput("f_chan", "Channel (SCN1A/SCN2A)", chan_ch),
        selectInput("f_sub", "Subtype", c("All", sort(unique(datasets$subtype)))),
        hr(), helpText("Organism then channel. Link column opens the paper.")),
      card(card_header("Master dataset table"), DTOutput("dt_datasets")),
      card(card_header("Row detail"), uiOutput("row_detail")),
      card(card_header("Gene map \u2014 mutation location per dataset"),
           helpText("Click a marker to open that dataset above."),
           plotOutput("gene_map", click = "map_click", height = "420px"),
           verbatimTextOutput("map_info")),
      )),

  nav_panel("t-SNE", icon = icon("braille"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        radioButtons("tsne_color", "Color by",
          c("Gene"="gene","Organism"="organism","Dataset"="dataset","Mutation"="mutation")),
        checkboxInput("tsne_nowt", "Hide controls (WT)", FALSE),
        hr(), helpText("Gene-level. Click a point to inspect.")),
      card(card_header("Sample embedding (gene-level)"),
           plotOutput("tsne_plot", click = "tsne_click", height = "560px")),
      card(card_header("Selected sample"), verbatimTextOutput("tsne_info")))),

  nav_panel("Gene DE", icon = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("de_org", "Organism", org_ch),
        selectInput("de_chan", "Channel (SCN1A/SCN2A)", chan_ch),
        selectizeInput("de_gene", "Search gene", choices = NULL, options = list(placeholder = "e.g. SCN2A")),
        sliderInput("de_padj", "padj", 0, 1, 0.05, step = 0.01),
        hr(), helpText("Organism -> channel -> gene.")),
      card(card_header(textOutput("de_title")), plotOutput("de_plot", height = "540px")),
      card(card_header("Values"), DTOutput("de_table")))),

  nav_panel("DET", icon = icon("wave-square"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("det_org", "Organism", org_ch),
        selectInput("det_chan", "Channel (SCN1A/SCN2A)", chan_ch),
        selectizeInput("det_gene", "Search gene", choices = NULL, options = list(placeholder = "e.g. SCN2A")),
        sliderInput("det_q", "qvalue", 0, 1, 0.05, step = 0.01),
        checkboxInput("tsne_det_nowt", "t-SNE: hide controls (WT)", FALSE),
        hr(), helpText("Organism -> channel -> gene. Sub-tab: transcript t-SNE.")),
      navset_tab(
        nav_panel("Bargraph",
          card(card_header(textOutput("det_title")), plotOutput("det_plot", height = "560px")),
          card(card_header("Transcripts"), DTOutput("det_table"))),
        nav_panel("t-SNE (DET)",
          card(card_header("Sample embedding (transcript-level)"),
               plotOutput("tsne_det_plot", click = "tsne_det_click", height = "560px")),
          card(card_header("Selected sample"), verbatimTextOutput("tsne_det_info")))))),

  nav_panel("Pathways", icon = icon("diagram-project"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        sliderInput("path_min", "Min. datasets sharing", 1, 6, 1),
        hr(), helpText("Pathways across datasets = shared biology.")),
      card(card_header("Shared pathways"), plotOutput("path_plot", height = "560px")),
      card(card_header("Pathway x dataset"), DTOutput("path_table"))))
)

server <- function(input, output, session) {

  ds_filt <- reactive({
    d <- datasets
    if (keep_val(input$f_org))    d <- d[d$organism == input$f_org, ]
    if (keep_val(input$f_chan))   d <- d[d$gene == input$f_chan, ]
    if (keep_val(input$f_sub))    d <- d[d$subtype == input$f_sub, ]
    d
  })
  output$dt_datasets <- renderDT({
    d <- ds_filt()
    d$Link <- ifelse(is.na(d$url), "", paste0("<a href='", d$url, "' target='_blank'>open</a>"))
    d %>% select(dataset_id, gene, organism, assay, subtype, n_samples, mechanism, Link) %>%
      datatable(selection = "single", rownames = FALSE, escape = FALSE,
                options = list(pageLength = 18, scrollX = TRUE, scrollY = "420px")) %>%
      formatStyle("gene", target = "row",
        backgroundColor = styleEqual(c("SCN1A","SCN2A"), c("#eaf3f7","#fdeeea")))
  })
  output$row_detail <- renderUI({
    sel <- input$dt_datasets_rows_selected
    if (is.null(sel)) return(em("Select a row above."))
    r <- ds_filt()[sel, ]
    link <- if (!is.na(r$url)) tags$a(href = r$url, target = "_blank", "Open paper/data") else "no link"
    tagList(h5(strong(r$dataset_id),
      if (r$status == "our_dataset") span(class="badge bg-primary", "our dataset")),
      tags$ul(tags$li(strong("Link: "), link),
        tags$li(strong("Paper: "), r$source_paper, " (", r$pmid_doi, ")"),
        tags$li(strong("Genotypes: "), r$genotypes),
        tags$li(strong("System: "), r$tissue_system),
        tags$li(strong("Accessions: "), tags$code(r$accession)),
        tags$li(strong("Mechanism: "), r$mechanism)))
  })
  output$dt_proteome <- renderDT({
    if (is.null(proteome)) return(datatable(data.frame(note="no proteome_master.tsv")))
    datatable(proteome, rownames = FALSE, options = list(pageLength = 5, scrollX = TRUE))})
  output$dt_leads <- renderDT({
    if (is.null(leads)) return(datatable(data.frame(note="no leads_to_follow.tsv")))
    datatable(leads, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))})

  # ---- gene map ----
  output$gene_map <- renderPlot({
    genes <- c("SCN1A","SCN2A"); yb <- c(SCN1A = 2, SCN2A = 1)
    ex <- gene_model; mk <- gene_mut
    p <- ggplot() + xlim(-0.05, 1.20) + ylim(0.3, 3.2)
    for (g in genes) {
      y <- yb[[g]]; e <- ex[ex$gene == g, ]
      p <- p + annotate("segment", x = 0, xend = 1, y = y, yend = y, color = "grey70", linewidth = 0.8) +
        annotate("rect", xmin = e$x0, xmax = pmax(e$x1, e$x0 + 0.004), ymin = y - 0.08, ymax = y + 0.08,
                 fill = gene_pal[[g]], color = NA) +
        annotate("text", x = -0.02, y = y + 0.22, label = g, hjust = 0, fontface = "bold", size = 5)
      gw <- mk[mk$gene == g & mk$kind == "genewide", ]
      if (nrow(gw) > 0) p <- p + annotate("segment", x = 0, xend = 1, y = y - 0.20, yend = y - 0.20,
                 color = gene_pal[[g]], linewidth = 2, alpha = 0.25)
      pts <- mk[mk$gene == g & mk$kind %in% c("point","promoter"), ]
      if (nrow(pts) > 0) {
        pts$xx <- ifelse(is.na(pts$frac), 0, pts$frac)
        p <- p + annotate("point", x = pts$xx, y = y + 0.08, shape = 25, size = 3.2,
                 fill = "black", color = "black") +
          annotate("text", x = pts$xx, y = y + 0.38, label = pts$label, angle = 25,
                 hjust = 0, size = 3.2)
      }
    }
    p + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(40,60,20,20))
  })
  output$map_info <- renderText({
    cl <- input$map_click
    if (is.null(cl)) return("Click a mutation marker (triangle) to identify the dataset.")
    yb <- c(SCN1A = 2, SCN2A = 1)
    g <- if (cl$y > 1.5) "SCN1A" else "SCN2A"
    pts <- gene_mut[gene_mut$gene == g & gene_mut$kind %in% c("point","promoter"), ]
    if (nrow(pts) == 0) return("No point markers on this gene.")
    pts$xx <- ifelse(is.na(pts$frac), 0, pts$frac)
    i <- which.min(abs(pts$xx - cl$x))
    if (abs(pts$xx[i] - cl$x) > 0.06) return("No marker near the click.")
    paste0("Dataset: ", pts$dataset[i], "\nMutation: ", pts$label[i],
           "\nGene: ", g, "\nNote: ", pts$note[i])
  })

  tsne_view <- reactive({ df <- tsne_df; if (isTRUE(input$tsne_nowt)) df <- df %>% filter(mutation != "WT"); df })
  output$tsne_plot <- renderPlot({
    cc <- input$tsne_color %||% "gene"
    ggplot(tsne_view(), aes(TSNE1, TSNE2, color = .data[[cc]])) +
      geom_point(size = 3, alpha = 0.85) +
      { if (cc == "gene") scale_color_manual(values = gene_pal) else NULL } +
      labs(color = tools::toTitleCase(cc), x = "t-SNE 1", y = "t-SNE 2") +
      theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank())
  })
  output$tsne_info <- renderText({
    cl <- input$tsne_click; if (is.null(cl)) return("Click a point.")
    np <- nearPoints(tsne_view(), cl, xvar="TSNE1", yvar="TSNE2", threshold=20, maxpoints=1)
    if (nrow(np)==0) return("No point nearby.")
    paste0("Sample:   ", np$sample, "\nDataset:  ", np$dataset, "\nGene:     ", np$gene,
           "\nOrganism: ", np$organism, "\nMutation: ", np$mutation)
  })

  de_ids <- reactive({
    ids <- datasets_all$dataset_id
    if (keep_val(input$de_org))  ids <- ids[datasets_all$organism == input$de_org]
    if (keep_val(input$de_chan)) ids <- intersect(ids, datasets_all$dataset_id[datasets_all$gene == input$de_chan])
    ids })
  de_genes <- reactive({ sort(unique(de_df$gene[de_df$dataset %in% de_ids()])) })
  observeEvent(de_genes(), { g <- de_genes()
    updateSelectizeInput(session, "de_gene", choices = g,
      selected = if ("SCN2A" %in% g) "SCN2A" else g[1], server = TRUE) })
  de_sub <- reactive({ req(input$de_gene)
    de_df %>% filter(gene == input$de_gene, dataset %in% de_ids()) %>%
      mutate(sig = padj < input$de_padj,
             label = if ("contrast" %in% names(de_df)) paste(dataset, contrast, sep=" - ") else dataset) })
  output$de_title <- renderText(sprintf("%s - log2FC across datasets", input$de_gene %||% ""))
  output$de_plot <- renderPlot({ d <- de_sub(); validate(need(nrow(d)>0, "No data for this combination."))
    d$label <- factor(d$label, levels = d$label[order(d$log2FoldChange)])
    ggplot(d, aes(log2FoldChange, label, fill = log2FoldChange, alpha = sig)) +
      geom_col() + geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_gradient2(low = "#207394", mid = "grey85", high = "#CD5241") +
      scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.35), guide = "none") +
      labs(x = "log2 fold change", y = NULL, caption = "Solid = passes padj") +
      theme_minimal(base_size = 13) })
  output$de_table <- renderDT({ de_sub() %>% arrange(padj) %>%
      select(any_of(c("dataset","contrast","log2FoldChange","padj"))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10)) })

  det_ids <- reactive({
    ids <- datasets_all$dataset_id
    if (keep_val(input$det_org))  ids <- ids[datasets_all$organism == input$det_org]
    if (keep_val(input$det_chan)) ids <- intersect(ids, datasets_all$dataset_id[datasets_all$gene == input$det_chan])
    ids })
  det_genes <- reactive({ sort(unique(det_df$gene_name[det_df$dataset %in% det_ids()])) })
  observeEvent(det_genes(), { g <- det_genes()
    updateSelectizeInput(session, "det_gene", choices = g,
      selected = if ("SCN2A" %in% g) "SCN2A" else g[1], server = TRUE) })
  det_sub <- reactive({ req(input$det_gene)
    det_df %>% filter(gene_name == input$det_gene, dataset %in% det_ids()) %>%
      mutate(sig = qvalue < input$det_q, label = paste(transcript_name, dataset, contrast, sep=" - ")) })
  output$det_title <- renderText(sprintf("%s - transcripts (DET)", input$det_gene %||% ""))
  output$det_plot <- renderPlot({ d <- det_sub(); validate(need(nrow(d)>0, "No transcripts for this combination."))
    d <- d %>% arrange(desc(abs(log2FC))) %>% head(40)
    d$label <- factor(d$label, levels = d$label[order(d$log2FC)])
    ggplot(d, aes(log2FC, label, fill = log2FC, alpha = sig)) +
      geom_col() + geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_gradient2(low = "#207394", mid = "grey85", high = "#CD5241") +
      scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.35), guide = "none") +
      labs(x = "log2FC (transcript)", y = NULL, caption = "top 40 by |log2FC|; solid passes qvalue") +
      theme_minimal(base_size = 11) })
  output$det_table <- renderDT({ det_sub() %>% arrange(qvalue) %>%
      select(any_of(c("transcript_name","transcript_id","dataset","contrast","log2FC","qvalue"))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })

  tsne_det_view <- reactive({ df <- tsne_det_df
    if (keep_val(input$det_org)) df <- df %>% filter(organism == input$det_org)
    if (isTRUE(input$tsne_det_nowt)) df <- df %>% filter(mutation != "WT"); df })
  output$tsne_det_plot <- renderPlot({ df <- tsne_det_view(); validate(need(nrow(df)>0, "No samples."))
    ggplot(df, aes(TSNE1, TSNE2, color = mutation)) + geom_point(size = 3, alpha = 0.85) +
      labs(color = "Mutation", x = "t-SNE 1", y = "t-SNE 2", subtitle = "transcript-level") +
      theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank()) })
  output$tsne_det_info <- renderText({ cl <- input$tsne_det_click; if (is.null(cl)) return("Click a point.")
    np <- nearPoints(tsne_det_view(), cl, xvar="TSNE1", yvar="TSNE2", threshold=20, maxpoints=1)
    if (nrow(np)==0) return("No point nearby.")
    paste0("Sample:   ", np$sample, "\nDataset:  ", np$dataset, "\nGene:     ", np$gene,
           "\nOrganism: ", np$organism, "\nMutation: ", np$mutation) })

  path_shared <- reactive({ paths_df %>% group_by(pathway) %>%
      summarise(n_datasets = n_distinct(dataset),
                datasets = paste(sort(unique(dataset)), collapse=", "), .groups="drop") %>%
      filter(n_datasets >= input$path_min) %>% arrange(desc(n_datasets)) })
  output$path_plot <- renderPlot({ ps <- path_shared(); validate(need(nrow(ps)>0, "No pathways.")); ps <- head(ps, 40)
    ps$pathway <- factor(ps$pathway, levels = ps$pathway[order(ps$n_datasets)])
    ggplot(ps, aes(n_datasets, pathway, fill = n_datasets)) + geom_col() +
      scale_fill_gradient(low = "#a9cfdf", high = "#207394") +
      labs(x = "number of datasets", y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none") })
  output$path_table <- renderDT(datatable(path_shared(), rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)))
}

shinyApp(ui, server)
