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
gene_pal <- c(SCN1A = "#207394", SCN2A = "#CD5241")

# filtro seguro: ignora se input NULL/vazio/All
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

# mapa dataset -> canal (gene do canal estudado) p/ cascata nas buscas
chan_of <- setNames(datasets_all$gene, datasets_all$dataset_id)

org_ch  <- c("All", sort(unique(datasets$organism)))
chan_ch <- c("All","SCN1A","SCN2A")

ui <- page_navbar(
  title = "SCN1A / SCN2A Atlas",
  theme = bs_theme(version = 5, primary = "#207394"),

  nav_panel("Datasets", icon = icon("table"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("f_org", "Organism", org_ch),
        selectInput("f_chan", "Canal (SCN1A/SCN2A)", chan_ch),
        selectInput("f_sub", "Subtype", c("All", sort(unique(datasets$subtype)))),
        selectInput("f_status", "Status", c("All", sort(unique(datasets$status)))),
        hr(), helpText("Organismo -> canal. Coluna Link abre o paper.")),
      card(card_header("Master dataset table"), DTOutput("dt_datasets")),
      card(card_header("Detalhe da linha"), uiOutput("row_detail")),
      navset_tab(nav_panel("Proteome", DTOutput("dt_proteome")),
                 nav_panel("Leads", DTOutput("dt_leads"))))),

  nav_panel("t-SNE", icon = icon("braille"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        radioButtons("tsne_color", "Colorir por",
          c("Gene"="gene","Organism"="organism","Dataset"="dataset","Mutation"="mutation")),
        checkboxInput("tsne_nowt", "Esconder controles (WT)", FALSE),
        hr(), helpText("Gene-level. Clique num ponto p/ inspecionar.")),
      card(card_header("Sample embedding (gene-level)"),
           plotOutput("tsne_plot", click = "tsne_click", height = "560px")),
      card(card_header("Amostra clicada"), verbatimTextOutput("tsne_info")))),

  nav_panel("Gene DE", icon = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("de_org", "Organism", org_ch),
        selectInput("de_chan", "Canal (SCN1A/SCN2A)", chan_ch),
        selectizeInput("de_gene", "Gene buscado", choices = NULL, options = list(placeholder = "ex. SCN2A")),
        sliderInput("de_padj", "padj", 0, 1, 0.05, step = 0.01),
        hr(), helpText("Organismo -> canal -> gene.")),
      card(card_header(textOutput("de_title")), plotOutput("de_plot", height = "540px")),
      card(card_header("Valores"), DTOutput("de_table")))),

  nav_panel("DET", icon = icon("wave-square"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("det_org", "Organism", org_ch),
        selectInput("det_chan", "Canal (SCN1A/SCN2A)", chan_ch),
        selectizeInput("det_gene", "Gene buscado", choices = NULL, options = list(placeholder = "ex. SCN2A")),
        sliderInput("det_q", "qvalue", 0, 1, 0.05, step = 0.01),
        checkboxInput("tsne_det_nowt", "t-SNE: esconder controles (WT)", FALSE),
        hr(), helpText("Organismo -> canal -> gene. Sub-aba: t-SNE por transcrito.")),
      navset_tab(
        nav_panel("Bargraph",
          card(card_header(textOutput("det_title")), plotOutput("det_plot", height = "560px")),
          card(card_header("Transcritos"), DTOutput("det_table"))),
        nav_panel("t-SNE (DET)",
          card(card_header("Sample embedding (transcript-level)"),
               plotOutput("tsne_det_plot", click = "tsne_det_click", height = "560px")),
          card(card_header("Amostra clicada"), verbatimTextOutput("tsne_det_info")))))),

  nav_panel("Pathways", icon = icon("diagram-project"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        sliderInput("path_min", "Min. datasets compartilhando", 1, 6, 1),
        hr(), helpText("Vias em varios datasets = biologia compartilhada.")),
      card(card_header("Vias compartilhadas"), plotOutput("path_plot", height = "560px")),
      card(card_header("Via x dataset"), DTOutput("path_table"))))
)

server <- function(input, output, session) {

  ds_filt <- reactive({
    d <- datasets
    if (keep_val(input$f_org))    d <- d[d$organism == input$f_org, ]
    if (keep_val(input$f_chan))   d <- d[d$gene == input$f_chan, ]
    if (keep_val(input$f_sub))    d <- d[d$subtype == input$f_sub, ]
    if (keep_val(input$f_status)) d <- d[d$status == input$f_status, ]
    d
  })

  output$dt_datasets <- renderDT({
    d <- ds_filt()
    d$Link <- ifelse(is.na(d$url), "",
      paste0("<a href='", d$url, "' target='_blank'>abrir</a>"))
    d %>% select(dataset_id, gene, organism, assay, subtype, n_samples, mechanism, status, Link) %>%
      datatable(selection = "single", rownames = FALSE, escape = FALSE,
                options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatStyle("gene", target = "row",
        backgroundColor = styleEqual(c("SCN1A","SCN2A"), c("#eaf3f7","#fdeeea")))
  })

  output$row_detail <- renderUI({
    sel <- input$dt_datasets_rows_selected
    if (is.null(sel)) return(em("Selecione uma linha acima."))
    r <- ds_filt()[sel, ]
    link <- if (!is.na(r$url)) tags$a(href = r$url, target = "_blank", "Abrir paper/dados") else "sem link"
    tagList(
      h5(strong(r$dataset_id),
         if (r$status == "our_dataset") span(class="badge bg-primary", "our dataset")),
      tags$ul(
        tags$li(strong("Link: "), link),
        tags$li(strong("Paper: "), r$source_paper, " (", r$pmid_doi, ")"),
        tags$li(strong("Genotipos: "), r$genotypes),
        tags$li(strong("Sistema: "), r$tissue_system),
        tags$li(strong("Acessos: "), tags$code(r$accession)),
        tags$li(strong("Mecanismo: "), r$mechanism),
        tags$li(strong("Notas: "), r$notes)))
  })

  output$dt_proteome <- renderDT({
    if (is.null(proteome)) return(datatable(data.frame(nota="sem proteome_master.tsv")))
    datatable(proteome, rownames = FALSE, options = list(pageLength = 5, scrollX = TRUE))
  })
  output$dt_leads <- renderDT({
    if (is.null(leads)) return(datatable(data.frame(nota="sem leads_to_follow.tsv")))
    datatable(leads, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
  })

  # t-SNE gene-level (com filtro WT)
  tsne_view <- reactive({
    df <- tsne_df
    if (isTRUE(input$tsne_nowt)) df <- df %>% filter(mutation != "WT")
    df
  })
  output$tsne_plot <- renderPlot({
    cc <- input$tsne_color %||% "gene"
    ggplot(tsne_view(), aes(TSNE1, TSNE2, color = .data[[cc]])) +
      geom_point(size = 3, alpha = 0.85) +
      { if (cc == "gene") scale_color_manual(values = gene_pal) else NULL } +
      labs(color = tools::toTitleCase(cc), x = "t-SNE 1", y = "t-SNE 2") +
      theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank())
  })
  output$tsne_info <- renderText({
    cl <- input$tsne_click
    if (is.null(cl)) return("Clique num ponto no grafico.")
    np <- nearPoints(tsne_view(), cl, xvar = "TSNE1", yvar = "TSNE2", threshold = 20, maxpoints = 1)
    if (nrow(np) == 0) return("Nenhum ponto perto.")
    paste0("Sample:   ", np$sample, "\nDataset:  ", np$dataset,
           "\nGene:     ", np$gene, "\nOrganism: ", np$organism,
           "\nMutation: ", np$mutation)
  })

  # ---- Gene DE: organismo -> canal -> gene ----
  de_ids <- reactive({
    ids <- datasets_all$dataset_id
    if (keep_val(input$de_org))  ids <- ids[datasets_all$organism == input$de_org]
    if (keep_val(input$de_chan)) ids <- intersect(ids, datasets_all$dataset_id[datasets_all$gene == input$de_chan])
    ids
  })
  de_genes <- reactive({ sort(unique(de_df$gene[de_df$dataset %in% de_ids()])) })
  observeEvent(de_genes(), {
    g <- de_genes()
    updateSelectizeInput(session, "de_gene", choices = g,
      selected = if ("SCN2A" %in% g) "SCN2A" else g[1], server = TRUE)
  })
  de_sub <- reactive({
    req(input$de_gene)
    de_df %>% filter(gene == input$de_gene, dataset %in% de_ids()) %>%
      mutate(sig = padj < input$de_padj,
             label = if ("contrast" %in% names(de_df)) paste(dataset, contrast, sep=" - ") else dataset)
  })
  output$de_title <- renderText(sprintf("%s - logFC across datasets", input$de_gene %||% ""))
  output$de_plot <- renderPlot({
    d <- de_sub(); validate(need(nrow(d) > 0, "Sem dados p/ esta combinacao."))
    d$label <- factor(d$label, levels = d$label[order(d$log2FoldChange)])
    ggplot(d, aes(log2FoldChange, label, fill = log2FoldChange, alpha = sig)) +
      geom_col() + geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_gradient2(low = "#207394", mid = "grey85", high = "#CD5241") +
      scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.35), guide = "none") +
      labs(x = "log2 fold change", y = NULL, caption = "Solido = passa no padj") +
      theme_minimal(base_size = 13)
  })
  output$de_table <- renderDT({
    de_sub() %>% arrange(padj) %>% select(any_of(c("dataset","contrast","log2FoldChange","padj"))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  # ---- DET: organismo -> canal -> gene ----
  det_ids <- reactive({
    ids <- datasets_all$dataset_id
    if (keep_val(input$det_org))  ids <- ids[datasets_all$organism == input$det_org]
    if (keep_val(input$det_chan)) ids <- intersect(ids, datasets_all$dataset_id[datasets_all$gene == input$det_chan])
    ids
  })
  det_genes <- reactive({ sort(unique(det_df$gene_name[det_df$dataset %in% det_ids()])) })
  observeEvent(det_genes(), {
    g <- det_genes()
    updateSelectizeInput(session, "det_gene", choices = g,
      selected = if ("SCN2A" %in% g) "SCN2A" else g[1], server = TRUE)
  })
  det_sub <- reactive({
    req(input$det_gene)
    det_df %>% filter(gene_name == input$det_gene, dataset %in% det_ids()) %>%
      mutate(sig = qvalue < input$det_q,
             label = paste(transcript_name, dataset, contrast, sep=" - "))
  })
  output$det_title <- renderText(sprintf("%s - transcritos (DET)", input$det_gene %||% ""))
  output$det_plot <- renderPlot({
    d <- det_sub(); validate(need(nrow(d) > 0, "Sem transcritos p/ esta combinacao."))
    d <- d %>% arrange(desc(abs(log2FC))) %>% head(40)
    d$label <- factor(d$label, levels = d$label[order(d$log2FC)])
    ggplot(d, aes(log2FC, label, fill = log2FC, alpha = sig)) +
      geom_col() + geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_gradient2(low = "#207394", mid = "grey85", high = "#CD5241") +
      scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.35), guide = "none") +
      labs(x = "log2FC (transcrito)", y = NULL, caption = "top 40 por |log2FC|; solido passa no qvalue") +
      theme_minimal(base_size = 11)
  })
  output$det_table <- renderDT({
    det_sub() %>% arrange(qvalue) %>%
      select(any_of(c("transcript_name","transcript_id","dataset","contrast","log2FC","qvalue"))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- t-SNE de DET (com filtro organismo + WT) ----
  tsne_det_view <- reactive({
    df <- tsne_det_df
    if (keep_val(input$det_org)) df <- df %>% filter(organism == input$det_org)
    if (isTRUE(input$tsne_det_nowt)) df <- df %>% filter(mutation != "WT")
    df
  })
  output$tsne_det_plot <- renderPlot({
    df <- tsne_det_view(); validate(need(nrow(df) > 0, "Sem amostras."))
    ggplot(df, aes(TSNE1, TSNE2, color = mutation)) +
      geom_point(size = 3, alpha = 0.85) +
      labs(color = "Mutation", x = "t-SNE 1", y = "t-SNE 2", subtitle = "transcript-level") +
      theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank())
  })
  output$tsne_det_info <- renderText({
    cl <- input$tsne_det_click
    if (is.null(cl)) return("Clique num ponto no grafico.")
    np <- nearPoints(tsne_det_view(), cl, xvar = "TSNE1", yvar = "TSNE2", threshold = 20, maxpoints = 1)
    if (nrow(np) == 0) return("Nenhum ponto perto.")
    paste0("Sample:   ", np$sample, "\nDataset:  ", np$dataset,
           "\nGene:     ", np$gene, "\nOrganism: ", np$organism,
           "\nMutation: ", np$mutation)
  })

  # ---- Pathways ----
  path_shared <- reactive({
    paths_df %>% group_by(pathway) %>%
      summarise(n_datasets = n_distinct(dataset),
                datasets = paste(sort(unique(dataset)), collapse=", "), .groups="drop") %>%
      filter(n_datasets >= input$path_min) %>% arrange(desc(n_datasets))
  })
  output$path_plot <- renderPlot({
    ps <- path_shared(); validate(need(nrow(ps) > 0, "Sem vias.")); ps <- head(ps, 40)
    ps$pathway <- factor(ps$pathway, levels = ps$pathway[order(ps$n_datasets)])
    ggplot(ps, aes(n_datasets, pathway, fill = n_datasets)) +
      geom_col() + scale_fill_gradient(low = "#a9cfdf", high = "#207394") +
      labs(x = "numero de datasets", y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$path_table <- renderDT(
    datatable(path_shared(), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)))
}

shinyApp(ui, server)
