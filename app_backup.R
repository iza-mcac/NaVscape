library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

DATA_DIR <- if (dir.exists("data")) "data" else if (dir.exists(file.path("app","data"))) file.path("app","data") else "."
p  <- function(f) file.path(DATA_DIR, f)
rt <- function(f) suppressMessages(readr::read_tsv(f, show_col_types = FALSE))

load_or_demo <- function(fname, demo_fn) {
  fp <- p(fname)
  if (file.exists(fp)) list(data = rt(fp), demo = FALSE) else list(data = demo_fn(), demo = TRUE)
}

datasets_all <- rt(p("datasets_master.tsv"))
datasets <- datasets_all %>% filter(assay %in% c("RNA-seq","ATAC-seq"))
proteome <- if (file.exists(p("proteome_master.tsv"))) rt(p("proteome_master.tsv")) else NULL
leads    <- if (file.exists(p("leads_to_follow.tsv"))) rt(p("leads_to_follow.tsv")) else NULL

tsne_df  <- rt(p("tsne_coords.tsv"))
de_df    <- rt(p("de_all_datasets.tsv"))
paths_df <- rt(p("pathways.tsv"))

gene_pal <- c(SCN1A = "#207394", SCN2A = "#CD5241")

ui <- page_navbar(
  title = "SCN1A / SCN2A Atlas",
  theme = bs_theme(version = 5, primary = "#207394"),

  nav_panel("Datasets", icon = icon("table"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectInput("f_gene", "Gene", c("All", sort(unique(datasets$gene)))),
        selectInput("f_org", "Organism", c("All", sort(unique(datasets$organism)))),
        selectInput("f_sub", "Subtype", c("All", sort(unique(datasets$subtype)))),
        selectInput("f_status", "Status", c("All", sort(unique(datasets$status)))),
        hr(), helpText("Clique numa linha p/ detalhes.")),
      card(card_header("Master dataset table"), DTOutput("dt_datasets")),
      card(card_header("Detalhe da linha"), uiOutput("row_detail")),
      navset_tab(
        nav_panel("Proteome", DTOutput("dt_proteome")),
        nav_panel("Leads", DTOutput("dt_leads")))
    )),

  nav_panel("t-SNE", icon = icon("braille"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        radioButtons("tsne_color", "Colorir por",
          c("Gene"="gene","Organism"="organism","Dataset"="dataset","Mutation"="mutation")),
        hr(), helpText("Clique num ponto p/ inspecionar.")),
      card(card_header("Sample embedding (dataset-corrected)"),
           plotOutput("tsne_plot", click = "tsne_click", height = "560px")),
      card(card_header("Amostra clicada"), verbatimTextOutput("tsne_info"))
    )),

  nav_panel("Gene DE", icon = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        selectizeInput("de_gene", "Gene", choices = NULL,
                       options = list(placeholder = "ex. SCN2A")),
        sliderInput("de_padj", "padj", 0, 1, 0.05, step = 0.01),
        hr(), helpText("logFC do gene em cada dataset/contraste.")),
      card(card_header(textOutput("de_title")), plotOutput("de_plot", height = "540px")),
      card(card_header("Valores"), DTOutput("de_table"))
    )),

  nav_panel("Pathways", icon = icon("diagram-project"),
    layout_sidebar(
      sidebar = sidebar(width = 260,
        sliderInput("path_min", "Min. datasets compartilhando", 1, 6, 1),
        hr(), helpText("Vias em varios datasets = biologia compartilhada.")),
      card(card_header("Vias compartilhadas"),
           plotOutput("path_plot", height = "560px")),
      card(card_header("Via x dataset"), DTOutput("path_table"))
    ))
)

server <- function(input, output, session) {

  ds_filt <- reactive({
    d <- datasets
    if (input$f_gene   != "All") d <- d[d$gene == input$f_gene, ]
    if (input$f_org    != "All") d <- d[d$organism == input$f_org, ]
    if (input$f_sub    != "All") d <- d[d$subtype == input$f_sub, ]
    if (input$f_status != "All") d <- d[d$status == input$f_status, ]
    d
  })

  output$dt_datasets <- renderDT({
    ds_filt() %>%
      select(dataset_id, gene, organism, assay, subtype, n_samples, mechanism, status, source_paper) %>%
      datatable(selection = "single", rownames = FALSE,
                options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatStyle("gene", target = "row",
        backgroundColor = styleEqual(c("SCN1A","SCN2A"), c("#eaf3f7","#fdeeea")))
  })

  output$row_detail <- renderUI({
    sel <- input$dt_datasets_rows_selected
    if (is.null(sel)) return(em("Selecione uma linha acima."))
    r <- ds_filt()[sel, ]
    tagList(
      h5(strong(r$dataset_id),
         if (r$status == "our_dataset") span(class="badge bg-primary", "our dataset")),
      tags$ul(
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

  output$tsne_plot <- renderPlot({
    cc <- input$tsne_color
    ggplot(tsne_df, aes(TSNE1, TSNE2, color = .data[[cc]])) +
      geom_point(size = 3, alpha = 0.85) +
      { if (cc == "gene") scale_color_manual(values = gene_pal) else NULL } +
      labs(color = tools::toTitleCase(cc), x = "t-SNE 1", y = "t-SNE 2") +
      theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank())
  })

  output$tsne_info <- renderText({
    cl <- input$tsne_click
    if (is.null(cl)) return("Clique num ponto no grafico.")
    np <- nearPoints(tsne_df, cl, xvar = "TSNE1", yvar = "TSNE2", threshold = 20, maxpoints = 1)
    if (nrow(np) == 0) return("Nenhum ponto perto - tente de novo.")
    paste0("Sample:   ", np$sample,
           "\nDataset:  ", np$dataset,
           "\nGene:     ", np$gene,
           "\nOrganism: ", np$organism,
           "\nMutation: ", np$mutation)
  })

  updateSelectizeInput(session, "de_gene",
    choices = sort(unique(de_df$gene)),
    selected = if ("SCN2A" %in% de_df$gene) "SCN2A" else sort(unique(de_df$gene))[1],
    server = TRUE)

  de_sub <- reactive({
    req(input$de_gene)
    de_df %>% filter(gene == input$de_gene) %>%
      mutate(sig = padj < input$de_padj,
             label = if ("contrast" %in% names(de_df)) paste(dataset, contrast, sep = " - ") else dataset)
  })

  output$de_title <- renderText(sprintf("%s - logFC across datasets", input$de_gene %||% ""))

  output$de_plot <- renderPlot({
    d <- de_sub(); req(nrow(d) > 0)
    d$label <- factor(d$label, levels = d$label[order(d$log2FoldChange)])
    ggplot(d, aes(log2FoldChange, label, fill = log2FoldChange, alpha = sig)) +
      geom_col() + geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_gradient2(low = "#207394", mid = "grey85", high = "#CD5241") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
      labs(x = "log2 fold change (mut vs ctrl)", y = NULL, caption = "Solido = passa no padj") +
      theme_minimal(base_size = 13)
  })

  output$de_table <- renderDT({
    de_sub() %>% arrange(padj) %>%
      select(any_of(c("dataset","contrast","log2FoldChange","padj"))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  path_shared <- reactive({
    paths_df %>%
      group_by(pathway) %>%
      summarise(n_datasets = n_distinct(dataset),
                datasets = paste(sort(unique(dataset)), collapse = ", "), .groups = "drop") %>%
      filter(n_datasets >= input$path_min) %>%
      arrange(desc(n_datasets))
  })

  output$path_plot <- renderPlot({
    ps <- path_shared(); req(nrow(ps) > 0)
    ps <- head(ps, 40)
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
