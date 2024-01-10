### Libraries ------------------------------------------------------------------
library(shiny)
library(shinythemes)
library(shinyjs) 
library(shinyWidgets)
library(mailtoR) #[!]
library(scales) #[!]
library(dplyr)
library(ggplot2)

### Read core script and create global variables -------------------------------
source('auxiliary.R')

phenos <- readr::read_csv("ukbb_table.csv")[["english_name"]]

gPCs <- 1:40
covar.ids <- readr::read_table('covariate_ids.txt', col_names=FALSE)$X1
pheno.names <- readr::read_table('phenos_names.txt', col_names=FALSE)$X1

lenient_cutoff_choices <- c(1e-8,1e-7,1e-6)
stringent_cutoff_choices <- c(1e-10)

metric_names_df <- data.frame(english_name=c("Pearson r",
                                             "Spearman \u03c1",
                                             "Cosine Similarity",
                                             "Prevalence at Top 10%ile",
                                             "Prevalence at Top 1%ile",
                                             "Odds Ratio at Top 10%ile",
                                             "Odds Ratio at Top 1%ile",
                                             "Percentile-Prevalence Pearson r",
                                             "Percentile-Prevalence Spearman \u03c1",
                                             "Percentile-Average Phenotype Pearson r",
                                             "Percentile-Average Phenotype Spearman \u03c1"),
                              col_name=c("PEARSON_CORR","SPEARMAN_RHO","COSINE_SIM",
                                         "TOP10PCT_PREV","TOP1PCT_PREV","TOP10PCT_OR",
                                         "TOP1PCT_OR","PERCENTILE_AVE_PREV_SLOPE",
                                         "PERCENTILE_AVE_PREV_SPEARMAN",
                                         "PERCENTILE_AVE_PHENO_SLOPE",
                                         "PERCENTILE_AVE_PHENO_SPEARMAN"))

### Define UI ------------------------------------------------------------------
about_page <- tabPanel(
  title = "Home",
  setBackgroundColor(
    color = c("#F5FFFA", "#E6E6FA"),
    gradient = "linear",#E6E6FA"
    direction = "bottom"
  ),
  tagList(
    lapply(1:6, function(i) br())),
  fluidRow(
    column(width=3),
    column(width=6,
           
           div(style = "text-align:center;,
      padding-left: 4px; padding-right: 4px; padding-top: 4px;
      padding-bottom: 4px;padding:4px;
      dborder-color:black;background-color: white;",
               h1("Polygenic Scores and Stratification Through the Lens of Principal Components", align='center'),
               br(),
               div(
                 style = "text-align:center;",
                 actionButton("select_phenotype", "Select Phenotype", icon("paper-plane"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 actionButton("learn_more", "Learn More", icon("circle-info"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 br()),
               br(), width=8,align="center"),
    )),
  tagList(
    lapply(1:3, function(i) br())),
  h4("Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song. (2024+) 'Highly parameterized polygenic scores tend to overfit to population stratification via random effects.' Submitted.", align='center'),
  h4(paste('Last updated:',Sys.Date()), align='center'),
  tagList(
    lapply(1:4, function(i) br())),
  fluidRow(
    column(2),
    column(2),
    column(1),
    column(1,
           div(  align="center",  
                 tags$a(href='https://people.eecs.berkeley.edu/~yss/group.html',
                        tags$img(src="UC_Berkeley_Seal_80px.png", height='70', width='70'),
                        ''))),
    column(1,
           tags$a(href='https://github.com/songlab-cal',
                  tags$img(src="github-mark.png", height='70', width='70')
           ))
  ),
  tagList(
    lapply(1:3, function(i) br())))


## About page 
source('Summary_about_page.R', local = TRUE)

theme_flatly <- shinytheme("flatly")

# Gene_page <- tabPanel(
#   
#   title='Gene',value="Gene",
#   titlePanel("Genetic Priority Score"),
#   theme=theme_flatly
# )

Phenotype_page <- tabPanel(
  
  title='Phenotype', value="Phenotype",
  titlePanel("Phenotype Stratification"),
  sidebarLayout(
    sidebarPanel(
      # Title
      # Select Phenotype     
      h3("Phenotype Overview"),
      fluidRow(
        column(12,selectInput(inputId = "pheno_selector",
                             label = "Pick a Phenotype",
                             choices = phenos))
      ),
      helpText("Phenotype examples: Albumin, Standing height, Cholesterol"),
      fluidRow(
        column(12,selectInput(inputId = "pheno_version_selector",
                              label = "Pick Phenotype Version",
                              choices = c("Original",
                                          "Inverse Rank Normal Transformed")))
      ),
      h4("UK Biobank (UKB) Metadata"),
      uiOutput("phenoInfoSummary"),
      uiOutput("UKBB_Reference"),
      h4("Polygenic Scores (PGSs)"),
      uiOutput("PGSSummary"),
      uiOutput("zeroPGS")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stratification Profile",
                 fluidRow(column(6,h4("1. Cosine Similarity Profile"),
                                 plotOutput(outputId = "phenoCossimPlot")),
                          column(6,h4("2. Pearson Correlation Profile"),
                                 plotOutput(outputId = "phenoCorrPlot"))
                          ),
                 tags$li("The profiles above show the distribution of correlation and similarity measures of the phenotype and genetic PCs (gPCs)."),
                 tags$li("For example, the bar above gPC_1 depicts the Pearson correlation and cosine similarity between the phenotype and gPC1."),
                 tags$li("The gPC with the largest magnitude is coloured red.")),
        tabPanel("Other Statistics",
                 h4("Stratification Profile Summary and Prediction Inflation"),
                 column(12, align="center", 
                        tableOutput("phenoSummary"),
                        ),
                 actionButton("PhenoStratShowBtn", "Show Details"), # [TESTING]
                 actionButton("PhenoStratHideBtn", "Hide Details"), # [TESTING]
                 shinyjs::hidden(shiny::tags$ul(id = "pheno_strat_details",
                                       shiny::tags$li(p(shiny::tags$b("Cosine Similarity Evenness"), 
                                                        "computes the Shannon entropy of the 40 unsigned and normalized cosine similarities. Namely, if $(x_1,\\ldots,x_{40})$ denote the cosine similarities with the 40 gPCs as shown in the Stratification Profile tab, then let $y_i=|x_i|/(|x_1|+\\ldots+|x_{40}|)$ be the unsigned and normalized cosine similarity with the $i$th gPC, so that $y_1+\\ldots+y_{40}=1$. The Shannon entropy is $\\sum_{i=1}^{40}y_i\\log (1/y_i)$, which captures the evenness of the distribution of cosine similarities. Larger values imply greater evenness (Max $=\\log(40)\\approx 3.69$; Min = $0$).")), 
                                       shiny::tags$li(p("Similarly,",
                                                        shiny::tags$b("Pearson r Evenness"),
                                                        "computes the Shannon entropy of the 40 unsigned and normalized Pearson correlations. Larger values imply greater evenness (Max $=\\log(40)\\approx 3.69$; Min = $0$).")),
                                       shiny::tags$li(p(shiny::tags$b("Mean Incremental R\u00B2"), 
                                                        "reports the average incremental R\u00B2, across 200 random polygenic scores (rPGSs) constructed by arbitrarily choosing 10% of all autosomal variants and independently assigning standard Gaussian effects to them (i.e., $\\beta_j\\sim N(0,1)$ for variant $j$ chosen). For each rPGS, R\u00B2 was computed on a linear model including the rPGS, alongside age, sex and the top 20 gPCs. This is subtracted by the baseline model including just the latter covariates, to obtain the Incremental R\u00B2.")),
                                       shiny::tags$li(p(shiny::tags$b("Incremental R\u00B2 p-value "), 
                                                        "reports the significance of the Mean Incremental R\u00B2 described previously. Significance is obtained by permuting the phenotypes randomly across the individuals, and repeating the same procedure of computing average incremental R\u00B2, for 100 times. The empirical p-value reports the fraction of permuted Mean Incremental R\u00B2's that beat the observed ($(\\#\\text{beat} + 1)/101$). A smaller p-value indicates greater significance of the observed Mean Incremental R\u00B2.")) 
                                       
                                       
                 ))
                 )
      )                                    
    )),
)

PGS_page <- tabPanel(
  
  title='PGS', value="PGS",
  titlePanel("PGS Stratification"),
  tabsetPanel(
    tabPanel("Quick Overview",
             p(uiOutput("pgsSentence")),
             p(uiOutput("PGSmsg")),
             p(HTML("To select another phenotype, go back to <b>Phenotype</b> tab.")),
             fluidRow(
               column(5,
                      h4(uiOutput("pgsChromDistSentence")),
                      plotOutput(outputId = "pgsStratPlot")),
               column(7,h4(uiOutput("pgsStratSentence")),
                      tableOutput("pgsStratSummary"))
             )),
    tabPanel("Polygenic Architecture and Performance Sensitivity",
             p(uiOutput("ifnoPGS")))
  )
)

ui <- fluidPage(
  withMathJax(),
  shinyjs::useShinyjs(),
  tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$']]}
                });
                </script>
                ")),
  navbarPage(
    id='tabset',
    title = 
      div(
        div(
          img(src="dna_pca_3.png",
              height = "45px",width = "45px",
              style = "position: relative; margin:-20px -20px; display:right-align;")),''),
    position = 'fixed-top',
    windowTitle="Principal Components and Effect Randomness",
    header=tags$style(type="text/css", "body {padding-top: 70px;},
  
  "),
    
    theme=shinytheme('cosmo'),
    about_page,
    Phenotype_page,
    PGS_page,
    Aboutinfo_page
  ))

server <- function(input, output,session) {
  ## Homepage selectors
  observeEvent(input$select_phenotype, {
    updateTabsetPanel(session, "tabset", selected = "Phenotype")
  })
  
  observeEvent(input$learn_more, {
    updateTabsetPanel(session, "tabset", selected = "About")
  })
  
  ############## Phenotype Tab ##############
  # Quick Information about Phenotype ------------------------------------------
  summarizedPheno <- reactive({
    req(input$pheno_selector)
    summarizePheno(name=input$pheno_selector)
  })
  
  output$phenoInfoSummary <- renderUI({
    category <- req(try(summarizedPheno()[["Category"]]))
    field <- summarizedPheno()[["Field"]]
    line_1 <- paste0("Category: ", category)
    line_2 <- paste0("UKB Field: ", field)
    HTML(paste(line_1, line_2, sep = '<br/>'))
  })
  
  # URL to UKBB Data Portal
  output$UKBB_Reference <- renderUI(a(href=summarizedPheno()[["URL"]], 
                                      "Reference", 
                                      target="_blank"))
  
  # Quick Information about PGSs associated with Phenotype ---------------------
  # This enables selection of PGS matched to a selected phenotype
  observeEvent(input$pheno_selector, {
    message("Table event observed -- phenotype selected")
    no_avail_pgs <- (ukbb_table %>% 
                       subset(english_name==input$pheno_selector))[["no_avail_pgs"]]
    if (no_avail_pgs == 0) {
      choices <- NA
    } else if (no_avail_pgs == 1) {
      choices <- "Clumping and thresholding (lenient)"
    } else if (no_avail_pgs == 2) {
      choices <- c("Clumping and thresholding (lenient)",
                   "Clumping and thresholding (stringent)")
    }
    updateSelectInput(session, "pgs_selector", choices = choices)
  })
  
  output$PGSSummary <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==1) {
      line_1 <- paste0("There is ", no_avail_pgs, " available for phenotype selected.")
      line_2 <- "View <b>PGS</b> tab to explore metrics for the PGS."
      line_3 <- "."
      HTML(paste(line_1,line_2,sep='<br/>'))
    } else if (no_avail_pgs > 1) { 
      line_1 <- paste0("There are ", no_avail_pgs, " available for phenotype selected.")
      line_2 <- "Choose a PGS below and explore metrics under <b>PGS</b> tab."
      HTML(paste(line_1,line_2,sep='<br/>'))
    } else {
      HTML("No PGSs generated for phenotype selected.")
    }
  })
  
  output$zeroPGS <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==0) {
      return()
    } else {
      fluidRow(
        column(12,selectInput(inputId = "pgs_selector", 
                             label = "",
                             choices = ""))
      )
    }
  })
  
  # Phenotype stratification profile -------------------------------------------
  output$phenoCossimPlot <- renderPlot(req(try(plot(getPhenoCossimPlot(x=input$pheno_selector,
                                                                       y=input$pheno_version_selector) + 
                                                      theme(text=element_text(family='DejaVu Sans'))
  ))))
  output$phenoCorrPlot <- renderPlot(req(try(plot(getPhenoCorrPlot(x=input$pheno_selector,
                                                                   y=input$pheno_version_selector) + 
                                                    theme(text=element_text(family='DejaVu Sans'))))))
  
  # Phenotype stratification and rPGS statistics -------------------------------
  phenoStats <- reactive({getPhenoStats(x=input$pheno_selector,
                                        y=input$pheno_version_selector)})
  output$phenoSummary <- renderTable(req(try(phenoStats()$TABLE)), 
                                     digits = 4,
                                     striped = TRUE,
                                     bordered = TRUE)
  output$phenoSentence <- renderUI(req(try(phenoStats()$SENTENCE)))
  
  # Enable hiding and showing of phenotype stratification statistics details
  observeEvent(input$PhenoStratShowBtn,
               {show("pheno_strat_details",asis=TRUE)})
  observeEvent(input$PhenoStratHideBtn,
               {hide("pheno_strat_details",asis=TRUE)})
  
  ############## PGS Tab ##############
  # PGS stratification profile and chromosome plot -----------------------------
  output$PGSmsg <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==0) {
      HTML(summarizedPheno()[["Zero_PGS_Sentence"]])
    } else {
      HTML("We analyze population stratification of the PGS on its training cohort of 288,728 individuals of European descent.")
    }
  })
  # PGS Stratification and Quick Overview Data 
  pgsStratStats <- reactive({getPGSStratStats(x=input$pheno_selector,
                                              y=input$pgs_selector)})
  
  output$pgsStratSummary <- renderTable(req(try(pgsStratStats()$TABLE)),
                                        digits = 4,
                                        hover = TRUE,
                                        bordered = TRUE)
  output$pgsSentence <- renderUI(req(try(pgsStratStats()$SENTENCE)))
  output$pgsStratPlot <- renderPlot(req(try(pgsStratStats()$PLOT)))
  output$pgsChromDistSentence <- renderUI(req(try(pgsStratStats()$CHROM_DIST_SENTENCE)))
  output$pgsStratSentence <- renderUI(req(try(pgsStratStats()$PGS_STRAT_SENTENCE)))
  
  # PGS polygenic architecture and performance sensitivity
  output$ifnoPGS <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==0) {
      HTML(summarizedPheno()[["Zero_PGS_Sentence"]])
    } else {
      HTML("We analyze the perturbed-fixed architecture of the PGS, and the PGS sensitivity on a held-out test cohort of 68,931 individuals of European descent.")
    }
  })
  
  ## Gene Page
  
  ##Evidence tables
  
  ### DOE Evidence   
  
  ##Evidence tables
  
  ## Phenotype Page
  
  ## Evidence phenotype
  
  ## Evidence tables
  
  
}

shinyApp(ui, server)
