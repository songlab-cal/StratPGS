### Libraries -----------------------------------------------------------
library(shiny)
library(shinythemes)
library(shinyjs) 
library(shinyWidgets)
library(mailtoR) #[!]
library(scales) #[!]
library(dplyr)
library(ggplot2)

### Read core script and create global variables ------------------------
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

### Define UI -------------------------------------------------------------
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
                 #actionButton("Search_by_phenotype", "Search by Indication",class="btn btn-info"),
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
                 p("The profiles show the distribution of correlation and similarity measures of the phenotype and genetic PCs (gPCs).")),
        tabPanel("Other Statistics")
      )                                    
    )),
)

PGS_page <- tabPanel(
  
  title='PGS', value="PGS",
  titlePanel("PGS Stratification")
)

ui <- fluidPage(
  withMathJax(),
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
              height = "45px",width = "45px",style = "position: relative; margin:-20px -20px; display:right-align;")),''),
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
  
  ## Phenotype Tab
  # Quick Information about Phenotype --------------------------------------------
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
  
  # Quick Information about PGSs associated with Phenotype ------------------------
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
  
  # Phenotype stratification profile ------------------------------------------
  output$phenoCossimPlot <- renderPlot(req(try(plot(getPhenoCossimPlot(x=input$pheno_selector,
                                                                       y=input$pheno_version_selector) + 
                                                      theme(text=element_text(family='DejaVu Sans'))
  ))))
  output$phenoCorrPlot <- renderPlot(req(try(plot(getPhenoCorrPlot(x=input$pheno_selector,
                                                                   y=input$pheno_version_selector) + 
                                                    theme(text=element_text(family='DejaVu Sans'))))))
  
  ## Gene Page
  
  ##Evidence tables
  
  ### DOE Evidence   
  
  ##Evidence tables
  
  ## Phenotype Page
  
  ## Evidence phenotype
  
  ## Evidence tables
  
  
}

shinyApp(ui, server)
