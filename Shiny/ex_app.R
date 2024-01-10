#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://urldefense.com/v3/__http://shiny.rstudio.com/__;!!HrbR-XT-OQ!XTLJTajDBjHW_UuUobS4sK4S3CRNVMaFJ-BCZ9r5D5J6N5W9fHD-aMs9XH_slXMvOf3NXHuQqAS0BJuXLg$ 
#
# This version dated: 12/20/23
# 

library(shiny)
library(shinythemes)
library(shinyjs) # [TESTING]
library(dplyr)
library(ggplot2)
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

# Define UI 
ui <- fluidPage(
  useShinyjs(), # [TESTING]
  theme = shinytheme("cosmo"),
  withMathJax(),
  tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$']]}
                });
                </script>
                ")),
  tags$head(
    tags$style(HTML("
      .selectize-input {
        height: 35px;
        width: 400px;
        font-size: 12pt;
        padding-top: 5px;
      }
    "))
  ),
  # Application title
  titlePanel("Polygenic Scores and Stratification Through the Lens of Principal Components"),
  p("Created by: Alan Aw"),
  p("To learn more about our work, check out our preprint: Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song. (2024+) 'Highly parameterized polygenic scores tend to overfit to population stratification via random effects.' (in preparation). We recommend browsing this app on a Desktop for optimal user experience."),
  # Sidebar with user input selection 
  sidebarLayout(
    # SIDE BAR PANEL
    sidebarPanel(
      HTML('<center><img src="visualization.png" width="300"></center>'),
      HTML('<center>High-dimensional statistics meets human genetics.</center>'),
      h3("Pick a Phenotype"),
      fluidRow(
        column(7,selectInput(inputId = "pheno_selector",
                             label = "",
                             choices = phenos))
      ),
      h3("Summary"),
      h4("UK Biobank (UKB) Phenotype Metadata"),
      uiOutput("phenoInfoSummary"),
      uiOutput("UKBB_Reference"),
      p("Explore stratification of phenotype under ", shiny::tags$b("Phenotype Stratification"), " tab"),
      h4("Polygenic Scores (PGSs)"),
      uiOutput("PGSSummary"),
      uiOutput("zeroPGS"),
      h3("Acknowledgements"),
      p("Our work is funded by both NIH (R35-GM134922) and the Chan Zuckerberg Initiative Foundation (CZF2019-002449). Our work is approved by the UKB under application number 33751. The graphic above was generated with the help of DALL·E 3."),
      h4("Feedback and Troubleshooting"),
      p("We welcome any feedback on improving the accessibility of this dashboard. Please send all feedback to alan(dot)aw(at)pennmedicine(dot)upenn(dot)edu.")
    ),
    # MAIN PANEL
    mainPanel(
      tabsetPanel(
        # Phenotype Stratification
        tabPanel("Phenotype Stratification",
                 p(""),
                 p("We analyze population stratification of the phenotype on a cohort of 487,296 individuals, mostly of self-identified European ancestry."),
                 fluidRow(
                   column(12,h3("Pick Phenotype Version"),
                          # Input: Select dataset (in the future, build upload functionality)
                          selectInput(inputId = "pheno_version_selector",
                                      label = "",
                                      choices = c("Original", 
                                                  "Inverse Rank Normal Transformed")))
                 ),
                 h3("Stratification Information"),
                 p("Below summarizes stratification bias and relationships with the top 40 genetic PCs (gPCs) for the selected phenotype and version choice."),
                 fluidRow(column(6,h4("1. Cosine Similarity Profile"),
                                 plotOutput(outputId = "phenoCossimPlot")),
                          column(6,h4("2. Pearson Correlation Profile"),
                                 plotOutput(outputId = "phenoCorrPlot")),
                          h4(uiOutput("phenoSentence")),
                          tableOutput("phenoSummary"),
                          actionButton("PhenoStratShowBtn", "Show Details"), # [TESTING]
                          actionButton("PhenoStratHideBtn", "Hide Details"), # [TESTING]
                          hidden(shiny::tags$ul(id = "pheno_strat_details",
                                                shiny::tags$li(p("The", shiny::tags$b("Cosine Similarity Profile"), 
                                                                 "and", shiny::tags$b("Pearson Correlation Profile"),
                                                                 "summarize empirical distributions of cosine similarities and Pearson correlations between the phenotype vector and the gPC vectors. The gPC with the largest magnitude is coloured red.")), 
                                                shiny::tags$li(p(shiny::tags$b("Cosine Similarity Evenness"), 
                                                                 "computes the Shannon entropy of the 40 unsigned and normalized cosine similarities. Namely, if $(x_1,\\ldots,x_{40})$ denote the cosine similarities with the 40 gPCs as shown in the left plot above, then let $y_i=|x_i|/(|x_1|+\\ldots+|x_{40}|)$ be the unsigned and normalized cosine similarity with the $i$th gPC, so that $y_1+\\ldots+y_{40}=1$. The Shannon entropy is $\\sum_{i=1}^{40}y_i\\log (1/y_i)$, which captures the evenness of the distribution of cosine similarities. Larger values imply greater evenness (Max $=\\log(40)\\approx 3.69$; Min = $0$).")), 
                                                shiny::tags$li(p("Similarly,",
                                                                 shiny::tags$b("Pearson r Evenness"),
                                                                 "computes the Shannon entropy of the 40 unsigned and normalized Pearson correlations. Larger values imply greater evenness (Max $=\\log(40)\\approx 3.69$; Min = $0$).")),
                                                shiny::tags$li(p(shiny::tags$b("Mean Incremental R\u00B2"), 
                                                                 "reports the average incremental R\u00B2, across 200 random polygenic scores (rPGSs) constructed by arbitrarily choosing 10% of all autosomal variants and independently assigning standard Gaussian effects to them (i.e., $\\beta_j\\sim N(0,1)$ for variant $j$ chosen). For each rPGS, R\u00B2 was computed on a linear model including the rPGS, alongside age, sex and the top 20 gPCs. This is subtracted by the baseline model including just the latter covariates, to obtain the Incremental R\u00B2.")),
                                                shiny::tags$li(p(shiny::tags$b("Incremental R\u00B2 p-value "), 
                                                                 "reports the significance of the Mean Incremental R\u00B2 described previously. Significance is obtained by permuting the phenotypes randomly across the individuals, and repeating the same procedure of computing average incremental R\u00B2, for 100 times. The empirical p-value reports the fraction of permuted Mean Incremental R\u00B2's that beat the observed. A smaller p-value indicates greater significance of the observed Mean Incremental R\u00B2.")) 
                                                
                                                
                          ))),
        ),
        
        # Panel: PGS Stratification
        tabPanel("PGS Stratification",
                 p(""),
                 h3("Stratification and Quick Overview"),
                 p(uiOutput("PGSmsg")),
                 p(uiOutput("pgsSentence")),
                 fluidRow(
                   id="strat_diagnosis",
                   column(5,
                          h4(uiOutput("pgsChromDistSentence")),
                          plotOutput(outputId = "pgsStratPlot")),
                   column(7,h4(uiOutput("pgsStratSentence")),
                          tableOutput("pgsStratSummary"))
                 ),
                 h3("Performance Sensitivity and Polygenic Architecture"),
                 p(uiOutput("ifnoPGS")),
                 uiOutput("PGScutoffselect"),
                 h4(uiOutput("pgsPerturbFixedHeader")),
                 uiOutput("pgsPerturbFixedSentence"),
                 tableOutput("pgsPerturbFixedSummary"),
                 uiOutput("PGSmetricselect"),
                 h4(uiOutput("pgsSensitivityHeader")),
                 uiOutput("pgsSensitivitySentence"),
                 fluidRow(
                   id="sensitivity_diagnosis",
                   column(7,
                          plotOutput(outputId = "pgsSensitivityPlot")),
                   column(5,
                          tableOutput("pgsSensitivitySummary"))
                 ),
        ),
        # Panel: Exogenous Covariates (age, sex, etc.)
        tabPanel("Exogenous Covariates (Beta)",
                 # Let user choose Covariate
                 fluidRow(
                   p("We have computed the distributions of correlations and cosine similarities between an exogenous covariate and genetic PCs."),
                   p("These computations do not require phenotype input; rather they are an exploratory data-analytic device to help the geneticist identify unusual stratification of a non-genetic covariate along the genetic PC axes."),
                   h3("Pick a Covariate"),
                   # Input: Select dataset (in the future, build upload functionality)
                   selectInput(inputId = "covar",
                               label = "",
                               choices = covar.ids)
                 ),
                 # Show table of statistics
                 fluidRow(
                   tabPanel(title = "Relevance to Genetic PCs",
                            p("Below summarizes statistical relationships between the selected covariate and the top 40 genetic PCs of the cohort."),
                            fluidRow(column(6,h4("1. Cosine Similarity Profile"),
                                            plotOutput(outputId = "cossimPlot")),
                                     column(6,h4("2. Pearson Correlation Profile"),
                                            plotOutput(outputId = "corrPlot")
                                            #actionButton("PhenoStratShowBtn", "Show Details"), # [TESTING]
                                            #actionButton("PhenoStratHideBtn", "Hide Details"), # [TESTING]
                                            )),
                            
                   )
                )
        ), 
        # Summary of Approach
        tabPanel("Methodology",
                 h3("Overview"),
                 p("The core thesis of our work is that polygenic scores (PGSs) using too many variants are overparameterized, with the overparameterization capturing population structure."),
                 p("Our thesis is justified by the following arguments."),
                 tags$li("Theoretically, PGSs can be viewed as linear combinations of PCs, or linear combinations of singular vectors with skewed weights that capture genetic structure."),
                 tags$li("This implies that choosing completely random variant effects for a PGS leads to ''better than random'' performance, if the phenotype also happens to be stratified. We observe this in UKB data."),
                 tags$li("Empirically, permuting or flipping the signs of some variants --- specifically variants with high GWAS p-value --- in a PGS does not weaken the performance in a held-out cohort."),
                 tags$li("The extent to which PGS performance is weakened by the above procedure is driven by stratification of PGS by PCs in the discovery cohort."),
                 h3("Computational Experiments"),
                 p("Our experiments involving generating random variant effects, as well as perturbing effects of existing PGSs, are summarized in the graphic below."),
                 HTML('<center><img src="outline_graphic.jpg" width="800"></center>'),
                 h3("Exogenous Covariate Stratification"),
                 p("Our tab exploring stratification of exogenous covariates by PCs is work in progress.")
        )
      )
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  # Enable selection of PGS matched to a selected phenotype
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
  
  # Enable hiding and showing of phenotype stratification statistics details
  observeEvent(input$PhenoStratShowBtn,
               {show("pheno_strat_details")})
  observeEvent(input$PhenoStratHideBtn,
               {hide("pheno_strat_details")})
  
  # Information about Phenotype (1)
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
  
  output$PGSSummary <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==1) {
      line_1 <- paste0("There is ", no_avail_pgs, " available for phenotype selected.")
      line_2 <- "View <b>PGS Stratification</b> tab to explore metrics for the PGS."
      line_3 <- "."
      HTML(paste(line_1,line_2,sep='<br/>'))
    } else if (no_avail_pgs > 1) { 
      line_1 <- paste0("There are ", no_avail_pgs, " available for phenotype selected.")
      line_2 <- "Choose a PGS below and explore metrics under <b>PGS Stratification</b> tab."
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
        column(7,selectInput(inputId = "pgs_selector", 
                             label = "", 
                             choices = ""))
      )
    }
  })
  
  # PGS cutoff selector
  output$PGScutoffselect <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]] 
    if (no_avail_pgs==0) {
      return(NULL)
    } else {
      pgs_type <- req(try(input$pgs_selector))
      if (pgs_type=="Clumping and thresholding (lenient)") {
        fluidRow(
          column(5,selectInput(inputId = "pgs_cutoff_selector", 
                               label = "Select perturbation cutoff", 
                               choices = c(1e-6,1e-7,1e-8)))
        )
      } else if (pgs_type=="Clumping and thresholding (stringent)") {
        fluidRow(
          column(5,selectInput(inputId = "pgs_cutoff_selector", 
                               label = "Select perturbation cutoff", 
                               choices = c(1e-10)))
        )
      }
    }
  })
  
  # PGS performance metric selector
  # PGS cutoff selector
  output$PGSmetricselect <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]] 
    if (no_avail_pgs==0) {
      return(NULL)
    } else {
      fluidRow(
        column(5,selectInput(inputId = "pgs_metric_selector", 
                             label = "Select performance metric", 
                             choices = metric_names_df$english_name))
      )
    }
  })
  # Information about Phenotype (2)
  # URL to UKBB Data Portal
  output$UKBB_Reference <- renderUI(a(href=summarizedPheno()[["URL"]], 
                                          "Reference", 
                                          target="_blank"))
  
  # Information about PGS (1)
  output$PGSmsg <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==0) {
      HTML("No PGS data for chosen phenotype.")
    } else {
      HTML("We analyze population stratification of the PGS on its training cohort of 288,728 individuals of European descent.")
    }
  })
  
  output$ifnoPGS <- renderUI({
    no_avail_pgs <- summarizedPheno()[["No_Avail_PGS"]]
    if (no_avail_pgs==0) {
      HTML("Not Applicable for chosen phenotype.")
    } else {
      HTML("We analyze sensitivity of the PGS on a held-out test cohort of 68,931 individuals of European descent.")
    }
  })
  
  # Plots for covariate
  output$cossimPlot <- renderPlot(req(try(plot(getCossimPlot(x=input$covar) + theme(text=element_text(family='DejaVu Sans'))))))
  output$corrPlot <- renderPlot(req(try(plot(getCorrPlot(x=input$covar) + theme(text=element_text(family='DejaVu Sans'))))))
  
  # Plots for phenotype
  output$phenoCossimPlot <- renderPlot(req(try(plot(getPhenoCossimPlot(x=input$pheno_selector,
                                                                       y=input$pheno_version_selector) + 
                                                      theme(text=element_text(family='DejaVu Sans'))
                                                    ))))
  output$phenoCorrPlot <- renderPlot(req(try(plot(getPhenoCorrPlot(x=input$pheno_selector,
                                                                   y=input$pheno_version_selector) + 
                                                    theme(text=element_text(family='DejaVu Sans'))))))
  
  # Table for phenotype
  phenoStats <- reactive({getPhenoStats(x=input$pheno_selector,
                                        y=input$pheno_version_selector)})
  output$phenoSummary <- renderTable(req(try(phenoStats()$TABLE)), 
                                     digits = 4,
                                     bordered = TRUE)
  output$phenoSentence <- renderUI(req(try(phenoStats()$SENTENCE)))
  
  # PGS Stratification and Quick Overview Data 
  pgsStratStats <- reactive({getPGSStratStats(x=input$pheno_selector,
                                         y=input$pgs_selector)})
  
  #output$pgsStratSummary <- DT::renderDataTable(req(try(pgsStratStats()$TABLE))) #[!] Try regular table for now
  output$pgsStratSummary <- renderTable(req(try(pgsStratStats()$TABLE)),
                                        digits = 4,
                                        bordered = TRUE)
  output$pgsSentence <- renderUI(req(try(pgsStratStats()$SENTENCE)))
  output$pgsStratPlot <- renderPlot(req(try(pgsStratStats()$PLOT)))
  output$pgsChromDistSentence <- renderUI(req(try(pgsStratStats()$CHROM_DIST_SENTENCE)))
  output$pgsStratSentence <- renderUI(req(try(pgsStratStats()$PGS_STRAT_SENTENCE)))
  
  # PGS Perturb-Fixed Architecture Data
  pgsPerturbFixed <- reactive({getPGSStats1(x=input$pheno_selector,
                                  y=input$pgs_selector,
                                  z=input$pgs_cutoff_selector)})
  output$pgsPerturbFixedSentence <- renderUI(req(try(pgsPerturbFixed()$SENTENCE)))
  output$pgsPerturbFixedSummary <- renderTable(req(try(pgsPerturbFixed()$TABLE)),
                                            digits=4,bordered=TRUE)
  output$pgsPerturbFixedHeader <- renderUI(req(try(pgsPerturbFixed()$PERTURB_FIX_HEADER)))
  
  # PGS Sensitivity Data
  pgsSensitivity <- reactive({getPGSStats2(x=input$pheno_selector,
                                           y=input$pgs_selector,
                                           z=input$pgs_cutoff_selector,
                                           w=input$pgs_metric_selector)})
  output$pgsSensitivitySummary <- renderTable(req(try(pgsSensitivity()$TABLE)),
                                              bordered=TRUE)
  output$pgsSensitivityPlot <- renderPlot(req(try(pgsSensitivity()$PLOT)))
  output$pgsSensitivityHeader <- renderUI(req(try(pgsSensitivity()$HEADER)))
  output$pgsSensitivitySentence <- renderUI(req(try(HTML(pgsSensitivity()$SENTENCE))))
}



# Run the application 
shinyApp(ui = ui, 
         server = server)