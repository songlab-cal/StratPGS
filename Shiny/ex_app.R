#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://urldefense.com/v3/__http://shiny.rstudio.com/__;!!HrbR-XT-OQ!XTLJTajDBjHW_UuUobS4sK4S3CRNVMaFJ-BCZ9r5D5J6N5W9fHD-aMs9XH_slXMvOf3NXHuQqAS0BJuXLg$ 
#
# This version dated: 11/14/23
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
  p("To learn more about our work, check out our preprint: Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song. (2023+) 'The pitfalls of overparameterized polygenic scores and new diagnostics to combat overfitting.' (in preparation). We recommend browsing this app on a Desktop for optimal user experience."),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    # SIDE BAR PANEL
    sidebarPanel(
      HTML('<center><img src="visualization.png" width="300"></center>'),
      HTML('<center>High-dimensional statistics meets human genetics, in the style of Володи́мирч Та́тлін.</center>'),
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
      uiOutput("zeroPGS")
      
      # fluidRow(
      #   column(7,selectInput(inputId = "pgs_selector", 
      #                        label = "Polygenic Score (PGS)", 
      #                        choices = ""))
      # )
      
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
                 p(uiOutput("PGSmsg")),
                 h4("1. Quick Overview"),
                 p(uiOutput("pgsSentence")),
                 plotOutput(outputId = "pgsStratPlot"),
                 h4("2. Stratification Statistics"),
                 tableOutput("pgsStratSummary") # [!] To do: Fix the table (do a 3 x 4?)
                 # [!] Possibility: combine the table and the plot in one row
                 # [!]              add second row that allows selection of performance metric 
                 #                  (view plot and sensitivity)
                 #DT::dataTableOutput("pgsStratSummary")
        ),
        # Panel: Exogenous Covariates (age, sex, etc.)
        tabPanel("Exogenous Covariates (Beta)",
                 # Let user choose Covariate
                 fluidRow(
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
                            fluidRow(column(6,h3("1. Cosine Similarity"),
                                            plotOutput(outputId = "cossimPlot"),
                                            h4("SNP Information (dbSNP)"),
                                            fluidRow(column(3,uiOutput(outputId = "PICS2TopDBSNP")),
                                                     column(3,uiOutput(outputId = "PICS2StableDBSNP")))),
                                     column(6,h3("2. Pearson Correlation"),
                                            plotOutput(outputId = "corrPlot"),
                                            actionButton("PhenoStratShowBtn", "Show Details"), # [TESTING]
                                            actionButton("PhenoStratHideBtn", "Hide Details"), # [TESTING]
                                            hidden(shiny::tags$ul(id = "stab_stats",
                                              shiny::tags$li(p("The", shiny::tags$b("No. Slices Reporting Variant with Positive Probability"),
                                                               "measures how often the stable variant has a posterior probability larger than zero, when performing fine-mapping across all user-defined slices (five in total for this study).")), 
                                              shiny::tags$li(p("Based on the quantity measured above,", 
                                                               shiny::tags$b("Max Slice-Slice F_ST"), 
                                                               "computes the maximum",
                                                               HTML(paste0(shiny::tags$i('F'),shiny::tags$sub('ST'))),
                                                               "(a measure of genetic differentiaton; see",
                                                               shiny::tags$a("Holsinger and Weir, 2009", href="https://urldefense.com/v3/__https://www.nature.com/articles/nrg2611__;!!HrbR-XT-OQ!XTLJTajDBjHW_UuUobS4sK4S3CRNVMaFJ-BCZ9r5D5J6N5W9fHD-aMs9XH_slXMvOf3NXHuQqAQl-Ko_rg$ "),
                                                               ") between any pair of slices for which the variant has a positive posterior probability. Larger values imply greater genetic differentiation.")), 
                                              shiny::tags$li(p("Similarly,",
                                                               shiny::tags$b("Max Slice-Slice AF Difference"),
                                                               "computes the maximum",
                                                               shiny::tags$i('mean allele frequency difference'),
                                                               "(a measure of genetic heterogeneity) between any pair of slices for which the variant has a positive posterior probability. Larger values imply greater genetic heterogeneity."))
                                            )))),
                            
                   )
                )
        ), 
        # Summary of Approach
        tabPanel("Methodology",
                 h3("Cohort"),
                 p("We use 502,506 individuals from the UK Biobank."),
                 h3("Relative Performance"),
                 p("To be added.")
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
  
  # Enable hiding and revealing of functional plots (Genome-wide tab) # [TESTING]
  observeEvent(input$funcAnnotPlotShowBtn,
               {show("funcAnnotPlot")})
  observeEvent(input$funcAnnotPlotHideBtn,
               {hide("funcAnnotPlot")})
  
  observeEvent(input$funcDensityPlotShowBtn,
               {show("funcDensityPlot")})
  observeEvent(input$funcDensityPlotHideBtn,
               {hide("funcDensityPlot")})
  
  # Enable hiding and showing of stability statistics details (Single gene tab) # [TESTING]
  observeEvent(input$PhenoStratShowBtn,
               {show("pheno_strat_details")})
  observeEvent(input$PhenoStratHideBtn,
               {hide("pheno_strat_details")})
  
  # Enable hiding and showing of gene expression table (Side panel) # [TESTING]
  observeEvent(input$geneExpressionShowBtn,
               {show("geneExpressionSumTab")})
  observeEvent(input$geneExpressionHideBtn,
               {hide("geneExpressionSumTab")})
  
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
      line_2 <- "Select PGS below to explore metrics."
      HTML(paste(line_1,line_2,sep='<br/>'))
    } else if (no_avail_pgs > 1) { 
      line_1 <- paste0("There are ", no_avail_pgs, " available for phenotype selected.")
      line_2 <- "Choose a PGS below to explore metrics."
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
                                     rownames = TRUE)
  output$phenoSentence <- renderUI(req(try(phenoStats()$SENTENCE)))
  
  # PGS Data 
  pgsStats <- reactive({getPGSStratStats(x=input$pheno_selector,
                                         y=input$pgs_selector)})
  
  #output$pgsStratSummary <- DT::renderDataTable(req(try(pgsStats()$TABLE))) #[!] Try regular table for now
  output$pgsStratSummary <- renderTable(req(try(pgsStats()$TABLE)))
  output$pgsSentence <- renderUI(req(try(pgsStats()$SENTENCE)))
  output$pgsStratPlot <- renderPlot(req(try(pgsStats()$PLOT)))
  
  # Gene Expression metadata [!]
  output$geneExpressionPlot <- renderPlot(req(try(summarizeGeneExp(chr=input$table_selector, 
                                                                   input.gene.name=input$gene_selector)[['PLOT']])))
  output$geneExpressionSumTab <- renderTable(req(try(summarizeGeneExp(chr=input$table_selector, 
                                                               input.gene.name=input$gene_selector)[['SUM_TABLE']])), 
                                             digits = 4, rownames = TRUE)
  output$GTEx <-renderUI(a(href=paste0('https://urldefense.com/v3/__https://www.gtexportal.org/home/gene/__;!!HrbR-XT-OQ!XTLJTajDBjHW_UuUobS4sK4S3CRNVMaFJ-BCZ9r5D5J6N5W9fHD-aMs9XH_slXMvOf3NXHuQqATL1azpKA$ ', 
                                       stringr::str_split(input$gene_selector,"[.]")[[1]][1]),
                           "GTEx",target="_blank"))
  output$Ensembl <- renderUI(a(href=paste0('https://urldefense.com/v3/__http://grch37.ensembl.org/Human/Search/Results?q=__;!!HrbR-XT-OQ!XTLJTajDBjHW_UuUobS4sK4S3CRNVMaFJ-BCZ9r5D5J6N5W9fHD-aMs9XH_slXMvOf3NXHuQqARtvQyDNw$ ', 
                                           stringr::str_split(input$gene_selector,"[.]")[[1]][1],
                                           ';site=ensembl;facet_species=Human'), 
                               "Ensembl", target="_blank"))
  
  
  # PICS comparison tables and dbSNP URLs
  PICSResult <- reactive({
    summarizePICSResult(chr=input$table_selector,
                        input.gene.name=input$gene_selector,
                        ps=input$potential_set)})
  
  output$PICS2TopVSStable <- renderTable(req(try(PICSResult()$TOP_VS_STABLE)), rownames = TRUE)
  
  output$TopStableMatch <- renderUI({
    req(try(top_variant <- PICSResult()$TOP_VS_STABLE["rsID", "Top"]))
    req(try(stable_variant <- PICSResult()$TOP_VS_STABLE["rsID", "Stable"]))
    if (top_variant == stable_variant){
      h4("Top and stable variants", strong("match"), paste0("(", top_variant, ")"))
    } else {
      h4("Top and stable variants are", strong("different"), paste0("(", top_variant, " vs ", stable_variant, ")"))
    }
  })
  
  # Create a paragraph-sentence summarizing user-selected gene
  summarizedGeneExp <- reactive({
    req(input$table_selector, input$gene_selector)
    summarizeGeneExp(chr=input$table_selector, input.gene.name=input$gene_selector)
  })
  
  output$GeneInfoSummary <- renderUI({
    gene_symbol <- summarizedGeneExp()[['GENE_SYMBOL']]
    req(try(ensembl_gene <- summarizedGeneExp()[['ENSEMBL_GENE']]))
    biotype <- summarizedGeneExp()[['BIOTYPE']]
    if (is.na(gene_symbol)) {
      p("The selected gene is ",
        strong(ensembl_gene), 
        ", with no gene name available on Biomart. Click on links below for more information.")
    } else {
      p("The selected gene is ",
        strong(ensembl_gene),
        " (Gene Name: ",
        strong(gene_symbol),
        "). Its biotype classification is: ",
        strong(biotype),
        ". Click on links below for more information.")
    }
  })
  
  output$PICS2StableStats <- renderTable(req(try(PICSResult()$STABLE_STATS)), 
                                         rownames = TRUE,
                                         digits = 5)
  output$PICS2TopDBSNP <- renderUI(a(href=req(try(PICSResult()$TOP_DBSNP)), 
                                     "Top Variant", 
                                     target="_blank"))
  output$PICS2StableDBSNP <- renderUI(a(href=req(try(PICSResult()$STAB_DBSNP)), 
                                        "Stable Variant", 
                                        target="_blank"))
  
  # PICS functional annotation 
  output$FuncAnnotTable <- renderTable(req(try(getFuncAnnotTable(chr=input$table_selector,
                                                                 input.gene.name=input$gene_selector, 
                                                                 ps=input$potential_set, 
                                                                 func.class=input$func_annot_class,
                                                                 func.annot=input$func_annot_selector))), 
                                       rownames = TRUE,
                                       digits = 5)
  output$FuncAnnotDescription <- renderText(func_annot_descriptions[[input$func_annot_class]][input$func_annot_selector])
  output$FuncAnnotReference <- renderUI(a(href=func_annot_links[[input$func_annot_class]][input$func_annot_selector], 
                                          "Reference", 
                                          target="_blank"))
  
  # Genome-wide pairwise plots of functional annotation
  output$funcDensityPlot <- renderPlot({
    # If functional annotation input has not refreshed, pick the first annotation in the correct list
    func.class=input$gw_func_annot_class
    func.annot=input$gw_func_annot_selector
    ps=input$gw_ps
    func_annot_list <- readLines(paste0('all-func-annots/', func.class, '_annot_names.txt'))
    req(func.annot %in% func_annot_list)
    # Draw density plot
    getDensityPlot(func.class, func.annot, ps)})
  
  output$funcAnnotPlot <- renderPlot({
    # If functional annotation input has not refreshed, pick the first annotation in the correct list
    func.class=input$gw_func_annot_class
    func.annot=input$gw_func_annot_selector
    ps=input$gw_ps
    func_annot_list <- readLines(paste0('all-func-annots/', func.class, '_annot_names.txt'))
    req(func.annot %in% func_annot_list)
    # Draw pairwise plot
    getPairwisePlot(func.class, func.annot, ps)})
  
  output$funcAnnotDescription <- renderText(func_annot_descriptions[[input$gw_func_annot_class]][input$gw_func_annot_selector])
  output$funcAnnotReference <- renderUI(a(href=func_annot_links[[input$gw_func_annot_class]][input$gw_func_annot_selector], 
                                          "Reference", 
                                          target="_blank"))
  
  # Explain Enformer log-root scale
  output$EnformerFormula <- renderUI({
    withMathJax(
      sprintf("Note: For Enformer annotations, original perturbations are transformed by the function \\(f(x) = \\log(x^\\frac{1}{4} + 1)\\) (log-root scale) to improve visualization.")
    )
  })
  
  # Explain enrichment scores for Ensembl (Single gene tab)
  output$EnsemblNote <- renderUI({
    req(try(annot <- input$func_annot_selector))
    if (annot %in% c('CTCF Binding Site', 
                     'Enhancer', 'Open chromatin', 
                     'Promoter Flanking Region', 
                     'Promoter', 'TF binding site')) {
      p("Note: Enrichments are calculated by allocating a score of 1 or 0 for each input variant based on its inclusion in a given annotated region, and then taking the ratio of the fine-mapped variant (top or stable) and the mean score.")
    }
  })
  
  # Explain enrichment scores for Ensembl (Genome-wide tab)
  output$GWEnsemblNote <- renderUI({
    req(try(annot <- input$gw_func_annot_selector))
    if (annot %in% c('CTCF Binding Site', 
                     'Enhancer', 'Open chromatin', 
                     'Promoter Flanking Region', 
                     'Promoter', 'TF binding site')) {
      p("Note: Enrichments are calculated by allocating a score of 1 or 0 for each input variant based on its inclusion in a given annotated region, and then taking the ratio of the fine-mapped variant (top or stable) and the mean score.")
    }
  })
  
  # Graphic summarizing PICS2
  output$PICS2Graphic <- renderImage({
    list(
      src = file.path("www/pics2_tikz.jpg"),
      #contentType = "image/jpeg",
      width = 3085/4,
      height = 1568/4
    )
  }, deleteFile = FALSE)
}



# Run the application 
shinyApp(ui = ui, 
         server = server)