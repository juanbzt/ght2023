#################################################################################
#                                 APP
################################################################################

# This script allows you to merge the differential expression data from all datasets analysed into a unique dataframe
# and performes a search in it.

# Users must have the following packages installed. To install them, use the install.packages("") function.
library(stringr)
library(plyr)
library(shiny)
library(reactable)
library(ggplot2)
library(shinyWidgets)
library(reactablefmtr)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(survival)
library(survminer)

# IMPORT THE RESULTS OF THE DIFFERENTIAL EXPRESSION ANALYSIS PERFORMED ON DIFFERENT DATASETS (THESE FILES ARE ON THE "DATA" FOLDER)
GSE3933_tumor_vs_no_tumor <- read.table("data/GSE3933_tumor_vs_no_tumor.txt", header = TRUE)
GSE6956_tumor_vs_adjacent <- read.table("data/GSE6956_tumor_vs_adjacent.txt", header = TRUE)
GSE16560_Lethal_vs_Indolent <- read.table("data/GSE16560_Lethal_vs_Indolent.txt", header = TRUE)
GSE25136_Recurrent_vs_Non_recurrent <- read.table("data/GSE25136_Recurrent_vs_Non_recurrent.txt", header = TRUE)
GSE34312_tumor_vs_adjacent <- read.table("data/GSE34312_tumor_vs_adjacent.txt", header = TRUE)
GSE35988_metastatic_vs_benign <- read.table("data/GSE35988_metastatic_vs_benign.txt", header = TRUE)
GSE35988_metastatic_vs_PCa <- read.table("data/GSE35988_metastatic_vs_PCa.txt", header = TRUE)
GSE35988_PCa_vs_Benign <- read.table("data/GSE35988_PCa_vs_Benign.txt", header = TRUE)
GSE46602_PCa_vs_Benign <- read.table("data/GSE46602_PCa_vs_Benign.txt", header = TRUE)

# LOAD SURVIVAL AND EXPRESSION DATA FROM SBONER DATASET (THESE FILES ARE ON THE "DATA" FOLDER)
Sboner_survival <- read.table("data/Sboner_survival.txt", header = TRUE)
Sboner_expression <- read.table("data/Sboner_expression.txt", header = TRUE)

##### CREATE A DATAFRAME CONTAINING INFORMATION FROM ALL THE DATASETS
# Create a dataframe by binding the results from the output of differential expression analyses
{
  df <- data.frame()
  df <- rbind.data.frame(df, GSE3933_tumor_vs_no_tumor, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE6956_tumor_vs_adjacent, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE16560_Lethal_vs_Indolent, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE25136_Recurrent_vs_Non_recurrent, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE34312_tumor_vs_adjacent, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE35988_metastatic_vs_benign, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE35988_metastatic_vs_PCa, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE35988_PCa_vs_Benign, stringsAsFactors = F)
  df <- rbind.data.frame(df, GSE46602_PCa_vs_Benign, stringsAsFactors = F)
  df$Gene.title = NULL # Remove the Gene.Title column for table's aesthetic
  df <- df %>% 
    mutate_if(is.numeric, round, digits = 3) # Format of numeric data
  }


##### CREATE FUNCTIONS TO ADD VISUAL COMPONENTS TO RESULT TABLES
# Create function to visualize p-value significance in reactables
PValue_significance <- function(value){
  if (!is.numeric(value)) {
    return(value)
  } else if (is.na(value)) {
    return(value)
  } else if (value <= 0.05) {
    tagList(value,
            div(style = list(display = "inline-block", marginLeft = "8px", color = "#059305"), icon("check-circle", lib = "font-awesome")))
  } else {
    tagList(value,
            div(style = list(display = "inline-block", marginLeft = "8px", color = "#d3d3d3"), icon("minus-circle", lib = "font-awesome")))
  } }


# Create function to visualize gene expression bars in reactables
bar_chart_pos_neg <- function(value, max_value = 1, height = "16px",
                              neg_fill = "#005ab5", pos_fill = "#dc3220") {
  neg_chart <- div(style = list(flex = "1 1 0"))
  pos_chart <- div(style = list(flex = "1 1 0"))
  width <- paste0(abs(value / 5) * 100, "%")
  if (!is.numeric(value)) {
    return(value)
  } else if (is.na(value)) {
    return(value)
  } else if (value <= 0) {
    bar <- div(style = list(marginLeft = "8px", background = neg_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"), value, bar)
    neg_chart <- tagAppendChild(neg_chart, chart)
  } else {
    bar <- div(style = list(marginRight = "8px", background = pos_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center"), bar, value)
    pos_chart <- tagAppendChild(pos_chart, chart)
  }
  
  div(style = list(display = "flex"), neg_chart, pos_chart)
}

###############################
    # USER INTERFACE #
###############################

# Define UI for application
ui <- fluidPage(
  
   # App title
   fluidRow(
     style="margin-left: auto; margin-right: auto;",
     class="container",
     column(12,
         HTML('
          <div style="display: flex; justify-content: center; align-items: center;">
              <div>
                <img src="ADN.png" style="height: 150px; object-fit: contain; margin-right: 20px; margin-top: 15px;"/>
              </div>  
              <div>
                  <h1>Gene Hunter</h1>
              </div>
              <div>
                <img src="ADN.png" style="height: 150px; object-fit: contain; margin-left: 20px; margin-top: 15px;"/>
              </div>
            </div>
        ')
      )
    ),
   
   
   # CONFIGURE SIDE PANEL
   sidebarLayout(column(2,
     wellPanel( h5(strong("Enter your gene!")),
                textInput(inputId = "wanted_gene", label = NULL, placeholder = "Gene symbol"), 
               actionButton(inputId = "go", label = "Go!", icon("dna", lib = "font-awesome")),
               br(),
               br(),
               h6("Select the datasets you wish to include in the search"),
               pickerInput("Datasets",
                           choices = c("Sboner", "Ashida", "Sun", "Wallace", "Lapointe", "Mortensen", "Grasso (1)", "Grasso (2)", "Grasso (3)"),
                           options = list('actions-box' = TRUE),
                           selected = c("Sboner", "Ashida", "Sun", "Wallace", "Lapointe", "Mortensen", "Grasso (1)", "Grasso (2)", "Grasso (3)"),
                           multiple = TRUE,
                           inline = FALSE)),
     hr(),
     actionLink ("labweb", "Laboratory of Cancer & Inflammation - University of Buenos Aires", icon ("flask", lib = "font-awesome"), 
                   onclick = "window.open ('http://cancer.qb.fcen.uba.ar/', '_blank')")),
     mainPanel(width = 10,
       navbarPage("Gene Hunter",
     
      # CONFIGURE EACH TAB
       # Principal Tab
       tabPanel(icon("home", lib = "font-awesome"),
                fluidRow(column( 7,
                                 h3("Welcome!"), 
                                 helpText("This is the main page of the", strong("Gene Hunter app of the Laboratory of Cancer & Inflammation"), "(University of Buenos Aires). It is a tool that 
                                              allows you to look for a gene in several datasets at once.",
                                          br(),
                                          downloadButton("Bottomright", outputId = "Download", label = "Download table")
                          )),
                          column(5,
                                 h6 (style="color:grey", "Note: Grasso's dataset presents 3 different comparisons. You can choose which comparisons to include in the analysis:"),
                                 h6 (style="color:grey", "(1) PCa vs. Benign"),
                                 h6 (style="color:grey", "(2) Metastatic vs. Benign"),
                                 h6 (style="color:grey", "(3) Metastatic vs. PCa")
                                 )),
                 hr(),
                 
                 # Output format
                 fluidRow(
                 style="margin-left: auto; margin-right: auto;",
                 column(12,
                 br(),
                 reactableOutput("results"))),
                hr()
                 ),
       
        
        # Strict Search tab
        tabPanel("Strict Search",
                 fluidRow(column( 7,
                                  h3("Strict Search"),
                                  helpText("This section allows you to perform a stricter search: when entering the name of a gene,
                                           the finder will show only results that appear on the datasets",
                                           strong("exactly"),"as the entered pattern.",
                                           p("When is it helpful? It is very helpful when the gene name you are looking for appears in other genes' names,", 
                                             em(" for example: AR, ERG, and other short names", "."))),
                                  downloadButton(outputId = "Download2", label = "Download strict table")),
                          column(5,
                                 h6 (style="color:grey", "Note: Grasso's dataset presents 3 different comparisons. You can choose which comparisons to include in the analysis:"),
                                 h6 (style="color:grey", "(1) PCa vs. Benign"),
                                 h6 (style="color:grey", "(2) Metastatic vs. Benign"),
                                 h6 (style="color:grey", "(3) Metastatic vs. PCa")
                          )),
                          hr(),
                 
                 # Output Format
                 fluidRow(
                   style="margin-left: auto; margin-right: auto;",
                   column(12,
                          br(),
                          reactableOutput("results_strict"))),
                 hr()
                 ),
                 
        
       # Plots tab: Using a navbarMenu creates a dropdown menu
       navbarMenu ("Plots",
                   tabPanel("LogFC vs Author",
                            fluidRow(column(7,
                             h3("Plot: LogFC vs Author"),

                             helpText ("This section allows you to generate different plots from the data available in this app. 
                                       In this case, you are plotting 'Author' on the x axis and 'LogFC' on the y axis. The P-value for each probe in each comparison is shown."),
                             downloadButton(outputId = "downloadplot", label = "Download plot")),
                             column(5, 
                                    h6 (style="color:grey", "Note: Grasso's dataset presents 3 different comparisons. You can choose which comparisons to include in the analysis:"),
                                    h6 (style="color:grey", "(1) PCa vs. Benign"),
                                    h6 (style="color:grey", "(2) Metastatic vs. Benign"),
                                    h6 (style="color:grey", "(3) Metastatic vs. PCa"))),
                            hr (),
                           
                             # Output Format
                             fluidRow(
                               style="margin-left: auto; margin-right: auto;",
                               column(12,
                                      br(),
                                      plotOutput("results_plot"))),
                             hr()),
                   
                   tabPanel ("Volcano Plot",
                            fluidRow(column(7,
                                   h3("Volcano Plot"),
                                   helpText ("This section allows you to generate different plots from the data available 
                                             in this app. In this case, you are plotting 'logFC' on the x axis and 'P Value' 
                                             on the y axis."),
                                   downloadButton(outputId = "downloadplot2", label = "Download plot")),
                                      column(5, 
                                            h6 (style="color:grey", "Note: Grasso's dataset presents 3 different comparisons. You can choose which comparisons to include in the analysis:"),
                                            h6 (style="color:grey", "(1) PCa vs. Benign"),
                                            h6 (style="color:grey", "(2) Metastatic vs. Benign"),
                                            h6 (style="color:grey", "(3) Metastatic vs. PCa"),
                             hr (),
                             helpText ("Please define the cutoff values you want"),
                             textInput(inputId = "wanted_logfc", label = "Log Fold Change cutoff:",
                                       placeholder = "Choose a cutoff value for Log Fold Change", value = 1.5),
                             textInput(inputId = "wanted_padj", label = "Adj.P.Value cutoff:",
                                       placeholder = "Choose a cutoff value for P-Adj", value = 0.05))),
                            hr(),
                             
                             # Output Format
                             fluidRow(
                               style="margin-left: auto; margin-right: auto;",
                               column(12,
                                      br(),
                                      plotOutput("results_plot2")))),
                   
                   tabPanel ("Survival",
                             fluidRow(column(7, 
                                             h3("Kaplan-Meier plot"), 
                                             helpText("This section allows you to generate different plots from the data available
                                             in this app. In this case, you are plotting Kaplan-Meier plots using data from the", 
                                             strong("Sboner et al."), "dataset."),
                                             downloadButton(outputId = "downloadplot3", label = "Download plot"),
                                             downloadButton(outputId = "downloadSurvival", label = "Download Survival table"))),
                             hr (),
          
                             # Output Format
                             fluidRow(
                               style="margin-left: auto; margin-right: auto;",
                               column(12,
                                      br(),
                                      plotOutput("results_plot3"),
                                      reactableOutput("Survival"))),
                             hr()),

                  tabPanel ("Boxplot: Probes per Comparison",
                            fluidRow(column(7,
                                             h3("Boxplot: Probes per Comparison"),
                                             helpText("This section allows you to generate different plots from the data available 
                                             in this app. In this case, you are plotting the probes available for each comparison."),
                                              downloadButton(outputId = "downloadplot4", label = "Download plot")),
                            column(5, 
                                   h6 (style="color:grey", "Note: Grasso's dataset presents 3 different comparisons. You can choose which comparisons to include in the analysis:"),
                                   h6 (style="color:grey", "(1) PCa vs. Benign"),
                                   h6 (style="color:grey", "(2) Metastatic vs. Benign"),
                                   h6 (style="color:grey", "(3) Metastatic vs. PCa"))),
                            hr(),
                            
                            # Output Format
                             fluidRow (
                               style="margin-left: auto; margin-right: auto;",
                               column(12,
                                      br(),
                                      plotOutput("results_plot4"))),
                             hr())),

              
        # DATASETS INFORMATION
        tabPanel("Datasets",
                 h3("Datasets"),
                 helpText("In this section you can find information about the datasets uploaded in our Gene Finder, 
                          as well as the links to their respective papers."),
                 hr(),
                 
                 # First row of 3 datasets
                 fluidRow(
                   style="margin-left: auto; margin-right: auto;",
                   #SBoner dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Sboner"),
                               style = "color: midnightblue;"),
                            h5( a("GSE16560", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16560")),
                            h5("Compares: Lethal vs Indolent",
                               p("Samples: 281"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-3-8"))
                        )
                   ),
                   #Ashida dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Ashida"),
                               style = "color: midnightblue;"),
                            h5( a("GSE34312", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34312")),
                            h5("Compares: Tumor vs Adjacent",
                               p("Samples: 20"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://www.ncbi.nlm.nih.gov/pubmed/22275508"))
                          )
                   ),
                   #Sun dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Sun"),
                               style = "color: midnightblue;"),
                            h5( a("GSE25136", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136")),
                            h5("Compares: Recurrent vs Non recurrent",
                               p("Samples: 79"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"),
                                 href ="https://www.ncbi.nlm.nih.gov/pubmed/19343730"))
                            )
                          )
                   ),
                 
                 # Second row of 3 datasets
                 fluidRow(
                   style="margin-left: auto; margin-right: auto;",
                   #Wallace dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Wallace"),
                               style = "color: midnightblue;"),
                            h5( a("GSE6956", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956")),
                            h5("Compares: Tumor vs Adjacent",
                               p("Samples: 89"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://www.ncbi.nlm.nih.gov/pubmed/18245496"))
                          )
                   ),
                   #Lapointe dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Lapointe"),
                               style = "color: midnightblue;"),
                            h5( a("GSE3933", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3933")),
                            h5("Compares: Tumor vs No tumor",
                               p("Samples: 112"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC321763/"))
                            )
                          ),
                   #Grasso dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Grasso"),
                               style = "color: midnightblue;"),
                            h5( a("GSE35988", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35988")),
                            h5("Compares: Tumor, PCa y Benign     ",
                               p("Samples: 244"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://www.ncbi.nlm.nih.gov/pubmed/22722839"))
                          )
                   )
                   ),
                 
                 # Third row 
                 fluidRow(
                   style="margin-left: auto; margin-right: auto;",
                   #Mortensen dataset
                   column(4,
                          style="display: flex; justify-content: center; align-items: center",
                          wellPanel(
                            h4(strong("Mortensen"),
                               style = "color: midnightblue;"),
                            h5( a("GSE46602", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46602")),
                            h5("Compares: PCa vs Benign",
                               p("Samples: 50"),
                               a(icon("file-alt", 
                                      lib = "font-awesome", "fa-3x"), 
                                 href ="https://www.ncbi.nlm.nih.gov/pubmed/26522007"))
                            )
                          )
                   )
                 )
       )
       ),
   position = "left",
   fluid = TRUE)
)


#####################
     # SERVER #
#####################

server <- function(input, output, session) {

##Subset dataframe by selected authors
     subsetted <- reactive(
       df[df$Author %in% input$Datasets, ])
  
  
##### GENERAL SEARCH
  GoSearch <- eventReactive(input$go, {
              input$wanted_gene
             })
  
  # Set the search result in a reactive function so it can be called later from the table and from the download button.
  YourGene <- reactive({

    # Search for a gene in the dataframe. This search is approximated.
    search <- subsetted()[which(agrepl(GoSearch(), subsetted()$Gene.symbol, 
                          ignore.case = TRUE, 
                          fixed = TRUE,
                          max.distance = list(cost = 0, all = 0, insertions= 0, deletions = 0, substitutions = 0))),]
  })
  
  # RESULTS TABLE
  output$results <- renderReactable({
   
    # Create a table using reactable
    reactable(YourGene(), 
              searchable = F, 
              striped = TRUE, 
              groupBy = "Gene.symbol", 
              fullWidth = TRUE, 
              highlight = TRUE,
              resizable = TRUE,
              filterable = TRUE,
              columns = list(
                adj.P.Val = colDef(cell= function(value){
                  PValue_significance(value)}),
                P.Value = colDef(cell= function(value){
                  PValue_significance(value)}),
                logFC = colDef(cell = function(value) {
                  bar_chart_pos_neg(value)},
                  align = "center",
                  minWidth = 400)))
  }) 
  
  # DOWNLOAD
  output$Download <- downloadHandler(filename = function() {paste(GoSearch(), ".txt", sep = "")}, 
                                     content = function(file) {write.table(YourGene(), file, row.names=T, sep="\t")})
  
  
#### STRICT SEARCH

  ## Set the search result in a reactive function so it can be called later from the table/plots and from the download button
  YourGene2 <- reactive({  
    # Search for a gene in the dataframe. This is a strict search.
    search_2 <- subsetted()[which(subsetted()[["Gene.symbol"]] == GoSearch()),]
  })
  
  
   # RESULT TABLE
  output$results_strict <- renderReactable({
    # Create a table using reactable
    reactable(YourGene2(), 
              searchable = F, 
              striped = TRUE, 
              groupBy = "Comparison", 
              fullWidth = TRUE, 
              highlight = TRUE,
              resizable = TRUE,
              filterable = TRUE,
              columns = list(
                adj.P.Val = colDef(cell= function(value){
                  PValue_significance(value)}),
                P.Value = colDef(cell= function(value){
                  PValue_significance(value)}),
              logFC = colDef(cell = function(value) {
                bar_chart_pos_neg(value)},
                             align = "center",
                             minWidth = 400))
    )
  }) 

  # DOWNLOAD
  output$Download2 <- downloadHandler(filename = function() {paste(GoSearch(), " - Strict.txt", sep = "")}, 
                                     content = function(file) {write.table(YourGene2(), file, row.names=T, sep="\t")})
  
  
#### PLOTS
  #AUTHOR VS LOG FC PLOT
  # Create a plot inside a reactive consumer using ggplot2 package.
  YourPlot = reactive({
    ggplot (YourGene2(), aes (x=Author, y=logFC, fill=ID)) + 
      geom_col(
        mapping= NULL,
        data = NULL,
        position = position_dodge2(preserve = 'single'),
        width = 1,
        na.rm = FALSE,
        show.legend = NA,
        inherit.aes = TRUE) +
      ggtitle (GoSearch()) + 
      theme_bw() +
      theme(plot.title = element_text(face="bold", size = 15, hjust= 0.5),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 13),
            axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_text(data = YourGene2(), 
                aes(label = adj.P.Val, 
                    y = (abs(logFC) + 0.1)*sign(logFC),
                    x=Author),
                position = position_dodge(width = 1), 
                size=3)
      })
  
  
  #VOLCANO PLOT
  # Create a plot inside a reactive consumer using EnhancedVolcano package.
   YourPlot2 <- reactive({
    EnhancedVolcano(toptable = YourGene2(),
                    lab = rownames(YourGene2()),
                    x='logFC',
                    y='adj.P.Val',
                    title= GoSearch(),
                    subtitle = 15,
                    caption = paste("Log Fold Change cutoff = ", input$wanted_logfc, "; Adj.P.Value cutoff = ", input$wanted_padj),
                    xlab=bquote(~Log~'Fold Change'),
                    ylab=bquote(~-Log~"(P-Adj)"),
                    selectLab = "",
                    # Set the cut-off values chosen by the user
                    pCutoff= as.numeric(input$wanted_padj),
                    FCcutoff= as.numeric(input$wanted_logfc),
                    pointSize=3.0,
                    labSize=3.0,
                    labCol = "black",
                    boxedLabels = TRUE,
                    labhjust = 1.0,
                    legendLabels =c("NS","LogFC","P-Adj","P-Adj & LogFC"),
                    legendPosition="bottom") +
       theme_bw() +
       theme(plot.title = element_text(face="bold", size = 15, hjust= 0.5),
             axis.title = element_text(size = 13),
             axis.text = element_text(size = 13),
             legend.title = element_blank())
  })
  

   # SURVIVAL
   ## Set the search result in a reactive function so it can be called later from the table and from the download button
   YourGene3 <- reactive({  
     # Search for a gene in the dataframe. This is a strict search.
     search_3 <- df[which(df[["Gene.symbol"]] == GoSearch()),]
   })
   # Look for the ID associated to the genes name
   ID_interest <- reactive({
     unlist(subset(YourGene3(), Author == "Sboner", select = c("ID")))
   })
   
   #Create a table with survival data for the gene of interest. In this table there is a column that assigns each patient a "Low" or 
   # "High" status, acording to the expression of the gene of interest.
  Sur_table <- reactive({
    Sboner_survival %>%
      mutate(Gene_expression = Sboner_expression[,ID_interest()],
             Levels = case_when(Gene_expression <= quantile(Gene_expression, 0.25) ~ "LOW",
                                Gene_expression >= quantile(Gene_expression, 0.75) ~ "HIGH")) %>%
      mutate(Levels = factor(Levels, ordered = TRUE, levels = c("LOW", "HIGH")))
    })
    
   Sur_analysis <- reactive({
     table <- Sur_table()
     fit <- survfit(Surv(table$fup_month, table$Survival) ~ table$Levels, data = Sur_table())
     fit
   })
   
   
   #KAPLAN-MEIER PLOT
   # Create a plot inside a reactive consumer using survminer package.
   YourPlot3 <- reactive({
     Survival <- Sur_analysis()
     plot(Survival, 
          col = c("#E7B800", "#2E9FDF"), 
          conf.int = F, 
          mark.time = TRUE, 
          lwd = 2, 
          xmax = 200,
          ylab = "Survival",
          xlab = "Time (months)",
          main = paste("Survival:", GoSearch()))
     legend(0, 0.2, c("LOW", "HIGH"), lty = 1, col = c("#E7B800", "#2E9FDF"))
     })
  
   
   #BOXPLOT: MEAN OF PROBES PER COMPARISON
   # Create a plot inside a reactive consumer using ggplot2 package.
   YourPlot4 <- reactive({
     ggplot(YourGene2(), aes(y=logFC, x=Comparison)) +
       geom_boxplot(alpha = 0.2, notch=FALSE, fill="aliceblue", color= "azure4") +
       geom_point(aes(color=ID), shape=18, size=4, position = "jitter") +
       ggtitle (paste(GoSearch(), ": Probes per comparison")) +
       theme_bw() +
       theme(plot.title = element_text(face="bold", size = 15, hjust= 0.5),
             axis.title = element_text(size = 13),
             axis.text = element_text(size = 13),
             axis.text.x = element_text(angle = 45, hjust = 1))
   })   


  # PLOTS OUTPUTS
   output$results_plot <- renderPlot({
    YourPlot()})
   
   output$results_plot2 <- renderPlot({
     YourPlot2()})

   output$results_plot3 <- renderPlot({
     YourPlot3()})
   
   # This table contains the survival data for the gene of interest
   output$Survival <- renderReactable({
     reactable(Sur_table(),
               striped = TRUE,
               searchable = F, 
               fullWidth = TRUE, 
               highlight = TRUE,
               resizable = TRUE,
               filterable = F,
               columns = list(
                 Age = colDef(align = "center"),
                 fup_month = colDef(name = "Follow up months",
                                    align = "center"),
                 gleason = colDef(name = "Gleason Score",
                                  align = "center"),
                 gleason_major = colDef(name = "Major Gleason",
                                        align = "center"),
                 gleason_minor = colDef(name = "Minor Gleason",
                                        align = "center"),
                 Survival = colDef(align = "center"),
                 Gene_expression = colDef(name = "Gene expression",
                                          align = "center"),
                 Levels = colDef(align = "center")
               ))  })

   output$results_plot4 <- renderPlot({
     YourPlot4()})
   
   
  # DOWNLOAD PLOTS
  output$downloadplot <- downloadHandler(filename = function() {paste(GoSearch(), " - LogFC vs Author.png", sep = "")}, 
                                         content = function(file) {
                                           device <- function (..., width, height){
                                             grDevices::png (..., width = width, height = height,
                                                             res = 300, units = "in")
                                           }
                                           ggsave (file, plot=YourPlot() , device=device)
                                         })
  
  output$downloadplot2 <- downloadHandler(filename = function() {paste(GoSearch(), " - Volcano Plot.png", sep = "")}, 
                                         content = function(file) {
                                           device <- function (..., width, height){
                                             grDevices::png (..., width = width, height = height,
                                                             res = 300, units = "in")
                                           }
                                           ggsave (file, plot=YourPlot2() , device=device)
                                         })
  
  output$downloadplot3 <- downloadHandler(filename = function() {paste(GoSearch(), " - Kaplan-Meier Plot.png", sep = "")}, 
                                          content = function(file) {
                                            png(file)
                                            plot(Sur_analysis(), 
                                                 col = c("#E7B800", "#2E9FDF"), 
                                                 conf.int = F, 
                                                 mark.time = TRUE, 
                                                 lwd = 2, 
                                                 xmax = 200,
                                                 ylab = "Survival",
                                                 xlab = "Time (months)",
                                                 main = paste("Survival:", GoSearch()))
                                            legend("bottomleft", c("Low", "HIGH"), lty = 1, col = c("#E7B800", "#2E9FDF"))
                                            dev.off()
                                          })
  
  output$downloadSurvival <- downloadHandler(filename = function() {paste(GoSearch(), " - Survival.txt", sep = "")}, 
                                             content = function(file) {write.table(Sur_table(), file, row.names = TRUE, sep="\t")})

  output$downloadplot4 <- downloadHandler(filename = function() {paste(GoSearch(), " - Probes.png", sep = "")}, 
                                          content = function(file) {
                                            device <- function (..., width, height){
                                              grDevices::png (..., width = width, height = height,
                                                              res = 300, units = "in")
                                            }
                                            ggsave (file, plot=YourPlot4() , device=device)
                                          })
  } 


#### Run the application ####
shinyApp(ui = ui, server = server)
