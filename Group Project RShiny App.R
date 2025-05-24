# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

# Check if shiny is installed; if not, install it
if (!require("shiny")) {
  install.packages("shiny")
}

# Check if ggplot2 is installed; if not, install it
if (!require("ggplot2")) {
  install.packages("ggplot2")
}

# Check if gridExtra is installed; if not, install it
if (!require("gridExtra")) {
  install.packages("gridExtra")
}

# Check if treeio is installed; if not, install it
if (!require("treeio")) {
  install.packages("treeio")
}

# Check if ggtree is installed; if not, install it
if (!require("ggtree")) {
  install.packages("ggtree")
}

# Check if bslib is installed; if not, install it
if (!require("bslib")) {
  install.packages("bslib")
}

# Check if DT is installed; if not, install it
if (!require("DT")) {
  install.packages("DT")
}

library(shiny)
library(ggplot2)
library(treeio)
library(ggtree)
library(gridExtra)
library(bslib)
library(DT)

# First, set your working directory to source file location.

# Define UI
ui <- page_fluid(
  titlePanel("Evolution and Likelihood of Antibiotic Resistance"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("This application is best viewed in full-screen."),
      br(),
      br(),
      helpText("If you have run this application before, please remove any previous output files from your working directory folder."),
      br(),
      br(),
      tooltip(fileInput("genomeFile", "Upload your file of genomes:", accept = c(".txt", ".fasta", ".fa", ".fna")),"FASTA format required. File limit 2 GB.",placement="right"),
      tooltip(fileInput("geneFile","Upload your antibiotic resistance gene sequence file:",accept = c(".txt", ".fasta", ".fa", ".fna")),"FASTA format required. File limit 2 KB.",placement="right"),
      tooltip(sliderInput("percentIdentity", "Select the percent identity threshold:", min = 0, max = 1, value = 0.25),"Percent identity is the nucleotide similarity between two regions. Only genomes with a gene region with a nucleotide similarity at or above this threshold, compared to the input gene sequence, will be included in the calculations.", placement = "right"),
      tooltip(sliderInput("preservationThreshold","Select the nucleotide  preservation threshold:", min = 0, max = 1, value = 0.25), "Nucleotide preservation is the percentage of nucleotides in an identified gene region that match the input gene sequence. Genomes with a nucleotide preservation value below this threshold will be considered more likely to have antibiotic resistance, and vice versa.", placement="right"),
      helpText("Note: if no tree is created based off your initial parameters, try decreasing the percent identity metric."),
      br(),
      br(),
      actionButton("runGoCode", "Run Genome Processing")
    ),
    
    mainPanel(tabsetPanel(
      #create a nucleotide tree tab
      tabPanel(
        title = "Nucleotide Evolutionary Tree",
        plotOutput("nucleotidePlot")
      ),
      #create a protein tree tab
      tabPanel(
        title = "Protein Evolutionary Tree",
        plotOutput("proteinPlot")
      ),
      #create a likelihood table tab
      tabPanel(
        title = "Antibiotic Resistance Likelihood Data",
        dataTableOutput("likelihoodPlot")
      ),
      #create a preservation table tab
      tabPanel(
        title = "Protein Preservation Data",
        dataTableOutput("preservationPlot")
      )
    )
    )
    
  )
)

server <- function(input, output) {
  # Increase max file limit.
  options(shiny.maxRequestSize=2*1024*1024^2)
  
  observeEvent(input$runGoCode, {
    # Determine the file path or URL to use
    genomeFilePath <- NULL
    geneFilePath <- NULL
    percentIdentityValue <- NULL
    
    # Require all inputs
    if (!is.null(input$genomeFile)) {
      genomeFilePath <- input$genomeFile$datapath
    } 
    
    if (!is.null(input$geneFile)) {
      geneFilePath <- input$geneFile$datapath
    }
    
    if (!is.null(input$percentIdentity)) {
      percentIdentityValue <- input$percentIdentity
    }
    
    if (!is.null(input$preservationThreshold)) {
      preservationThresholdValue <- input$preservationThreshold
    }
    
    req(genomeFilePath)  # Ensure a genome file path is provided
    req(geneFilePath) # Ensure a gene file path is provided
    req(percentIdentityValue) # Ensure a percent identity is provided
    req(preservationThresholdValue) # Ensure a preservation threshold is provided
    
    # Create a progress indicator
    progress <- shiny::Progress$new()
    progress$set(message = "Computing data", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Update progress indicator
    progress$inc(1/5, detail = paste("Compiling program"))
    
    # Compile the Go program
    compile_result <- system("go build -o P4S-term-project.exe main.go", intern = TRUE)
    print(compile_result)  # For debugging, to see compile output
    
    # Update progress indicator
    progress$inc(1/5, detail = paste("Running program"))
    
    # Run the compiled Go program
    run_result <- system(paste("./P4S-term-project", genomeFilePath, geneFilePath, percentIdentityValue, preservationThresholdValue), intern = TRUE)
    print(run_result)  # For debugging, to see runtime output
    
    # Update progress indicator
    progress$inc(1/5, detail = paste("Reading in the trees"))
    
    # Read and plot the nucleotide tree data
    tree <- read.nhx("NucleotideNewickOutput.txt")
    p <- ggtree(tree) + geom_text(aes(label=E), hjust=-.3)
    
    # Read and plot the protein tree data
    tree2 <- read.nhx("ProteinNewickOutput.txt")
    p2 <- ggtree(tree2) + geom_text(aes(label=E), hjust=-.3) + geom_label(aes(x=branch, label=S), fill='lightgreen')
    
    # Update progress indicator
    progress$inc(1/5, detail = paste("Reading in the preservation and likelihood data"))
    
    # Read the likelihood data from the CSV file
    likelihoodData <- read.csv("NucleotideLikelihoodTable.csv", header = FALSE)
    colnames(likelihoodData) <- c("Genome Name", "Likelihood")

    # Read the likelihood data from the CSV file
    preservationData <- read.csv("ProteinPreservationTable.csv", header = FALSE)
    colnames(preservationData) <- c("Genome Name", "Preservation")
    
    # Update progress indicator
    progress$inc(1/5, detail = paste("Plotting data"))

    # Render the plot in the Shiny app
    output$nucleotidePlot <- renderPlot({
      msaplot(p,fasta="NucleotideFastaOutput.txt",offset = 4.5, height=.75)
    }, height = 22*nrow(preservationData))
    output$proteinPlot <- renderPlot({
      msaplot(p2,fasta="ProteinFastaOutput.txt",offset = 4.5, height=.75)
    }, height = 22*nrow(preservationData))
    output$likelihoodPlot <- renderDT(likelihoodData)
    output$preservationPlot <- renderDT(preservationData)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)