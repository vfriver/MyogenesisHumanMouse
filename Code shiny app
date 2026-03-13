# ==============================================================================
# Gene Expression Explorer: Human vs. Mouse Myogenesis
# ==============================================================================
# LOCAL EXECUTION INSTRUCTIONS:
#
# 1. PREREQUISITES:
#    Ensure you have R and RStudio installed. 
#    Install dependencies by running the following in your console:
#    install.packages(c("shiny", "DT", "plotly", "tidyr", "dplyr", "shinythemes"))
#
# 2. FILE ORGANIZATION:
#    Place this script ('app.R') and the data file ('tableHM.csv') 
#    in the same folder on your computer.
#
# 3. LAUNCH:
#    - Open this file in RStudio.
#    - Click the 'Run App' button at the top right of the editor panel.
#    - Or run: shiny::runApp() in the R console.
#
# 4. INTERACTION:
#    Search for a gene in the table and click its row to update the plot.
# ==============================================================================

# Your Shiny code starts here...

tableHM <- read.csv("/home//log2TPM_combined.csv")  # Adjust the file name 
#change the name of columns, if requiered
colnames(tableHM)[1] <- "GENE"
colnames(tableHM)
#####
library(shiny)
library(DT)
library(plotly)
library(tidyr)
library(dplyr)
library(shinythemes)

# Logic: Human (Cols 2-5), Mouse (Cols 6-9)
# Ensure tableHM is loaded in your environment

ui <- fluidPage(
  theme = shinytheme("flatly"), 
  titlePanel("Gene expression explorer: Human vs. Mouse myogenesis"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Instructions"),
      tags$ul(
        tags$li("Search for a gene in the table."),
        tags$li("Select the row to generate the Violin Plot.")
      ),
      hr(),
      plotlyOutput("genePlot", height = "450px")
    ),
    
    mainPanel(
      width = 8,
      h3("Gene Dataset"),
      DTOutput("tableGenes")
    )
  )
)

server <- function(input, output) {
  
  output$tableGenes <- renderDT({
    datatable(tableHM, 
              selection = 'single', 
              rownames = FALSE,
              options = list(
                pageLength = 12,
                searchHighlight = TRUE,
                language = list(search = "Search Gene:"), 
                dom = 'ftip' 
              ))
  })
  
  output$genePlot <- renderPlotly({
    s <- input$tableGenes_rows_selected
    
    if (length(s) == 0) {
      return(plot_ly() %>% 
               add_annotations(text = "Please select a gene\nfrom the table", 
                               showarrow = F, font = list(size = 16, color = "grey")))
    }
    
    gen_data <- tableHM[s, ]
    gene_name <- gen_data$GENE 
    
    df_plot <- data.frame(
      Value = as.numeric(gen_data[1, 2:9]),
      Species = c(rep("Human (LHCN-M2)", 4), rep("Mouse (C2C12)", 4))
    )
    
    plot_ly(df_plot, 
            x = ~Species, 
            y = ~Value, 
            split = ~Species, 
            type = 'violin',
            # --- INTERNAL BOXPLOT SETTINGS ---
            box = list(visible = T, width = 0.15), 
            meanline = list(visible = T),
            # --- POINTS DISABLED ---
            points = FALSE, 
            # ------------------------
            line = list(color = 'black', width = 1.5),
            color = ~Species,
            colors = c("#2c3e50", "#e67e22")) %>%
      layout(
        title = list(text = paste0("<b>Gene: ", gene_name, "</b>"), y = 0.95),
        yaxis = list(title = "Normalized expression", zeroline = F),
        xaxis = list(title = ""),
        showlegend = FALSE,
        margin = list(t = 60)
      )
  })
}

shinyApp(ui = ui, server = server)
