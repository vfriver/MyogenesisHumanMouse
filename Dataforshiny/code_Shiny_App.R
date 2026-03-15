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
library(shiny)
library(DT)
library(plotly)
library(tidyr)
library(dplyr)
library(shinythemes)

tableHM <- read.csv("/home/log2TPM_combined.csv")  # Adjust the file name 
#change the name of columns, if requiered
colnames(tableHM)[1] <- "GENE"
colnames(tableHM)
# Logic: Human (Cols 2-5), Mouse (Cols 6-9)
# Ensure tableHM is loaded in your environment
if(!"mean_human" %in% colnames(tableHM)){
  tableHM <- tableHM %>%
    mutate(
      mean_human = rowMeans(.[, 2:5], na.rm = TRUE),
      mean_mouse = rowMeans(.[, 6:9], na.rm = TRUE)
    )
}

# USER INTERFACE (UI) 
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Gene Expression Explorer: Human vs. Mouse Myogenesis"),
  
  sidebarLayout(
    # LEFT COLUMN: DATA TABLE
    sidebarPanel(
      width = 4,
      h4("Gene Dataset"),
      p("Select a row to update the plots on the right."),
      DTOutput("tableGenes")
    ),
    
  
    mainPanel(
      width = 8,
    
      div(
        h4("Global Correlation (Mean Expression)"),
        plotlyOutput("corPlot", height = "350px")
      ),
      hr(),
      div(
        h4("Distribution per Species"),
        plotlyOutput("genePlot", height = "350px")
      )
    )
  )
)

# --- (SERVER) ---
server <- function(input, output, session) {
  

  output$tableGenes <- renderDT({
    datatable(tableHM[, c("GENE", "mean_human", "mean_mouse")], 
              selection = 'single', 
              rownames = FALSE,
              options = list(pageLength = 15, dom = 'ftp', scrollX = TRUE))
  })
  
  # Red point
  output$corPlot <- renderPlotly({
    s <- input$tableGenes_rows_selected
    
    p <- plot_ly(tableHM, x = ~mean_human, y = ~mean_mouse, type = 'scatter', mode = 'markers',
                 marker = list(color = 'rgba(200, 200, 200, 0.4)', size = 6),
                 text = ~GENE, hoverinfo = 'text', name = "All Genes") %>%
      layout(xaxis = list(title = "Average log2(TPM+1), Human"),
             yaxis = list(title = "Average log2(TPM+1), Mouse"),
             showlegend = FALSE,
             margin = list(t = 30))
    
    # Add point
    if (length(s)) {
      selected_gene <- tableHM[s, ]
      p <- p %>% add_markers(x = selected_gene$mean_human, 
                             y = selected_gene$mean_mouse,
                             marker = list(color = 'red', size = 12, 
                                           line = list(color = 'black', width = 1)),
                             name = "Selected",
                             inherit = FALSE)
    }
    p
  })
  
  # Vlnplot
  output$genePlot <- renderPlotly({
    s <- input$tableGenes_rows_selected
    
    if (length(s) == 0) {
      return(plot_ly() %>% 
               add_annotations(text = "Please select a gene from the table", 
                               showarrow = F, font = list(size = 16, color = "grey")))
    }
    
    gen_data <- tableHM[s, ]
    df_plot <- data.frame(
      Value = as.numeric(gen_data[1, 2:9]),
      Species = c(rep("Human", 4), rep("Mouse", 4))
    )
    
    plot_ly(df_plot, x = ~Species, y = ~Value, split = ~Species, type = 'violin',
            box = list(visible = T),
            meanline = list(visible = T),
            color = ~Species,
            colors = c("#2c3e50", "#e67e22")) %>%
      layout(title = list(text = paste0("<b>Gene: ", gen_data$GENE, "</b>")),
             yaxis = list(title = "Normalized expression"),
             xaxis = list(title = ""),
             showlegend = FALSE)
  })
}

shinyApp(ui = ui, server = server)
