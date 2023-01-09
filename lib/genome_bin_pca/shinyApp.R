
library('shiny')
library('tidyverse')

load("pca_data.RData")

ui <- fluidPage(

  selectInput("bin_size", "Genome Bin Size:",
              c("100"=100,
                "1,000"=1000,
                "10,000"=10000)),
  selectInput("species", "Species:",
              #c("Listeria monocytogenes"="s__Listeria_monocytogenes",
                #"Escherichia coli"="s__Escherichia_coli",
                #"Salmonella enterica"="s__Salmonella_enterica")),
              c(
                "Listeria monocytogenes"="Lmonocytogenes_4b",
                "Escherichia coli"="EcoliO157EDL933",
                "Salmonella enterica"="Senterica_ATCC14028"
                )),
              
  fluidRow(
    column(12, plotOutput("graph"))
  ),

)


server <- function(input, output, session) {

  output$table <- renderTable(df_list[[paste0(input$bin_size, input$species)]] %>% dplyr::select(date,type,PC1,PC2))
  output$graph <- renderPlot({
    df_plot <- df_list[[paste0(input$bin_size, input$species)]]
    ggplot(df_plot, aes(x=PC1,y=PC2)) + geom_point(aes(pch=type, color=date), size=3) + theme_bw(16)}, res=100)
  
  #df_plot_reactive <- reactive({df_plot[[as.numeric(input$bin_size)]]})
  
}

shinyApp(ui, server)