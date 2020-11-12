library(shiny)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)

load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/all_methods.RData')
source(file = '/home/sbp29/R/Projects/adaptivefitness/scripts/SDPG/plot_growth_curves.R')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("The Famous 28 Liquid Growth Results"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "orfs", 
                  label = "1. Select Group ORF/s",  # Give the input a label to be displayed in the app
                  choices = unique(all_results$orf_name[all_results$orf_name != 'BF_control']),
                  selected = unique(all_results$orf_name[all_results$orf_name != 'BF_control']), multiple = T),
      selectInput(inputId = "condition", 
                  label = "2. Select Condition/s",  # Give the input a label to be displayed in the app
                  choices = unique(all_results$condition),
                  selected = unique(all_results$condition), multiple = T),
      selectInput(inputId = "replicate",
                  label = "3. Select Replicate/s",  # Give the input a label to be displayed in the app
                  choices = unique(all_results$exp_rep),
                  selected = unique(all_results$exp_rep), multiple = T),
      radioButtons(inputId = "gc",
                   label = "4. Growth Curves",
                   choices = list("Yes" = "Y",
                                  "No" = "N"),
                   selected = "N"),
      sliderInput(inputId = "n_col",
                  label = "5. Modify Growth Curve Plot",
                  min = 1, max = 5,
                  value = 1),
      
      submitButton(text = "Plot", icon = NULL, width = NULL),
      
      width = 3),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("growthPlot")
      # br(),
      # plotOutput("boxPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$growthPlot <- renderPlot({
    plot_growth_curves(input$condition, input$replicate, input$orfs, input$gc, input$n_col)
  }, height = 900, width = 1200
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

