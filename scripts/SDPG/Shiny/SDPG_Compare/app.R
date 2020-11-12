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

load('/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/pair_comp.RDATA')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/SDPG/make_comp_plot.R')

load('/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/results.RDATA')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/SDPG/make_box_plot.R')

##### CLEAN DATA
lidres <- lidres[lidres$orf_name != 'YHR021W-A' &
                   lidres$orf_name != 'BOR' &
                   lidres$orf_name != 'REF' &
                   !is.na(lidres$orf_name),]
# Removing plates with bad colony grids
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'GLU' & lidres$rep != 'R3' & lidres$hours >= 20),]
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'CAS' & lidres$rep == 'R2' & lidres$hours >= 20),]

# Saturation information
lidres$time <- NA
lidres$time[lidres$density == 1536 & lidres$hours == 48 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 1536 & lidres$hours %in% c(72, 78) & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 24 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 36 & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[is.na(lidres$time)] <- "Other"

lidres <- lidres[lidres$time == 'Saturated',]

# Factorize the Experimental Arms
lidres$arm <- factor(lidres$arm, levels = c('GLU','CAS','SDA'))
# Factorize the ORFs
lidres$orf_name <- factor(lidres$orf_name, levels = sort(unique(lidres$orf_name)))

# Removing outliers
for (d in unique(lidres$density)) {
  for (a in unique(lidres$arm[lidres$density == d])) {
    for (r in unique(lidres$rep[lidres$density == d &
                                lidres$arm == a])) {
      for (h in unique(lidres$hours[lidres$density == d &
                                    lidres$arm == a &
                                    lidres$rep == r])) {
        for (o in unique(lidres$orf_name[lidres$density == d &
                                         lidres$arm == a &
                                         lidres$rep == r &
                                         lidres$hours == h])) {
          temp <- lidres$fitness[lidres$density == d &
                                   lidres$arm == a &
                                   lidres$rep == r &
                                   lidres$hours == h &
                                   lidres$orf_name == o]
          m <- median(temp, na.rm = T)
          madev <- mad(temp)
          ul <- m + 3*madev
          ll <- m - 3*madev

          lidres$fitness[lidres$density == d &
                           lidres$arm == a &
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness > ul] <- NA
          lidres$fitness[lidres$density == d &
                           lidres$arm == a &
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness < ll] <- NA


        }
        cont_median <- median(lidres$average[lidres$density == d &
                                               lidres$arm == a &
                                               lidres$rep == r &
                                               lidres$hours == h &
                                               lidres$orf_name == 'BF_control'], na.rm = T)
        lidres$ccs[lidres$density == d &
                     lidres$arm == a &
                     lidres$rep == r &
                     lidres$hours == h] <- lidres$cs_median[lidres$density == d &
                                                              lidres$arm == a &
                                                              lidres$rep == r &
                                                              lidres$hours == h] * cont_median
      }
    }
  }
}

# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("The Famous 28 Pair-wise Comparison"),
   
   sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "g1_orf_name", 
                    label = "1. Select Group 1 ORF/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g1_orf_name[pair_comp$g1_orf_name != 'REF']),
                    selected = "BF_control", multiple = T),
        selectInput(inputId = "g1_density", 
                    label = "2. Select Group 1 Density/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g1_density),
                    selected = "1536", multiple = T),
        selectInput(inputId = "g1_arm", 
                    label = "3. Select Group 1 Condition/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g1_arm),
                    selected = "GLU", multiple = T),
        selectInput(inputId = "g1_rep",
                    label = "4. Select Group 1 Replicate/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g1_rep),
                    selected = "R1", multiple = T),
        
        selectInput(inputId = "g2_orf_name", 
                    label = "5. Select Group 2 ORF/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g2_orf_name[pair_comp$g2_orf_name != 'REF']),
                    selected = "BF_control", multiple = T),
        selectInput(inputId = "g2_density", 
                    label = "6. Select Group 2 Density/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g2_density),
                    selected = "1536", multiple = T),
        selectInput(inputId = "g2_arm", 
                    label = "7. Select Group 2 Condition/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g2_arm),
                    selected = "GLU", multiple = T),
        selectInput(inputId = "g2_rep",
                    label = "8. Select Group 2 Replicate/s",  # Give the input a label to be displayed in the app
                    choices = unique(pair_comp$g2_rep),
                    selected = "R1", multiple = T),
        
        radioButtons(inputId = "mat_type",
                     label = "9. Matrix Type",
                     choices = list("Wilcoxon P-value" = "pvalue",
                                    "Cliff's Delta Estimate" = "cliffdelta",
                                    "Effect Size Magnitude" = "magnitude",
                                    "Phenotype" = "phenotype"),
                     selected = "pvalue"),
        submitButton(text = "Plot", icon = NULL, width = NULL),
        # actionButton("plot_comp_plot", "Pair-wise Fitness Comparison"),
        # actionButton("plot_box_plot", "Fitness Box Plot"),
        width = 3),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("compPlot")
        # br(),
        # plotOutput("boxPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$compPlot <- renderPlot({
    ggpubr::ggarrange(make_comp_plot(pair_comp, input$mat_type,
                      input$g1_orf_name, input$g1_density, input$g1_arm, input$g1_rep,
                      input$g2_orf_name, input$g2_density, input$g2_arm, input$g2_rep),
                      make_box_plot(lidres,
                                    input$g1_orf_name, input$g1_density, input$g1_arm, input$g1_rep,
                                    input$g2_orf_name, input$g2_density, input$g2_arm, input$g2_rep),
                      nrow = 2, heights = c(4,1))
     }, height = 1200, width = 1000
     )
  
  # output$boxPlot <- renderPlot({
  #                  make_box_plot(lidres,
  #                                input$g1_orf_name, input$g1_density, input$g1_arm, input$g1_rep,
  #                                input$g2_orf_name, input$g2_density, input$g2_arm, input$g2_rep)
  #   }, height = 200, width = 1000
  #   )

  # observeEvent(input$plot_comp_plot,
  #              output$compPlot <- renderPlot({
  #                make_comp_plot(pair_comp, input$mat_type,
  #                   input$g1_orf_name, input$g1_density, input$g1_arm, input$g1_rep,
  #                   input$g2_orf_name, input$g2_density, input$g2_arm, input$g2_rep)
  #  }, height = 800, width = 1000))
  # 
  # observeEvent(input$plot_box_plot,
  #              output$boxPlot <- renderPlot({
  #                make_box_plot(lidres,
  #                               input$g1_orf_name, input$g1_density, input$g1_arm, input$g1_rep,
  #                               input$g2_orf_name, input$g2_density, input$g2_arm, input$g2_rep)
  #              }, height = 200, width = 1000))
  
}

# Run the application 
shinyApp(ui = ui, server = server)

