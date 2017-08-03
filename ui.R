
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyTree)

shinyUI(fluidPage(
  titlePanel('Pima gene expression analysis'),
  sidebarLayout(
    sidebarPanel(
      flowLayout(
        tags$h5("Select values to show:"),
        uiOutput( "stats"),
        tags$h5("Select tissue, trait, comparison or search it below:"),
        shinyTree(outputId = "tree", 
                  checkbox = TRUE,
                  search = T)
      )
    ),
    mainPanel(
      # verbatimTextOutput("selTxt")
      DT::dataTableOutput('mytable1' )
    )
  )
))
