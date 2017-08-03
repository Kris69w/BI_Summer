
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(shinyTree)

# read in data
myDf <- readRDS("DEGTableWithAllComps.rds")
# read in column names
listOfColNames <- colnames(myDf)[2:ncol(myDf)]
stats <- c("adj.P.Val","rstat")

shinyServer(function(input, output, session) {
  output$stats <- renderUI({
    checkboxGroupInput(inputId = "stats", 
                       label = NA, 
                       choices = as.list(stats), 
                       selected = stats, 
                       inline = T) 
  })
  
  output$tree <- renderTree({
    # Tissue
    tissue <- as.character(unique(lapply(X = strsplit(x = listOfColNames, split = "_", fixed = T),FUN = "[[",3)))
    tree <- list()
    for (myTissue in tissue) {
      # Get samples from this tissue
      currTissueSamples <- listOfColNames[grep(pattern = paste("_",myTissue,"_",sep = ""), x = listOfColNames, ignore.case = F, fixed = T)]
      # get available biopsies for this tissue
      currTissueClasses <- as.character(unique(lapply(X = strsplit(x = currTissueSamples, split = "_", fixed = T),FUN = "[[",2)))
      currTissueClasses <- paste("B",currTissueClasses, sep = "_")
      subClass <- list()
      for (myClass in currTissueClasses) {
        # get samples from this biopsy
        currClassSamples <- currTissueSamples[grep(pattern = paste(myClass,"_",sep = ""), x = currTissueSamples, ignore.case = F, fixed = T)]
        # Get comps from this class
        currClassComps <- as.character(unique(lapply(X = strsplit(x = currClassSamples, split = "_", fixed = T),FUN = function(x) {
          last <- length(x) - 1
          x <- paste(x[4:last],collapse = "_")
          return(x) })))
        subComps <- list()
        for (myComp in sort(currClassComps)) {
          subComps[[myComp]] <- paste(myClass,myTissue,myComp, sep = "_")
        }
        subClass[[myClass]] <- subComps
      }
      tree[[myTissue]] <- subClass
    } 
    return(tree)
  })
  output$mytable1 <- DT::renderDataTable({
    #output$selTxt <- renderText({
    myDf[,1] = paste("<a href=https://www.ncbi.nlm.nih.gov/gene?term=",as.character(myDf[,1]),"%5BSymbol%5D%20AND%209606%5Btaxid%5D&cmd=DetailsSearch>",as.character(myDf[,1]),"</a>", sep = "")
    if (is.null(input$tree)) {
      DT::datatable(myDf[,c(colnames(myDf)[1]),drop = FALSE], escape = FALSE)   
    } else {
      elements <- get_selected(input$tree, format = "slices")
      colNames <- unlist(lapply(X = elements ,FUN = function(x){
        comp <- names(x[[1]][[1]][1])
        if (!(is.null(comp))) {
          class <- names(x[[1]][1])
          tissue <- names(x[1])
          colName_pre <- paste(class,tissue,comp, sep = "_")
          colName <- paste(colName_pre, input$stats, sep = "_")
          return(colName)
        }
      }))
      DT::datatable(data = myDf[,c(colnames(myDf)[1],colNames),drop = FALSE], 
                    rownames = F,
                    escape = FALSE,
                    #filter = "top",
                    # extensions = 'Buttons',
                    options = list(
                      pageLength = 10,
                      lengthMenu = c(10,100,500, 1000, 2000,5000),
                      # dom = 'Bfrtip',
                      # buttons = c('copy', 'csv', 'excel'),
                      initComplete = JS(
                        "function(settings, json) {",
                        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                        "}")
                    )
      )    
    }
  })
})


