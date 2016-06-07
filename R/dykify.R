
dykify <- function(datafile = "intersections-subset"){
  library(shiny)
  library(RColorBrewer)

  assign("datafile", datafile, envir = .GlobalEnv)
  data( list = datafile )

  ui <- dykifier.ui()
  shinyApp(ui, dykifier.server)
}
