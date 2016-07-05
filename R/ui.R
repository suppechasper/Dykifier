
dykifier.ui <- function(){


fluidPage(
    fluidRow(
      column(3, 
        plotOutput( "dykes", brush=brushOpts("dykes_brush", resetOnNew=TRUE), 
                    width=250, height=500),
         downloadButton('downloadData', 'Download phantom dikes') 
      ),
      column(5, 
        plotOutput( "dykes.zoom", width=500, height=500)
      ),
      column(4, 
        sliderInput("d.ext", "Distance:", min=0, max=100, step = 0.5, value=3),
        sliderInput("a.ext", "Angle:", min=0, max=1, step = 0.001, value=0.1),
        actionButton("update", label = "Update Probabilities"),
        hr(),
        plotOutput("histogram", height=180, click = "h_click"),
        plotOutput("lhistogram", height=240),
        sliderInput("nbins", "#bins:", min=10, max=200, step = 1, value=100)
      )
    )
)

}


