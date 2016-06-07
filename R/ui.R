
dykifier.ui <- function(){


fluidPage(
    fluidRow(
      column(3, 
        plotOutput( "dykes", brush=brushOpts("dykes_brush", resetOnNew=TRUE), 
                    width=250, height=500) 
      ),
      column(5, 
        plotOutput( "dykes.zoom", width=500, height=500)
      ),
      column(4, 
        sliderInput("d.ext", "Distance:", min=0, max=2000, step = 10, value=1000),
        sliderInput("a.ext", "Angle:", min=0, max=1, step = 0.001, value=0.1),
        plotOutput("histogram"),
        sliderInput("p.thres", "Probability:", min=0, max=1, step = 0.01, value=0.2),
        actionButton("update", label = "Apply Changes")
      )
    )
)

}


