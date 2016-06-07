



dykifier.server <- function(input,output,session){

  X <- c()
  for(i in 1:length(segments) ){
    X = rbind(X,segments[[i]])
  }

  type <- rep(3, nrow(X) )
  type[ rtype == "quartz diorite" ] <- 1
  type[ rtype == "basalt" ] <- 2

  tcolors <- brewer.pal("Set1", n=3)


  ranges <- reactiveValues(x = NULL, y = NULL)


  intersection.probs <- reactive({
    input$update

    d.ext <- isolate( input$d.ext )
    a.ext <- isolate( input$a.ext )
    dn <- 1 / ( sqrt(2*pi)*d.ext)
    an <- 1 / ( sqrt(2*pi)*a.ext)
    
    dv <- 2*d.ext^2
    av <- 2*a.ext^2
    
    withProgress(message = 'Updating probabilities', value = 0, {
    
    p <- matrix(NA, ncol=ncol(angles1), nrow=nrow(angles1) )
    for(i in 1:nrow(angles1) ){
      #incProgress( 1 / nrow(angles1) )
      pd <- dn * exp( -  dists[i,]^2 / dv )
      pa1 <- av * exp( - (angles1[i, ] - 1)^2 / av )
      pa2 <- av * exp( - (angles2[i, ] - 1)^2 / av )
      p[i, ] <- pd*pa1*pa2
    }
    
    })
    p <- p / max(p, na.rm=TRUE)
    #updateSliderInput(session, "p.thres", value = 0.5*pm, max = pm, step=pm/100)
    porder <-  t(apply(-p,1,order))
    list(p=p, porder=porder)
  })

  intersection.selected <- reactive({
    input$update

    p <- intersection.probs()
    x <- list()
    index <- 1
    withProgress(message = 'Updating selection', value = 0, {
    for(i in 1:nrow(p$p)){
      incProgress( 1 / nrow(p$p) )
      if( p$p[i, p$porder[i,1] ] >= isolate(input$p.thres) ){
        x[[index]] = rbind( coords[i, ], 
                            coords[ KNN$nn.ind[i, p$porder[i, 1] ], ]
                          )
        index <- index +1
      } 
    }
    })
    x
  })

  ## Dykes plot ##
  output$dykes <- renderPlot({
    plot(X, pch=".", asp=1, xlab="x", ylab="y")
    for( i in 1:length(segments) ){
      lines( segments[[i]],  col=tcolors[type[i]] )
    }
     
    for(segs in intersection.selected() ){
      lines( segs, col="gray" )
    } 
   # for(i in 1:nrow(dirs) ){
   #    l1 <- rbind(  3 * sdevs[i,1]*dirs[i,] + means[i,] ,
   #           -3 * sdevs[i,1]*dirs[i,] + means[i,] )
   #     lines(l1, col="#FF990075" )
   # }

  })


  observe({
    brush <- input$dykes_brush
    if (!is.null(brush)) {
      ranges$x <<- c(brush$xmin, brush$xmax)
      ranges$y <<- c(brush$ymin, brush$ymax)

    } else {
      ranges$x <<- NULL
      ranges$y <<- NULL
    }
  })

  output$dykes.zoom <- renderPlot({
    plot(X, pch=".", asp=1, xlim = ranges$x, ylim = ranges$y, xlab="x", ylab="y")
    for( i in 1:length(segments) ){
      lines( segments[[i]],  col=tcolors[type[i]] )
    }
    #for(i in 1:nrow(dirs) ){
    #   l1 <- rbind(  3 * sdevs[i,1]*dirs[i,] + means[i,] ,
    #          -3 * sdevs[i,1]*dirs[i,] + means[i,] )
    #    lines(l1, col="#FF990075" )
    #}
    for(segs in intersection.selected()  ){
      lines( segs, col="gray" )
    }


  })

  output$histogram <- renderPlot({
    p <- intersection.probs()
    hist(p$p[1:nrow(p$p) + nrow(p$p)*( p$porder[, 1]-1  )], 100, 
         main="", xlab="Probability")
    abline(v=input$p.thres)
  }) 

}



