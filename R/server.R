



dykifier.server <- function(input,output,session){
  library("inline")
  threshold <- reactiveValues()
  threshold$value <- 0.02

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
    a.ext <- isolate( input$a.ext )
    an <- 1 / ( sqrt(2*pi)*a.ext)
    
    av <- 2*a.ext^2
    
    withProgress(message = 'Updating probabilities', value = 0, {
  
     

    k <- ncol(angles1)
    na <- nrow(angles1)
    p <- matrix(NA, ncol=3, nrow = k*na )

    index <- 1
    for(i in 1:na ){
      incProgress( 0.5 / na )
      si = i %% length(segments)
      if(si==0){
        si = length(segments)
      }
      d.ext <- sdevs[si, 1] * isolate( input$d.ext )
      dv <- 2*d.ext^2
      dn <- 1 / ( sqrt(2*pi)*d.ext)
      pd <- dn * exp( -  dists[i,]^2 / dv )
      pa1 <- av * exp( - (angles1[i, ] - 1)^2 / av )
      pa2 <- av * exp( - (angles2[i, ] - 1)^2 / av )
      #p[i, ] <- pd*pa1*pa2
      p[index:(index+k-1), ] <- cbind(pd*pa1*pa2, rep(i, k), KNN$nn.index[i, ] )

      index <- index+k
    }
    
    })

    porder <- order(p[,1], decreasing=TRUE)
    segs <-  matrix(NA, ncol=5, nrow=na )

    cpp_segs_src <- '
      Rcpp::IntegerVector porder(po);
      Rcpp::NumericMatrix p(probs);
      Rcpp::NumericMatrix segs(hsegs);
      Rcpp::NumericMatrix coords(x);
      int np = porder.size();
      int na = segs.nrow();
      bool *used = new bool[na];
      for(int i=0; i<na; i++){
        used[i] = false;
      }
      int index = 0;
    
      for(int i=0; i<np; i++){
        int i1 = p(porder[i]-1, 1)-1;
        if( used[i1] ){
          continue;
        }
        int i2 = p(porder[i]-1, 2)-1;
        if( used[i2] ){
          continue;
        }
      
        segs(index, 0) = p( porder[i]-1, 0);
        segs(index, 1) = coords(i1, 0);
        segs(index, 2) = coords(i1, 1);
        segs(index, 3) = coords(i2, 0);
        segs(index, 4) = coords(i2, 1);
        used[i1] = true;
        used[i2] = true;
        index++; 
     }
     delete[] used;
     return segs;
   '
   cpp_segs <- cxxfunction(
                   signature(po="integer", probs="matrix", 
                             hsegs="matrix", x="matrix" ), 
                   cpp_segs_src, plugin="Rcpp")



   segs <- cpp_segs(porder, p, segs, coords)
   
   # used <- rep(FALSE, na )
   # index <- 1
   # for(i in 1:length(porder) ){
   #   incProgress( 0.5 / (k*na) )
   #   i1 <- p[porder[i], 2]
   #   if( used[i1] ){
   #     next
   #   }
   #   i2 <- p[porder[i], 3]
   #   if( used[i2] ){
   #     next
   #   }
   #   segs[index, ] <- c( p[ porder[i], 1], coords[i1, ], coords[i2, ] )
   #   used[i1] = TRUE
   #   used[i2] = TRUE
   #   index <- index+1
  #  }

  #  })

    segs <- segs[complete.cases(segs), ]
    segs[,1] <- segs[,1] / max(segs[,1])
    colnames(segs) <- c("Probability", "x1", "y1", "x2", "y2")
    segs
  })



  intersection.selected <- reactive({

    segs <- intersection.probs()
    segs <- segs[ which( segs[,1] >= threshold$value ) , ]

    segs
  })



  ## Dykes plot ##
  output$dykes <- renderPlot({
    plot(X, pch=".", asp=1, xlab="x", ylab="y")
    for( i in 1:length(segments) ){
      lines( segments[[i]],  col=tcolors[type[i]] )
    }
    hsegs <- intersection.selected()
    segments( hsegs[,2], hsegs[, 3], hsegs[,4], hsegs[,5], col="gray") 
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

    } 
  })



  output$dykes.zoom <- renderPlot({
    plot(X, pch=".", asp=1, xlim = ranges$x, ylim = ranges$y, xlab="x", ylab="y")
    for( i in 1:length(segments) ){
      lines( segments[[i]],  col=tcolors[type[i]], lwd=3 )
    }
    #for(i in 1:nrow(dirs) ){
    #   l1 <- rbind(  3 * sdevs[i,1]*dirs[i,] + means[i,] ,
    #          -3 * sdevs[i,1]*dirs[i,] + means[i,] )
    #    lines(l1, col="#FF990075" )
    #}
    hsegs <- intersection.selected()
    segments( hsegs[,2], hsegs[, 3], hsegs[,4], hsegs[,5], col="gray", lwd=2) 


  })



  output$histogram <- renderPlot({
    par(mar=c(4,4,1,1) )
    hsegs <- intersection.probs()
    hist(hsegs[,1], 200, main="", xlab="Probability")
    abline(v=threshold$value)
  })

  ob.hclick <- observeEvent(input$h_click, {
    if( !is.null( input$h_click) ){
      ev <- input$h_click
      threshold$value <- input$h_click$x
    }
  })

  output$lhistogram <- renderPlot({
    par(mar=c(4,4,1,1) )
    hsegs <- intersection.selected()
    xd <- hsegs[,2] - hsegs[,4]
    yd <- hsegs[,3] - hsegs[,5]
    ls <- sqrt(xd^2 + yd^2)
    hist(ls, input$nbins, main="", xlab="Lengths")
  })  


  output$downloadData <- downloadHandler( "phantom-dikes.csv",
    content = function(file) {
      write.csv(intersection.selected(),  file)
    }
  )

}



