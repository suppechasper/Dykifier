



dykifier.server <- function(input,output,session){
  library("inline")
  threshold <- reactiveValues()
  threshold$pvalue <- 0.02
  threshold$lvalue <- -0.01

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
  
     

    ns <- length(segments)
    k <- ncol(angles1)
    na <- nrow(angles1)
    p <- matrix( NA, ncol=5, nrow = k*na )

    index <- 1
    for(i in 1:na ){
      incProgress( 0.5 / na )
      si = i %% length(segments)
      if(si==0){
        si = length(segments)
      }
      d.ext <- sdevs[si, 1] * isolate( input$d.ext )
      dv  <- 2*d.ext^2
      dn  <- 1 / ( sqrt(2*pi)*d.ext )
      pd  <- dn * exp( -  dists[i,]^2 / dv )
      pa1 <- av * exp( - (angles1[i, ] - 1)^2 / av )
      pa2 <- av * exp( - (angles2[i, ] - 1)^2 / av )
      #p[i, ] <- pd*pa1*pa2
      p[index:(index+k-1), ] <- cbind(pd*pa1*pa2, rep(i, k),
                                                  KNN$nn.index[i, ],
                                                  rep(i, k) %% ns, 
                                                  KNN$nn.index[i, ]  %% ns )

      index <- index+k
    }
    p[ p[, 4] == 0, 4 ] = ns
    p[ p[, 5] == 0, 5 ] = ns

    })

    porder <- order(p[,1], decreasing=TRUE)
    segs <-  matrix(NA, ncol=7, nrow=na )

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
        segs(index, 5) = p( porder[i]-1, 3 );
        segs(index, 6) = p( porder[i]-1, 4 );        
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
    colnames(segs) <- c("Probability", "x1", "y1", "x2", "y2", "sid1", "sid2")
    segs
  })



  intersection.selected <- reactive({

    segs <- intersection.probs()
    segs <- segs[ which( segs[,1] >= threshold$pvalue ) , ]

    segs
  })




  intersection.selected.connected <- reactive({

    segs <- intersection.selected()
    id.segs <- rep( NA, length(segments) )
    id.ph <- rep( NA, nrow(segs) )
    did <- 0
    for( i in 1:nrow(segs) ){
      cid <- NA
      if( !is.na( id.segs[ segs[i,6] ]) ){
        cid = id.segs[ segs[i,6] ]
      }
      if( !is.na( id.segs[ segs[i,7] ]) ){
        if( !is.na(cid) ){
          #replace with single id
          eid <- id.segs[ segs[i,7] ]
          ind <- which(id.segs == eid)
          id.segs[ind] = cid
          ind <- which(id.ph == eid)
          if(length(ind) > 0 ){
            id.ph[ind] = cid
          }
        }
        else{
          cid = id.segs[ segs[i,7] ]
        }
      }
      if( is.na(cid) ){
        did <- did + 1 
        cid <- did
      }
      id.ph[i] = cid
      id.segs[ segs[i,6] ] = cid
      id.segs[ segs[i,7] ] = cid
    }
    ndykes = length( unique(id.segs) )
    
    for( i in 1:length(id.segs) ){
      if( is.na( id.segs[i] ) ){
        did = did+1
        id.segs[i] = did
      }
    }

    dlengths = rep(0, did+1)
    for(i in 1:did){
      ind <- which(id.ph == i)
      if(length(ind) > 0 ){
        xd <- (segs[ind,2] - segs[ind,4])^2
        yd <- (segs[ind,3] - segs[ind,5])^2
        dlengths[i] <- sqrt( sum(xd+yd) )
      }
      ind <- which(id.segs == i)
      if(length(ind) > 0 ){
        dlengths[i] <- dlengths[i] + sum(slength[ind] )
      }
    }
    
    list( ndykes = ndykes, nsegs = length( unique(id.segs) ) - ndykes, id.segs = id.segs, 
          id.ph = id.ph, dlengths = dlengths )
  })




  ## Dykes plot ##
  output$dykes <- renderPlot({
    plot(X, pch=".", asp=1, xlab="x", ylab="y")
    
    csegs <- intersection.selected.connected()

    for( i in 1:length(segments) ){
      if( csegs$dlengths[ csegs$id.segs[i] ] > threshold$lvalue ){ 
        lines( segments[[i]],  col=tcolors[type[i]], lwd=3 )
      }
    }

    hsegs <- intersection.selected()
    for( i in 1:nrow(hsegs) ){
      if( csegs$dlengths[ csegs$id.ph[i] ] > threshold$lvalue){
        segments( hsegs[i,2], hsegs[i, 3], hsegs[i,4], hsegs[i,5], col="gray", lwd=2) 
      }
    }

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
    
    csegs <- intersection.selected.connected()
    for( i in 1:length(segments) ){
      if( csegs$dlengths[csegs$id.segs[i]] > threshold$lvalue ){ 
        lines( segments[[i]],  col=tcolors[type[i]], lwd=3 )
      }
    }

    hsegs <- intersection.selected()
    for( i in 1:nrow(hsegs) ){
      if(csegs$dlengths[csegs$id.ph[i] ] > threshold$lvalue ){
        segments( hsegs[i,2], hsegs[i, 3], hsegs[i,4], hsegs[i,5], col="gray", lwd=2) 
      }
    }

  })



  output$histogram <- renderPlot({
    par(mar=c(4,4,1,1) )
    hsegs <- intersection.probs()
    hist(hsegs[,1], 200, main="", xlab="Probability")
    abline(v=threshold$pvalue)
  })

  ob.hclick <- observeEvent(input$h_click, {
    if( !is.null( input$h_click) ){
      ev <- input$h_click
      threshold$pvalue <- input$h_click$x
    }
  })

  output$lhistogram <- renderPlot({
    par(mar=c(4,4,1,1) )
    hsegs <- intersection.selected()
    xd <- hsegs[,2] - hsegs[,4]
    yd <- hsegs[,3] - hsegs[,5]
    ls <- sqrt(xd^2 + yd^2)
    hist(ls, input$nlbins, main="", xlab="Phantom Segment Lengths")
  })  

  output$dhistogram <- renderPlot({
    par(mar=c(4,4,1,1) )
    dsegs <- intersection.selected.connected()
    hist(dsegs$dlengths, input$ndbins, main="", xlab="Dike Lengths")
    abline(v=threshold$lvalue)
  })

  ob.dhclick <- observeEvent(input$dh_click, {
    if( !is.null( input$dh_click) ){
      ev <- input$dh_click
      threshold$lvalue <- input$dh_click$x
    }
  })  

  output$downloadData <- downloadHandler( "phantom-dikes.csv",
    content = function(file) {
      write.csv(intersection.selected(),  file)
    }
  )

  output$ndykes.info <- renderText({ 
      cdykes <- intersection.selected.connected()
      paste("Number of connected dykes: ", cdykes$ndykes )
  })
    
  output$nsegs.info <- renderText({ 
      cdykes <- intersection.selected.connected()
      paste("Number of disconnected segments left: ", cdykes$nsegs )
  })


}



