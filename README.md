# Dykie Hallucinator
Shiny app for hallucinating dykes.


## Installation ##

1. Install R
2. Start R
3. In R run the commands:
    1. install.packages( c("shiny", "RColorBrewer", "inline", "devtools") )
    2. library("devtools")
    3. devtools::install_github( "suppechasper/Dykifier" )
     
    
## Running ##

1. Start R.
2. In R run the commands:
   1. library(Dykifier)
   2. dykify() for a subset, or dykify("intersections") for the whole data set


The plot on the left shows an overview. By dragging an area a zoomed in version
is shown in the middle plot. The colors are quartz diorite (red), basalt
(green), other (blue) and the gray lines are the hallucinated dykes.

The sliders adjust how the probability of connecting dyke segments is computed.
For angle a smaller value means that more emphasis is put on the segments
connecting in a straight line, larger values allow more kinks. For the distance
a smaller value puts emphasis on connecting only dykes that are close by, a
larger value will still prefer smaller distance but penalizes larger distances
less.  The probability slider is a threshold for excluding connections that are
smaller than the selected value. The plot shows the distributions of
probability values for connecting segments.
After changing the values the Apply Updates button needs to be clicked to
effect the changes.

