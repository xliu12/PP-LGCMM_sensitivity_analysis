
# separate plot for originally 
#significant & non-significant mediation effect

library(mvtnorm)
library(tidyverse)
options(dplyr.print_max = 1e9)
options(tibble.width = Inf)

# load("~/Google Drive/Labs/JohnnyLab/BSAgcm/fixC_srML_Z.N200.RData")

# significant mediation in orginal analysis ----
sensiplot.sig <- function( plotdat_all, ifsig_original=TRUE, dire_bpath,
                           xaxis="rCIm", yaxis="rCIy",facetgrid="rCSm ~ rCSy"){
  
  plotdat=plotdat_all %>%
    filter(!is.na(range.lo))
  
  colnames(plotdat)[colnames(plotdat)==xaxis]='xaxis'
  sensiplot = ggplot( data = plotdat ) +
    ylab(yaxis) + xlab(xaxis) + xlim(-1,1) + ylim(-1,1) +
    facet_grid( as.formula(facetgrid), labeller = label_both ) 
  sensiplot = sensiplot+
    scale_fill_manual(guide = "legend",name="After accounting for potential confounder C",
                      values = c("Mediation remains significant but point estimate of the b path in the opposite direction"="yellow2"
                                 ,"Mediation becomes non-significant (point estimate of the b path in the opposite direction)"="orange"
                                 ,"Mediation becomes non-significant (point estimate of the b path in the same direction)"="red"
                      ) 
    ) +
    scale_linetype_manual(name="After accounting for potential confounder C",guide = "legend"
                          ,values=c("Point estimate of the b path = 0"=2
                                    ,"z-statistic of testing the b path = critical value"=3
                                    ,"Range limits"=1)
    ) + 
    theme_bw() +
    theme( axis.title = element_text( size = 13 )
           , legend.text = element_text( size = 13 )
           , legend.title = element_text( size = 13 )
           , legend.position = "bottom"
           , legend.box = "vertical"
           , legend.direction = "vertical"
           , legend.key.width = unit(1.5,"cm")
           , legend.key = element_rect(fill = NULL)
           , strip.text = element_text(size = 13)) 
  sensiplot
  # direction ribbons ----
  # ribbons (run first, so the lines would not be covered)
  if( length(which(!is.na(plotdat$bound.losigb)))!=0  ){
    if( !dire_bpath ){ 
      sensiplot = sensiplot +
        # sigb right ribbon
        geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (bound.losigb<zerob)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation remains significant but point estimate of the b path in the opposite direction') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (-1>bound.losigb) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        # geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (-1>zerob) & (-1<bound.losigb) ), colour=NA, alpha=0.7,
        #             aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (-1>zerob) & (-1>bound.losigb)), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) 
      if( nrow(filter(plotdat,  (xaxis >= noconfounding) & (-1>zerob) & (-1<bound.losigb) ))>0 ){
        sensiplot = sensiplot +
          geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (-1>zerob) & (-1<bound.losigb) ), colour=NA, alpha=0.7,
                      aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) 
      }
      sensiplot = sensiplot +  
        # sigb left ribbon
        geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (bound.upsigb>zerob)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation remains significant but point estimate of the b path in the opposite direction') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (bound.upsigb>1) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        # geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (zerob>1) & (bound.upsigb<1) ), colour=NA, alpha=0.7,
        #             aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (zerob>1) & (bound.upsigb>1) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') )
      if( nrow( filter(plotdat,  (xaxis < noconfounding) & (zerob>1) & (bound.upsigb<1) ) )>0 ){
        sensiplot = sensiplot + 
          geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (zerob>1) & (bound.upsigb<1) ), colour=NA, alpha=0.7,
                      aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') )
      }
      sensiplot
    }
    
    if( dire_bpath ){ 
      sensiplot = sensiplot +
        # sigb left ribbon
        geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding) & (bound.losigb<zerob)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation remains significant but point estimate of the b path in the opposite direction') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding) & (-1>bound.losigb) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        # geom_ribbon(data = filter(plotdat,  (xaxis >= noconfounding) & (-1>zerob) & (-1<bound.losigb) ), colour=NA, alpha=0.7,
        #             aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding) & (-1>zerob) & (-1>bound.losigb)), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) 
      if( nrow(filter(plotdat,  (xaxis <= noconfounding) & (-1>zerob) & (-1<bound.losigb) ))>0 ){
        sensiplot = sensiplot +
          geom_ribbon(data = filter(plotdat,  (xaxis <= noconfounding) & (-1>zerob) & (-1<bound.losigb) ), colour=NA, alpha=0.7,
                      aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) 
      }
      sensiplot = sensiplot +  
        # sigb right ribbon
        geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) & (bound.upsigb>zerob)  ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation remains significant but point estimate of the b path in the opposite direction') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) & (bound.upsigb>1) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=zerob, fill='Mediation becomes non-significant (point estimate of the b path in the opposite direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=zerob, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        # geom_ribbon(data = filter(plotdat,  (xaxis < noconfounding) & (zerob>1) & (bound.upsigb<1) ), colour=NA, alpha=0.7,
        #             aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) & (zerob>1) & (bound.upsigb>1) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.losigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') )
      if( nrow( filter(plotdat,  (xaxis > noconfounding) & (zerob>1) & (bound.upsigb<1) ) )>0 ){
        sensiplot = sensiplot + 
          geom_ribbon(data = filter(plotdat,  (xaxis > noconfounding) & (zerob>1) & (bound.upsigb<1) ), colour=NA, alpha=0.7,
                      aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='Mediation becomes non-significant (point estimate of the b path in the same direction)') )
      }
      sensiplot
    }
  }
  
  # lines ----
  sensiplot = sensiplot +
    # range limits
    geom_line( aes(x=xaxis, y=range.lo, linetype="Range limits"),size=0.5 ) +
    geom_line( aes(x=xaxis, y=range.up, linetype="Range limits"),size=0.5 ) +
    # zerob
    geom_line( data = filter(plotdat,  (xaxis < noconfounding)),
               aes(x=xaxis, y=zerob, linetype='Point estimate of the b path = 0'),size=0.5 ) +
    geom_line( data = filter(plotdat,  (xaxis >= noconfounding)),
               aes(x=xaxis, y=zerob, linetype='Point estimate of the b path = 0'),size=0.5 ) +
    # sigb left lines
    geom_line(data = filter(plotdat,  (xaxis < noconfounding)   ), # # & (bound.losigb<bound.upsigb) cases with na zerob: bound.lo=bound.up
              aes(x=xaxis, y=bound.losigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +
    geom_line(data = filter(plotdat,  (xaxis < noconfounding) ),
              aes(x=xaxis, y=bound.upsigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +
    # sig right lines
    geom_line(data = filter(plotdat,  (xaxis >= noconfounding)  ),
              aes(x=xaxis, y=bound.losigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +  
    geom_line(data = filter(plotdat,  (xaxis >= noconfounding)   ),#& (bound.upsigb>bound.losigb)
              aes(x=xaxis, y=bound.upsigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) 
  
  sensiplot
  return(sensiplot)
}

# non-significant mediation in orginal analysis ----
sensiplot.nonsig <- function( plotdat_all, ifsig_original=FALSE, dire_bpath,
                              xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy" ){
  # only display results within the range limits
  plotdat=filter(plotdat_all, !is.na(range.lo)) 
  colnames(plotdat)[colnames(plotdat)==xaxis]="xaxis"
  
  sensiplot = ggplot( data = plotdat ) +
    ylab(yaxis) + xlab(xaxis) + xlim(-1,1) + ylim(-1,1) +
    facet_grid( as.formula(facetgrid), labeller = label_both ) 
  
  sensiplot = sensiplot+
    scale_fill_manual(guide = "legend",name="Shaded area",
                      values = c("the b path becomes significant"="green4")) +
    scale_linetype_manual(name="Line",guide = "legend"
                          ,values=c("Point estimate of the b path = 0"=2 # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
                                    ,"z-statistic of testing the b path = critical value"=3
                                    ,"Range limits"=1)
    ) + 
    theme_bw() +
    theme( axis.title = element_text( size = 13 )
           , legend.text = element_text( size = 13 )
           , legend.title = element_text( size = 13 )
           , legend.position = "bottom"
           , legend.box = "vertical"
           , legend.direction = "vertical"
           , legend.key.width = unit(1.5,"cm")
           , legend.key = element_rect(fill = NULL)
           , strip.text = element_text(size = 13)) 
  sensiplot
  
  if( length(which(!is.na(plotdat$bound.losigb)))!=0 ){
    # ribbons----
    if(!dire_bpath){
      sensiplot= sensiplot+
        # sigb left
        geom_ribbon(data = filter(plotdat, (xaxis<noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat, (xaxis<noconfounding) & (-1>bound.losigb) ), alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=-1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis<noconfounding) & (-1<bound.losigb) ),colour=NA,alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='the b path becomes significant') ) +
        # sigb left
        geom_ribbon(data = filter(plotdat, (xaxis>=noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymin=bound.upsigb, ymax=1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat, (xaxis>=noconfounding) & (1<bound.upsigb) ), alpha=0.7,
                    aes(x=xaxis, ymin=bound.losigb, ymax=1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis>=noconfounding) & (1>bound.upsigb) ),colour=NA,alpha=0.7,
                    aes(x=xaxis, ymin=-1, ymax=bound.losigb, fill='the b path becomes significant') ) 
      sensiplot
      
    }
    
    if(dire_bpath){
      sensiplot= sensiplot+
        # sigb right
        geom_ribbon(data = filter(plotdat, (xaxis>noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymax=bound.losigb, ymin=-1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat, (xaxis>noconfounding) & (-1>bound.losigb) ), alpha=0.7,
                    aes(x=xaxis, ymax=bound.upsigb, ymin=-1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis>noconfounding) & (-1<bound.losigb) ),colour=NA,alpha=0.7,
                    aes(x=xaxis, ymax=1, ymin=bound.upsigb, fill='the b path becomes significant') ) +
        # sigb right
        geom_ribbon(data = filter(plotdat, (xaxis<=noconfounding) ), colour=NA, alpha=0.7,
                    aes(x=xaxis, ymin=bound.upsigb, ymax=1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat, (xaxis<=noconfounding) & (1<bound.upsigb) ), alpha=0.7,
                    aes(x=xaxis, ymin=bound.losigb, ymax=1, fill='the b path becomes significant') ) +
        geom_ribbon(data = filter(plotdat,  (xaxis<=noconfounding) & (1>bound.upsigb) ),colour=NA,alpha=0.7,
                    aes(x=xaxis, ymin=-1, ymax=bound.losigb, fill='the b path becomes significant') ) 
      sensiplot
      
    }
  }
  
  # lines----
  sensiplot = sensiplot +
    # range limits
    geom_line( aes(x=xaxis, y=range.lo, linetype="Range limits"),size=0.5 ) +
    geom_line( aes(x=xaxis, y=range.up, linetype="Range limits"),size=0.5 ) +
    # zerob
    geom_line( aes(x=xaxis, y=zerob, linetype='Point estimate of the b path = 0'),size=0.5 ) +
    # sigb lo left
    geom_line(data = filter(plotdat,  (xaxis<noconfounding) ),#(bound.losigb<zerob)&
              aes(x=xaxis, y=bound.losigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +
    geom_line(data = filter(plotdat,  (xaxis<noconfounding) ),#(bound.losigb<zerob)&
              aes(x=xaxis, y=bound.upsigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +
    # sigb right
    geom_line(data = filter(plotdat,  (xaxis>=noconfounding) ),
              aes(x=xaxis, y=bound.losigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 ) +
    geom_line(data = filter(plotdat,  (xaxis>=noconfounding) ),
              aes(x=xaxis, y=bound.upsigb, linetype='z-statistic of testing the b path = critical value'),size=0.5 )
  
  sensiplot
  
  return(sensiplot)
  
}







