minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){

	  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

    major.ticks <- pretty(lims,n=5)
    if(missing(mn)) mn <- min(major.ticks)
      if(missing(mx)) mx <- max(major.ticks)

      major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

        labels <- sapply(major.ticks,function(i)
			             as.expression(bquote(10^ .(i)))
				               )
        axis(ax,at=major.ticks,labels=labels,...)

	  n <- n+2
	  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
	    minors <- minors[-c(1,n)]

	    minor.ticks = c(outer(minors,major.ticks,`+`))
	      minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


	      axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}
