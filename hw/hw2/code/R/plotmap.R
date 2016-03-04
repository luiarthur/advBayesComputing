library(LatticeKrig) # quilt.plot

plotmap <- function(y,s,main.plot='',bks=quantile(y,c(.1,.9)),margin=0,
                    col.map=colorRampPalette(c('dark blue','grey90','dark red'))(length(y)),
                    ylim.map = range(s[,2])+c(-margin,margin),
                    xlim.map = range(s[,1])+c(-margin,margin),...) {
  y <- ifelse(y>bks[2],bks[2],y)
  y <- ifelse(y<bks[1],bks[1],y)
  quilt.plot(s[,1],s[,2],y,
             fg='grey90',bty='n',main=main.plot,
             ylim=ylim.map,
             xlim=xlim.map,
             breaks=seq(bks[1],bks[2],len=length(y)+1),
             col= col.map,...)
}

plotmap3d <- function(x,y,z,val,clim.map=NULL,col.map=NULL,...) {
  library(plot3D)
  n <- length(x)
  if (is.null(col.map)) col.map <- colorRampPalette(c("darkblue","lightblue","grey85","pink","darkred"))(n)

  if (is.null(clim.map)) {
    clim.map[1] <- quantile(val,.975)
    clim.map[2] <- quantile(val,.025)
  }

  val.map <- val
  val.map <- ifelse(val.map > clim.map[2], clim.map[2], val.map)
  val.map <- ifelse(val.map < clim.map[1], clim.map[1], val.map)

  scatter3D(x,y,z,ticktype="detailed",clim=clim.map,col=col.map,colvar=val.map,...)
}
