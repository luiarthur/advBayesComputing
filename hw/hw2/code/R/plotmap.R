library(LatticeKrig) # quilt.plot

plotmap <- function(y,s,main.plot='',bks=quantile(y,c(.1,.9)),margin=0,
                    col.map=colorRampPalette(c('dark blue','grey90','dark red'))(length(y)),
                    ylim.map = range(s[,2])+c(-margin,margin),
                    xlim.map = range(s[,1])+c(-margin,margin)) {
  quilt.plot(s[,2],s[,1],y,
             fg='grey90',bty='n',main=main.plot,
             ylim=ylim.map,
             xlim=xlim.map,
             breaks=seq(bks[1],bks[2],len=length(y)+1),
             col= col.map)
}
