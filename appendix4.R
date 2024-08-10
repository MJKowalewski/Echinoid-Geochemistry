###############################
# Appendix 4 - Principal Component Analysis
# custom functions
# written by M. Kowalewski (kowalewski@ufl.edu)
# Last updated: August 10, 2024
###############################

# pca.main.F function====
pca.main.F <- function(x, scale=F, pcs=c(1,2), ...) {
  y <- prcomp(x, scale.=scale)
  eval <- y$sdev^2
  p.eval <- round(100*eval/sum(eval),2)
  names(eval) <- paste0('PC',1:ncol(x))
  names(p.eval) <- paste0('PC',1:ncol(x))
  plot(y$x[,pcs[1]], y$x[,pcs[2]],
       xlab=paste0('PC', pcs[1], ' [', p.eval[pcs[1]],'%]'),
       ylab=paste0('PC',pcs[2], ' [',p.eval[pcs[2]],'%]'), ...)
  list(evals=eval, perc.evals=p.eval, evectors=y$rotation,
              scores=y$x, scale=paste('scale =',scale), pcs=pcs)
}

# pca.vec.F function====
pca.vec.F <- function(x, scale=F, pcs=c(1,2), eload=F, 
                      x.axis.label=T, lab.offset=1.25, ...) {
  y <- prcomp(x, scale.=scale)
  eval <- round(100*y$sdev^2/sum(y$sdev^2),2)
    plot(0, 0, type='n', xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), las=1, xpd=NA,
         ylab='', xlab='', axes=F, ...)
     axis(4, las=1, cex.axis=0.7, at=seq(-1,1,0.5), labels=seq(-1,1,0.5))
     mtext(side=4, line=3, cex=0.7, paste0('PC', pcs[2], ' [',eval[pcs[2]],']'))
     box()
     if (x.axis.label) mtext(side=3, line=3, paste0('PC', pcs[1], '[',eval[pcs[1]],']'), cex=0.7)
     lines(cos(seq(0, 2 * pi, length = 200)), sin(seq(0, 2 * pi, length = 200)),
          lwd=3, col='skyblue')
     lines(sqrt(0.5) * cos(seq(0, 2 * pi, length = 200)),
          sqrt(0.5) * sin(seq(0, 2 * pi, length = 200)), lwd=1, col='skyblue')
     if (!eload) {
      for (j in 1:ncol(x)) { 
        arrows(0, 0, cor(y$x[,pcs[1]], x[,j]), cor(y$x[,pcs[2]], x[,j]),
               length=0.05)
        points(lab.offset*cor(y$x[,pcs[1]], x[,j]), lab.offset*cor(y$x[,pcs[2]], x[,j]),
               pch=16, cex=2)
        text(lab.offset*cor(y$x[,pcs[1]], x[,j]), lab.offset*cor(y$x[,pcs[2]], x[,j]),
             colnames(x)[j], cex=0.7, col='white')
        axis(side=3, cex.axis=0.7, at=seq(-1,1,0.5), labels=seq(-1,1,0.5))
      }
    }
    if (eload) {
      arrows(0, 0, y$rotation[,pcs[1]], y$rotation[,pcs[2]], length=0.05)
      points(lab.offset*y$rotation[,pcs[1]], lab.offset*y$rotation[,pcs[2]], cex=2,
             pch=16)
      text(lab.offset*y$rotation[,pcs[1]], lab.offset*y$rotation[,pcs[2]], cex=0.6,
           colnames(x), col='white')
      
    }
  }  

# pca.scree.F function====
pca.scree.F <- function(x, scale=F, barcol='red', 
                        times=100, cfint=c(0.005, 0.995), ...) {
  y <- prcomp(x, scale.=scale)
  eval <- round(100*y$sdev^2/sum(y$sdev^2),2)
  out1 <- NULL
  for(i in 1:times) {
    ry <- x[sample(1:nrow(x), replace=T),]
    if (prod(apply(ry, 2, sd)) > 0) {
      revals <- prcomp(ry, scale.=scale)$sdev^2
      reval <- 100*revals/sum(revals)
      out1 <- rbind(out1, reval)
    }
  }
  plot(1:ncol(x), eval, type='n', xlab='PCs', ylab='',
     axes=F, ylim=c(0, max(out1)))
  axis(1, cex.axis=0.7,
     at=1:ncol(x), labels=1:ncol(x))
  axis(4, las=1, cex.axis=0.7)
  mtext(side=4, line=3, '% variance', cex=0.7)
  abline(h=0.7*mean(eval), lwd=1, col='skyblue3', xpd=T)
  points(1:ncol(x), eval, type='o', ...)
  for (i in 1:ncol(x)) points(c(i,i),
                            c(quantile(out1[,i], prob=cfint[1]), 
                              quantile(out1[,i], prob=cfint[2])),
                            type='l', col=barcol, lwd=1.5)
  mtext(side=3, line=-1, paste('iter =', times), cex=0.6, adj=0.05)
  box()
}

pca.vec.F2 <- function(x, scale=F, pcs=c(1,2), eload=F, 
                      x.axis.label=T, lab.offset=1.25, ...) {
  y <- prcomp(x, scale.=scale)
  eval <- round(100*y$sdev^2/sum(y$sdev^2),2)
  plot(0, 0, type='n', xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), las=1, xpd=NA,
       ylab='', xlab='', axes=F, ...)
  axis(2, las=1, cex.axis=0.7, at=seq(-1,1,0.5), labels=seq(-1,1,0.5))
  axis(1, las=1, cex.axis=0.7, at=seq(-1,1,0.5), labels=seq(-1,1,0.5))
  mtext(side=2, line=3, cex=0.7, paste0('PC', pcs[2], ' [',eval[pcs[2]],'%]'))
  box()
  if (x.axis.label) mtext(side=1, line=3, paste0('PC', pcs[1], ' [',eval[pcs[1]],'%]'), cex=0.7)
  lines(cos(seq(0, 2 * pi, length = 200)), sin(seq(0, 2 * pi, length = 200)),
        lwd=3, col='skyblue')
  lines(sqrt(0.5) * cos(seq(0, 2 * pi, length = 200)),
        sqrt(0.5) * sin(seq(0, 2 * pi, length = 200)), lwd=1, col='skyblue')

  if (!eload) {
    for (j in 1:ncol(x)) { 
      arrows(0, 0, cor(y$x[,pcs[1]], x[,j]), cor(y$x[,pcs[2]], x[,j]),
             length=0.05)
      text(lab.offset*cor(y$x[,pcs[1]], x[,j]), lab.offset*cor(y$x[,pcs[2]], x[,j]),
           colnames(x)[j], cex=0.6, col='black')
    }
  }
  if (eload) {
    arrows(0, 0, y$rotation[,pcs[1]], y$rotation[,pcs[2]], length=0.05)
    text(lab.offset*y$rotation[,pcs[1]], lab.offset*y$rotation[,pcs[2]], cex=0.6,
         colnames(x), col='black')
    
  }
}  


