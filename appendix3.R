###
# Appendix 3 - R Script (Gorzelak et al.)
# written by M. Kowalewski (kowalewski@ufl.edu)
# Last updated: August 10, 2024
# salinity data in Appendix 1 from: 
# https://hub.arcgis.com/maps/0a053328e78e40d18f95eb1922df4f3a/about
# temperature data in Appendix 1 from:
# https://www.ncei.noaa.gov/access/gulf-of-mexico-climate/gulf-data.html
# coastal coordinates for map figures (Appendix 5) downloaded from:
# https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset
###

# libraries and functions====
options(scipen = 100)
library(ppcor)
library(asbio)
library(DescTools)
source(file='appendix4.R')
mycol <- c('skyblue1', 'blue3', 'red1','gray40') # colors for taxa (alphabetically by genus)
pchtax <- 21:24
pchreg <- c(21, 22, 24)
mystats <- function(x) c(median(x), mean(x), sd(x), length(x))
medboot <- function(x, times=99, prob=c(0.5, 0.025, 0.975, 0.25, 0.75)) {
  out1 <- vector(mode='numeric', length=times)
  for (i in 1:times) {
    y <- sample(x, replace=T)
    out1[i] <- median(y)
  }
  out1 <- c(out1, median(x))
  out1 <- out1 - (median(out1) - median(x))
  out2 <- as.vector(quantile(out1, prob=prob))
  names(out2) <- prob
  return(out2)
}

# Datasets ====
bulk <- read.csv('appendix1.csv', stringsAsFactors=T) # bulk analyses
nn <- read.csv('appendix2.csv', stringsAsFactors=T) # nano analyses
map2 <- read.csv('appendix5.csv', header=F, na.strings=c(NA, '.', ''))

# group Encope aberrans and  Encope michelini as Encope spp.
bulk$species <- as.character(bulk$species)
bulk$species[bulk$species == 'Encope aberrans'] <- 'Encope spp.'
bulk$species[bulk$species == 'Encope michelini'] <- 'Encope spp.'
bulk$species <- as.factor(bulk$species)
bulk$Pb <- bulk$Pb/1000 # convert from micromolar to milimolar

# change assignment and salinity estimates for the outer Cedar Key
# stations to 'Northern Gulf'. Not clear which variant is most appropriate
# no major conclusions changed depending on station assignment decision
# salinity station: 21FLSEAS_WQX-37SEAS105  31.5000000 31.5000 31.00000 32.00000 31.000 32.000 28.8703 -82.7790 1.6000000
repl.sal <- c(31.5000000, 31.5000, 31.00000, 32.00000, 31.000, 32.000)
bulk$Region[which(bulk$Station == 'CK-2')] <- 'Northern Gulf'
bulk[which(bulk$Station == 'CK-2'), c("sal_mean", "sal_med", "sal_q25", "sal_q75", "sal_min", "sal_max")] <- repl.sal
# bulk$Region[which(bulk$Station == 'CK-5')] <- 'Northern Gulf'
# bulk[which(bulk$Station == 'CK-5'), c("sal_mean", "sal_med", "sal_q25", "sal_q75", "sal_min", "sal_max")] <- repl.sal

# Table 1 Summary of Sites ====
sites <- as.factor(paste(bulk$Region, paste0('Site-', bulk$Station)))
taxtable <- tapply(bulk$species, sites, table)
taxregion <- as.data.frame(tapply(as.factor(bulk$species), bulk$Region, table))
taxregion <- tapply(as.factor(bulk$species), bulk$Region, table)
taxregion2 <- NULL
for (i in 1:length(taxregion)) taxregion2 <- rbind(taxregion2, taxregion[[i]])
rownames(taxregion2) <- names(taxregion)
taxregion2

taxtable2 <- NULL
for (i in 1:length(taxtable)) taxtable2 <- rbind(taxtable2, taxtable[[i]])
table1 <- cbind(latitude=tapply(bulk$Latitude, sites, mean),
      longitude=tapply(bulk$Longitude, sites, mean),
      depth=tapply(bulk$Depth, sites, mean), taxtable2)
table1
write.csv(table1, 'Table 1 Site Summary.csv')

# Maps ====
#------------------------ COASTLINE MAP ------------------------
# coastal coordinates downloaded from:
# https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset


plgs2 <- which(is.na(map2[2,]))

#- FIGURE 1 - Map of all sites====
# color scheme
reg.col <- c('seagreen4', 'coral', 'black')
wt.col <- 'white' # water color
land.col <- 'gray70' # color of land area
coast.col <- 'gray20' # color of coastline
old.stud <- 'gray20'

pdf("FIG 1 MAP_SITES.pdf", height=6.1, width=6.8)
plot.lim <- t(map2[,2:5]) # full map
plot(plot.lim, type='n', xlab='', ylab='', las=1, axes=F)
polygon(t(map2[,2:5]), col=wt.col, border=NA)
for (i in 2:length(plgs2)) {
  if (i < length(plgs2)) {
    a1 <- plgs2[i]+1
    a2 <- plgs2[i+1]-1
    polygon(t(map2[,a1:a2]), col=land.col, border=coast.col)
  }
}
points(bulk$Longitude, bulk$Latitude, pch=21, 
       bg=adjustcolor(reg.col[bulk$Region], .2),
       col=reg.col[bulk$Region], cex=0.8)
text(-82.92, 28.85, 'CK', cex=0.8, col=reg.col[1])
text(-80.85, 24.45, 'FK', cex=0.8, col=reg.col[2])
text(-83.7, 29.2, 'NG', cex=0.8, col=reg.col[3])
text(-82.5, 30, 'Florida', cex=0.8)
text(-84.8, 27, 'Gulf of Mexico', srt=0, cex=0.8)
points(c(-80, -80), c(30, 30.9), type='l', lwd=2.5)
points(-80, 30.8, pch=17, cex=1.5)
points(-80, 30.4, pch=15, cex=2.5, col=wt.col)
points(-80, 30.4, pch='N', cex=1)
points(c(-88, -87), c(24.5, 24.5), type='l', lwd=2.5)
text(-87.5, 24.3, '100 km', cex=0.8)
axis(1, at=seq(-88,-80,2), labels=c(expression(-88~degree),
     expression(-86~degree), expression(-84~degree),
     expression(-82~degree), expression(-80~degree)))
axis(2, las=1, at=24:31, labels=c(expression(24~degree), expression(25~degree),
                                  expression(26~degree), expression(27~degree),
                                  expression(28~degree), expression(29~degree),
                                  expression(30~degree), expression(31~degree)))
mtext(side=1, line=3, expression(Longitude~degree))
mtext(side=2, line=3, expression(Latitude~degree))
box()
dev.off()

# Figure 2 salinity and temperature plots====
pdf('FIG 2 temperature and salinity by region.pdf', width=5, height=5)
plot(bulk$meanT, bulk$sal_med, cex.axis=0.7, type='n', ylim=c(25,38), xlim=c(16,32),
     ylab='salinity (\u2030)', las=1, xlab=bquote('temperature '(degree~'C')))
points(bulk$meanT, bulk$sal_med, pch=21, col=reg.col[bulk$Region], bg=reg.col[bulk$Region])
for (i in 1:nrow(bulk)) points(c(bulk$winter[i], bulk$summer[i]), rep(bulk$sal_med[i], 2),
                               lwd=0.8, type='l', col=reg.col[bulk$Region][i])
for (i in 1:nrow(bulk)) points(rep(bulk$meanT[i], 2), c(bulk$sal_q25[i], bulk$sal_q75[i]),
       lwd=0.8, type='l', col=reg.col[bulk$Region][i])
for (i in 1:nrow(bulk)) {
  trng <- range(c(bulk$summer[i], bulk$fall[i], bulk$spring[i], bulk$winter[i]))
  points(c(trng[1], trng[2]), rep(bulk$sal_med[i], 2),
                               lwd=0.8, type='l', col=reg.col[bulk$Region][i])
}
points(bulk$winter, bulk$sal_med, pch=16, col=reg.col[bulk$Region], cex=0.5)
points(bulk$summer, bulk$sal_med, pch=16, col=reg.col[bulk$Region], cex=0.5)
points(bulk$spring, bulk$sal_med, pch=16, col=reg.col[bulk$Region], cex=0.5)
points(bulk$fall, bulk$sal_med, pch=16, col=reg.col[bulk$Region], cex=0.5)
points(bulk$meanT, bulk$sal_med, pch=21, col=reg.col[bulk$Region], bg=reg.col[bulk$Region])
legend('topleft', pt.bg=reg.col, cex=0.7, pt.cex=0.9, text.col=reg.col, pch=21, levels(bulk$Region))
dev.off()

# NOTE! - salinity and temperature are strongly correlated at the regional level (effective n=3)
cor.test(bulk$sal_med, bulk$meanT)
cor.test(bulk$sal_med, bulk$Depth)
cor.test(bulk$meanT, bulk$Depth)

# BULK ANALYSES ====
# correlation and partial correlations bulk====
# analyses restricted to most relevant elements, elements with values beyond detection limits removed
elratios <- list(Mg=bulk$Mg, Sr=bulk$Sr, Ba=bulk$Ba, S=bulk$S, Li=bulk$Li,
                 Pb=bulk$Pb, Zn=bulk$Zn, Mn=bulk$Mn)
outcorrep <- NULL
for (i in 1:length(elratios)) {
tempcor <- cor.test(bulk$meanT, elratios[[i]])
salcor <- cor.test(bulk$sal_med, elratios[[i]])
depthcor <- cor.test(bulk$Depth, elratios[[i]])
spcors <- spcor(cbind(elratios[[i]], temp=bulk$meanT, sal=bulk$sal_mean, 
                      depth=bulk$Depth))
pres <- cbind(rt=tempcor$estimate, pt=tempcor$p.value,
              rs=salcor$estimate, ps=salcor$p.value,
              rd=depthcor$estimate, pd=depthcor$p.value,
          srt=spcors[[1]][2,1], spt=spcors[[2]][2,1], 
          srs=spcors[[1]][3,1], sps=spcors[[2]][3,1],
          srd=spcors[[1]][4,1], spd=spcors[[2]][4,1])
rownames(pres) <- names(elratios)[i]
outcorrep <- rbind(outcorrep, pres)
}
round(outcorrep, 5)
write.csv(round(outcorrep, 5), 'Table 2 correlations.csv')


# body size====
size.cor <- NULL
for(i in 1:length(elratios)) {
spcorbs <- spcor(na.omit(cbind(elratios[[i]], size=bulk$size, temp=bulk$meanT, sal=bulk$sal_mean,
      depth=bulk$Depth)))
size.cor <-  rbind(size.cor, c(rsp=spcorbs[[1]][1,2], p=spcorbs[[2]][1,2]))
}
rownames(size.cor) <- names(elratios)
size.cor

# salinity versus element ratios plots====
pdf('Figure 3 bulk elements vs salinity plots.pdf', height=9, width=6)
tempar <- par(mfrow=c(length(elratios)/2,2), oma=c(5,1,4,0), mar=c(0,5,1,2))
for (i in 1:length(elratios)) {
pcors <- pcor(cbind(temp=bulk$meanT, sal=bulk$sal_mean, elratios[[i]]))
plot(bulk$sal_med, elratios[[i]], pch=pchreg[bulk$Region], cex=0.9,
     las=1, xlab='', col=mycol[bulk$species],
     bg=adjustcolor(mycol[bulk$species],0.1), axes=F,
     ylab='')
mtext(side=2, line=4, paste0(names(elratios)[i],'/Ca [mmol/mol]'), cex=0.8)
mtext(side=3, line=-1.2, adj=0.05, cex=0.7,
      paste('r =', round(outcorrep[i,3],3)))
mtext(side=3, line=-2, adj=0.05, cex=0.7,
      paste('r* =', round(outcorrep[i,9],3)))
if (i == 1) legend(28, 205, pch=NA, pt.bg=adjustcolor(mycol,0.1),
                   text.col=mycol, col=mycol,
                   levels(bulk$species), xpd=NA, bty='n', x.intersp=0.5)
if (i == 2) legend(28, 3, pch=pchreg, pt.bg=adjustcolor('gray40',0.3),
                   text.col='gray30', col='gray40',
                   levels(bulk$Region), xpd=NA, bty='n', x.intersp=0.5)
# if (i == length(elratios)-1) mtext(side=1, line=3, bquote('temperature [C'*degree*']'))
if (i > length(elratios)-2) mtext(side=1, line=3, 'salinity [\u2030]', cex=0.8)
if (i > length(elratios)-2) axis(1)
axis(2, las=1)
box()
mtext(side=3, line=-1.5, adj=0.5, LETTERS[i], cex=0.8)
}
par(tempar)
dev.off()

# pca analysis Figure 4====
select.vars <- c(18:25)
select.pch <- pchtax
scaling <- T
bubbles <- 0.5 + 0.3*(bulk$sal_med-28)
for (k in 2:3) {
  mypcs <- c(1,k)
  if (k == 2) pdf('Figure 4.pdf', width=6.2, height=5)
  if (k == 3) pdf('Figure PC2 vs PC3 not included.pdf', width=6.2, height=5)
  lay.m <- rbind(c(1,1,1,2), c(1,1,1,3), c(1,1,1,4))
  layout(lay.m)
  tempar <- par(mar=c(0,0,0,0), oma=c(5,5,5,5), xpd=NA)
  outpca <- pca.main.F(bulk[,select.vars], scale=scaling, cex=bubbles, las=1,
                       pch=select.pch[bulk$species], pcs=mypcs,
                       col=reg.col[bulk$Region],
                       bg=adjustcolor(reg.col[bulk$Region], 0.1))
  min.x <- min(outpca$scores[,mypcs[1]])
  max.y <- max(outpca$scores[,mypcs[2]]) + 0.25*diff(range(outpca$scores[,mypcs[2]]))
  legend(min.x, max.y, pch=21, pt.bg=adjustcolor(reg.col, 0.1), col=reg.col,
         levels(bulk$Region), text.col=reg.col, bty='n', xpd=NA)
  legend(min.x + abs(min.x), max.y, pch=select.pch, 
         levels(bulk$species), bty='n', xpd=NA)
  mtext(side=3, line=-1.5, cex=0.9, adj=0.97, 'A')
  pca.vec.F(bulk[,select.vars], pcs=mypcs, scale=scaling, eload=F, lab.offset=1.05)
  mtext(side=3, line=-1.5, cex=0.9, adj=0.97, 'B')
  pca.vec.F(bulk[,select.vars], pcs=mypcs, scale=scaling, eload=T, x.axis.label=F, lab.offset=1.05)
  mtext(side=3, line=-1.5, cex=0.9, adj=0.97, 'C')
  pca.scree.F(bulk[,select.vars], scale=scaling, times=10000)
  mtext(side=3, line=-1.5, cex=0.9, adj=0.97, 'D')
  par(tempar)
  dev.off()
}
outpca <- pca.main.F(bulk[,select.vars], scale=scaling)

# pc scores versus external variables====
cor.test(bulk$size, outpca$scores[,1], na.action=na.omit)
cor.test(bulk$size, outpca$scores[,2], na.action=na.omit)
cor.test(bulk$sal_mean, outpca$scores[,1], na.action=na.omit)
cor.test(bulk$sal_mean, outpca$scores[,2], na.action=na.omit)
cor.test(bulk$meanT, outpca$scores[,1], na.action=na.omit)
cor.test(bulk$meanT, outpca$scores[,2], na.action=na.omit)
cor.test(bulk$Depth, outpca$scores[,1], na.action=na.omit)
cor.test(bulk$Depth, outpca$scores[,2], na.action=na.omit)

# Mg vs. S plot Figure 5====
pdf('Figure 5.pdf', width=6, height=6)
plot(bulk$Mg, bulk$S, pch=pchtax[bulk$species], col=reg.col[bulk$Region],
     bg=adjustcolor(reg.col[bulk$Region],0.1),
     cex=1.2, las=1, xlab='Mg/Ca [mmol/mol]', ylab='S/Ca [mmol/mol]', ylim=c(6.5,13))
legend('bottomright', pch=pchtax, levels(bulk$species), bty='n', cex=0.8)
legend('topleft', pch=21, pt.bg=adjustcolor(reg.col,0.1), bty='n',
       levels(bulk$Region), text.col=reg.col, col=reg.col, cex=0.8)
dev.off()

# NANO-SCALE ANALYSES====

corcf <- round(cor(cbind(mgca=nn$Mg.Ca, H=nn$H)),2)
nnL <- nn[nn$taxon=='Leodia sexiesperforata',]
nnE <- nn[nn$taxon=='Encope michelini',]
corcfL <- round(cor(cbind(mgca=nnL$Mg.Ca, H=nnL$H)),2)
corcfE <- round(cor(cbind(mgca=nnE$Mg.Ca, H=nnE$H)),2)

# NANO correlation tests====
# FIGURE 6 NANO 1====
bytaxon <- list(all=nn, Encope=nn[nn$taxon=='Encope michelini',],
                Leodia=nn[nn$taxon=='Leodia sexiesperforata',])
pdf('figure 6.pdf', width=7, height=5)
tempar <- par(mfrow=c(2,3), mar=c(2,3,1,0), oma=c(3,3,0,1))
for (i in 1:3) {
plot(bytaxon[[i]]$Mg.Ca ~ bytaxon[[i]]$zone, ylim=range(nn$Mg.Ca),
     ylab='', xlab='', las=1)
  outt <- wilcox.test(bytaxon[[i]]$Mg.Ca ~ bytaxon[[i]]$zone)
  if (outt$p.value < 0.001) outt$p.value <- round(outt$p.value, 4)
  if (outt$p.value > 0.001) outt$p.value <- round(outt$p.value, 3)
  if (outt$p.value < 0.0001) outt$p.value <- '< 0.0001'
  mtext(side=3, line=-1.5, xpd=NA, cex=0.7, adj=0.05,
      bquote(italic('p')==.(outt$p.value)))
  mtext(side=3, line=-1.5, adj=0.98, LETTERS[i], cex=0.8)
  if (i == 1) mtext(side=2, line=3.5, xpd=NA, 'Mg/Ca [mmol/mol]')
}
for (i in 1:3) {
plot(bytaxon[[i]]$H ~ bytaxon[[i]]$zone, ylim=range(nn$H),
     ylab='H', xlab='', las=1)
  outt <- wilcox.test(bytaxon[[i]]$H ~ bytaxon[[i]]$zone)
  if (outt$p.value < 0.001) outt$p.value <- round(outt$p.value, 4)
  if (outt$p.value > 0.001) outt$p.value <- round(outt$p.value, 3)
  if (outt$p.value < 0.0001) outt$p.value <- '< 0.0001'
  mtext(side=3, line=-1.5, xpd=NA, cex=0.7, adj=0.05,
      bquote(italic('p')==.(outt$p.value)))
  mtext(side=3, line=-1.5, adj=0.98, LETTERS[i+3], cex=0.8)
  if (i == 1) mtext(side=2, line=3.5, xpd=NA, 'H')
  if (i == 1) mtext(side=1, line=3, 'pooled data')
  if (i == 2) mtext(side=1, line=3, 'Encope michelini', font=3)
  if (i == 3) mtext(side=1, line=3, 'Leodia sexiesperforata', font=3)
}
par(tempar)
dev.off()

# FIGURE NANO 7====
pdf('Figure 7 RESULTS bivariate plots of nano variables.pdf', width=7, height=5)
plot(nn$Mg.Ca, nn$H, pch=pchtax[2:3][nn$taxon], cex=1.2, ylab='H', xlab='Mg/Ca', las=1,
     bg=adjustcolor(mycol[2:3][nn$taxon],.3), col=mycol[2:3][nn$taxon])
legend('bottomright', pch=pchtax[2:3], levels(nn$taxon), cex=0.8, 
       pt.bg=adjustcolor(mycol[2:3],.3), col=mycol[2:3], bty='n', text.font=3)
mtext(side=3, line=-1.5, adj=0.03, paste('r (all) =', corcf[1,2]), cex=0.6)
mtext(side=3, line=-2.5, adj=0.03, paste('r (L.s) =', corcfL[1,2]), cex=0.6)
mtext(side=3, line=-3.5, adj=0.03, paste('r (E.m) =', corcfE[1,2]), cex=0.6)
points(nn$Mg.Ca[nn$Region == 'Northern Gulf'], nn$H[nn$Region == 'Northern Gulf'],
       pch=16, cex=0.6)
dev.off()

#stat descriptors and wilcoxon test====
tapply(nn$Mg.Ca, nn$zone, mystats)
wilcox.test(nn$Mg.Ca ~ nn$zone)
tapply(nn$H, nn$zone, mystats)
wilcox.test(nn$H ~ nn$zone)
wilcox.test(nn$Mg.Ca ~ nn$taxon)
wilcox.test(nn$H ~ nn$taxon)

# figure 8 and table 3
byzoneMg <- matrix(unlist(tapply(nn$Mg.Ca, nn$zone, mystats)),2,4, byrow=T)
bytaxonMg <- matrix(unlist(tapply(nn$Mg.Ca, nn$taxon, mystats)),2,4, byrow=T)
bytypeMg <- as.data.frame(tapply(nn$Mg.Ca, list(nn$taxon, nn$zone), mystats))
bytypeMg <- matrix(unlist(bytypeMg), 4, 4, byrow=T)
byzoneH <- matrix(unlist(tapply(nn$H, nn$zone, mystats)),2,4, byrow=T)
bytaxonH <- matrix(unlist(tapply(nn$H, nn$taxon, mystats)),2,4, byrow=T)
bytypeH <- as.data.frame(tapply(nn$H, list(nn$taxon, nn$zone), mystats))
bytypeH <- matrix(unlist(bytypeH), 4, 4, byrow=T)
pooled <- rbind(byzoneMg, bytaxonMg, bytypeMg, byzoneH, bytaxonH, bytypeH)
colnames(pooled) <- c('median', 'mean', 'std.dev', 'n')
Groups <- rep(c("Inner (all data)", "Outer (all data)",
            "E. michelini (all data)",
            "L. sexiesperforata (all data)",
            "E. michelini (inner)",
            "L. sexiesperforata (inner)",
            "E. michelini (outer)",
            "L. sexiesperforata (outer)"),2)
Variable <- c(rep('Mg/Ca', 8), rep('H', 8))
meds <- data.frame(Groups, Variable, pooled)
write.csv(meds, 'table3.csv')

# Figure Nano 8====

medlist <- split(meds, meds$Variable)
my.lwd <- 1
my.arr <- 0.15
q25 <- c(tapply(nn$Mg.Ca, paste(nn$zone, nn$taxon), quantile, prob=0.25),
         tapply(nn$H, paste(nn$zone, nn$taxon), quantile, prob=0.25))
q75 <- c(tapply(nn$Mg.Ca, paste(nn$zone, nn$taxon), quantile, prob=0.75),
         tapply(nn$H, paste(nn$zone, nn$taxon), quantile, prob=0.75))
q50 <- c(tapply(nn$Mg.Ca, paste(nn$zone, nn$taxon), quantile, prob=0.5),
         tapply(nn$H, paste(nn$zone, nn$taxon), quantile, prob=0.5))
mCIMg <- tapply(nn$Mg.Ca, paste(nn$zone, nn$taxon), MedianCI,
                method='boot', R=9999, conf.level=0.95)
mCIMg <- tapply(nn$Mg.Ca, paste(nn$zone, nn$taxon), MedianCI,
                conf.level=0.95)
mCIH <- tapply(nn$H, paste(nn$zone, nn$taxon), MedianCI,
               method='boot', R=9999, conf.level=0.95)
pdf('Figure 8 RESULTS taxon-zone comparisons for nano variables.pdf', width=6, height=5)
plot(medlist[[2]]$median, medlist[[1]]$median, type='n',
     las=1, xlab='Mg/Ca [mmol/mol]', ylab='H',
     xlim=c(0.14, 0.22), ylim=c(3.5, 4.8))
arrows(1.008*medlist[[2]]$median[5], 1.01*medlist[[1]]$median[5], 
       0.992*medlist[[2]]$median[7], 0.99*medlist[[1]]$median[7],
       lwd=my.lwd, length=my.arr, col='gray50')
arrows(1.005*medlist[[2]]$median[6], 1.01*medlist[[1]]$median[6], 
       0.995*medlist[[2]]$median[8], 0.99*medlist[[1]]$median[8],
       lwd=my.lwd, length=my.arr, col='gray50')
text(0.157, 4.3, 'E. michelini', pos=1, font=3, cex=0.8)
text(0.21, 4.2, 'L. sexiesperforata', pos=1, font=3, cex=0.8)
text(0.195, 4.75, 'outer', pos=1, cex=0.8)
text(0.17, 4, 'inner', pos=1, cex=0.8)

for (i in 1:4) {
  points(cbind(mCIMg[[i]][2:3], mCIH[[i]][1]), type='l',
         lwd=1, col='gray50')
  points(cbind(mCIMg[[i]][1], mCIH[[i]][2:3]), type='l',
         lwd=1, col='gray50')
  points(cbind(mCIMg[[i]][1], mCIH[[i]][1]), pch=16,
         lwd=1, col='gray50')
}
for (i in 1:4) {
  points(c(q25[i], q75[i]), rep(q50[i+4], 2), lwd=2, type='l')
  points(rep(q50[i], 2), c(q25[i+4], q75[i+4]), lwd=2, type='l')
}
points(medlist[[2]]$median[5:8], medlist[[1]]$median[5:8], pch=21,
       bg='white')
dev.off()

# table 4====
corlist <- list(nn, nnE, nnL)
namelist <- c('all', 'E. michelini', 'L.sexiesperforata')
mymethod <- 'pearson' 
OUT1 <- NULL
for (i in 1:3) {
  MgH <- cor.test(corlist[[i]]$Mg.Ca, corlist[[i]]$H, method=mymethod)
  MgHt <- cbind(r=MgH$estimate, "lower limit"=MgH$conf.int[1],
                "upper limit"=MgH$conf.int[2], t=MgH$parameter, p=MgH$p.value)
  print(MgHt)
  finalout <- rbind(MgHt)
  print(finalout)
  rownames(finalout) <- paste(namelist[i], c('Mg vs H'))
  OUT1 <- rbind(OUT1, finalout)
}
OUT1
write.csv(OUT1, 'table4.csv')
