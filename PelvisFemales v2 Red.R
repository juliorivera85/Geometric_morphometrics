#load packages needed
library(PCDimension)
library(geiger)
library(mvMORPH)
library(geomorph)
library(car)
library(ips)
library(scatterplot3d)
library(phytools) 

# set the directory to read in the *.nts files
setwd("C:/Users/emart120/ASU Dropbox/Martins Lab/Research/CT Project/CT People/Logan Kenny (CT Scans)/analysis/nts/25 landmark females (1 per species)")
setwd("~/Library/CloudStorage/Dropbox/CT Project/CT People/Logan Kenny (CT Scans)/analysis/nts/25 landmark females (1 per species)")
filelistALL <- list.files(pattern=".nts")
datALL <- readmulti.nts(filelistALL)

##########################################################################
### Run a Generalized Procrustes Analysis
###   This creates a 3D description of shape from the raw coordinate data
##########################################################################

Y.gpa <- gpagen(datALL, Proj=TRUE, ProcD=TRUE, curves=NULL, surface=NULL)
summary(Y.gpa)
plot(Y.gpa)
y <- two.d.array(Y.gpa$coords)

for (i in 1:length(attributes(Y.gpa$coords)$dimnames[[3]])) {
  temp <- strsplit(attributes(Y.gpa$coords)$dimnames[[3]][i], "_")
  attributes(Y.gpa$coords)$dimnames[[3]][i] <- paste("Sceloporus_",temp[[1]][2], sep="")
}

attributes(Y.gpa$coords)$dimnames[[3]]

attributes(Y.gpa$Csize)$names <- attributes(Y.gpa$coords)$dimnames[[3]]

##creating a geomorph dataframe for downstream analyses
gfd <- geomorph.data.frame(Y.gpa)

## Checking outliers for possible landmarking problems
outliers <- plotOutliers(Y.gpa$coords) 
M <- mshape(Y.gpa$coords)   

##########################################################################
######## identify comparison groups ##################
##########################################################################

setwd("C:/Users/emart120/ASU Dropbox/Martins Lab/Research/CT Project/CT People/Logan Kenny (CT Scans)/analysis/")
setwd("~/Library/CloudStorage/Dropbox/CT Project/CT People/Logan Kenny (CT Scans)/analysis/")
group <- read.csv("2DScelData25.csv")

### get data only for species for which we have CT scans 
filelist2 <- gsub(".nts", "", attributes(Y.gpa$Csize)$names)
vars <- group$species %in% filelist2
group2 <- group$species[vars]
group2 <- as.data.frame(group2)
names(group2) <- "species" 
group2 <- group

## S. taeniocnemis is likely to be viviparous
#group2$parity[72] <- parity[72]<-"v"
#group2$parity[73] <- parity[73]<-"v"
group2$parity[74] <- parity[72]<-"v"
group2$parity[75] <- parity[73]<-"v"

parity <- group2$parity <- as.factor(group2$parity)
species <- group2$species <- as.character(group2$species)
arb <- group2$arb <- as.factor(group2$arb)
size <- group2$SVL
catsize <- group2$catsize <- as.factor(group2$catsize)

## associating the categorical data with species it belongs to
names(parity) <- group2$species; parity
names(size) <- group2$species; size
names(catsize) <- group2$species; catsize
names(arb) <- group2$species; arb

gfd$parity <- parity 
gfd$svl <- size
gfd$catsize <- catsize
gfd$arb <- arb

gfd$species <- attributes(gfd$Csize)$name

## calculate morphological diparity between groups
morphol.disparity(f1=gfd$coords ~ 1, groups= ~ parity, data=gfd, iter=999) 
morphol.disparity(f1=gfd$coords ~ 1, groups= ~ arb, data=gfd, iter=999) 

#############################
######## Make Figure 4
# egg=black, live=green, terrestrial = brown/oraange, arboreal = purple
#############################

#PANCOVA Plotting of parity & height (Fig. 4A)
groupV<- subset(group,group$parity=="v")
head(groupV)
groupO<- subset(group,group$parity=="o")
head(groupO)
plot(groupO$SVL, groupO$H, type="n", xlab="SVL (mm)", ylab="Pelvis height (mm)", cex.lab = 1.5,
   ylim = c(0,6), xlim=c(40,110))

#create shading & plot points
#newx.v <- seq(0, max(groupV$SVL), length.out=100)
newx.v <- seq(-15, 120, length.out=100)
pred.H.v = predict(lm(groupV$H ~ groupV$SVL))
#newy.v <- seq(1, max(pred.H.v), length.out=110)
newy.v <- seq(0, 10.2, length.out=100)
y.hi.v <- newy.v + 0.28
y.lo.v <- newy.v - 0.28
polygon(c(newx.v, rev(newx.v)), c(y.lo.v, rev(y.hi.v)), col="lightgreen", border=NA, density = 100)
points(groupV$SVL,groupV$H,pch=16,col="darkolivegreen3", cex=1.25)
abline(lm(groupV$H ~ groupV$SVL),col="darkolivegreen3", lwd=2)

#create shading & plot points
#newx <- seq(0, max(groupO$SVL), length.out=100)
newx <- seq(-15, 120, length.out=100)
pred.H = predict(lm(groupO$H ~ groupO$SVL))
newy <- seq(-1.5, 12.3, length.out=100)
#newy <- seq(0, max(pred.H), length.out=100)
y.hi <- newy + 0.28
y.lo <- newy - 0.28
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="darkgray", border=NA, density = 100)
points(groupO$SVL,groupO$H,pch=16,col="black", cex=1.25)
abline(lm(groupO$H~groupO$SVL),col="black", lwd=2)

#######################################

#PANCOVA Plotting of parity & length
groupV<- subset(group,group$parity=="v")
head(groupV)
groupO<- subset(group,group$parity=="o")
head(groupO)
plot(groupO$SVL,groupO$L, type="n", xlab="SVL (mm)", ylab="Pelvis length (mm)", cex.lab = 1.5,
   ylim = c(4,20), xlim=c(40,110))

#create shading & plot points
#newx <- seq(0, max(groupV$SVL), length.out=100)
newx <- seq(-17, 120, length.out=100)
pred.L = predict(lm(groupV$L ~ groupV$SVL))
newy <- seq(-4, 20, length.out=100)
y.hi <- newy + 0.84
y.lo <- newy - 0.84
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="lightgreen", border=NA, density = 80)
points(groupV$SVL,groupV$L,pch=16,col="darkolivegreen3", cex=1.25)
abline(lm(groupV$L~groupV$SVL),col="darkolivegreen3", lwd=2)

#create shading & plot points
#newx.t <- seq(0, max(groupT$SVL), length.out=100)
newx.t <- seq(0, 120, length.out=100)
pred.L.t = predict(lm(groupO$L ~ groupO$SVL))
newy.t <- seq(-2, 22, length.out=100)
#newy.t <- seq(0, max(pred.H.t), length.out=100)
y.hi.t <- newy.t + 0.84
y.lo.t <- newy.t - 0.84
polygon(c(newx.t, rev(newx.t)), c(y.lo.t, rev(y.hi.t)), col="darkgray", border=NA, density = 80)
points(groupO$SVL,groupO$L,pch=16,col="black", cex=1.25)
abline(lm(groupO$L ~ groupO$SVL),col="black", lwd=2)

##########################################

#PANCOVA Plotting of parity & width !!!!!!!!!!!!!!!!!!!!!!!!!!!!
groupV<- subset(group,group$parity=="v")
head(groupV)
groupO<- subset(group,group$parity=="o")
head(groupO)
plot(groupO$SVL,groupO$L, type="n", xlab="SVL (mm)", ylab="Pelvis length (mm)", cex.lab = 1.5,
   ylim = c(4,20), xlim=c(40,110))

#create shading & plot points
#newx <- seq(0, max(groupV$SVL), length.out=100)
newx <- seq(-17, 120, length.out=100)
pred.L = predict(lm(groupV$L ~ groupV$SVL))
newy <- seq(-4, 20, length.out=100)
y.hi <- newy + 0.84
y.lo <- newy - 0.84
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="lightgreen", border=NA, density = 80)
points(groupV$SVL,groupV$L,pch=16,col="darkolivegreen3", cex=1.25)
abline(lm(groupV$L~groupV$SVL),col="darkolivegreen3", lwd=2)

#create shading & plot points
#newx.t <- seq(0, max(groupT$SVL), length.out=100)
newx.t <- seq(0, 120, length.out=100)
pred.L.t = predict(lm(groupO$L ~ groupO$SVL))
newy.t <- seq(-2, 22, length.out=100)
#newy.t <- seq(0, max(pred.H.t), length.out=100)
y.hi.t <- newy.t + 0.84
y.lo.t <- newy.t - 0.84
polygon(c(newx.t, rev(newx.t)), c(y.lo.t, rev(y.hi.t)), col="darkgray", border=NA, density = 80)
points(groupO$SVL,groupO$L,pch=16,col="black", cex=1.25)
abline(lm(groupO$L ~ groupO$SVL),col="black", lwd=2)

##########################################
##########################################
##########################################

#PANCOVA Plotting of arboreality & height (Fig. 4C)
groupT<- subset(group,group$arb==0)
head(groupT)
groupA<- subset(group,group$arb==1)
head(groupA)
plot(groupA$SVL,groupA$H, type="n", xlab="SVL (mm)", ylab="Pelvis height (mm)", cex.lab = 1.5,
   ylim = c(4,13), xlim=c(40,110))

#create shading & plot points
#newx <- seq(0, max(groupA$SVL), length.out=100)
newx <- seq(1, 120, length.out=100)
pred.H = predict(lm(groupA$H ~ groupA$SVL))
newy <- seq(-0.7, 12.5, length.out=100)
#newy <- seq(0, max(pred.H), length.out=100)
y.hi <- newy + 0.68
y.lo <- newy - 0.68
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="olivedrab3", border=NA, density = 80)
points(groupA$SVL,groupA$H,pch=16,col="darkgreen", cex=1.25)
abline(lm(groupA$H~groupA$SVL),col="darkgreen", lwd=2)

#create shading & plot points
#newx.t <- seq(0, max(groupT$SVL), length.out=100)
newx.t <- seq(0, 120, length.out=100)
pred.H.t = predict(lm(groupT$H ~ groupT$SVL))
newy.t <- seq(1.8, 10.5, length.out=100)
#newy.t <- seq(0, max(pred.H.t), length.out=100)
y.hi.t <- newy.t + 0.68
y.lo.t <- newy.t - 0.68
polygon(c(newx.t, rev(newx.t)), c(y.lo.t, rev(y.hi.t)), col="orange2", border=NA, density = 80)
points(groupT$SVL,groupT$H,pch=16,col="brown", cex=1.25)
abline(lm(groupT$H ~ groupT$SVL),col="brown", lwd=2)

##########################################


#PANCOVA Plotting of arboreality & length  (Fig. 4B)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
groupT<- subset(group,group$arb==0)
head(groupT)
groupA<- subset(group,group$arb==1)
head(groupA)
plot(groupA$SVL,groupA$W, type="n", xlab="SVL (mm)", ylab="Pelvis width (mm)", cex.lab = 1.5,
   ylim = c(4,20), xlim=c(40,110))

#create shading & plot points
#newx <- seq(0, max(groupA$SVL), length.out=100)
newx <- seq(10, 120, length.out=100)
pred.H = predict(lm(groupA$W ~ groupA$SVL))
newy <- seq(1, 17.5, length.out=100)
#newy <- seq(0, max(pred.H), length.out=100)
y.hi <- newy + 1.03
y.lo <- newy - 1.03
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="olivedrab3", border=NA, density = 80)
points(groupA$SVL,groupA$W,pch=16,col="darkgreen", cex=1.25)
abline(lm(groupA$W~groupA$SVL),col="darkgreen", lwd=2)

#create shading & plot points
#newx.t <- seq(0, max(groupT$SVL), length.out=100)
newx.t <- seq(0, 120, length.out=100)
pred.H.t = predict(lm(groupT$W ~ groupT$SVL))
newy.t <- seq(-3, 19.8, length.out=100)
#newy.t <- seq(0, max(pred.H.t), length.out=100)
y.hi.t <- newy.t + 0.35
y.lo.t <- newy.t - 0.35
polygon(c(newx.t, rev(newx.t)), c(y.lo.t, rev(y.hi.t)), col="orange2", border=NA, density = 80)
points(groupT$SVL,groupT$W,pch=16,col="brown", cex=1.25)
abline(lm(groupT$W ~ groupT$SVL),col="brown", lwd=2)


##########################################

#PANCOVA Plotting of arboreality & width  (Fig. 4B)
groupT<- subset(group,group$arb==0)
head(groupT)
groupA<- subset(group,group$arb==1)
head(groupA)
plot(groupA$SVL,groupA$W, type="n", xlab="SVL (mm)", ylab="Pelvis width (mm)", cex.lab = 1.5,
   ylim = c(4,20), xlim=c(40,110))

#create shading & plot points
#newx <- seq(0, max(groupA$SVL), length.out=100)
newx <- seq(10, 120, length.out=100)
pred.H = predict(lm(groupA$W ~ groupA$SVL))
newy <- seq(1, 17.5, length.out=100)
#newy <- seq(0, max(pred.H), length.out=100)
y.hi <- newy + 1.03
y.lo <- newy - 1.03
polygon(c(newx, rev(newx)), c(y.lo, rev(y.hi)), col="olivedrab3", border=NA, density = 80)
points(groupA$SVL,groupA$W,pch=16,col="darkgreen", cex=1.25)
abline(lm(groupA$W~groupA$SVL),col="darkgreen", lwd=2)

#create shading & plot points
#newx.t <- seq(0, max(groupT$SVL), length.out=100)
newx.t <- seq(0, 120, length.out=100)
pred.H.t = predict(lm(groupT$W ~ groupT$SVL))
newy.t <- seq(-3, 19.8, length.out=100)
#newy.t <- seq(0, max(pred.H.t), length.out=100)
y.hi.t <- newy.t + 0.35
y.lo.t <- newy.t - 0.35
polygon(c(newx.t, rev(newx.t)), c(y.lo.t, rev(y.hi.t)), col="orange2", border=NA, density = 80)
points(groupT$SVL,groupT$W,pch=16,col="brown", cex=1.25)
abline(lm(groupT$W ~ groupT$SVL),col="brown", lwd=2)


##########################################



##########################################################################
####### non-phylo analyses 
##########################################################################

#############################
####### PCA  
### see https://rdrr.io/cran/geomorph/f/vignettes/geomorph.PCA.Rmd for tutorial
#############################

par(mfrow=c(1,1)) # 1 graph at a time

pca.lands <- gm.prcomp(A=Y.gpa$coords) 
summary(pca.lands) ## % variance explained is here
plot(pca.lands,pch=21,cex=1.5,bg="black")

### plot of min and max differences in PC1
plotRefToTarget(pca.lands$shapes$shapes.comp1$min,
  pca.lands$shapes$shapes.comp1$max, method="vector")

### plot of min and max differences in PC2
plotRefToTarget(pca.lands$shapes$shapes.comp2$min,
  pca.lands$shapes$shapes.comp2$max,method="vector")

### plot of min and max differences in PC3
plotRefToTarget(pca.lands$shapes$shapes.comp3$min,
  pca.lands$shapes$shapes.comp3$max,method="vector")

### plot of min and max differences in PC4
plotRefToTarget(pca.lands$shapes$shapes.comp4$min,
  pca.lands$shapes$shapes.comp4$max,method="vector")

### plot of min and max differences in PC5
plotRefToTarget(pca.lands$shapes$shapes.comp5$min,
  pca.lands$shapes$shapes.comp5$max,method="vector")
  
  
#PCA with tree
pca.lands.phy <- gm.prcomp(A=Y.gpa$coords, phy=tree) 
summary(pca.lands.phy) ## % variance explained is here
plot(pca.lands.phy, pch=21, cex=1, bg="black", phylo=TRUE)


#######################################
#### Figure 2. Parity and PC Graph ####
#######################################
par(mfrow=c(1,1)) # 1 graph at a time
par(mar=c(5, 5, 3, 3)) # add margins
colors <- c("green", "black") # mark parity
plot(pca.lands, pch=22, cex = 1.5, bg = gfd$parity, 
     ylab="PC2 (14.1%)", xlab="PC1 (21.1%)", main = "Pelvis Shape") 
legend("topright", c("live-bearers", "egg-layers"), pch=22, pt.bg = unique(colors))
text(pca.lands$x[,1:2],rownames(pca.lands$x[,1:2]),pos=3,cex=.5)

##### procrustes regression for shape, basically an ANOVA 
#parity and habitat
all.aov <- procD.lm(gfd$coords ~ species + gfd$parity * gfd$arb * gfd$svl, iter=999) 
all.aov <- procD.lm(gfd$coords ~ gfd$parity * gfd$arb * gfd$svl, iter=999) 
summary(all.aov)

##########################
##########################
####### PHYLOGENY ########
##########################
##########################

## this code reads in the tree and matches the species on the tree with the species
## we have skull data for and drops the ones we don't need.
tree <- read.nexus("scelop_timecalib_2016.nex")
par(mfrow=c(1,1)) # 1 graph at a time
plot(tree)
tree$tip.label[[82]]<-"Sceloporus_prezygus"  # changing serrifer to prezygus
tree$tip.label[[116]]<-"Sceloporus_taeniocnemis"  # changing taeniocnemis_mvz4213 to taeniocnemis
tree$tip.label[[117]]<-"Sceloporus_acanthinus"  # changing acanthinus_anmo1932 to acanthinus
tree$tip.label[[125]]<-"Sceloporus_cowlesi"  # changing cowelsi to cowlesi
tree$tip.label[[112]]<-"Sceloporus_formosus"  # changing formosus_anmo1248 to formosus
tree$tip.label[[66]]<-"Sceloporus_magister"  # changing magister_magister_rom14488 to magister
tree$tip.label[[128]]<-"Sceloporus_occidentalis"  # changing occidentalis_uwbm6281 to occidentalis
tree$tip.label[[72]]<-"Sceloporus_orcutti"  # changing orcutti_uwbm7654 to orcutti
tree$tip.label[[42]]<-"Sceloporus_pyrocephalus"  # changing pyrocephalus_utar53473 to pyrocephalus
#tree$tip.label[[106]]<-"Sceloporus_scitulus"  # changing formosus_scitulus to scitulus
tree$tip.label[[38]]<-"Sceloporus_siniferus"  # changing siniferus_mvz236299 to siniferus
tree$tip.label[[91]]<-"Sceloporus_torquatus"  # changing torquatus_uwbm6600 to torquatus
tree$tip.label[[63]]<-"Sceloporus_zosteromus"  # changing zosteromus_ADG74 to zosteromus
tree$tip.label[[53]]<-"Sceloporus_scalaris"  # changing scalaris_rwb06247 to scalaris)

tree.tps <- tree$tip
#find species that we don't have data for and remove the tips
To.Drop <-  tree.tps[!tree.tps %in% group2$species]
tree2 <- drop.tip(tree, To.Drop)
plotTree(tree2)
edgelabels()
#find species we have shape data for but not in the tree
group2$species[!group2$species %in% tree.tps]

# add cozumelae, sister to variabilis
node <- which(tree2$tip.label=="Sceloporus_variabilis"); node
plot(tree2)
edgelabels()  ### find edge number below variabilis = 6
tmp<-tree2$edge.length[[6]]; tmp<-tmp/2; tmp  ## place it half way down
tree3 <- bind.tip(tree2, tip.label="Sceloporus_cozumelae", where=node, position=tmp, edge.length=tmp)
plotTree(tree3)
edgelabels()

# add hartwegi as sister to taeniocnemis
node <- which(tree3$tip.label=="Sceloporus_taeniocnemis"); node
plot(tree3); edgelabels(cex=0.3)  ### find edge number below taeniocnemis = 87
tmp<-tree3$edge.length[[87]]; tmp<-tmp/2; tmp  ## place it half way down
tree4 <- bind.tip(tree3, tip.label="Sceloporus_hartwegi", where=node, position=tmp, edge.length=tmp)
plotTree(tree4)
edgelabels()

# add utiformis as sister to grandeveus
node <- which(tree4$tip.label=="Sceloporus_grandaevus"); node
plot(tree4); edgelabels()  ### find edge number below grandaevus = 15
tmp<-tree4$edge.length[[15]]; tmp<-tmp/2; tmp  ## place it half way down
tree5 <- bind.tip(tree4, tip.label="Sceloporus_utiformis", where=node, position=tmp, edge.length=tmp)
plotTree(tree5)
edgelabels()

# add tristichus as sister to cowlesi
node <- which(tree5$tip.label=="Sceloporus_cowlesi"); node
plot(tree5); edgelabels()  ### find edge number below cowlesi = 89
tmp<-tree5$edge.length[[102]]; tmp<-tmp/2; tmp  ## place it half way down
tree6 <- bind.tip(tree5, tip.label="Sceloporus_tristichus", where=node, position=tmp, edge.length=tmp)
plotTree(tree6)

#write.tree(tree6, file="tree6")

## checking names to make sure data matches tree
names(Y.gpa$Csize) <- gsub(pattern = ".nts_1", replacement = "", x = names(Y.gpa$Csize))
dimnames(Y.gpa$coords)[[3]]<-gsub(pattern = ".nts_1", replacement = "", x = dimnames(Y.gpa$coords)[[3]])
l1<-sort(names(Y.gpa$Csize))
l2<-sort(tree6$tip.label)
identical(l1,l2)
setdiff(l1,l2)


####################################################
# 2D distances
####################################################

#PLhypo = 1 and 7; Used for SLOUCH Length
#PLnohypo = 1 and 6; Used for SLOUCH Length (to check results)
#PLilium1 = 20 to 22
#PLilium2 = 21 to 23
#PLischiumT1 = 14 to 16
#PLischiumT2 = 15 to 17
#PLischiumLnohypo = 3 to 6
#PLischiumLhypo = 3 to 7

#PWpubis = 10 to 11
#PWback = 18 to 19
#PWilium = 22 to 23; Used for SLOUCH Width
#PWischium = 16 to 17

#PHdvf1 = 16 to 20
#PHdvf2 = 17 to 21
##PHdvf = average; Used for SLOUCH Height
#PHes1 = 10 to 24
#PHes2 = 11 to 25

#PGAT12.1 = 13 to 15
#PGAT12.2 = 12 to 14
#PGAT13.1 = 3 to 13
#PGAT13.2 = 3 to 12
#PGAT2.1 = 3 to 15
#PGAT2.2 = 3 to 14
#PGAT34 = 2 to 3
#PGAT4.1 = 2 to 12
#PGAT4.2 = 2 to 13

lmks <- matrix(c(1,7, 1,6, 20,22, 21,23, 14,16, 15,17, 3,6, 3,7, 10,11, 18,19, 22,23, 
          16,17, 16,20, 17,21, 10,24, 11,25, 13,15, 12,14, 3,13, 3,12, 3,15, 3,14, 
          2,3, 2,12, 2,13, 16,18), ncol=2, byrow=TRUE, 
          dimnames = list(c("PLhypo", "PLnohypo", "PLilium1", "PLilium2", "PLischiumT1", "PLischiumT2",
            "PLischiumLnohypo", "PLischiumLhypo", "PWpubis", "PWback", "PWilium", "PWischium",
            "PHdvf1","PHdvf2","PHes1","PHes2","PGAT12.1","PGAT12.2","PGAT13.1","PGAT13.2",
            "PGAT2.1","PGAT2.2","PGAT34","PGAT4.1","PGAT4.2","Harb"), c("start", "end")))
#TwoD <- interlmkdist(dat, lmks)
TwoD <- interlmkdist(datALL, lmks)
TwoD <- data.frame(TwoD)

rownames(TwoD)<- c("Sceloporus_acanthinus", "Sceloporus_adleri", "Sceloporus_aeneus", "Sceloporus_angustus", "Sceloporus_arenicolus", "Sceloporus_bicanthalis",  "Sceloporus_carinatus", "Sceloporus_chrysostictus", "Sceloporus_clarkii", "Sceloporus_consobrinus", "Sceloporus_couchii", "Sceloporus_cowlesi", "Sceloporus_cozumelae", "Sceloporus_cyanogenys", "Sceloporus_dugesii", "Sceloporus_edwardtaylori", "Sceloporus_formosus", "Sceloporus_formosus_scitulus", "Sceloporus_graciosus", "Sceloporus_grammicus", "Sceloporus_grandaevus", "Sceloporus_hartwegi", "Sceloporus_horridus", "Sceloporus_hunsakeri", "Sceloporus_internasalis", "Sceloporus_jalapae", "Sceloporus_jarrovii", "Sceloporus_licki", "Sceloporus_magister", "Sceloporus_malachiticus", "Sceloporus_megalepidurus", "Sceloporus_melanorhinus", "Sceloporus_merriami", "Sceloporus_minor", "Sceloporus_mucronatus", "Sceloporus_nelsoni", "Sceloporus_occidentalis", "Sceloporus_olivaceus", "Sceloporus_orcutti", "Sceloporus_parvus", "Sceloporus_poinsettii", "Sceloporus_prezygus", "Sceloporus_pyrocephalus", "Sceloporus_scalaris", "Sceloporus_siniferus", "Sceloporus_spinosus", "Sceloporus_squamosus", "Sceloporus_taeniocnemis", "Sceloporus_teapensis", "Sceloporus_tristichus", "Sceloporus_undulatus", "Sceloporus_utiformis", "Sceloporus_variabilis", "Sceloporus_virgatus", "Sceloporus_woodi", "Sceloporus_zosteromus")

#Lengths
gfd$PLhypo <- TwoD$PLhypo

#Widths
gfd$PWilium <- TwoD$PWilium

#Heights
gfd$PHdvf <- ((TwoD$PHdvf1 + TwoD$PHdvf2)/2)

##########################################
### 2D Adjustments for PC Graphs Below ###
##########################################

L <- gfd$PLhypo
W <- gfd$PWilium
H <- gfd$PHdvf
names(L) <- names(W) <- names(H) <- gfd$species

Clade1 = c("Sceloporus_grammicus", "Sceloporus_megalepidurus", "Sceloporus_jarrovii", "Sceloporus_cyanogenys", 
            "Sceloporus_poinsettii", "Sceloporus_prezygus", "Sceloporus_minor", "Sceloporus_dugesii", 
            "Sceloporus_mucronatus")

Clade2 <- c("Sceloporus_internasalis", "Sceloporus_acanthinus", "Sceloporus_hartwegi", 
            "Sceloporus_taeniocnemis", "Sceloporus_malachiticus", "Sceloporus_formosus", 
            "Sceloporus_adleri", "Sceloporus_formosus scitulus")

Clade3 <- "Sceloporus_bicanthalis"

Eggs <- names(L[!(names(L) %in% Clade1)])
Eggs <- Eggs[!(Eggs %in% Clade2)]
Eggs <- Eggs[!(Eggs %in% Clade3)]

mean(L[which(names(L) %in% Clade1)])
mean(W[which(names(W) %in% Clade1)])
mean(H[which(names(H) %in% Clade1)])
mean(L[which(names(L) %in% Clade2)])
mean(W[which(names(W) %in% Clade2)])
mean(H[which(names(H) %in% Clade2)])
mean(L[which(names(L) %in% Eggs)])
mean(W[which(names(W) %in% Eggs)])
mean(H[which(names(H) %in% Eggs)])

sd(L[which(names(L) %in% Clade1)])/sqrt(length(Clade1))
sd(W[which(names(W) %in% Clade1)])/sqrt(length(Clade1))
sd(H[which(names(H) %in% Clade1)])/sqrt(length(Clade1))
sd(L[which(names(L) %in% Clade2)])/sqrt(length(Clade2))
sd(W[which(names(W) %in% Clade2)])/sqrt(length(Clade2))
sd(H[which(names(H) %in% Clade2)])/sqrt(length(Clade2))
sd(L[which(names(L) %in% Eggs)])/sqrt(length(Eggs))
sd(W[which(names(W) %in% Eggs)])/sqrt(length(Eggs))
sd(H[which(names(H) %in% Eggs)])/sqrt(length(Eggs))


Ha<-(TwoD$Harb/size)

sL <- (L/size)
sW <- (W/size)
sH <- (H/size)
mean(sL[which(names(sL) %in% Clade1)])
mean(sW[which(names(sW) %in% Clade1)])
mean(sH[which(names(sH) %in% Clade1)])
mean(sL[which(names(sL) %in% Clade2)])
mean(sW[which(names(sW) %in% Clade2)])
mean(sH[which(names(sH) %in% Clade2)])
mean(sL[which(names(sL) %in% Eggs)])
mean(sW[which(names(sW) %in% Eggs)])
mean(sH[which(names(sH) %in% Eggs)])

sd(sL[which(names(sL) %in% Clade1)])/sqrt(length(Clade1))
sd(sW[which(names(sW) %in% Clade1)])/sqrt(length(Clade1))
sd(sH[which(names(sH) %in% Clade1)])/sqrt(length(Clade1))
sd(sL[which(names(sL) %in% Clade2)])/sqrt(length(Clade2))
sd(sW[which(names(sW) %in% Clade2)])/sqrt(length(Clade2))
sd(sH[which(names(sH) %in% Clade2)])/sqrt(length(Clade2))
sd(sL[which(names(sL) %in% Eggs)])/sqrt(length(Eggs))
sd(sW[which(names(sW) %in% Eggs)])/sqrt(length(Eggs))
sd(sH[which(names(sH) %in% Eggs)])/sqrt(length(Eggs))

cH <- (H/gfd$Csize)
cW <- (W/gfd$Csize)
cL <- (L/gfd$Csize)

##########################################
########## 2D Tests with Parity ##########
##########################################

t.test(sL ~ parity)
t.test(sW ~ parity)
t.test(Ha ~ arb)

###############################################
############### PC 1 to 4 and Body Size & LWH ###############
###############################################
library(ape)
library(geiger)
library(nlme)

tree <- read.tree("tree6")

# set up data frame for generalized analysis
Nsp <- length(tree$tip.label)
x <- size
#x <- sL
#x <- sH
#x <- sW

y <- pca.lands$x[,1]
#y <- pca.lands$x[,2]
#y <- pca.lands$x[,3]
#y <- pca.lands$x[,4]
x2 <- parity
xy <- data.frame(x, y, x2, row.names=attributes(size)$names)
xy <- xy[tree$tip.label, ]
  # puts the data in the same order as the phylogeny

plot(xy$y ~ xy$x)
cor.test(xy$y,xy$x)
  bm.1<-corBrownian(phy=tree)
  bm.gls<-gls(y~x,correlation=bm.1, data=xy)
  summary(bm.gls)
  V<-corMatrix(Initialize(bm.1,xy)); 
  a.Y <- matrix(1,1,Nsp) %*% solve(V) %*% xy$y/sum(solve(V)); 
  a.X <- matrix(1,1,Nsp) %*% solve(V) %*% xy$x/sum(solve(V)); 
  FICr <- (xy$y-a.Y) %*% solve(V) %*% (xy$x-a.X)/sqrt(((xy$y-a.Y) %*% 
    solve(V) %*% (xy$y-a.Y))*((xy$x-a.X) %*% solve(V) %*% (xy$x-a.X))); 
  FICr
  ou.2 <- corMartins(1, phy = tree, fixed = FALSE)
  ou.gls <- gls(y~x,correlation = ou.2, data = xy)
  summary(ou.gls)
  PGLSrA <- unname(unlist(ou.gls$modelStruct))
  PGLSrA
  V<-corMatrix(Initialize(ou.2,xy)); 
  a.Y <- matrix(1,1,Nsp) %*% solve(V) %*% xy$y/sum(solve(V));  
  a.X <- matrix(1,1,Nsp) %*% solve(V) %*% xy$x/sum(solve(V)); 
  PGLSr <-(xy$y-a.Y) %*% solve(V) %*% (xy$x-a.X)/sqrt(((xy$y-a.Y) %*% 
    solve(V) %*% (xy$y-a.Y))*((xy$x-a.X) %*% solve(V) %*% (xy$x-a.X))); 
  PGLSr


####################
##PGLS on procrustus shape and not PC axes
###################

gfd$coords #procrustus coordinates

#size
size.pgls <- procD.pgls(coords ~ Csize, phy=tree, data=gfd)
anova(size.pgls)
summary(size.pgls)

#Length
L.pgls <- procD.pgls(coords ~ gfd$PLhypo, phy=tree, data=gfd)
anova(L.pgls)
summary(L.pgls)

#Width
W.pgls <- procD.pgls(coords ~ gfd$PWilium, phy=tree, data=gfd)
anova(W.pgls)
summary(W.pgls)

#Heigth
H.pgls <- procD.pgls(coords ~ gfd$PHdvf, phy=tree, data=gfd)
anova(H.pgls)
summary(H.pgls)


####################
##PGLS on procrustus shape and parity
###################

parity.pgls <- procD.pgls(coords ~ parity, phy=tree, data=gfd)
anova(parity.pgls)
summary(parity.pgls)

####################
##PGLS on procrustus shape and habitat
###################

hab.pgls <- procD.pgls(coords ~ arb, phy=tree, data=gfd)
anova(hab.pgls)
summary(hab.pgls)


####################
##PGLS on procrustus shape and parity and hbiatat
###################

habpar.pgls <- procD.pgls(coords ~ arb * parity, phy=tree, data=gfd)
anova(habpar.pgls)
summary(habpar.pgls)



phylosig(tree, size, method="lambda", test = TRUE)
phylosig(tree, size, method="K", test = TRUE)
phylosig(tree, sL, method="lambda", test = TRUE)
phylosig(tree, sL, method="K", test = TRUE)
phylosig(tree, sW, method="lambda", test = TRUE)
phylosig(tree, sW, method="K", test = TRUE)
phylosig(tree, sH, method="lambda", test = TRUE)
phylosig(tree, sH, method="K", test = TRUE)


##############################################################################################
# ancestral reconstructions - Figure 3
##############################################################################################

#####For Continuous Variable######
#height
fitML.cont <-ace(sH, tree, type="continuous", method="ML")
anc.rec <- contMap(tree, sH, fsize=.75, type="phylogram")

#length
fitML.cont <-ace(sL, tree, type="continuous", method="ML")
anc.rec <- contMap(tree, sL, fsize=.75, type="phylogram")

#width
fitML.cont <-ace(sW, tree, type="continuous", method="ML")
anc.rec <- contMap(tree, sL, fsize=.75, type="phylogram")

##############################################################################################
# SLOUCH analysis for head dimension, PCs, and body size or color category
##############################################################################################

setwd("C:/Users/emart120/ASU Dropbox/Martins Lab/Research/CT Project/CT People/Logan Kenny (CT Scans)/analysis/")

require(ape)
require(ouch)
require(car)
require(xtable)
require(ggplot2)
require(gridExtra)
require(phytools)
require(slouch)
require(phangorn)

tree2 <- read.tree("tree6")
tree2$edge.length<-
   tree2$edge.length/max(nodeHeights(tree2)[,2])*1  ### rescale to total length of 1.0
plot(tree2)
add.scale.bar(tree2)

sceldat <- read.csv("2DScelData.csv")
rownames(sceldat) <- sceldat[,1]

######################
###### Arboreal ######
######################
#setting node labels as regimes where: 
#0=not arboreal
#1=arboreal
#tree2$node.label <- as.factor(c(0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

#tree2 <- paintBranches(tree2, edge=c(7,9,22,23,29,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,54,55,65,66,84,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111), state=1, anc.state=0)

#plot(tree2)

#x <- sceldat$arb
#y <- sceldat$SVL
##################################

######################
####### Parity #######
######################
#setting node labels as regimes where: 
#0=oviparous
#1=viviparous
tree2$node.label <- as.factor(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0))

tree2 <- paintBranches(tree2, edge=c(17,26,27,28,29,30,31,32,33,34,39,40,41,42,43,44,45,46,85,86,87,88,89,90,91,92,98,99,100,101,102,103,104), state=1, anc.state=0)

plot(tree2)

x <- sceldat$par
y <- sceldat$SVL
#######################################

#Pelvis, size corrected dimensions; appears to be better explained by parity
y.L <- (sceldat$sL)
y.W <- (sceldat$sW)
y.H <- (sceldat$sH)

# puts the data in the same order as the phylogeny
xy <- data.frame(x, y, row.names=row.names(sceldat))
xy <- xy[tree2$tip.label, ]

xy.L <- data.frame(x, y.L, row.names=row.names(sceldat))
xy.L <- xy.L[tree2$tip.label, ]

xy.W <- data.frame(x, y.W, row.names=row.names(sceldat))
xy.W <- xy.W[tree2$tip.label, ]

xy.H <- data.frame(x, y.H, row.names=row.names(sceldat))
xy.H <- xy.H[tree2$tip.label, ]

tipdata <- xy$y
names(tipdata) <- tree2$tip

tipdata.L <- xy.L$y.L
names(tipdata.L) <- tree2$tip

tipdata.W <- xy.W$y.W
names(tipdata.W) <- tree2$tip

tipdata.H <- xy.H$y.H
names(tipdata.H) <- tree2$tip

# Make factors for SLOUCH analysis
xy$x <- as.factor(xy$x)
slfact <- xy$x
names(slfact) <- tree2$tip

####Creating column for Optima####
sceldat$OU1<-"OU1"
OU1 <- as.factor(sceldat$OU1)
names(OU1) <- tree2$tip


#########
## PELVIS LENGTH
##########
## Brownian motion 
#For Parity: Diffusion variance (warning: outside grid), for this, change the first seq # so it includes the decimal variance (ie. the variance is 0.000255, so the seq(0.0001,...) will encompass it all) 
bm.model.L <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.L, hillclimb=TRUE)
summary(bm.model.L)
logLik(bm.model.L)
hillclimbplot(bm.model.L)

## OU1 
slmodel.OU1.L <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 100),
                      vy_values = seq(0.001, 50, length.out = 100),
                      response=tipdata.L, fixed.fact=OU1, hillclimb=TRUE)	
regimeplot(slmodel.OU1.L)			 
summary(slmodel.OU1.L)
logLik(slmodel.OU1.L)
hillclimbplot(slmodel.OU1.L)

## OU multiple
slmodel.L <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.L, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.L)			 
summary(slmodel.L)
logLik(slmodel.L)
hillclimbplot(slmodel.L)

########
## PELVIS WIDTH
########
## Brownian motion 
bm.model.W <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.W, hillclimb=TRUE) #, lower=0.0005)
summary(bm.model.W)
logLik(bm.model.W)
hillclimbplot(bm.model.W)

## OU1 
slmodel.OU1.W <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.W, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.W)			 
summary(slmodel.OU1.W)
logLik(slmodel.OU1.W)
hillclimbplot(slmodel.OU1.W)					  

## OU multiple
slmodel.W <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.W, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.W)			 
summary(slmodel.W)
logLik(slmodel.W)
hillclimbplot(slmodel.W)

########
## PELVIS HEIGHT
########
## Brownian motion 
bm.model.H <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.H, hillclimb=TRUE)
summary(bm.model.H)
plot(bm.model.H)
logLik(bm.model.H)
hillclimbplot(bm.model.H)

## OU1
#All lizards in clade are evolving towards the same optima#
slmodel.OU1.H <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.H, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.H)			 
plot(slmodel.OU1.H)
summary(slmodel.OU1.H)
logLik(slmodel.OU1.H)
hillclimbplot(slmodel.OU1.H)				

## OU multiple
#Multiple optima (ie. one for egg laying and one for live bearing)
slmodel.H <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.H, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.H)			 
plot(slmodel.H)
summary(slmodel.H)
logLik(slmodel.H)
hillclimbplot(slmodel.H)


##############################################################################################
### PANCOVA ###
##############################################################################################

setwd("C:/Users/emart120/ASU Dropbox/Martins Lab/Research/CT Project/CT People/Logan Kenny (CT Scans)/analysis/")

## Load libraries ##
library(ctv)
library(ape)
library(geiger)
library(nlme)
library(phangorn)
library(phytools)

### INPUTS #############################################################################  
########################################################################################  

# Load data and sort it alphabetically (set working directory first)
Data<-read.csv(file='2DScelData.csv',header=TRUE)
Data<-Data[order(Data$species),]

# OTUs names
Names<-Data[,1]

# Continuous response variable (run 1 at a time; update column number for pelvis)
Trait<-Data[,8]                # pelvisL
#Trait<-Data[,9]               # pelvisW
#Trait<-Data[,10]               # pelvisH


# Categorical explanatory variable (the thing we are interested in, locomotion vs. parity, update column number for parity and run, then update for locomotion then run)
Factor<-Data[,4]                # Parity

# Continuous explanatory variable (update SVL column)
Covariate<-Data[,6]            # Raw SVL


# Load and visualize phylogeny
tree<-read.tree("tree6")
plot(tree,cex=0.4)
nodelabels(cex=0.4)

# Confirm name match; super important, will tell if tree or data is messed up (must match)
sum(Names == sort(tree$tip.label))

## Specify evolutionary regime ##
# The next step assigns the clades that should be mapped in the apomorphic character (ovi is the plesiomorphic trait and vivi is apomorphic trait because it was nearly involved; and then look at locomotion)
# state of the factor in a matrix, where each row corresponds to a group of species 
# that have inherited the derived state after a single change in the phylogeny.
# populate FastClades for aphomorphic trait; so all of the individuals who are viviparous will be entered; order doesn't matter, but do it by clade, so more closely related species should be next to eachother (3, enter, 3, enter, and put "xxx" if there is not a group of 3 in the clade; group by colors from phylogeny article)

FastClades=t(array(c('Sceloporus_grammicus','Sceloporus_megalepidurus','Sceloporus_jarrovii','Sceloporus_cyanogenys','Sceloporus_poinsettii','Sceloporus_prezygus','Sceloporus_minor','Sceloporus_dugesii','Sceloporus_mucronatus',
                     'Sceloporus_internasalis','Sceloporus_acanthinus','Sceloporus_hartwegi','Sceloporus_taeniocnemis','Sceloporus_formosus','Sceloporus_adleri','Sceloporus_formosus_scitulus','xxx','xxx',
                     'Sceloporus_bicanthalis','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx'),dim=c(9,3)))


##########################
### ANALYSES #############  
##########################

#changes in pelvis dimension are independent from phylogenetic (evolutionary) history
# Non-phylogenetic ANCOVA
source("TIPS.R")

#dimension of pelvis are result of a stocastic process (random evolution)
# PANCOVA under Brownian motion
source("BM.R")

#some signal (along phylogeny) constraining evolution, species more closely related will look more similar; intermediate between the two above
# PANCOVA accounting for phylogenetic signal
source("PS.R")



##########################################################################
##### Clustering
##########################################################################

# MCLUSTDA ##########################################################
library(dplyr)
library(mclust)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(psych)

set.seed(12213)
setwd("C:/Users/emart120/ASU Dropbox/Martins Lab/Research/CT Project/CT People/Logan Kenny (CT Scans)/analysis/")
Data <- read.csv("2DScelData25.csv")

XZ <- Data %>% select('L', 'W', 'H', 'SVL', 'arb', 'par', 'sp')

XYZ <- XZ %>% mutate_at(c('L', 'W', 'H', 'SVL'), ~(scale(.) %>% as.vector))
XYZ

#use 60% of dataset as training set and 40% as test set
sample <- sample(c(TRUE, FALSE), nrow(XYZ), replace=TRUE, prob=c(0.6,0.4))
X.train <- XYZ[sample, 1:4]
X.test <- XYZ[!sample, 1:4]
Class1.train <- XYZ[sample, 5] ### arboreality
Class1.test <- XYZ[!sample, 5]
clPairs(X.train, Class1.train)

Class2.train <- XYZ[sample,6] ### parity
Class2.test <- XYZ[!sample,6]

clPairs(X.train, Class1.train)
clPairs(X.train, Class2.train)

BIC <- mclustBIC(X.train)
plot(BIC)
summary(BIC)
mod1 <- Mclust(X.train, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")
table(Class1.train, mod1$classification)
plot(mod1, what = "uncertainty")

mod2 <- MclustDA(X.train, Class1.train, modelType = "EDDA")
summary(mod2)

mod3 <- MclustDA(X.train, Class1.train)
summary(mod3)


### TRAINING ###

MclustDA_Arb <- MclustDA(X.train, Class1.train) 
MclustDA_Par <- MclustDA(X.train, Class2.train) 

summary(MclustDA_Arb, parameters=TRUE)
plot(MclustDA_Arb, what = "classification", symbols=c(21, 22))

summary(MclustDA_Par, parameters=TRUE)
plot(MclustDA_Par, what = "classification", symbols=c(21, 22))

### TESTING ###

pca.pelvis <- prcomp(XYZ[1:4], scale=TRUE)
var <- get_pca_var(pca.pelvis)
var$contrib
fviz_pca_var(pca.pelvis, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

summary(MclustDA_Arb, newdata = X.test, newclass = Class1.test)
plot(MclustDA_Arb, what = "classification", newdata = X.test)
fviz_pca_biplot(pca.pelvis, label="var", repel=TRUE, habillage = XYZ$arb)

summary(MclustDA_Par, newdata = X.test, newclass = Class2.test)
plot(MclustDA_Par, what = "classification", newdata = X.test)
fviz_pca_biplot(pca.pelvis, label="var", repel=TRUE, habillage = XYZ$par)


################################# varimax ###########################
ncomp=2
pca.pelvis.rotated <- psych::principal(XZ[,-5:-7], rotate="varimax", nfactors=ncomp, scores=TRUE)
#pca.pelvis.rotated$scores  # Scores returned by principal()
#pca.pelvis.rotated$values  # eigenvalues
pca.pelvis.rotated$loadings

biplot(pca.pelvis.rotated, XZ[,5], col = c("blue", "red"), main = "PCA with pelvis dimensions")














