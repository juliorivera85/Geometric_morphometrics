##### this is the code for digitizing and analyzing skull scans
install.packages(c("PCDimension", "geiger", "mvMORPH", "car", "geomorph", "ips"))

setwd("C:/Users/emart120/Dropbox (ASU)/MartinsLab_old/Scelop17_20ASU/CT Scan project (mostly pre-2022, JR & EM only)/analysis")

setwd("~/Library/CloudStorage/Dropbox/CT Scan project/analysis")

#load packages needed
library(PCDimension)
library(geiger)
library(mvMORPH)
require(car)
library(geomorph)
library(ips)
library(scatterplot3d)
library(phytools) 

###########################################
###########################################
###### START OF ANALYSIS
###########################################
###########################################

# the file name should match name of phylogeny and on spreadsheet
# this is for males
setwd("male/")
#this is for females
#setwd("female/")

filelist <- list.files(pattern=".nts")
dat <- readmulti.nts(filelist)
setwd("../")

#making groups from spreadsheet
group <- read.csv("scelinfo.csv")  #### species	belly	parity	size(S,M,L)
filelist2 <- gsub(".nts", "", filelist)
vars <- group$species %in% filelist2
group2 <- group$species[vars]
group2 <- as.data.frame(group2)
names(group2) <- "species" 
group2 <- merge(group, group2)

belly <- group2$belly <- as.factor(group2$belly)
parity <- group2$parity <- as.factor(group2$parity)
species <- group2$species <- as.character(group2$species)
size <- group2$size <- as.factor(group2$size)

dimnames(dat)[[3]] <- gsub(pattern = ".nts_1", replacement = "", x = dimnames(dat)[[3]])
## associating the categorical data with species it belongs to
names(belly) <- dimnames(dat)[[3]]
belly
names(parity) <- dimnames(dat)[[3]]
parity
names(size) <- dimnames(dat)[[3]]
size

## Start of Generalized Procrustes Analysis
## superimposition of the raw coordinate data
Y.gpa <- gpagen(dat, Proj=TRUE, ProcD=TRUE, curves=NULL, surface=NULL)
summary(Y.gpa)
plot(Y.gpa)
y <- two.d.array(Y.gpa$coords)

# Many objects are generated
attributes(Y.gpa) 
# 3D array of Procrustes coordinates
Y.gpa$coords 
# Vector of centroid sizes
Y.gpa$Csize 
# The number of GPA iterations until convergence was found
Y.gpa$iter
# Variance-covariance matrix among landmark coordinates
Y.gpa$points.VCV
# Variances of landmark points
Y.gpa$points.var
# The consensus (mean) configuration; akin to using mshape
Y.gpa$consensus
# Data frame with an n x (pk) matrix of Procrustes residuals and centroid size
Y.gpa$data
# Final convergence criterion value
Y.gpa$Q

#females
#YF.gpa <- gpagen(dat, Proj=TRUE, ProcD=TRUE, curves=NULL, surface=NULL)
#summary(YF.gpa)
#plot(YF.gpa)
#yf <- two.d.array(YF.gpa$coords)

##creating a geomorph dataframe for downstream analyses
gfd <- geomorph.data.frame(Y.gpa)
#gfdf <- geomorph.data.frame(YF.gpa)

gfd$belly <- belly
gfd$parity <- parity 
gfd$svl <- size

#gfdf$belly <- belly
#gfdf$parity <- parity 
#gfdf$svl <- size

## Checking outliers for possible digitization problems
outliers <- plotOutliers(Y.gpa$coords) 
M <- mshape(Y.gpa$coords)   
## Specimens falling above the upper quartile are plotted in red and their 
## address returned, for inspection by plotRefToTarget
# Example (for the first outlier)
plotRefToTarget(M,Y.gpa$coords[,,outliers[1]], method="TPS", label = T)

# Example (for the second outlier)
plotRefToTarget(M,Y.gpa$coords[,,outliers[2]], method="vector", label = T)
		
# Example (for the second outlier)
plotRefToTarget(M,Y.gpa$coords[,,outliers[3]], method="vector", label = T)
		
##This code allows to you see the deformation of all the species
plotRefToTarget(M,Y.gpa$coords[,,4], method="vector", label = F, verbose=T)


###########################
#### new stuff allometry test
###########################

allom <- procD.lm(Y.gpa$coords ~ log(Y.gpa$Csize), iter=999, R2=TRUE)
summary(allom)

#test and plot allometry
allom_plot <- plot(allom, type="regression", predictor=log(Y.gpa$Csize), reg.type="RegScore")

pdf(file="skull_regression.pdf")
nice_plot <- plot(allom_plot$RegScore ~ log(Y.gpa$Csize), xlab="Log(centroid size)", ylab="Shape regression score", pch=21, cex=1.5, bg="gray")
abline(lm(allom_plot$RegScore ~ log(Y.gpa$Csize)), lwd=2)
dev.off()

#predicted shape
preds <- shape.predictor(allom$GM$fitted, x=allom_plot$PredLine,
						predmin = min(allom_plot$PredLine),
						predmax = max(allom_plot$PredLine))
plotRefToTarget(preds$predmin, preds$predmax, method="vector")

#PCA
PCA.plot <- gm.prcomp(Y.gpa$coords)

pdf(file="size_PCA_12.pdf")
plot(PCA.plot$x[,1], PCA.plot$x[,2], xlab="PC1", ylab="PC2", pch=21, cex=Y.gpa$Csize/50, bg="gray")
dev.off()
pdf(file="size_PCA_34.pdf")
plot(PCA.plot$x[,3], PCA.plot$x[,4], xlab="PC3", ylab="PC4", pch=21, cex=Y.gpa$Csize/50, bg="gray")
dev.off()

############################
####### getting the phylogeny (skip it if you read in phylogeny)
############################

## this code reads in the tree and matches the species on the tree with the species
## we have skull data for and drops the ones we don't need.
#tree <- read.nexus("scelop_timecalib_2016.nex")

#tree$tip.label[[82]]<-"Sceloporus_prezygus"  # changing serrifer to prezygus
#tree.tps <- tree$tip
#find species that we don't have data for and remove the tips
#To.Drop <-  tree.tps[!tree.tps %in% group2$species]
#tree2 <- drop.tip(tree, To.Drop)

#find species we have shape data for but not in the tree
#group2$species[!group2$species %in% tree.tps]

# add cozumelae, sister to variabilis
#node <- which(tree2$tip.label=="Sceloporus_variabilis"); node
#plot(tree2)
#edgelabels()  ### find edge number below variabilis = 6
#tmp<-tree2$edge.length[[6]]; tmp<-tmp/2; tmp  ## place it half way down
#tree3 <- bind.tip(tree2, tip.label="Sceloporus_cozumelae", where=node, position=tmp, edge.length=tmp)
#plot(tree3)

# add hartwegi as sister to taeneocnemis
#node <- which(tree3$tip.label=="Sceloporus_taeniocnemis_mvz4213"); node
#plot(tree3); edgelabels()  ### find edge number below taeniocnemis = 91
#tmp<-tree3$edge.length[[91]]; tmp<-tmp/2; tmp  ## place it half way down
#tree4 <- bind.tip(tree3, tip.label="Sceloporus_hartwegi", where=node, position=tmp, edge.length=tmp)
#plot(tree4)

# add utiformis as sister to grandeveus
#node <- which(tree4$tip.label=="Sceloporus_grandaevus"); node
#plot(tree4); edgelabels()  ### find edge number below grandaevus = 15
#tmp<-tree4$edge.length[[15]]; tmp<-tmp/2; tmp  ## place it half way down
#tree5 <- bind.tip(tree4, tip.label="Sceloporus_utiformis", where=node, position=tmp, edge.length=tmp)
#plot(tree5)

#write.tree(tree5, file="tree5")

tree5 <- read.tree("tree5")

## checking names to make sure data matches tree
names(Y.gpa$Csize) <- gsub(pattern = ".nts_1", replacement = "", x = names(Y.gpa$Csize))
dimnames(Y.gpa$coords)[[3]]<-gsub(pattern = ".nts_1", replacement = "", x = dimnames(Y.gpa$coords)[[3]])
l1<-sort(names(Y.gpa$Csize))
l2<-sort(tree5$tip.label)
identical(l1,l2)
setdiff(l1,l2)

### for females
names(YF.gpa$Csize) <- gsub(pattern = ".nts_1", replacement = "", x = names(YF.gpa$Csize))
dimnames(YF.gpa$coords)[[3]]<-gsub(pattern = ".nts_1", replacement = "", x = dimnames(YF.gpa$coords)[[3]])
tree.tps <- tree5$tip
To.Drop <-  tree.tps[!tree.tps %in% l1]
tree6 <- drop.tip(tree5, To.Drop)
l1<-sort(names(YF.gpa$Csize))
l2<-sort(tree6$tip.label)
identical(l1,l2)
setdiff(l1,l2)

############################
####### non-phylo analyses 
############################
####### PCA  
### see https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.PCA.html for tutorial
par(mfrow=c(1,1)) # 1 graph at a time
# pca.lands <- plotTangentSpace(A=Y.gpa$coords, label=TRUE)## deprecated
pca.lands <- gm.prcomp(A=Y.gpa$coords) ## new version
#pca.lands <- gm.prcomp(A=Y.gpa$coords, phy=tree2) ## new version
summary(pca.lands) ## % variance explained is here

#par(mar=c(2, 2, 2, 2))
colors <- c("blue", "white") # mark belly color
colors <- colors[as.numeric(gfd$belly)]
plot(pca.lands, pch=22, cex = 1.5, bg = colors) 
#plot(pca.lands, pch=22, cex = 1.5, bg = gfd$parity) 
#plot(pca.lands, pch=22, cex = 1.5, bg = gfd$svl) 
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (34%)", pos = 4, font = 2)
text(0, 0.95*par()$usr[4], labels = "PC2 (9%)", pos = 4, font = 2)
legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))

 ### Visualize shape variation using picknplot.shape Because picknplot requires 
 ### user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(plot(pca.lands))

sum(pca.lands$sdev^2)
barplot(pca.lands$sdev^2/sum(pca.lands$sdev^2))

# Using broken stick model
brokenStick(1:ncol(pca.lands$x),ncol(pca.lands$x))
sum(brokenStick(1:ncol(pca.lands$x),ncol(pca.lands$x)))
# How many significant PCs?
bsDimension(summary(pca.lands)$PC.summary[2,])
#the first 4 PC axis are significant

plot(pca.lands$x[,1:2],pch=16,xlab="PC1",ylab="PC2")		
text(pca.lands$x[,1:2],rownames(pca.lands$x[,1:2]),pos=4,cex=.5)

scatterplot3d(pca.lands$x[,1:3], color="blue", pch=16)

PCA.scores.scel <- pca.lands$x[,1:4]  ### first 4 PCs only
write.csv(PCA.scores.scel, file="scel_pca_score.csv")

par(mfrow=c(1,1)) # 1 graph at a time
consensus <- apply(Y.gpa$coords, c(1,2), mean)
plot(consensus,asp=1, type="n")

for(i in 1:length(Y.gpa$Csize)) 
  {points(Y.gpa$coords[,,i])}
points(consensus, col="red", cex=2, pch=20)


##### procrustes regression for shape, basically an ANOVA for Goodall's F-test
##### if you want cont variables, same code for X should be cont.
belly.aov <- procD.lm(gfd$coords ~ belly, iter=999) # not significant
summary(belly.aov)


size.aov <- procD.lm(gfd$coords ~ size, iter=999) # significant
summary(size.aov)
pairwise(size.aov)
summary(pairwise(size.aov, groups=size))

size.aov.phy <- procD.pgls(gfd$coords ~ size, phy=tree5, iter=999) # significant
summary(size.aov.phy)
summary(pairwise(size.aov.phy, groups=size))

belsize.aov <- procD.lm(gfd$coords ~ belly + size, iter=999) # significant
summary(belsize.aov)

# Visualizing mean shapes from group as factor (based on MANOVA) i.e. blue vs white
# PCA using prcomp
PCAgroup.m <- prcomp(belly.aov$fitted)
# plots position of group means (based on color) 
plot(PCAgroup.m$x, asp=1, pch=19) 
# Mean shapes (according to groups)
GroupMeans.m <- unique(round(PCAgroup.m$x,4))
# Note that previous line picks up 3 PCs (check if useful enough)
GroupMeans.m 
# Prediction (using 4 PCs, and assuming two groups): 
# 14 is the # of landmarks (customize)
# 69 is the number of dimensions (3D data)
predGroup.m <- shape.predictor(arrayspecs(belly.aov$fitted, 69,3), 
                              x= PCAgroup.m$x[,1:4], Intercept = FALSE, 
                              Group1 = GroupMeans.m[1,1:4], 
                              Group2 = GroupMeans.m[2,1:4])
## Shape differences using consensus as reference and "TPS" as plotting method 
## (for 3D data a method different to "TPS" might be better, such as "vector" or "points"; see below):
## a thin-plate spline deformation grid is generated (better for 2D data). 
# Shape of group 1 (note that M corresponds to the general consensus) 
plotRefToTarget(M, predGroup.m$Group1) 
# Shape of group 1 magnified by 3X (to visualize differences better)
plotRefToTarget(M, predGroup.m$Group1, mag=5) 
# Shape of group 2
plotRefToTarget(M, predGroup.m$Group2)
# Shape of group 2 magnified by 3X (to visualize differences better)
plotRefToTarget(M, predGroup.m$Group2, mag=5)


## Differences between groups: Another visualization procedure
# For groups shown to be significantly different with procD.lm, first 
# average the data by groups, calculate group means and plot those means:
## calculate mean shape for each group
meansGroup.m <- aggregate(two.d.array(gfd$coords) ~ belly, FUN=mean)
# First row is first group
meansGroup.m[1,1]
# make mean vectors as matrix for group 1
grp1.mn.m <- matrix(as.numeric(meansGroup.m[1,-1]), ncol=3, byrow=T)
# Second row is another group
meansGroup.m[2,1]
# make mean vectors as matrix for second group
grp2.mn.m <- matrix(as.numeric(meansGroup.m[2,-1]), ncol=3, byrow=T)

## Shape of groups: Should be similar to the MANOVA predictions
# For group 1
plotRefToTarget(M, grp1.mn.m, method="vector") 
# For group 1 magnifying differences by 5 (to visualize differences better)
plotRefToTarget(M, grp1.mn.m, method="vector", mag=5) 
# For group 2
plotRefToTarget(M, grp2.mn.m, method="vector") 
# For group 2 magnifying differences by 5 (to visualize differences better)
plotRefToTarget(M, grp2.mn.m, method="vector", mag=5) 

parity.aov <- procD.lm(y ~ parity, data=gfd, iter=999) #not significant
summary(parity.aov)
#### diagnostic plots, including plotOutliers (check historial)
plot(parity.aov, type = "diagnostics", outliers = TRUE)
#### PC plot rotated to major axis of fitted values
plot(parity.aov, type = "PC", pch = 19, col = "blue")


###########################
###########################
# Modularity test
###########################
###########################

#Y.gpa$coords contains procrustes coordinates
#we will see if back of the head, front, and jaw are modular

#first we will creat a categorical list of the points and if they belong to
#A=front skull, B=back skull, C=jaw
#part.gp <- c("A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","A","A","A","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C")
#part.gp <- as.factor(part.gp)

#part2.gp<- c(rep("A",13), rep("B",21), rep("A",3), rep("C",32))
## as Julio initially defined in 3 parts, front/back/jaw
## Sanger et al (2011) found elongate anoles had tripartite skulls, but defined
##   the 3 parts as snout, orbital region, and braincase
part2.gp<-c(rep("0",69))
part2.gp[1:13]  <- "A" ; part2.gp[35:37] <- "A" 
part2.gp[14:34] <- "B" ; part2.gp[38:69] <- "C" 
part2.gp <- as.factor(part2.gp); part2.gp

# "Anolis model". Anterior vs posterior with line in mid-orbit (Sanger et al. 2011)
## landmarks: 1-10 front; 11-34 back, 35-37 front, 38-65 back, 66-67 front, 68-69 back
anolis.gp<-c(rep("0",69))
anolis.gp[1:10]  <- "A" ; anolis.gp[35:37] <- "A" 
anolis.gp[11:34] <- "B" ; anolis.gp[38:65] <- "B" 
anolis.gp[66:67] <- "A" ; anolis.gp[68:69] <- "B" 
anolis.gp <- as.factor(anolis.gp); anolis.gp

# "Mammalian model". front includes outer sides of back skull (Sanger et al. 2011)
## reflects developmental origins
## landmarks: 1-10 front; 11-30 back, 31-43 front, 44-45 back, 46-51 front, 52-59 back, 60-67 front, 68-69 back
## 31-43, 46-51, 60-67 count as "front"
mammal.gp<-c(rep("0",69))
mammal.gp[1:10]  <- "A" ; mammal.gp[31:43] <- "A" 
mammal.gp[11:30] <- "B" ; mammal.gp[44:45] <- "B" 
mammal.gp[46:51] <- "A" ; mammal.gp[52:59] <- "B" 
mammal.gp[60:67] <- "A" ; mammal.gp[68:69] <- "B" 
mammal.gp <- as.factor(mammal.gp); mammal.gp


part2.mod <- modularity.test(Y.gpa$coords, part2.gp, iter=999)
anolis.mod <- modularity.test(Y.gpa$coords, anolis.gp, iter=999)
mammal.mod <- modularity.test(Y.gpa$coords, mammal.gp, iter=999)

model.Z <- compare.CR(part2.mod, anolis.mod, mammal.mod, CR.null = TRUE)
summary(model.Z)
## tripartite model fits best

plot(part2.mod); plot(anolis.mod); plot(mammal.mod)

### Females
#part2F.mod <- modularity.test(YF.gpa$coords, part2.gp, iter=999)
#anolisF.mod <- modularity.test(YF.gpa$coords, anolis.gp, iter=999)
#mammalF.mod <- modularity.test(YF.gpa$coords, mammal.gp, iter=999)
#modelF.Z <- compare.CR(part2F.mod, anolisF.mod, mammalF.mod, CR.null = TRUE)
#summary(modelF.Z)
#plot(part2F.mod); plot(anolisF.mod); plot(mammalF.mod)

## Phylo Integration
part2.int <- integration.test(Y.gpa$coords, part2.gp, A2=NULL)
anolis.int <- integration.test(Y.gpa$coords, anolis.gp, A2=NULL)
mammal.int <- integration.test(Y.gpa$coords, mammal.gp, A2=NULL)
summary(part2.int); summary(anolis.int); summary(mammal.int)
#plot(part2.int)

tree5 <- read.tree("tree5")
plot(tree5)
part2.phyint <- phylo.integration(Y.gpa$coords, part2.gp, phy=tree5, A2=NULL)
anolis.phyint <- phylo.integration(Y.gpa$coords, anolis.gp, phy=tree5, A2=NULL)
mammal.phyint <- phylo.integration(Y.gpa$coords, mammal.gp, phy=tree5, A2=NULL)
summary(part2.phyint); summary(anolis.phyint); summary(mammal.phyint)
plot.pls(part2.phyint) #exact same resilts as normal integration

#tree5 <- read.tree("tree5"); plot(tree5)
### need to remove taxa not in data set
#part2F.phyint <- phylo.integration(YF.gpa$coords, part.gp, phy=tree5, A2=NULL)
#anolisF.phyint <- phylo.integration(YF.gpa$coords, part3.gp, phy=tree5, A2=NULL)
#mammalF.phyint <- phylo.integration(YF.gpa$coords, part4.gp, phy=tree5, A2=NULL)
#summary(part2F.phyint); summary(anolisF.phyint); summary(mammalF.phyint)

## Get linear measure of head length, head width, and head height

setwd("male/")
filelist <- list.files(pattern=".nts")
dat <- readmulti.nts(filelist)
setwd("../")

#now, create a matrix telling R where to start and stop measuring from in terms of the landmarks.
lmks <- matrix(c(15,37,50,51,12,43), ncol=2, byrow=TRUE, dimnames = list(c("head L", "head W", "head H"),c("start", "end")))

#this command will automatically measure all the distance between landmarks and produce a data frame 
## with the lengths of each sengment
result <- interlmkdist(dat, lmks)
#write.csv(result[1], file="result1.csv")
#write.csv(result[2], file="result2.csv")
#write.csv(result[3], file="result3.csv")
lin.dat <- read.csv("linear mes.csv")

lin.dat$head.L <- lin.dat$head.L/lin.dat$SVL
lin.dat$head.W <- lin.dat$head.W/lin.dat$SVL
lin.dat$head.H <- lin.dat$head.H/lin.dat$SVL

#ventral patch and linear measures
aov.head.L.col <- aov(head.L ~ color, data=lin.dat)
summary(aov.head.L.col)
TukeyHSD(aov.head.L.col) #white have longer head

aov.head.W.col <- aov(head.W ~ color, data=lin.dat)
summary(aov.head.W.col)
TukeyHSD(aov.head.W.col) # no difference

aov.head.H.col <- aov(head.H ~ color, data=lin.dat)
summary(aov.head.H.col)
TukeyHSD(aov.head.H.col) # no difference


#size and linear measures
aov.head.L.size <- aov(head.L ~ size, data=lin.dat)
summary(aov.head.L.size)
TukeyHSD(aov.head.L.size) #small species have longer heads than L

aov.head.W.size <- aov(head.W ~ size, data=lin.dat)
summary(aov.head.W.size)
TukeyHSD(aov.head.W.size) # no difference

aov.head.H.size <- aov(head.H ~ size, data=lin.dat)
summary(aov.head.H.size)
TukeyHSD(aov.head.H.size) #no difference


###########################
###########################
# 3 types of shape PCAs
###########################
###########################

tree5 <- read.tree("tree5")
plot(tree5)

#par(mar=c(2, 2, 2, 2))
#colors <- c("blue", "white") 
#colors <- colors[as.numeric(gfd$belly)]  # mark belly color
colors <- c("red", "yellow", "green") 
colors <- colors[as.numeric(gfd$svl)]  # mark body size
#plot(pca.GLS, pch=22, cex = 1.5, bg = gfd$parity) 

par(mfrow=c(1,1)) # 1 graph at a time
pca.lands <- gm.prcomp(A=Y.gpa$coords) ## traditional PCA
summary(pca.lands)
PS.shape <- physignal(A=Y.gpa$coords, phy=tree5, iter=999)### phy sig in shape
summary(PS.shape) ### K-mult test, is identical for PCA, phylomorphospace, and phy-PCA
PS.size <- physignal(A=Y.gpa$Csize, phy=tree5, iter=999) ### phy sig in size
summary(PS.size) ### K-mult test, is identical for PCA, phylomorphospace, and phy-PCA

plot(pca.lands, pch=22, cex = 1.5, bg = colors) 
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (34%)", pos = 4, font = 2) ## TIPS
text(0, 0.95*par()$usr[4], labels = "PC2 (9%)", pos = 4, font = 2) ## TIPS
#legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

plot(pca.lands$x[,1], Y.gpa$Csize, pch=16, xlab="PC1", ylab="Centroid Size", 
          col=colors)
text(pca.lands$x[,1], Y.gpa$Csize, names(pca.lands$x[,1]), pos=3, cex=.5)
### S. formosus replicates = 17, 19
points(c(unname(pca.lands$x[,1][17]), unname(pca.lands$x[,1][19])), 
       c(unname(Y.gpa$Csize[17]), unname(Y.gpa$Csize[19])), 
               pch = 17, cex = 4, col = "red")
### S. smargadinus replicates = 47, 48
points(c(unname(pca.lands$x[,1][47]), unname(pca.lands$x[,1][48])), 
       c(unname(Y.gpa$Csize[47]), unname(Y.gpa$Csize[48])), 
               pch = 16, cex = 4, col = "green")
### S. grammicus replicates = 21, 22
points(c(unname(pca.lands$x[,1][21]), unname(pca.lands$x[,1][22])), 
       c(unname(Y.gpa$Csize[21]), unname(Y.gpa$Csize[22])), 
               pch = 15, cex = 4, col = "orange")



#fit <- procD.lm(coords ~ log(Csize), data=gfd, iter=0, print.progress = FALSE)
#plotAllometry(fit, size = gfd$Csize, logsz = TRUE, method = "PredLine", pch = 19)
msh <- mshape(Y.gpa$coords) ### global consensus 
plotRefToTarget(pca.lands$shapes$shapes.comp1$min, msh) ### minimum deformation grid, PC1
plotRefToTarget(msh, pca.lands$shapes$shapes.comp1$max) ### maximum deformation grid, PC1 
plotRefToTarget(pca.lands$shapes$shapes.comp1$min, pca.lands$shapes$shapes.comp1$max, method = "vector", mag = 2)
par(mfrow=c(1,2)) # 2 graphs in a row
plot(Y.gpa$Csize, pca.lands$x[,1])     ### PC1 on Centroid size
cor.test(Y.gpa$Csize, pca.lands$x[,1]) ### PC1 on Centroid size
plot(Y.gpa$Csize, pca.lands$x[,2])     ### PC2 on Centroid size
cor.test(Y.gpa$Csize, pca.lands$x[,2]) ### PC2 on Centroid size

pca.lands2 <- gm.prcomp(A=Y.gpa$coords, phy=tree5) ## phylomorphospace, tree is used 
### only in special plots (not these), analysis results are identical to traditional PCA
plot(pca.lands2, pch=22, cex = 1.5, bg = colors) 
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (34%)", pos = 4, font = 2) ## TIPS
text(0, 0.95*par()$usr[4], labels = "PC2 (9%)", pos = 4, font = 2) ## TIPS
#legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

#### From Collyer & Adams 2021: An obvious difference in the dispersion of values in
####  PC plots from PA and Phy-PCA might suggest factors other than phylogeny
####  account for trait variation. An obvious difference in the dispersion
####  of values in PC plots from PA and Phy-PCA might suggest factors
####  other than phylogeny account for trends in the multivariate data space.
####

par(mfrow=c(1,3)) # 3 graphs in a row
## in male skulls, PCA = phylomorphospace, which are both different from phy-PCA. Phy matters.
plot(pca.lands, phylo = TRUE, main = "PCA") ### traditional PCA
plot(pca.lands2, phylo = TRUE, main = "PCA.w.phylo") ### phylomorphospace
plot(pca.GLS, phylo = TRUE, main = "phylo PCA")      ### Revell phy-PCA

par(mfrow=c(1,2)) # 2 graphs in a row
### Phylogenetic PCA - PCA based on GLS-centering; Revell 2009
pca.GLS <- gm.prcomp(A=Y.gpa$coords, phy=tree5, GLS=TRUE) 
summary(pca.GLS)
plot(pca.GLS, pch=22, cex = 1.5, bg = colors) 
#plot(pca.GLS, phylo = TRUE, pch=22, cex = 1.5, bg = colors) # with species names and tree
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (42%)", pos = 4, font = 2) ## GLS
text(0, 0.95*par()$usr[4], labels = "PC2 (21%)", pos = 4, font = 2) ## GLS
#legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

#msh <- mshape(Y.gpa$coords)
plotRefToTarget(pca.GLS$shapes$shapes.comp1$min, msh)
plotRefToTarget(msh, pca.GLS$shapes$shapes.comp1$max)
plotRefToTarget(pca.GLS$shapes$shapes.comp1$min, pca.GLS$shapes$shapes.comp1$max, method = "vector", mag = 2)

par(mfrow=c(1,2)) # 2 graphs in a row
plot(Y.gpa$Csize, pca.GLS$x[,1])     ### PC1 on Centroid size
cor.test(Y.gpa$Csize, pca.GLS$x[,1]) ### PC1 on Centroid size
plot(Y.gpa$Csize, pca.GLS$x[,2])     ### PC2 on Centroid size
cor.test(Y.gpa$Csize, pca.GLS$x[,2]) ### PC2 on Centroid size


### Phylogenetic PCA - PCA based on GLS-centering and transformed; Revell 2009
###    PCA is independent of phylogeny
phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = tree5, GLS = TRUE, transform = TRUE) 
summary(phylo.tPCA)
#plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA with transformation", pch=22, cex = 1.5, bg = colors)
plot(phylo.tPCA, main = "phylo PCA with transformation", pch=22, cex = 1.5, bg = colors)
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (42%)", pos = 4, font = 2) ## GLS, transformed
text(0, 0.95*par()$usr[4], labels = "PC2 (21%)", pos = 4, font = 2) ## GLS, transformed
#legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

#msh <- mshape(Y.gpa$coords)
plotRefToTarget(phylo.tPCA$shapes$shapes.comp1$min, msh)
plotRefToTarget(msh, phylo.tPCA$shapes$shapes.comp1$max)
plotRefToTarget(phylo.tPCA$shapes$shapes.comp1$min, phylo.tPCA$shapes$shapes.comp1$max, method = "vector", mag = 2)
plot(phylo.tPCA, time.plot = TRUE, pch = 22, bg = c(rep("red", 5), rep("green", 4)), cex = 2, 
     phylo.par = list(edge.color = "grey60", edge.width = 1.5, tip.txt.cex = 0.75,
                      node.labels = F, anc.states = F))

### PaCA - Alignment of data to physlogenetic signal rather than axis of 
### greatest variation, like in PCA - see Collyer and Adams 2021.

par(mfrow=c(1,3)) # 3 graphs in a row

#par(mfrow=c(1,1)) # 1 graph at a time
 # OLS method (rotation of PCA)
 PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = tree5, align.to.phy = TRUE)
 summary(PaCA.ols)
 plot(PaCA.ols, main = "PaCA using OLS", pch=22, cex = 1.5, bg = colors, phylo=TRUE)
 text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (97%)", pos = 4, font = 2) ## data aligned to phy sig, OLS rotation
 text(0, 0.95*par()$usr[4], labels = "PC2 (21%)", pos = 4, font = 2) ## data aligned to phy sig, OLS rotation
 #legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
 legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

 # GLS method (rotation of Phylogenetic PCA)
 PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = tree5, align.to.phy = TRUE, GLS = TRUE)
 summary(PaCA.gls)
 plot(PaCA.gls, main = "PaCA using GLS", pch=22, cex = 1.5, bg = colors, phylo=TRUE)
 text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (71%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
 text(0, 0.95*par()$usr[4], labels = "PC2 (26%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
 #legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
 legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))
 
#msh <- mshape(Y.gpa$coords)
plotRefToTarget(PaCA.gls$shapes$shapes.comp1$min, msh)
plotRefToTarget(msh, PaCA.gls$shapes$shapes.comp1$max)
plotRefToTarget(PaCA.gls$shapes$shapes.comp1$min, PaCA.gls$shapes$shapes.comp1$max, method = "vector", mag = 2)

par(mfrow=c(1,2)) # 2 graphs in a row
plot(Y.gpa$Csize, PaCA.gls$x[,1])     ### PC1 on Centroid size
cor.test(Y.gpa$Csize, PaCA.gls$x[,1]) ### PC1 on Centroid size
plot(Y.gpa$Csize, PaCA.gls$x[,2])     ### PC2 on Centroid size
cor.test(Y.gpa$Csize, PaCA.gls$x[,2]) ### PC2 on Centroid size

 # GLS method (rotation of Phylogenetic PCA with transformed data)
 PaCA.gls2 <- gm.prcomp(Y.gpa$coords, phy = tree5, 
   align.to.phy = TRUE, GLS = TRUE, transform = TRUE)
 summary(PaCA.gls2)
 plot(PaCA.gls2, main = "PaCA using GLS and transformed projection", pch=22, cex = 1.5, bg = colors)
 text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (71%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
 text(0, 0.95*par()$usr[4], labels = "PC2 (26%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
 #legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
 legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))

par(mfrow=c(1,2)) # 2 graphs in a row; phyPCA vs PACA
## in male skulls, results are very similar for phyPCA and PACA-GLS, 
plot(pca.GLS, phylo = TRUE, main = "phylo PCA")      ### Revell phy-PCA
plot(PaCA.gls, main = "PaCA using GLS", pch=22, cex = 1.5, bg = colors, phylo=TRUE)
PS.shape <- physignal(A=Y.gpa$coords, phy=tree5, iter=999)### phy sig in shape
summary(PS.shape) ### K-mult test, is identical for PCA, phylomorphospace, and phy-PCA
PS.shape$K.by.p # Phylogenetic signal profile
PS.size <- physignal(A=Y.gpa$Csize, phy=tree5, iter=999) ### phy sig in size
summary(PS.size) ### K-mult test, is identical for PCA, phylomorphospace, and phy-PCA
PS.size$K.by.p # Phylogenetic signal profile

par(mfrow=c(1,3)) # 3 graphs in a row
plot(pca.lands,  main = "PCA") ### traditional PCA
plot(pca.GLS, phylo = TRUE, main = "phylo PCA")      ### Revell phy-PCA 2009
plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS") ### Collyer & Adams 2021


####################################################
###### summary PCA analyses
colors <- c("blue", "white") 
colors <- colors[as.numeric(gfd$belly)]  # mark belly color
#colors <- c("red", "yellow", "green") 
#colors <- colors[as.numeric(gfd$svl)]  # mark body size
#colors <- c("red", "gray", "green") 
#colors <- colors[as.numeric(gfd$parity)]  # mark parity

par(mfrow=c(1,3)) # 3 graphs in a row
pca.lands <- gm.prcomp(A=Y.gpa$coords) ## traditional PCA
summary(pca.lands)
plot(pca.lands, main = "Traditional PCA", pch=22, cex = 1.5, bg = colors) ### traditional PCA
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (34%)", pos = 4, font = 2) ## TIPS
text(0, 0.95*par()$usr[4], labels = "PC2 (9%)", pos = 4, font = 2) ## TIPS
legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("viviparous", "oviparous", "unknown"), pch=22, pt.bg = unique(colors))

pca.GLS <- gm.prcomp(A=Y.gpa$coords, phy=tree5, GLS=TRUE) 
summary(pca.GLS)
plot(pca.GLS, main = "Phy-PCA (Revell 2009)", pch=22, cex = 1.5, bg = colors) ### Revell phy-PCA 2009
#plot(pca.GLS, phylo = TRUE, pch=22, cex = 1.5, bg = colors) # with species names and tree
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (42%)", pos = 4, font = 2) ## GLS
text(0, 0.95*par()$usr[4], labels = "PC2 (21%)", pos = 4, font = 2) ## GLS
legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("viviparous", "oviparous", "unknown"), pch=22, pt.bg = unique(colors))

PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = tree5, align.to.phy = TRUE, GLS = TRUE)
summary(PaCA.gls)
plot(PaCA.gls, main = "PaCA (Collyer & Adams 2021)", pch=22, cex = 1.5, bg = colors) ### Collyer & Adams 2021
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (71%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
text(0, 0.95*par()$usr[4], labels = "PC2 (26%)", pos = 4, font = 2) ## data aligned to phy sig, GLS rotation
legend("topleft", c("blue", "white"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))
#legend("topleft", c("viviparous", "oviparous", "unknown"), pch=22, pt.bg = unique(colors))

####################################################
### As in previous versions
colors <- c("black", "red", "yellow")
colors <- colors[as.numeric(gfd$svl)]  # mark body size

Y.gpa <- gpagen(dat, Proj=TRUE, ProcD=TRUE, curves=NULL, surface=NULL)
pca.lands2 <- gm.prcomp(A=Y.gpa$coords, phy=tree5) 

plot(pca.lands2, pch=22, cex = 1.5, bg = colors, axis1=1, axis2=2, phylo=F, 
   phylo.par=list(tip.labels=TRUE, node.labels=F, anc.states=F)) 
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 (34%)", pos = 4, font = 2) ## TIPS
text(0, 0.95*par()$usr[4], labels = "PC2 (9%)", pos = 4, font = 2) ## TIPS
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = c("black", "red", "yellow"))

plot(pca.lands, pch=22, cex = 1.5, bg = colors, axis1=3, axis2=4) 
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC3 (7%)", pos = 4, font = 2) ## TIPS
text(0, 0.95*par()$usr[4], labels = "PC4 (6%)", pos = 4, font = 2) ## TIPS
legend("topleft", c("Large", "Medium", "Small"), pch=22, pt.bg = unique(colors))


####################################################
# Phylogenetic SIGNAL & RATES & Procrustes ANOVAs
####################################################

par(mfrow=c(2,2)) # 4 graphs 
PS.shape <- physignal(A=Y.gpa$coords, phy=tree5, iter=999)### phy sig in shape
summary(PS.shape) ### K-mult test, is identical for all PCAs. Diff for each axis
plot(PS.shape)
#plot(PS.shape$PACA, phylo = TRUE)
PS.shape$K.by.p # Phylogenetic signal profile
plot(PS.shape$K.by.p) # Phylogenetic signal profile

PS.size <- physignal(A=Y.gpa$Csize, phy=tree5, iter=999) ### phy sig in size
summary(PS.size) ### K-mult test, is identical for all PCAs
plot(PS.size)
#plot(PS.size$PACA, phylo = TRUE)
PS.size$K.by.p # Phylogenetic signal profile
plot(PS.size$K.by.p) # Phylogenetic signal profile

 ### 3D plot with a phylogeny and time on the z-axis
 plot(PCA.w.phylo, time.plot = TRUE)
 plot(PCA.w.phylo, time.plot = TRUE, bg = "red", 
    phylo.par = list(tip.labels = TRUE, 
 tip.txt.cex = 2, edge.color = "blue", edge.width = 2))

##comparing rates of shape evolution of phylogenies from groups we designated above
compare.evol.rates(phy=tree5, A=Y.gpa$coords, gp=belly,  iter=999) #no patterm (effect = 0.6731, p=0.262)
compare.evol.rates(phy=tree5, A=Y.gpa$coords, gp=parity, iter=999) #no pattern (effect = -0.6659, p=0.771)
compare.evol.rates(phy=tree5, A=Y.gpa$coords, gp=size, iter=999)   #no pattern (effect = -0.0698, p=0.583)

## calculate morphological diparity between groups
morphol.disparity(f1=gfd$coords ~1, groups= ~ belly, iter=999) 
   # marginal pattern (0.0065 blue, 0.0097 white, p=0.059)
morphol.disparity(f1=coords ~ 1, groups= ~ parity, data=gfd, iter=999) 
   #marginal pattern (0.008 o, 0.004 u, 0.005 v, p=0.055 o vs v)
morphol.disparity(f1=gfd$coords ~ 1, groups= ~ size, data=gfd, iter=999) 
   # no pattern (0.008 L, 0.006 M, 0.007 S, all p > 0.2)

#### Procrustes anova with effect sizes (Z) and centroid size as covariate
gfd$logSize <- log(gfd$Csize) 
Csize.aov <- procD.lm(coords ~ Csize, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(Csize.aov)
belly1.aov <- procD.lm(coords ~ belly*Csize, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(belly1.aov) ### Z = NS for belly color, 3.3 for Csize, no interaction
belly2.aov <- procD.lm(coords ~ belly, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(belly2.aov) ### Z = 2.3 for belly color
belly3.aov <- procD.pgls(coords ~ belly*Csize, phy = tree5, SS.type = "III", data = gfd, iter=999, print.progress = FALSE) 
summary(belly3.aov) ### Z = 2.8 for belly color
plot(belly3.aov) ## check assumptions
belly4.aov <- procD.pgls(coords ~ belly, phy = tree5, SS.type = "III", data = gfd, iter=999, print.progress = FALSE) 
summary(belly4.aov) ### Z = 2.8 for belly color
plot(belly4.aov) ## check assumptions

size1.aov <- procD.lm(coords ~ size*Csize, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(size1.aov)   ### Z = NS for SVL, Csize, interaction
size2.aov <- procD.lm(coords ~ size, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(size2.aov)   ### Z = 3.3 for SVL
size3.aov <- procD.pgls(coords ~ size*Csize, phy = tree5, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(size3.aov)   ### Z = 3.3 for SVL, 1.8 for Csize, 3.4 for interxn
size4.aov <- procD.pgls(coords ~ size, phy = tree5, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(size4.aov)   ### Z = 2.4 for SVL
size5.aov <- procD.pgls(coords ~ Csize, phy = tree5, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(size5.aov)   ### Z = 2.4 for SVL
plot(size3.aov) ## check assumptions
plot(size4.aov) ## check assumptions

parity1.aov <- procD.lm(coords ~ parity*Csize, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(parity1.aov) ### Z = NS for parity, 3.4 for Csize, no interaction
parity2.aov <- procD.lm(coords ~ parity, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(parity2.aov) ### Z = 1.8 for parity
parity3.aov <- procD.pgls(coords ~ parity*Csize, phy = tree5, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(parity3.aov) ### Z = NS for parity, 2.5 for Csize, no interaction
parity4.aov <- procD.pgls(coords ~ parity, phy = tree5, SS.type = "III", data = gfd, print.progress = FALSE) 
summary(parity4.aov) ### Z = NS for parity
plot(parity3.aov) ## check assumptions
plot(parity4.aov) ## check assumptions

#######################
######## PCM ###
#######################

## Creating named factors with the groupings that showed relevant trends in the morphospace
## make sure your groupings are OK
belly #belly patches
parity #parity of each species
tree5 #the pruned tree


## PHYLOGENETIC SIGNAL ##
# Multivariate test for shape
PS.shape <- physignal(A=Y.gpa$coords,phy=tree5, iter=999)
# Test summary
summary(PS.shape) 
# Histogram
plot(PS.shape) 
## if K is less than 1, phylogenetic signal is low



##########################
### ADAPTIVE LANDSCAPE ###
##########################

infile<-read.csv("skull_chi4.csv", header=TRUE) ### pc sizes (long vs short for pc1-4, length, width, height)
tree5 <- read.tree("tree5")
plot(tree5)

###### belly patch analysis only
# Assign same order of tip labels, to avoid missmatchings in mvMORPH
rownames(infile)<-infile$X
infile <- infile[order(match(row.names(infile),
                                         tree5$tip.label)),]
# Confirm shape variables and tree have the same species & order
shp.var <- subset(infile[,9:12])
shp.varPhyPCA <- subset(infile[,13:18])
shp.varPACA <- subset(infile[,19:20])
row.names(shp.var)==tree5$tip.label

### REGIMES ###
## Specify regimes. Example using command paintSubTree: 

#  blue vs white
tree6 <- paintSubTree(tree5, node=3, state="Lv2", anc.state="Lv1", stem = TRUE)
tree6 <- paintSubTree(tree6, node=68, state="Lv2", anc.state="Lv1", stem = TRUE)
tree6 <- paintSubTree(tree6, node=17, state="Lv2", anc.state="Lv1", stem = TRUE)
tree6 <- paintSubTree(tree6, node=56, state="Lv2", anc.state="Lv1", stem = TRUE)
tree6 <- paintSubTree(tree6, node=54, state="Lv2", anc.state="Lv1", stem = TRUE)
plot(tree6)

#shp.var <- shp.varPhyPCA  ### repeat with PhyPCA or PhyPACA
#shp.var <- shp.varPACA  ### repeat with PhyPCA or PhyPACA

## BM: RATE MATRICES WITH MULTIPLE MEANS UNDER FULL MODEL ##
# Equal rate matrix with unequal means
BM1M.belly<-mvBM(tree6, shp.var, scale.height = T, model="BM1", 
                param=list(smean=FALSE))
# Proportional rate matrices 
BMPM.belly<-mvBM(tree6, shp.var, scale.height = T, 
                param=list(constraint="proportional",smean=FALSE))
# Independent rate matrices
BMMM.belly<-mvBM(tree6, shp.var, model="BMM", scale.height = T, 
                param=list(smean=FALSE))

# Compare the fitted models
resultsBMM.belly <- list(BM1M.belly,BMPM.belly,BMMM.belly); resultsBMM.belly
# Akaike weigths under AICc
weightsBMM.belly <- aicw(resultsBMM.belly, aicc=TRUE); weightsBMM.belly

## OU: FULL REGIMES for ventral coloration ##
# Random root: Diagonal SIGMA & free ALPHA (estimate 2 sigmas and 1 alpha)
OUMRDF.belly <- mvOU(tree6, shp.var, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",sigma="diagonal"))
# Random root: Free SIGMA & diagonal ALPHA
OUMRFD.belly <- mvOU(tree6, shp.var, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",alpha="diagonal"))
# Random root: Diagonal SIGMA & ALPHA
OUMRDD.belly <- mvOU(tree6, shp.var, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",sigma="diagonal",alpha="diagonal"))

# Compare the fitted models
resultsOUM.belly <- list(OUMRDF.belly,OUMRFD.belly, OUMRDD.belly )
# Akaike weigths under AICc
weightsOUM.belly <- aicw(resultsOUM.belly, aicc=TRUE)

## all BM and OU results together 
OU.belly <- list(OUMRDF.belly, OUMRFD.belly, OUMRDD.belly, BM1M.belly, BMPM.belly, BMMM.belly)
# Akaike weigths under AICc
weightsOU.belly <- aicw(OU.belly, aicc=TRUE)

###### size category analysis only
# Obtain shape variables (using 4 PCs)
# Confirm shape variables and tree have the same species & order
shp.var2 <- subset(infile[,9:12])
row.names(shp.var)==tree5$tip.label
#shp.var2 <- shp.varPhyPCA  ### repeat with PhyPCA or PhyPACA
#shp.var2 <- shp.varPACA  ### repeat with PhyPCA or PhyPACA

### REGIMES ###
## Specify regimes. Example using command paintSubTree: 
#  Small, Medium, Large

tree7<-tree5
tree7$edge.length<-
   tree7$edge.length/max(nodeHeights(tree7)[,2])*1  ### rescale to total length of 1.0
plot(tree7)
add.scale.bar(tree7)
tree7 <- paintBranches(tree7, edge=c(41,53,6,52,28,42,40,35,26,33,29,14,55,13,43,1,51,2,108,109,110,111,112,113,97,100,102,103,104,105,106,46,90,73,98,60,61,62,63,64,86,87,88,89,90,94,101,96,85,3), state="medium")
tree7 <- paintBranches(tree7, edge=c(47,9,37,32,56,38,48,24,23,44,45,36,27,25,30,31,39,34,22,107,99,81,82,83,84,92,93,95), state="large")
plot(tree7)
nodelabels(cex=.5, adj=c(.5,2))

## BM: RATE MATRICES WITH MULTIPLE MEANS UNDER FULL MODEL ##
# Equal rate matrix with unequal means
BM1M.size<-mvBM(tree7, shp.var2, scale.height = T, model="BM1", 
                param=list(smean=FALSE))
# Proportional rate matrices 
BMPM.size<-mvBM(tree7, shp.var2, scale.height = T, 
                param=list(constraint="proportional",smean=FALSE))
# Independent rate matrices
BMMM.size<-mvBM(tree7, shp.var2, model="BMM", scale.height = T, 
                param=list(smean=FALSE))

# Compare the fitted models
resultsBMM.size <- list(BM1M.size,BMPM.size,BMMM.size)
# Akaike weigths under AICc
weightsBMM.size <- aicw(resultsBMM.size, aicc=TRUE)

## OU: FULL REGIMES ##
# Random root: Diagonal SIGMA & free ALPHA
OUMRDF.size <- mvOU(tree7, shp.var2, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",sigma="diagonal"))
# Random root: Free SIGMA & diagonal ALPHA
OUMRFD.size <- mvOU(tree7, shp.var2, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",alpha="diagonal"))
# Random root: Diagonal SIGMA & ALPHA
OUMRDD.size <- mvOU(tree7, shp.var2, model="OUM", scale.height = T, 
               param=list(vcv="randomRoot",sigma="diagonal",alpha="diagonal"))

# Compare the fitted models
resultsOUM.size <- list(OUMRDF.size,OUMRFD.size, OUMRDD.size )
# Akaike weigths under AICc
weightsOUM.size <- aicw(resultsOUM.size, aicc=TRUE)

## all results together 
OU.size <- list(OUMRDF.size,OUMRFD.size, OUMRDD.size,BM1M.size,BMPM.size,BMMM.size)
# Akaike weigths under AICc
weightsOU.size <- aicw(OU.size, aicc=TRUE)


########################### 
#
# SLOUCH analysis for head dimension, PCs, and body size or color category
#
###########################
setwd("C:/Users/emart120/Dropbox (ASU)/MartinsLab_old/Scelop17_20ASU/CT Scan project (mostly pre-2022, JR & EM only)/analysis")

require(ape)
require(ouch)
require(car)
require(xtable)
require(ggplot2)
require(gridExtra)
require(phytools)
require(slouch)
require(phangorn)

tree2 <- read.tree("tree5")
tree2$edge.length<-
   tree2$edge.length/max(nodeHeights(tree2)[,2])*1  ### rescale to total length of 1.0
plot(tree2)
add.scale.bar(tree2)

sceldat <- read.csv("skull_chi4.csv")
rownames(sceldat) <- sceldat[,1]

################
### for body size categories, use this:
################
#ancestral reconstruction small/medium/large from Rivera et al 2021
#setting node labels as regimes where 0=small, 1=medium, 2=large
tree2$node.label <- as.factor(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,2,2,2,1,1,1,1,1,1,1,2,2,1,2,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1))
plot(tree2)
add.scale.bar(tree2)
tree2 <- paintBranches(tree2, edge=c(41,53,6,52,28,42,40,35,26,33,29,14,55,13,43,1,51,2,108,109,110,111,112,113,97,100,102,103,104,105,106,46,90,73,98,60,61,62,63,64,86,87,88,89,90,94,101,96,85,3), state="medium")
tree2 <- paintBranches(tree2, edge=c(47,9,37,32,56,38,48,24,23,44,45,36,27,25,30,31,39,34,22,107,99,81,82,83,84,92,93,95), state="large")
plot(tree2)
nodelabels(cex=.5, adj=c(.5,2))

#body size categories
sceldat$size2 <- 0
sceldat$size2[sceldat$size=="M"] <- 1
sceldat$size2[sceldat$size=="L"] <- 2
x <- sceldat$size2
y <- sceldat$SVL
################

################
### for ventral patches, use this instead:
################
#ancestral reconstruction of the evolution of blue/white from Ossip-Drahos et al 2018
#setting node labels as regimes where 1=white, 0=blue
tree2$node.label <- as.factor(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
tree2 <- paintBranches(tree2, edge=c(3,10,11,12,68,17,56,54), state=1, anc.state=0)
plot(tree2)
#color categories
sceldat$color2 <- 0
sceldat$color2[sceldat$color=="W"] <- 1
x <- sceldat$color2
y <- sceldat$SVL
################

#head dimensions; stHL = head.L / SVL
y.L <- sceldat$stHL
y.W <- sceldat$stHW
y.H <- sceldat$stHH
y.PC1 <- sceldat$tradPCA1
y.PC2 <- sceldat$tradPCA2
y.PhC1 <- sceldat$phyPCA1
y.PhC2 <- sceldat$phyPCA2
y.PAC1 <- sceldat$PACA1
y.PAC2 <- sceldat$PACA2

# puts the data in the same order as the phylogeny
xy <- data.frame(x, y, row.names=row.names(sceldat))
xy <- xy[tree2$tip.label, ]

xy.L <- data.frame(x, y.L, row.names=row.names(sceldat))
xy.L <- xy.L[tree2$tip.label, ]

xy.W <- data.frame(x, y.W, row.names=row.names(sceldat))
xy.W <- xy.W[tree2$tip.label, ]

xy.H <- data.frame(x, y.H, row.names=row.names(sceldat))
xy.H <- xy.H[tree2$tip.label, ]

xy.PC1 <- data.frame(x, y.PC1, row.names=row.names(sceldat))
xy.PC1 <- xy.PC1[tree2$tip.label, ]

xy.PC2 <- data.frame(x, y.PC2, row.names=row.names(sceldat))
xy.PC2 <- xy.PC2[tree2$tip.label, ]

xy.PhC1 <- data.frame(x, y.PhC1, row.names=row.names(sceldat))
xy.PhC1 <- xy.PhC1[tree2$tip.label, ]

xy.PhC2 <- data.frame(x, y.PhC2, row.names=row.names(sceldat))
xy.PhC2 <- xy.PhC2[tree2$tip.label, ]

xy.PAC1 <- data.frame(x, y.PAC1, row.names=row.names(sceldat))
xy.PAC1 <- xy.PAC1[tree2$tip.label, ]

xy.PAC2 <- data.frame(x, y.PAC2, row.names=row.names(sceldat))
xy.PAC2 <- xy.PAC2[tree2$tip.label, ]

##Plot a phenogram showing y trait
tipdata <- xy$y
names(tipdata) <- tree2$tip

tipdata.L <- xy.L$y.L
names(tipdata.L) <- tree2$tip

tipdata.W <- xy.W$y.W
names(tipdata.W) <- tree2$tip

tipdata.H <- xy.H$y.H
names(tipdata.H) <- tree2$tip

tipdata.PC1 <- xy.PC1$y.PC1
names(tipdata.PC1) <- tree2$tip

tipdata.PC2 <- xy.PC2$y.PC2
names(tipdata.PC2) <- tree2$tip

tipdata.PhC1 <- xy.PhC1$y.PhC1
names(tipdata.PhC1) <- tree2$tip

tipdata.PhC2 <- xy.PhC2$y.PhC2
names(tipdata.PhC2) <- tree2$tip

tipdata.PAC1 <- xy.PAC1$y.PAC1
names(tipdata.PAC1) <- tree2$tip

tipdata.PAC2 <- xy.PAC2$y.PAC2
names(tipdata.PAC2) <- tree2$tip

# Make factors for SLOUCH analysis
xy$x <- as.factor(xy$x)
slfact <- xy$x
names(slfact) <- tree2$tip

sceldat$OU1<-"OU1"
OU1 <- as.factor(sceldat$OU1)
names(OU1) <- tree2$tip


phylosig(tree2, x=tipdata.L, method="lambda", test=TRUE)
phylosig(tree2, x=tipdata.W, method="lambda", test=TRUE)
phylosig(tree2, x=tipdata.H, method="lambda", test=TRUE)


#########
## HEAD LENGTH
##########



## Brownian motion 
bm.model.L <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.L, hillclimb=TRUE)
summary(bm.model.L)
logLik(bm.model.L)
hillclimbplot(bm.model.L)

## OU1 
slmodel.OU1.L <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
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
## HEAD WIDTH
########
## Brownian motion 
bm.model.W <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.W, hillclimb=TRUE, lower=0.0005)
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
## HEAD HEIGHT
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
slmodel.H <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.H, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.H)			 
plot(slmodel.H)
summary(slmodel.H)
logLik(slmodel.H)
hillclimbplot(slmodel.H)

########
## PC1 
########
## Brownian motion
bm.model.PC1 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PC1, hillclimb=TRUE)
summary(bm.model.PC1)
plot(bm.model.PC1)
logLik(bm.model.PC1)
hillclimbplot(bm.model.PC1)

## OU1
slmodel.OU1.PC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PC1, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PC1)			 
plot(slmodel.OU1.PC1)
summary(slmodel.OU1.PC1)
logLik(slmodel.OU1.PC1)
hillclimbplot(slmodel.OU1.PC1)				

## OU multiple
slmodel.PC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PC1, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PC1)			 
plot(slmodel.PC1)
summary(slmodel.PC1)
logLik(slmodel.PC1)
hillclimbplot(slmodel.PC1)

########
## PC2 
########
## Brownian motion 
bm.model.PC2 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PC2, hillclimb=TRUE)
summary(bm.model.PC2)
plot(bm.model.PC2)
logLik(bm.model.PC2)
hillclimbplot(bm.model.PC2)

## OU1 
slmodel.OU1.PC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PC2, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PC2)			 
plot(slmodel.OU1.PC2)
summary(slmodel.OU1.PC2)
logLik(slmodel.OU1.PC2)
hillclimbplot(slmodel.OU1.PC2)				

## OU multiple
slmodel.PC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PC2, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PC2)			 
plot(slmodel.PC2)
summary(slmodel.PC2)
logLik(slmodel.PC2)
hillclimbplot(slmodel.PC2)

########
## PhC1 
########
## Brownian motion 
bm.model.PhC1 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PhC1, hillclimb=TRUE)
summary(bm.model.PhC1)
plot(bm.model.PhC1)
logLik(bm.model.PhC1)
hillclimbplot(bm.model.PhC1)

## OU1 
slmodel.OU1.PhC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PhC1, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PhC1)			 
plot(slmodel.OU1.PhC1)
summary(slmodel.OU1.PhC1)
logLik(slmodel.OU1.PhC1)
hillclimbplot(slmodel.OU1.PhC1)				

## OU multiple
slmodel.PhC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PhC1, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PhC1)			 
plot(slmodel.PhC1)
summary(slmodel.PhC1)
logLik(slmodel.PhC1)
hillclimbplot(slmodel.PhC1)

########
## PhC2 
########
## Brownian motion 
bm.model.PhC2 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PhC2, hillclimb=TRUE)
summary(bm.model.PhC2)
plot(bm.model.PhC2)
logLik(bm.model.PhC2)
hillclimbplot(bm.model.PhC2)

## OU1 
slmodel.OU1.PhC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PhC2, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PhC2)			 
plot(slmodel.OU1.PhC2)
summary(slmodel.OU1.PhC2)
logLik(slmodel.OU1.PhC2)
hillclimbplot(slmodel.OU1.PhC2)	

## OU multiple
slmodel.PhC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PhC2, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PhC2)			 
plot(slmodel.PhC2)
summary(slmodel.PhC2)
logLik(slmodel.PhC2)
hillclimbplot(slmodel.PhC2)


########
## PAC1 
########
## Brownian motion 
bm.model.PAC1 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PAC1, hillclimb=TRUE)
summary(bm.model.PAC1)
plot(bm.model.PAC1)
logLik(bm.model.PAC1)
hillclimbplot(bm.model.PAC1)

## OU1 
slmodel.OU1.PAC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PAC1, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PAC1)			 
plot(slmodel.OU1.PAC1)
summary(slmodel.OU1.PAC1)
logLik(slmodel.OU1.PAC1)
hillclimbplot(slmodel.OU1.PAC1)	

## OU multiple
slmodel.PAC1 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PAC1, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PAC1)			 
plot(slmodel.PAC1)
summary(slmodel.PAC1)
logLik(slmodel.PAC1)
hillclimbplot(slmodel.PAC1)

########
## PAC2 
########
## Brownian motion 
bm.model.PAC2 <- brown.fit(tree2, species=tree2$tip.label, 
                      sigma2_y_values = seq(0.1, 10, length.out = 80), 
                      response=tipdata.PAC2, hillclimb=TRUE)
summary(bm.model.PAC2)
plot(bm.model.PAC2)
logLik(bm.model.PAC2)
hillclimbplot(bm.model.PAC2)

## OU1 
slmodel.OU1.PAC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PAC2, fixed.fact=OU1, hillclimb=TRUE)
regimeplot(slmodel.OU1.PAC2)			 
plot(slmodel.OU1.PAC2)
summary(slmodel.OU1.PAC2)
logLik(slmodel.OU1.PAC2)
hillclimbplot(slmodel.OU1.PAC2)	

## OU multiple
slmodel.PAC2 <- slouch.fit(phy=tree2, species=tree2$tip.label, 
                      hl_values = seq(0.001, 40, length.out = 50),
                      vy_values = seq(0.001, 50, length.out = 50),
                      response=tipdata.PAC2, fixed.fact=slfact, hillclimb=TRUE)
regimeplot(slmodel.PAC2)			 
plot(slmodel.PAC2)
summary(slmodel.PAC2)
logLik(slmodel.PAC2)
hillclimbplot(slmodel.PAC2)


#################################
#################################
# BROWNIE ANALYSIS OF HEAD DIMENSIONS 
#################################
#################################

require(phytools)
tree5 <- read.tree("tree5")

pcdat<-read.csv("skull_chi4.csv", header=TRUE) 
         ### pc sizes (long vs short for pc1-4, length, width, height)

headL <- pcdat$head.L; names(headL) <- pcdat$species
headW <- pcdat$head.W; names(headW) <- pcdat$species
headH <- pcdat$head.H; names(headH) <- pcdat$species
color <- pcdat$color; names(color)  <- pcdat$species
size <- pcdat$size; names(size) <- pcdat$species

### can't do PCAs because Brownie analysis is singular
#tPCA1 <- pcdat$tradPCA1; names(tPCA1) <- pcdat$species
#tPCA2 <- pcdat$tradPCA2; names(tPCA2) <- pcdat$species
#pPCA1 <- pcdat$phyPCA1; names(pPCA1) <- pcdat$species
#pPCA2 <- pcdat$phyPCA2; names(pPCA2) <- pcdat$species
#PACA1 <- pcdat$PACA1; names(PACA1) <- pcdat$species
#PACA2 <- pcdat$PACA2; names(PACA2) <- pcdat$species
#####

tree.col <- make.simmap(tree5, color, model="SYM")
tree.size <- make.simmap(tree5, size, model="SYM")

brownie.lite(tree.col, headL, maxit=3000)
brownie.lite(tree.col, headW, maxit=3000)
brownie.lite(tree.col, headH, maxit=3000)
### can't do PCAs because Brownie analysis is singular
#brownie.lite(tree.col, tPCA1, maxit=3000)
#brownie.lite(tree.col, tPCA2, maxit=3000)
#brownie.lite(tree.col, pPCA1, maxit=3000)
#brownie.lite(tree.col, pPCA2, maxit=3000)
#brownie.lite(tree.col, PACA1, maxit=3000)
#brownie.lite(tree.col, PACA2, maxit=3000)

brownie.lite(tree.size, headL, maxit=3000)
brownie.lite(tree.size, headW, maxit=3000)
brownie.lite(tree.size, headH, maxit=3000)



#################################
#################################
# PGLS 
#################################
#################################

pcdat<-read.csv("skull_chi4.csv", header=TRUE) 
         ### pc sizes (long vs short for pc1-4, length, width, height)
tree5 <- read.tree("tree5")
plot(tree5)
tree<-tree5

l1<-sort(pcdat$species)
l2<-sort(tree$tip.label)
identical(l1,l2)
setdiff(l1,l2)

#give row names to data frame
rownames(pcdat) <- pcdat$species
headL <- pcdat$head.L/pcdat$SVL
headW <- pcdat$head.W/pcdat$SVL
headH <- pcdat$head.H/pcdat$SVL
svl <- pcdat$SVL

cor.test(pcdat$head.L,pcdat$tradPCA1)
cor.test((pcdat$head.L/pcdat$SVL),pcdat$tradPCA1)

cor.test(pcdat$head.L,pcdat$tradPCA1)
cor.test(pcdat$head.L,pcdat$tradPCA2)
cor.test(pcdat$head.L,pcdat$tradPCA3)
cor.test(pcdat$head.L,pcdat$tradPCA4)
cor.test(pcdat$head.W,pcdat$tradPCA1)
cor.test(pcdat$head.W,pcdat$tradPCA2)
cor.test(pcdat$head.W,pcdat$tradPCA3)
cor.test(pcdat$head.W,pcdat$tradPCA4)
cor.test(pcdat$head.H,pcdat$tradPCA1)
cor.test(pcdat$head.H,pcdat$tradPCA2)
cor.test(pcdat$head.H,pcdat$tradPCA3)
cor.test(pcdat$head.H,pcdat$tradPCA4)

cor.test(pcdat$head.L,pcdat$tradPCA1)
cor.test(pcdat$head.W,pcdat$tradPCA1)
cor.test(pcdat$head.H,pcdat$tradPCA1)


cor.test(headL,pcdat$tradPCA1)
cor.test(headW,pcdat$tradPCA1)
cor.test(headH,pcdat$tradPCA1)


cor.test(pcdat$head.L,pcdat$phyPCA1)
cor.test(pcdat$head.W,pcdat$phyPCA1)
cor.test(pcdat$head.H,pcdat$phyPCA1)

cor.test(pcdat$head.L,pcdat$phyPCA2)
cor.test(pcdat$head.W,pcdat$phyPCA2)
cor.test(pcdat$head.H,pcdat$phyPCA2)

cor.test(pcdat$head.L,pcdat$phyPCA3)
cor.test(pcdat$head.W,pcdat$phyPCA3)
cor.test(pcdat$head.H,pcdat$phyPCA3)

cor.test(pcdat$head.L,pcdat$phyPCA4)
cor.test(pcdat$head.W,pcdat$phyPCA4)
cor.test(pcdat$head.H,pcdat$phyPCA4)

cor.test(pcdat$head.L,pcdat$phyPCA5)
cor.test(pcdat$head.W,pcdat$phyPCA5)
cor.test(pcdat$head.H,pcdat$phyPCA5)

cor.test(pcdat$head.L,pcdat$phyPCA6)
cor.test(pcdat$head.W,pcdat$phyPCA6)
cor.test(pcdat$head.H,pcdat$phyPCA6)

cor.test(pcdat$head.L,pcdat$PACA1)
cor.test(pcdat$head.W,pcdat$PACA1)
cor.test(pcdat$head.H,pcdat$PACA1)

cor.test(pcdat$head.L,pcdat$PACA2)
cor.test(pcdat$head.W,pcdat$PACA2)
cor.test(pcdat$head.H,pcdat$PACA2)

#phylogenetic independent contrast
names(svl) <- names(headH) <- names(headW) <- names(headL) <- rownames(pcdat)
pic.headL <- pic(headL, tree) # (headL was already divided by SVL)
pic.headW <- pic(headW, tree)
pic.headH <- pic(headH, tree)
pic.headL <- pic(pcdat$head.L, tree)  # (headL was not divided by SVL)
pic.headW <- pic(pcdat$head.W, tree)
pic.headH <- pic(pcdat$head.H, tree)
pic.svl <- pic(svl, tree)

#corelation test on PIC 
cor.test(pic.svl, pic.headL)
cor.test(pic.svl, pic.headW)
cor.test(pic.svl, pic.headH)

#give row names to data frame
rownames(pcdat) <- pcdat$species
PC1 <- pcdat$tradPCA1
PC2 <- pcdat$tradPCA2
PC3 <- pcdat$tradPCA3
PC4 <- pcdat$tradPCA4

PC1 <- pcdat$phyPCA1
PC2 <- pcdat$phyPCA2
PC3 <- pcdat$phyPCA3
PC4 <- pcdat$phyPCA4
PC5 <- pcdat$phyPCA5
PC6 <- pcdat$phyPCA6

PC1 <- pcdat$PACA1
PC2 <- pcdat$PACA2

#phylogenetic independent contrast
names(PC4) <- names(PC3) <- names(PC2) <- names(PC1) <- rownames(pcdat)
pic.PC1 <- pic(PC1, tree5)
pic.PC2 <- pic(PC2, tree5)
pic.PC3 <- pic(PC3, tree5)
pic.PC4 <- pic(PC4, tree5)
pic.PC5 <- pic(PC5, tree5)
pic.PC6 <- pic(PC6, tree5)

#corelation test on PIC
cor.test(pic.svl, pic.PC1)
cor.test(pic.svl, pic.PC2)
cor.test(pic.svl, pic.PC3)
cor.test(pic.svl, pic.PC4)

cor.test(pic.headL, pic.PC1)
cor.test(pic.headW, pic.PC1)
cor.test(pic.headH, pic.PC1)

cor.test(pic.headL, pic.PC2)
cor.test(pic.headW, pic.PC2)
cor.test(pic.headH, pic.PC2)

cor.test(pic.headL, pic.PC3)
cor.test(pic.headW, pic.PC3)
cor.test(pic.headH, pic.PC3)

cor.test(pic.headL, pic.PC4)
cor.test(pic.headW, pic.PC4)
cor.test(pic.headH, pic.PC4)

cor.test(pic.headL, pic.PC5)
cor.test(pic.headW, pic.PC5)
cor.test(pic.headH, pic.PC5)

cor.test(pic.headL, pic.PC6)
cor.test(pic.headW, pic.PC6)
cor.test(pic.headH, pic.PC6)

X <- cbind(svl, pcdat$head.W)
phylomorphospace(tree2,X[,1:2],xlab="trait 1",ylab="trait 2")
### acanthinus or angustus is the outlier? 

cor.test(pcdat$SVL, pcdat$head.L)
cor.test(pcdat$SVL, pcdat$head.W)
cor.test(pcdat$SVL, pcdat$head.H)
cor.test(pcdat$SVL, pcdat$tradPCA1)
cor.test(pcdat$SVL, pcdat$tradPCA2)
cor.test(pcdat$SVL, pcdat$phyPCA1)
cor.test(pcdat$SVL, pcdat$phyPCA2)
cor.test(pcdat$SVL, pcdat$PACA1)
cor.test(pcdat$SVL, pcdat$PACA2)

library(nlme)
Nsp <- length(tree2$tip.label)
y <- pcdat$SVL
x <- pcdat$head.L  #### vary this, raw measures
x <- headH  #### vary this, already divided by SVL
xy <- data.frame(x, y, row.names=row.names(pcdat))
xy <- xy[tree2$tip.label, ]
TIPS <- gls(y~x, data = xy)
#summary(TIPS)
TIPSr <- cor(xy$x,xy$y)
bm.1<-corBrownian(phy=tree)
bm.gls<-gls(y~x,correlation=bm.1, data=xy)
#summary(bm.gls)
V<-corMatrix(Initialize(bm.1,xy)); 
a.Y <- matrix(1,1,Nsp) %*% solve(V) %*% xy$y/sum(solve(V)); 
a.X <- matrix(1,1,Nsp) %*% solve(V) %*% xy$x/sum(solve(V)); 
FICr <- (xy$y-a.Y) %*% solve(V) %*% (xy$x-a.X)/sqrt(((xy$y-a.Y) %*% 
   solve(V) %*% (xy$y-a.Y))*((xy$x-a.X) %*% solve(V) %*% (xy$x-a.X))); 
ou.2 <- corMartins(1, phy = tree, fixed = FALSE)
ou.gls <- gls(y~x,correlation = ou.2, data = xy)
#summary(ou.gls)
V<-corMatrix(Initialize(ou.2,xy)); 
a.Y <- matrix(1,1,Nsp) %*% solve(V) %*% xy$y/sum(solve(V));  
a.X <- matrix(1,1,Nsp) %*% solve(V) %*% xy$x/sum(solve(V)); 
PGLSr <-(xy$y-a.Y) %*% solve(V) %*% (xy$x-a.X)/sqrt(((xy$y-a.Y) %*% 
   solve(V) %*% (xy$y-a.Y))*((xy$x-a.X) %*% solve(V) %*% (xy$x-a.X))); 
TIPSr; FICr[1]; PGLSr[1]
cor.test(x,y)


################################################################  
################################################################  
#### Divergent Sympatry analysis
################################################################  
################################################################  
# this next bit creates a csv file called sp_dist 
# that has all the info for the chi-square.
# no need to run if you already have the chi-square values
# and could start below where we begin to plot.
# If you don't have the chi-sqaure values, run the following: 
################################################################  
################################################################  

## first get the PC axes (3 types) and put in 
## three types

PCframe <- function(pca, CutOff) {  # default cutoff of 80%
  EigenSum<-k<-0
  repeat { 
     k<-k+1
     EigenSum<-EigenSum + pca$d[k] / sum(pca$d)
     if(EigenSum>=CutOff/100) {
       break}}
  PCAframe<-data.frame(pca$x[,1:k])
  return(PCAframe)}

pca.lands <- gm.prcomp(A=Y.gpa$coords) ## traditional PCA
summary(pca.lands)
pca.lands$d[1]/sum(pca.lands$d) # provides the Eigenvalue for the first principal component.
tradPCA <- PCframe(pca.lands,56)

pca.GLS <- gm.prcomp(A=Y.gpa$coords, phy=tree5, GLS=TRUE) 
summary(pca.GLS)
phyPCA <- PCframe(pca.GLS,80)

PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = tree5, align.to.phy = TRUE, GLS = TRUE)
summary(PaCA.gls)
PACA <- PCframe(PaCA.gls,90)

PCRes<-cbind(tradPCA, phyPCA, PACA)
tmp<- c("tradPCA1", "tradPCA2", "tradPCA3", "tradPCA4", 
         "phyPCA1", "phyPCA2", "phyPCA3", "phyPCA4", "phyPCA5", "phyPCA6", 
         "PACA1", "PACA2") 
colnames(PCRes) <- tmp

#write.csv(PCRes, file="skull_chi3.csv")
# For divergent sympatry analyses (skull_chi3.csv)
# manually drop acanthinus & angustus - only for divergent sympatry analyses
# add column for 3-letter species names
# use formula to identify long and short for each axis (above or below zero)
# add head LWH, divide each by SVL, sort by above vs below average

# For other analyses (skull_chi4.csv)
# add head LWH, divide each by SVL, sort by above vs below average
lin.dat <- read.csv("linear mes.csv")
temp<-cbind(lin.dat, PCRes)
temp$stHL <- temp$head.L / temp$SVL
temp$stHW <- temp$head.W / temp$SVL
temp$stHH <- temp$head.H / temp$SVL
temp$stHLc<- "short"; temp$stHLc[temp$stHL >  mean(temp$stHL)]<- "long"
temp$stHWc<- "short"; temp$stHWc[temp$stHW >  mean(temp$stHW)]<- "long"
temp$stHHc<- "short"; temp$stHHc[temp$stHH >  mean(temp$stHH)]<- "long"
write.csv(temp, file="skull_chi4.csv")

########################################################################
### start here after PCAs are run and data file is ready
########################################################################
setwd("C:/Users/emart120/Dropbox (ASU)/MartinsLab_old/Scelop17_20ASU/CT Scan project/analysis/divergent sympatry")
#infile<-read.table("outputsizes.txt", sep = " ", header=FALSE, skip = 10) ## species names and body size categories
infile<-read.csv("skull_chi3.csv", header=TRUE) ### pc sizes (long vs short for pc1-4, length, width, height)

setwd("./snout")

########################################################################
allfiles <-list.files("./")
#### this loop sets up a new data frame combining skull data (infile) and sympatric species info (allfiles)
#### change variable of interest in this loop
sp_dist<-list() #create empty list
for (i in (1:length(allfiles))) {
  species<-read.table(allfiles[i], sep ="\t", header=FALSE, skip=1) 
  #species[,2][is.na(data[,2])] = 0
  newspecies<-species[1,] # read in new species
  newspecies$species<-infile$species[i]  
  newspecies$temp<-infile$tPCA1[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$tPCA2[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$tPCA3[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$tPCA4[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA1[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA2[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA3[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA4[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA5[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$pPCA6[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$PACA1[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
#  newspecies$temp<-infile$PACA2[i]      ### change variable name here (tPCA1 - tPCA4, pPCA1 - pPCA6, PACA1, PACA2)
  sp_dist<-rbind(sp_dist, newspecies) 
  #sp_dist$species[i]<-allfiles[i]
}
sp_dist$V1<-as.numeric(gsub(",", "", sp_dist$V1)) ## distance (km)
sp_dist$V2<-as.numeric(gsub(",", "", sp_dist$V2)) ## total pairs of species
sp_dist$V3<-as.numeric(gsub(",", "", sp_dist$V3)) ## same species pairs
sp_dist$V4<-as.numeric(gsub(",", "", sp_dist$V4)) ## diff species pairs
sp_dist$V5<-as.numeric(gsub(",", "", sp_dist$V5)) ## number diff species                       
sp_dist$V7<-as.numeric(gsub(",", "", sp_dist$V7)) ## number diff size species
sp_dist$V8<-as.numeric(gsub(",", "", sp_dist$V8)) ## number diff size species pairs
chi2<-list(); sign<-list(); Nsp <- length(sp_dist$species)

########################################################################
NSmall <- length(sp_dist$species[sp_dist$temp == "short"]) ### change variable name here
NLarge <- length(sp_dist$species[sp_dist$temp == "long"])  ### change variable name here
#### Start for loop
for (i in (1:nrow(sp_dist))) {
  obs_dif<-sp_dist$V7[i]   ## number diff size species
  obs_total<-sp_dist$V5[i] ## number diff species
  obs_same<-(obs_total - obs_dif)
  if (sp_dist$temp[i] == "short") {       ### change variable name here
    exp_dif<-(NLarge/Nsp); exp_same<-(NSmall/Nsp) }
  if (sp_dist$temp[i] == "long") {        ### change variable name here
    exp_dif<-(NSmall/Nsp); exp_same<-(NLarge/Nsp) }
  newchi2 <- ((((obs_dif-(exp_dif*obs_total))^2)/(exp_dif*obs_total))+
              (((obs_same-(exp_same*obs_total))^2)/(exp_same*obs_total)))
  if (exp_dif*obs_total > obs_dif) {
    newchi2=newchi2*(-1)  }
  sp_dist$chi2[i]<-newchi2	}  # end for loop
########################################################################
sp2 <- sp_dist[-c(11)]
sp2$tPCA1_chi2 <- sp_dist$chi2
#sp2$tPCA2_chi2 <- sp_dist$chi2
#sp2$tPCA3_chi2 <- sp_dist$chi2
#sp2$tPCA4_chi2 <- sp_dist$chi2
#sp2$pPCA1_chi2 <- sp_dist$chi2
#sp2$pPCA2_chi2 <- sp_dist$chi2
#sp2$pPCA3_chi2 <- sp_dist$chi2
#sp2$pPCA4_chi2 <- sp_dist$chi2
#sp2$pPCA5_chi2 <- sp_dist$chi2
#sp2$pPCA6_chi2 <- sp_dist$chi2
#sp2$PACA1_chi2 <- sp_dist$chi2
#sp2$PACA2_chi2 <- sp_dist$chi2

#setwd("C:/Users/emart120/Dropbox (ASU)/MartinsLab_old/Scelop17_20ASU/CT Scan project/analysis/divergent sympatry")
#write.csv(sp2, file="chi-square_results2.csv")

### Need randomization test
###  select taxa at random from (1) genus as a whole, (2) phylogenetic relatives, (3) geographic neighbors
###  to be sympatric congeners. Then run chi-sq analysis




####################################
####################################
# start of plotting 
####################################
####################################

setwd("C:/Users/emart120/Dropbox (ASU)/MartinsLab_old/Scelop17_20ASU/CT Scan project/analysis/divergent sympatry")
sp_dist <- read.csv("chi-square_results2.csv")
infile <- read.csv("skull_chi3.csv", header=TRUE) ### pc sizes (long vs short for 3 types of PCs, length, width, height)
dat <- cbind(infile, sp_dist)

par(mfrow=c(2,2))

plot(dat$tradPCA1, dat$tPCA1_chi2, main=paste("Traditional PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-2, 15), las=2)
abline(0,0,col='red')

plot(dat$tradPCA2, dat$tPCA2_chi2, main=paste("Traditional PCA2", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-2, 15), las=2)
abline(0,0,col='red')

plot(dat$tradPCA3, dat$tPCA3_chi2, main=paste("Traditional PCA3", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-2, 15), las=2)
abline(0,0,col='red')

plot(dat$tradPCA4, dat$tPCA4_chi2, main=paste("Traditional PCA4", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-2, 15), las=2)
abline(0,0,col='red')

par(mfrow=c(1,2))

plot(dat$phyPCA1, dat$pPCA1_chi2, main=paste("Phylogenetic PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')

plot(dat$phyPCA2, dat$pPCA2_chi2, main=paste("Phylogenetic PCA2", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')

plot(dat$phyPCA3, dat$pPCA3_chi2, main=paste("Phylogenetic PCA3", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')

plot(dat$phyPCA4, dat$pPCA4_chi2, main=paste("Phylogenetic PCA4", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')

plot(dat$phyPCA5, dat$pPCA5_chi2, main=paste("Phylogenetic PCA5", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')

plot(dat$phyPCA6, dat$pPCA6_chi2, main=paste("Phylogenetic PCA6", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 20), las=2)
abline(0,0,col='red')


par(mfrow=c(1,2))

plot(dat$phyPACA1, dat$PACA1_chi2, main=paste("Phylogenetically aligned PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 50), las=2)
abline(0,0,col='red')

plot(dat$phyPACA2, dat$PACA2_chi2, main=paste("Phylogenetically aligned PCA2", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 50), las=2)
abline(0,0,col='red')


par(mfrow=c(1,3))
plot(dat$tradPCA1, dat$tPCA1_chi2, main=paste("Traditional PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 50), las=2)
abline(0,0,col='red')

plot(dat$phyPCA1, dat$pPCA1_chi2, main=paste("Phylogenetic PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 50), las=2)
abline(0,0,col='red')

plot(dat$phyPACA1, dat$PACA1_chi2, main=paste("Phylogenetically aligned PCA1", sep=" "), 
      ylab="Divergent Sympatry (chi-sq)", xlab="PCA scores", ylim = c(-5, 50), las=2)
abline(0,0,col='red')


########################################################################################
### Bayesian Clustering analyses (Cadena et al. 2017)
###    Use clustvarsel (Scrucca and Raftery 2004) to reduce dimensionality in ways that are
###    useful to group discrimination, without a priori information about groups 
###    (Raftery and Dean 2006; Maugis et al. 2009a, 2009b). Use mclust 5.0 (Scrucca et al. 2016) 
###    to fit a wide range of NMMs. At one extreme, NMMs assuming one morphological group. 
###    At the other extreme, assume many groups.
########## mclust 

library(mclust)
library(mvtnorm)
library(ellipse)
library(clustvarsel)

dat <- read.csv("skull_chi4.csv", header=TRUE) 
#keep <- c("head.L", "head.W", "head.H", "SVL", "tradPCA1", "tradPCA2", "tradPCA3", "tradPCA4")
#keep <- c("tradPCA1", "tradPCA2", "tradPCA3", "tradPCA4")
#keep <- c("phyPCA1", "phyPCA2", "phyPCA3", "phyPCA4", "phyPCA5", "phyPCA6")
#keep <- c("head.L", "head.W", "head.H", "SVL")
keep <- c("head.L", "head.W", "head.H")
X <- dat[keep]/dat$SVL
clPairs(X, dat$size)
class <- dat$size
table(class)

#perform logarithmic transformation of the traits of interest:
#morpho.data <- dat
#morpho.data.ln <- log(dat[,c(3:6)])
morpho.data.ln <- log(X[,c(1:3)])
#morpho.data.ln <- (X[,c(1:6)]) ### or not
#change column names and examine resulting data frame
#colnames(morpho.data.ln)
#colnames(morpho.data.ln) <- c("Ln (Head length)", "Ln (Head width)", "Ln (Head height)", "Ln (SVL)")
colnames(morpho.data.ln) <- c("Ln (Head length)", "Ln (Head width)", "Ln (Head height)")
#colnames(morpho.data.ln) <- c("TradPCA1", "TradPCA2", "TradPCA3", "TradPCA4")
#colnames(morpho.data.ln) <- c("PhyPCA1", "PhyPCA2", "PhyPCA3", "PhyPCA4", "PhyPCA5", "PhyPCA6")
dim(morpho.data.ln)
summary(morpho.data.ln)

### or not
# perform PCA 
morpho.data.ln.pca <- prcomp(morpho.data.ln, center = T, scale. = F) #PCA using the covariance matrix
attributes(morpho.data.ln.pca)
morpho.data.ln.pca$scale
morpho.data.ln.pca$center
summary(morpho.data.ln.pca)
summary(morpho.data.ln.pca$x)
morpho.data.ln.pca$rotation

#Gaussian mixture models
#backward variable selection 
#mclust.options() #check Mclust options
OptMc <- mclust.options() #save default
mclust.options(hcUse="VARS") #change default as needed
morpho.data.ln.varsel.back <- clustvarsel(morpho.data.ln, G=1:30, search=c("greedy"), direction = c("backward"))
#attributes(morpho.data.ln.varsel.back)
summary(morpho.data.ln.varsel.back)
morpho.data.ln.varsel.back$subset
morpho.data.ln.varsel.back$steps.info
morpho.data.ln.varsel.back$search
morpho.data.ln.varsel.back$direction


#forward variable selection using logarithmic transformation of the six traits of interest
#(Wing, Tail, Blength, Bdepth, Bwidth and Tarsus), and examine results:
#mclust.options() #check Mclust options
#OptMc <- mclust.options() #save default
mclust.options(hcUse="VARS") #change default as needed
morpho.data.ln.varsel.for <- clustvarsel(morpho.data.ln, G=1:30, search=c("greedy"), direction = c("forward"))
#attributes(morpho.data.ln.varsel.for)
summary(morpho.data.ln.varsel.for)
morpho.data.ln.varsel.for$subset
morpho.data.ln.varsel.for$steps.info
morpho.data.ln.varsel.for$search
morpho.data.ln.varsel.for$direction


#Run mclust analysis with variables chosen
Mcluster.morpho.data.ln.subset <- Mclust(morpho.data.ln[,c(1,3)], G=1:30)
summary(Mcluster.morpho.data.ln.subset)
attributes(Mcluster.morpho.data.ln.subset)
summary(Mcluster.morpho.data.ln.subset$data)
dim(Mcluster.morpho.data.ln.subset$data)


#extract BIC values for the best model conditional on the number of groups
BIC.Best.Model.Per.G <- apply(Mcluster.morpho.data.ln.subset$BIC, 1, max, na.rm=T)
#BIC.Best.Model.Per.G <- apply(Mcluster.morpho.data.ln.pca.subset$BIC, 1, max, na.rm=T)
#Mcluster.morpho.data.ln.pca.subset$BIC

BIC <- mclustBIC(X)
plot(BIC)

mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)

plot(mod1, what = "classification")

table(dat$size, mod1$classification)
plot(mod1, what = "uncertainty")

ICL <- mclustICL(X)
summary(ICL)
plot(ICL)

LRT <- mclustBootstrapLRT(X, modelName = "VVE")
LRT

X<- morpho.data.ln[-2]
X<- morpho.data.ln
### identify number of groups and best model with no prior information
par(mfrow=c(1,1))
BIC <- mclustBIC(X)
plot(BIC)
summary(BIC)

###  calculate likelihoods and uncertainties
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification") ## groups with no prior information
table(class, mod1$classification)   ## do these groups match "class"? 
plot(mod1, what = "uncertainty")
ICL <- mclustICL(X)
summary(ICL)
plot(ICL)
LRT <- mclustBootstrapLRT(X, modelName = "EEE")
LRT

### compare fit of null model (1 cluster) with inferred model (2 or 3 clusters)
mod1 <- Mclust(X, G = 1, x = BIC)  ## 1 cluster
summary(mod1, parameters = TRUE)
mod1 <- Mclust(X, G = 2, x = BIC)  ## 2 clusters
summary(mod1, parameters = TRUE)
mod1 <- Mclust(X, G = 3, x = BIC)  ## 3 clusters
summary(mod1, parameters = TRUE)

### extract which data points are in which cluster
results <- data.frame(dat, cluster=mod1$classification)
results$sp[results$cluster==1]
results$sp[results$cluster==2]
results$sp[results$cluster==3]

##### EM algorithm to confirm fit
(hc1 <- hc(X, modelName = "VVV", use = "SVD"))
BIC1 <- mclustBIC(X, initialization = list(hcPairs = hc1)) # ellipsoidal, varying volume, shape and orientation
summary(BIC1)
(hc2 <- hc(X, modelName = "VVV", use = "VARS"))
BIC2 <- mclustBIC(X, initialization = list(hcPairs = hc2)) ### no standardization of variables
summary(BIC2)
(hc3 <- hc(X, modelName = "EEE", use = "SVD"))
BIC3 <- mclustBIC(X, initialization = list(hcPairs = hc3)) ### ellipsoidal, equal volume, shape and orientation
summary(BIC3)
BIC <- mclustBICupdate(BIC1, BIC2, BIC3)
summary(BIC)
plot(BIC)

### clustering and discriminant analysis
mod2 <- MclustDA(X, class, modelType = "EDDA")  ### single component with same covariance structure
summary(mod2)
plot(mod2, what = "scatterplot")
plot(mod2, what = "classification")
mod3 <- MclustDA(X, class) ### multiple components with diff covariance structures
summary(mod3)
plot(mod3, what = "scatterplot")
plot(mod3, what = "classification")
cv <- cvMclustDA(mod2, nfold = 10)
str(cv)
unlist(cv[3:6])
cv <- cvMclustDA(mod3, nfold = 10)
str(cv)
unlist(cv[3:6])


### density
mod4 <- densityMclust(X)
summary(mod4)
plot(mod4, what = "BIC")
plot(mod4, what = "density", data = X, breaks = 15)
plot(mod4, what = "diagnostic", type = "cdf")
plot(mod4, what = "diagnostic", type = "qq")

mod5 <- densityMclust(X)
summary(mod5)
plot(mod5, what = "BIC")
plot(mod5, what = "density", type = "hdr", data = X, points.cex = 0.5)
plot(mod5, what = "density", type = "persp")

### bootstrap
boot1 <- MclustBootstrap(mod1, nboot = 999, type = "bs")
summary(boot1, what = "se")
summary(boot1, what = "ci")
plot(boot1, what = "pro")
plot(boot1, what = "mean")

boot4 <- MclustBootstrap(mod4, nboot = 999, type = "bs")
summary(boot4, what = "se")
summary(boot4, what = "ci")
plot(boot4, what = "pro")
plot(boot4, what = "mean")

### dimension reduction
mod1dr <- MclustDR(mod1)
summary(mod1dr)
plot(mod1dr, what = "pairs")
plot(mod1dr, what = "boundaries", ngrid = 200)
mod1dr <- MclustDR(mod1, lambda = 1)
summary(mod1dr)
plot(mod1dr, what = "scatterplot")
plot(mod1dr, what = "boundaries", ngrid = 200)

### classification
mod2dr <- MclustDR(mod2)
summary(mod2dr)
plot(mod2dr, what = "scatterplot")
plot(mod2dr, what = "boundaries", ngrid = 200)
mod3dr <- MclustDR(mod3)
summary(mod3dr)
plot(mod3dr, what = "scatterplot")
plot(mod3dr, what = "boundaries", ngrid = 200)
