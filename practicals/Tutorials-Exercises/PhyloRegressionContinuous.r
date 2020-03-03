#Analyses for Tutorial 04: PCMs

library(phytools)

#A quick simulation
mytree<- pbtree(n=50, scale=1) #one way to simulate a tree
  #note: in R, one has functions to simulate BD trees, random splits trees, etc. using different functions
plot(mytree)
X<-fastBM(tree=mytree) #simulates a continuous trait on phylogeny under Brownian motion
Y<-fastBM(tree=mytree)
cor(X,Y) #Not zero

#Now, what happens if we do this many times?
nsim=1000
trees<-pbtree(n=50, scale=1, nsim = nsim) #100 PB trees
X<-lapply(1:nsim, function(j) fastBM(tree=trees[[j]]))
Y<-lapply(1:nsim, function(j) fastBM(tree=trees[[j]]))
cor.raw<-unlist(lapply(1:nsim, function(j) cor(X[[j]],Y[[j]])))
hist(cor.raw, xlim=c(-1,1))  


#Now condition on phylogeny:
X.pic<-lapply(1:nsim, function(j) pic(x = X[[j]],phy = trees[[j]]))
Y.pic<-lapply(1:nsim, function(j) pic(x = Y[[j]],phy = trees[[j]]))
cor.pic<-unlist(lapply(1:nsim, function(j) cor(X.pic[[j]],Y.pic[[j]])))
hist(cor.pic, xlim=c(-1,1)) #mean near zero, but LOOK at variation!

library(ape)	
library(geiger)

#Here is a large time-dated molecular phylogeny (a chronogram):
ManderTree<-read.tree("TutorialData/Mander.tre",tree.names=T)
plot(ManderTree)

#Now read in data for species in the genus Hydromantes
Hyd.dat<-read.csv("TutorialData/HydromantesMns.csv", header=TRUE, row.names = 1)
Hyd.dat  #Notice we read in the first row as row.names. These MUST match the names in the phylogeny

#Now, prune the tree to match the data and vice-versa:
Hydro.new<-treedata(phy=ManderTree,data = Hyd.dat, warnings=FALSE)
plot(Hydro.new$phy)  #We have matched data to the tree!

#### 1: Phylogenetically Independent Contrasts
#Read in data, phylogeny, and match/prune them to one another
ManderTree<-read.tree("TutorialData/Mander.tre",tree.names=T)
   plot(ManderTree)
Mander.dat<-read.csv("TutorialData/PlethodonMns.csv", header=TRUE, row.names = 1)
Gps<-Mander.dat[,12]; names(Gps)<-row.names(Mander.dat)
Mander.dat<-Mander.dat[,-12]  #must remove prior to 'treedata', which does not like factors

Pleth.new<-treedata(phy=ManderTree,data = Mander.dat, warnings=FALSE)
   plot(Pleth.new$phy)  
SVL.new<-Pleth.new$data[,1] #grabs appropriate column from data matrix after treedata
HL.new<-Pleth.new$data[,3]
  match(names(Gps),row.names(Pleth.new$data))  #double check that the factor matches
   
#Phylogenetically-naive analyses
plot(SVL.new,HL.new)
cor.test(SVL.new,HL.new)  #correlation of two traits: pretty high
anova(lm(HL.new~SVL.new))

#Generate PICs and test while conditioning on phylogeny
X.pic<-pic(x = SVL.new, phy = Pleth.new$phy) 
Y.pic<-pic(x = HL.new, phy = Pleth.new$phy) 
cor.test(X.pic,Y.pic)
plot(X.pic,Y.pic)
anova(lm(Y.pic~X.pic+0))   #regression through the origin (b/c order of taxa for contrast is arbitrary [can 'spin' on node])
summary(lm(Y.pic~X.pic+0))   #coefficients of the model


#### 2: Phylogenetic Generalized Least Squares
#To perform PGLS in R, we must first estimate the phylogenetic covariance matrix V (C in matrix form):
V<-corBrownian(phy=Pleth.new$phy)
C <- vcv.phylo(phy = Pleth.new$phy)

#Now we run the analysis:
library(nlme)
bm.gls<-gls(HL.new~SVL.new,correlation = V, data=data.frame(SVL.new, HL.new))
summary(bm.gls)    
anova(bm.gls)


#### 3: Phylogenetic Transform
library(RRPP)
rdf<-rrpp.data.frame(SVL=as.matrix(SVL.new),HL=as.matrix(HL.new),Gps=as.matrix(Gps), C=C)
res.PhyT<-lm.rrpp(HL~SVL, data = rdf, Cov = C, print.progress = FALSE)
anova(res.PhyT)    #Same F-ratio as before (note: p-value found using permutation in this case)
coef(res.PhyT)

#### Phylogenetic ANOVA
#Phylogenetic anova
aov.phylo(SVL.new~Gps, phy = Pleth.new$phy)
anova(lm.rrpp(SVL~Gps, data = rdf, Cov = C, print.progress = FALSE))
anova(gls(SVL.new~Gps,correlation = corBrownian(phy=Pleth.new$phy), data=data.frame(SVL.new, Gps)))  #identical to Phylo.transform


####Phylogenetic Signal
#Phylogenetic Signal
phylosig(tree=Pleth.new$phy, x=SVL.new, method="K", test=T, nsim=1000)  #phytools
  res<-physignal(A=SVL.new,phy = Pleth.new$phy, print.progress = FALSE)   #geomorph
summary(res)  #the same, but the latter method (from geomorph) can accomodate multivariate data
plot(res)

####Additional Simulation Approaches

#Simulate correlated data
library(MASS)
R<-matrix(0.7, nrow=2,ncol=2); diag(R)<-1
dat.sim<-mvrnorm(n=10000,Sigma = R,mu = c(0,0))
plot(dat.sim, asp=1)

#Simulate BM correlated data on phylogeny
mytree<- pbtree(n=100, scale=1) #one way to simulate a tree
R<-matrix(0.7, nrow=2,ncol=2); diag(R)<-1
dat.BM<-sim.char(phy = mytree,par = R,nsim = 1)[,,1]
plot(dat.BM, asp=1)

##Practical Questions

#To apply what you have learned, please do the following:

#1: Simulate a phylogeny (or read one into R)
#2: Simultate two continuous traits with some known correlation between them
#3: Estimate the evolutionary association of the traits, conditioned on the phylogeny, using one of the methods you have learned

#Upload your code to the assignment center on Canvas 
