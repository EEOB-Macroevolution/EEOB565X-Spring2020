#Analyses for Tutorial 05: Discrete PCMs

## 1: Read data, tree, and prune/match one to the other
library(geiger) #also loads ape
library(phytools)

#Here is a large time-dated molecular phylogeny (a chronogram):
tree<-read.tree("TutorialData/tree.64.tre",tree.names=T)
plot(tree, show.tip.label = FALSE)
mydata<-read.csv('TutorialData/DiscreteData.csv', row.names=1, header=TRUE)
mydata[1:10,]

#Match data with tree
data.pruned<-treedata(phy=tree,data = mydata, warnings=FALSE)
tree<-data.pruned$phy
mydata<-data.pruned$data

#Plot data on tree
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,1],sort(unique(mydata[,1]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,1]))),piecol=c("green", "orange"),cex=.3, offset=.2)

#####1: Single trait analysis: transition rate comparison
library(corHMM)
   #Set up initial rate matrices
rmat.er<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ER")  #equal transition rates
rmat.ard<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ARD") #unequal transition rates
rmat.er   #all rates the same
rmat.ard  #all rates different 
   #Set up data
trt1<-cbind(row.names(mydata),mydata[,1])
trt2<-cbind(row.names(mydata),mydata[,3])
trt3<-cbind(row.names(mydata),mydata[,5])

#Fit models: TRAIT 1
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,1],sort(unique(mydata[,1]))),piecol=c("red", "black"),cex=.3, offset=0)
res1.er<-corHMM(tree,trt1,rate.cat=1,rate.mat=rmat.er,node.states="marginal") 
res1.er
res1.ard<-corHMM(tree,trt1,rate.cat=1,rate.mat=rmat.ard,node.states="marginal") 
   #compare models: logL and AIC
res1.er$loglik 
res1.ard$loglik  #identical in this case: no need for formal LRT 
res1.er$AICc
res1.ard$AICc    #Choose simpler model (lower AICc b/c fewer parameters and same logL)
res1.er$solution  #rate transition parameters

#Fit models: TRAIT 2
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
res2.er<-corHMM(tree,trt2,rate.cat=1,rate.mat=rmat.er,node.states="marginal") 
res2.ard<-corHMM(tree,trt2,rate.cat=1,rate.mat=rmat.ard,node.states="marginal") 
#compare models: logL and AIC
res2.er$loglik 
res2.ard$loglik  #identical in this case: no need for formal LRT 
res2.er$AICc
res2.ard$AICc    #Choose simpler model (lower AICc b/c fewer parameters and same logL)
res2.er$solution  #rate transition parameters

#Fit models: TRAIT 3
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,5],sort(unique(mydata[,5]))),piecol=c("red", "black"),cex=.3, offset=0)
res3.er<-corHMM(tree,trt3,rate.cat=1,rate.mat=rmat.er,node.states="marginal") 
res3.ard<-corHMM(tree,trt3,rate.cat=1,rate.mat=rmat.ard,node.states="marginal") 
#compare models: logL and AIC
res3.er$loglik 
res3.ard$loglik
LRT<- -2*(res3.er$loglik - res3.ard$loglik)  #LRT test
LRT
1-pchisq( LRT,df = 1)  #probability of LRT
res3.er$AICc
res3.ard$AICc    #ARD preferred
res3.ard$solution  #rate transition parameters
   #transition from 


#####2: Two trait analysis: Trait association
#Compare dependent (correlated) change model (ARD) with independent change model (ER)  #See Pagel 1994
trtset12<-cbind(row.names(mydata),mydata[,1:2])
trtset34<-cbind(row.names(mydata),mydata[,3:4])

#Traits 1&2
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,1],sort(unique(mydata[,1]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,2],sort(unique(mydata[,2]))),piecol=c("green", "orange"),cex=.3, offset=0.2)
disc.res12.er<-corDISC(tree,trtset12,ntraits=2,model="ER",node.states="marginal") 
disc.res12.ard<-corDISC(tree,trtset12,ntraits=2,model="ARD",node.states="marginal") 
disc.res12.er$loglik
disc.res12.ard$loglik
LRT<- -2*(disc.res12.er$loglik - disc.res12.ard$loglik)  #LRT test
1-pchisq( LRT,df = 4)  #probability of LRT
disc.res12.er$AIC
disc.res12.ard$AIC   #STRONG SUPPORT for correlated (dependent) evolution
   #Implemented in phytools
tr1<-mydata[,1]; names(tr1)<-row.names(mydata)
tr2<-mydata[,2]; names(tr2)<-row.names(mydata)
disc.res12.b<-fitPagel(tree,x=tr1,y=tr2)
disc.res12.b
  #NOTE: LOOK AT PHYLOGENY AND DATA HERE!  WE have strong support for correlated evolution, but both traits are co-distributed and clustered
    #essentially, 1 evolutionary shift could explain the pattern: the 'BiSSE' problem (need ancestral state estimation: next week)
 ##### ALWAYS PLOT YOUR DATA! DON'T JUST RUN STATISTICS (statistical anlaysis alone here leads to mis-interpretation)

#Traits 3&4
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,4],sort(unique(mydata[,4]))),piecol=c("green", "orange"),cex=.3, offset=0.2)
disc.res34.er<-corDISC(tree,trtset34,ntraits=2,model="ER",node.states="marginal") 
disc.res34.ard<-corDISC(tree,trtset34,ntraits=2,model="ARD",node.states="marginal") 
disc.res34.er$loglik
disc.res34.ard$loglik 
disc.res34.er$AIC
disc.res34.ard$AIC # VERY STRONG SUPPORT FOR CORRELATED EVOLUTION


#####3: Two trait analysis: Order to transitions matters (changes in trait 2 DEPEND on values of trait 1: Maddison 1990)
  #Must define the two models for comparison
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,4],sort(unique(mydata[,4]))),piecol=c("green", "orange"),cex=.3, offset=0.2)
eq.rat<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=2, nstates=2, model="ER")
eq.rat
dir.rat<-eq.rat; dir.rat[1,2]<-2; dir.rat[3,4]<-3  #e.g., q12<>q34: see Table in Pagel 1994
dir.rat
  #NOTE: which values are allowed to differ dependently depends upon the directional hypothesis under investigation!

disc.res34.er<-corDISC(tree,trtset34,ntraits=2,model="ER",node.states="marginal") 
disc.res34.dir<-corDISC(tree,trtset34,ntraits=2,model="ER",rate.mat=dir.rat,node.states="marginal") #run with directional hypothesis
disc.res34.er$loglik
disc.res34.dir$loglik 
disc.res34.er$AIC
disc.res34.dir$AIC # VERY STRONG SUPPORT FOR DIRECTED DEPENDENT EVOLUTION

