#Analyses for Tutorial 06: Ancestral State Estimation in R
  #Parts of this tutorial are based on that of L. Revell: http://www.phytools.org/eqg2015/asr.html

library(phytools)
library(geiger)

######Discrete Anc. State Estimation
## Read data & tree
tree<-read.tree("TutorialData/anole.gp.tre",tree.names=T)
group<-read.csv('TutorialData/anole.gp.csv', row.names=1, header=TRUE,colClasses=c('factor'))
gp<-as.factor(t(group)); names(gp)<-row.names(group)

#Plot 
plot(tree,type="fan")
cols<-setNames(palette()[1:length(unique(gp))],sort(unique(gp)))
tiplabels(pie=model.matrix(~gp-1),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

### A:Parsimony estimation
library(castor)
gp.int<-as.numeric(gp); names(gp.int)<-names(gp)  
anc.mp <- asr_max_parsimony(tree,gp.int)

#Plot
plot(tree, type="fan")
cols<-setNames(palette()[1:length(unique(gp.int))],sort(unique(gp.int)))
tiplabels(pie=model.matrix(~as.factor(gp.int)-1),piecol=cols,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=anc.mp$ancestral_likelihoods,piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

### B: ML estimation
anc.ML<-ace(gp,tree,model="ER",type="discrete")
  round(anc.ML$lik.anc[1:10,],3)  #show result: 
#PLOT  
plot(tree,type="fan")
cols<-setNames(palette()[1:length(unique(gp))],sort(unique(gp)))
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=anc.ML$lik.anc,piecol=cols,cex=0.3)
tiplabels(pie=model.matrix(~gp-1),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)


#MCMC Stochastic Character Mapping (SIMMAP) based on Huelsenbeck et al. 2003
# simulate single stochastic character map using empirical Bayes method
tree.smp1<-make.simmap(tree,gp,model="ER")
plot(tree.smp1,cols,type="fan")  #one run.  Not overly useful. 

#Must do many times
tree.smp<-make.simmap(tree,gp,model="ER",nsim=100)
anc.smp<-summary(tree.smp,plot=FALSE)
plot(anc.smp, type="fan", fsize=0.8,ftype="i")

rm(list = ls())

######Continuous Anc. State Estimation
#Data from Mahler et al. 2010. Evolution
tree<-read.tree("TutorialData/anole.svl.tre",tree.names=T)
svl<-read.csv('TutorialData/anole.svl.csv', row.names=1, header=TRUE)
svl<-as.matrix(svl)[,1]  #change to vector
anc.cont.ML<-fastAnc(tree,svl,vars=TRUE,CI=TRUE)
  #anc.cont.ML #output: estimate and 95% CI
anc.cont.ML$ace  #ancestral estimates
plot(tree)
nodelabels()


#PLOT as color map
tree.col<-contMap(tree,svl,plot=FALSE)  #runs Anc. St. Est. on branches of tree
plot(tree.col,type="fan")

#using 'ace' function in APE
anc.cont.ML2<-ace(x=svl,phy=tree, type="continuous", method="ML")
anc.cont.ML2$ace  #with APE
anc.cont.ML$ace   #with phytools: the same estimates

anc.cont.gls<-ace(x=svl,phy=tree, corStruct = corBrownian(1, tree), method="GLS")  #same as ML  (see Schluter et al. 1997)
anc.cont.gls$ace   #GLS: SCP. the same

rm(list = ls())
######################## Simulation with known nodes
tree<-pbtree(n=100,scale=1)
## simulate data with a trend
x<-fastBM(tree,internal=TRUE,mu=3)
phenogram(tree,x,ftype="off")

x.tip<-x[match(tree$tip.label,names(x))]
phenogram(tree,x.tip)  #traitgram under BM with no ancestral information

#estimate with no prior 
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
x<-x[tree$tip.label]
## let's see how bad we do if we ignore the trend
plot(a,fastAnc(tree,x),xlab="true values",
     ylab="estimated states under BM")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated without prior information")

## incorporate prior knowledge
pm<-setNames(c(1000,rep(0,tree$Nnode)),
             c("sig2",1:tree$Nnode+length(tree$tip.label)))
## the root & two randomly chosen nodes
nn<-as.character(c(length(tree$tip.label)+1,
                   sample(2:tree$Nnode+length(tree$tip.label),2)))
pm[nn]<-a[as.character(nn)]
## prior variance
pv<-setNames(c(1000^2,rep(1000,length(pm)-1)),names(pm))
pv[as.character(nn)]<-1e-100
## run MCMC
mcmc<-anc.Bayes(tree,x,ngen=100000,
                control=list(pr.mean=pm,pr.var=pv,
                             a=pm[as.character(length(tree$tip.label)+1)],
                             y=pm[as.character(2:tree$Nnode+length(tree$tip.label))]))

anc.est<-colMeans(mcmc$mcmc[201:1001,as.character(1:tree$Nnode+length(tree$tip.label))])

plot(a,anc.est,xlab="true values",
     ylab="estimated states using informative prior")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated using informative prior")
