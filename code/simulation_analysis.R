###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Analyse results of simulations and draw plots                         ##
##                                                                       ##
## James Liley and Chris Wallace, 2020                                   ##
## Correspondence: JL, james.liley@igmm.ed.ac.uk                         ##
##                                                                       ##
###########################################################################
#
#
# This R script will deterministically generate all plots and outputs used
#  in the paper above. Note that plots relating to the analysis of TWAS 
#  data are generated in a separate script.
#
# This script requires all requisite packages to be installed and the 
#  working directory to be set to the file one level up from where this 
#  script is saved (it should contain subdirectories ./code, ./data, 
#  ./simulations and ./outputs). All outputs are written to ./outputs.
#
# Simulations are generated using the R script ./code/run_simulation.R. 
#  Each run of the simulation is dependent on a random seed. A matrix of 
#  completed simulation outputs with associated random seeds, with FDR 
#  control level alpha equal to 0.1, is stored in ./data/simmat0.1.txt, 
#  and a matrix of simulation outputs with alpha=0.01 is stored in 
#  ./data/simmat0.01.txt
#  

###########################################################################
## Packages, directories and switches                                    ##
###########################################################################

# Packages

source("code/functions.R") #or library(cfdr)
library(mnormt)
library(mgcv)
library(pbivnorm)
library(MASS)
library(fields)
library(matrixStats)
library(latex2exp)

# Directories

# This variable gives the location of the working directory. It should 
#  contain a document README, and four folders  named ./code, ./data, 
#  ./simulations, and ./outputs. 
cfdr_dir="./"

# This variable gives the location to save outputs (plots, tables) to.
output_dir=paste0(cfdr_dir,"outputs/")



# Switches

# set to F to draw plots rather than saving them as PDFs
save_pdf=T




###########################################################################
## Read simulation data                                                  ##
###########################################################################

# Read tables. They do not have headers so we regenerate names.

# General simulations, with FDR control level 0.1. Parameters θ=(nhyp, n1p,
#  n1q, n1pq, sp, sq, distribution) are chosen from a fixed distribution. 
sim_gen_hfdr=read.table(paste0(cfdr_dir,"./data/sim_gen_high_fdr.txt"),fill=TRUE,header=F)

# Simulations as above in which the true number of associations (n1p+n1pq) 
#  is fixed at 0. Parameters θ'= θ\{n1p,n1pq} are chosen from a 
#  distribution satisfying θ'~ θ|n1p=n1pq=0
sim_gen_hfdr0=read.table(paste0(cfdr_dir,"./data/sim_gen_high_fdr_null.txt"),fill=TRUE,header=F)



# General simulations, with FDR control level 0.01
sim_gen_lfdr=read.table(paste0(cfdr_dir,"./data/sim_gen_low_fdr.txt"),fill=TRUE,header=F)

# Corresponding null simulations as above
sim_gen_lfdr0=read.table(paste0(cfdr_dir,"./data/sim_gen_low_fdr_null.txt"),fill=TRUE,header=F)



# Simulations with parameters coming from one of several fixed values, 
#  rather than randomly chosen
sim_fixed=read.table(paste0(cfdr_dir,"./data/sim_fixed.txt"),fill=TRUE,header=F)


# Simulations of only ECDF-based cFDR under various patterns of correlation 
#  between observations, with FDR control level 0.1
sim_cov=read.table(paste0(cfdr_dir,"./data/sim_cov.txt"),fill=TRUE,header=F)

# Corresponding null simulations as above
sim_cov0=read.table(paste0(cfdr_dir,"./data/sim_cov_null.txt"),fill=TRUE,header=F)


# Simulations with n1pq=0 (P and Q unrelated)
sim_ind=read.table(paste0(cfdr_dir,"./data/sim_unrelated.txt"),fill=TRUE,header=F)



# Simulations of only parametric cFDR, on the same datasets as sim_gen_hfdr, 
#  using a parametric rather than empirical estimate of Pr(Q<q|HP0) 
sim_paradj=read.table(paste0(cfdr_dir,"./data/sim_parametric_adjustment.txt"),fill=TRUE,header=F)






# Column names for table
nfold=3
vars=c("hit_p",as.vector(
  outer(as.vector(
    outer(as.vector(
      outer(as.vector(
        outer(0:1,1:4,function(x,y) paste0("_fdr",y,"_adj",x))),
        1:2,function(x,y) paste0(x,"_dist",y))),
      1:7,function(x,y) paste0("hit_cf",y,x))),
    c("",paste0("_fold",1:nfold)),function(x,y) paste0(x,y))))

nsv=c("seed","alpha","nhyp","dist1","dist2","n1p","n1q","n1pq","sp","sq",
  paste0("fit_",c("pi0","pi1","pi2","tau1","tau2","s1","s2","conv")),"pi0_null","sigma_null","cor_mode","rho")

nfdp=gsub("hit","fdp",vars)
ntdr=gsub("hit","tdr",vars)

nx=c(nsv, nfdp,ntdr)


colnames(sim_gen_hfdr)=nx; 
colnames(sim_gen_hfdr0)=nx; 

colnames(sim_gen_lfdr)=nx
colnames(sim_gen_lfdr0)=nx

colnames(sim_cov)=nx
colnames(sim_cov0)=nx

colnames(sim_fixed)=nx
colnames(sim_ind)=nx


colnames(sim_paradj)=nx




###########################################################################
## Some functions                                                        ##
###########################################################################

# This function is from https://stat.ethz.ch/pipermail/r-help/2007-June/133192.html,
#  implementing gatz95
weighted.var.se <- function(x, w, na.rm=FALSE) {
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}
vw=function(x,wmat) apply(wmat,1,function(w) weighted.var.se(x,w))

# This function plots a polygon around a pointwise standard error envelope
pse=function(x,y,s,alpha=0.95,env=TRUE, col="black",lty=1,lwd=1) {
  z=-qnorm((1-alpha)/2)
  cf=col2rgb(col)/256; cfill=rgb(cf[1],cf[2],cf[3],alpha=0.2)
  if (env) polygon(c(x,rev(x)),c(y + z*s,rev(y-z*s)),border=NA,col=cfill);
  lines(x,y,lty=lty,lwd=lwd,col=col)
}


# get mean and 95% CI
mci=function(x,alpha=0.95) {
  se=sd(x,na.rm=T)/sqrt(length(x))
  z=-qnorm((1-alpha)/2)
  c(mean(x,na.rm=T),mean(x,na.rm=T) - z*se,mean(x,na.rm=T) + z*se)
}


# Does one set of values (TDP or FDP) exceed another?
#  Formally: given parameter x, values y1,y2 (paired):
#  Within each eighth of x values, checks the following 
#  holds:
#   Wilcoxon test of ranks rejects identity of distribution at p<1e-3 (or pcut)
#   y1>y2 more often than y2>y1
#   
exceeds=function(x,y1,y2,div=8,pcut=5e-3) {
  w=which(!is.na(x+y1+y2)); x=x[w]; y1=y1[w]; y2=y2[w]
  xr=range(x)
  check=rep(0,div)
  cuts=xr[1] + (0:div)*(xr[2]-xr[1])/div
  cuts[div+1]=cuts[div+1]+1 # ensure largest x value is included
  for (i in 1:div) {
    xsub=which(x >= cuts[i] & x<cuts[i+1])
    ytest=(y1-y2)[xsub]
    mp=suppressWarnings(wilcox.test(ytest)$p.value)
    posdif=length(which(ytest> (1e-5)))-length(which(ytest< -(1e-5)))
    mdif=mean(ytest)
    if (mp<pcut & mdif>0 & posdif>0) check[i]=1
  }
  return(length(which(check>0))> 0.749*length(check))
}



###########################################################################
## Parameters                                                            ##
###########################################################################

sdx_general=0.15 # kernel width for smoothing, as proportion of x-axis range
lwd_general=2 # line width
normalise=F # set to TRUE to plot TDR relative to TDR of p-values, rather than absolute



###########################################################################
## FDR control using CDF methods                                         ##
###########################################################################
# Figure 2, supp fig 7


for (xdist in 1:2) {
  for (ia in 1:2) {
    
    
    atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
    if (ia==1) {
      alpha_cut=0.1
      attach(sim_gen_hfdr) 
    } else {
      alpha_cut=0.01
      attach(sim_gen_lfdr)
    }
    
    if (save_pdf) pdf(paste0(output_dir,"fdr_control_alpha",ia,"_dist",xdist,".pdf"),width=5,height=5)
    
    tp=fdp_p
    if (xdist==1) { # True null distribution
      t1=fdp_cf1_fdr1_adj0_dist1 # cFDR ECDF, original FDR control method
      t2=fdp_cf1_fdr2_adj0_dist1 # cFDR ECDF, naive new method
      t3=fdp_cf1_fdr3_adj0_dist1 # cFDR ECDF, leave-out-block
      t4=fdp_cf1_fdr4_adj0_dist1 # cFDR ECDF, leave-one-out
    } else { # Empirical estimate of null distribution
      t1=fdp_cf1_fdr1_adj0_dist2 
      t2=fdp_cf1_fdr2_adj0_dist2 
      t3=fdp_cf1_fdr3_adj0_dist2 
      t4=fdp_cf1_fdr4_adj0_dist2 
    }
    x=n1p+n1pq
    
    w=which(pmin(tp,t1,t2,t3,t4)>-1 & #txp-tx5 < 0.02 &
        is.finite(rowSums(cbind(tp,t1,t2,t3,t4))))
    txp=tp[w]; 
    tx1=t1[w]; tx2=t2[w]; tx3=t3[w]; tx4=t4[w];
    xx=x[w]
    
    sdx=sdx_general # gaussian smoothing
    mwd=lwd_general # line width
    
    sdx=sdx*(max(xx)-min(xx)) 
    
    x0=seq(min(xx),max(xx),length.out=100)
    wt=outer(x0,xx,function(x,y) dnorm(x-y,sd=sdx))
    wt=wt/rowSums(wt)
    
    xf= (wt %*% xx) # weighted x-values
    fp= (wt %*% txp); sep=sqrt(vw(txp,wt))
    f1= (wt %*% tx1); se1=sqrt(vw(tx1,wt))
    f2= (wt %*% tx2); se2=sqrt(vw(tx2,wt))
    f3= (wt %*% tx3); se3=sqrt(vw(tx3,wt))
    f4= (wt %*% tx4); se4=sqrt(vw(tx4,wt))
    
    
    # Null simulations
    if (ia==1) {
      nullp=sim_gen_hfdr0$fdp_p
      if (xdist==1) { # True null distribution
        null1=sim_gen_hfdr0$fdp_cf1_fdr1_adj0_dist1 
        null2=sim_gen_hfdr0$fdp_cf1_fdr2_adj0_dist1
        null3=sim_gen_hfdr0$fdp_cf1_fdr3_adj0_dist1
        null4=sim_gen_hfdr0$fdp_cf1_fdr4_adj0_dist1
      } else { # Empirical estimate of null distribution
        null1=sim_gen_hfdr0$fdp_cf1_fdr1_adj0_dist2 
        null2=sim_gen_hfdr0$fdp_cf1_fdr2_adj0_dist2 
        null3=sim_gen_hfdr0$fdp_cf1_fdr3_adj0_dist2 
        null4=sim_gen_hfdr0$fdp_cf1_fdr4_adj0_dist2 
      }
    } else { # if ia==2
      nullp=sim_gen_lfdr0$fdp_p
      if (xdist==1) { # True null distribution
        null1=sim_gen_lfdr0$fdp_cf1_fdr1_adj0_dist1 
        null2=sim_gen_lfdr0$fdp_cf1_fdr2_adj0_dist1
        null3=sim_gen_lfdr0$fdp_cf1_fdr3_adj0_dist1
        null4=sim_gen_lfdr0$fdp_cf1_fdr4_adj0_dist1
      } else { # Empirical estimate of null distribution
        null1=sim_gen_lfdr0$fdp_cf1_fdr1_adj0_dist2 
        null2=sim_gen_lfdr0$fdp_cf1_fdr2_adj0_dist2 
        null3=sim_gen_lfdr0$fdp_cf1_fdr3_adj0_dist2 
        null4=sim_gen_lfdr0$fdp_cf1_fdr4_adj0_dist2 
      }
    }
    
    # Confidence intervals at n1p+n1pq==0
    alpha_ci=0.95
    mnullp=mean(nullp); cnullp=binom.test(sum(nullp),length(nullp),alpha_cut)$conf.int[1:2]
    mnull1=mean(null1); cnull1=binom.test(sum(null1),length(null1),alpha_cut)$conf.int[1:2]
    mnull2=mean(null2); cnull2=binom.test(sum(null2),length(null2),alpha_cut)$conf.int[1:2]
    mnull3=mean(null3); cnull3=binom.test(sum(null3),length(null3),alpha_cut)$conf.int[1:2]
    mnull4=mean(null4); cnull4=binom.test(sum(null4),length(null4),alpha_cut)$conf.int[1:2]
    
    
    
    prange=c(0,max(c(fp,f1,f2,f3,f4,cnullp,cnull1,cnull2,cnull3,cnull4)));
    yx="FDR"
    plot(0,type="n",xlim=c(-50,max(xx)),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
    
    pse(xf,fp,sep,col="black",lty=1,lwd=mwd)
    pse(xf,f1,se1,col="red",lty=1,lwd=mwd) # cFDR original
    pse(xf,f2,se2,col="blue",lty=1,lwd=mwd) # cFDR naive
    pse(xf,f3,se3,col="black",lty=2,lwd=mwd) # cFDR leave-out-block
    pse(xf,f4,se4,col="black",lty=3,lwd=mwd) # cFDR leave-one-out
    
    nx=0.01*max(xx) #draw the lines for confidence intervals at pi0=0 this far apart
    lines(c(0,0)-2*nx,cnullp,col="black",lty=1,lwd=mwd); points(-2*nx,mnullp,pch=16,col="black")
    lines(c(0,0)-1*nx,cnull1,col="red",lty=1,lwd=mwd); points(-1*nx,mnull1,pch=16,col="red")
    lines(c(0,0)+0*nx,cnull2,col="blue",lty=1,lwd=mwd); points(0*nx,mnull2,pch=16,col="blue")
    lines(c(0,0)+1*nx,cnull3,col="black",lty=2,lwd=mwd); points(1*nx,mnull3,pch=16,col="black")
    lines(c(0,0)+2*nx,cnull4,col="black",lty=3,lwd=mwd); points(2*nx,mnull4,pch=16,col="black")
    
    abline(h=alpha_cut,col="gray")
    
    legend("bottomright",bg="white",
      c("P-val.",
        "Orig.",
        "Naive",
        "LOO",
        "LOB"),
      col=c("black","red","blue","black","black"),
      lty=c(1,1,1,2,3))
    
    # exceedances
    if (ia==1 & xdist==1) {
      print(paste("M4 exceeds M1: ",exceeds(xx,tx4,tx1))) # does FDR using method 4 exceed FDR using method 1?
      print(paste("M3 exceeds M1: ",exceeds(xx,tx3,tx1))) # does FDR using method 3 exceed FDR using method 1?
      print(paste("M2 exceeds M4: ",exceeds(xx,tx2,tx4))) # does FDR using method 2 exceed FDR using method 4?
      print(paste("M2 exceeds M3: ",exceeds(xx,tx2,tx3))) # does FDR using method 2 exceed FDR using method 4?
    }
    
    if (save_pdf) dev.off()
    
  }
}


###########################################################################
## Power of ECDF estimator                                               ##
###########################################################################
## Figure 3, supp fig 8

for (ia in 1:2) {
  
  atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
  if (ia==1) {
    alpha_cut=0.1
    attach(sim_gen_hfdr) 
  } else {
    alpha_cut=0.01
    attach(sim_gen_lfdr)
  }
  
  if (save_pdf) pdf(paste0(output_dir,"power1_alpha",ia,".pdf"),width=5,height=5)
  
  tp=tdr_p
  t1=tdr_cf1_fdr3_adj1_dist1 # cFDRn 
  t2=tdr_cf1_fdr3_adj0_dist1 # cFDR
  t3=tdr_cf1_fdr1_adj1_dist1 # cFDR, old method
  x=n1p+n1pq
  
  w=which(pmin(tp,t1,t2,t3) > -1)
  
  txp=tp[w]; 
  tx1=t1[w]; tx2=t2[w]; tx3=t3[w]; 
  xx=x[w]
  
  sdx=sdx_general # gaussian smoothing
  mwd=lwd_general # line width
  
  sdx=sdx*(max(xx)-min(xx))  
  
  x0=seq(min(xx),max(xx),length.out=100)
  wt=outer(x0,xx,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xf= (wt %*% xx);
  fp= (wt %*% txp); sep=sqrt(vw(txp,wt))
  f1= (wt %*% tx1); se1=sqrt(vw(tx1,wt))
  f2= (wt %*% tx2); se2=sqrt(vw(tx2,wt))
  f3= (wt %*% tx3); se3=sqrt(vw(tx3,wt))
  
  if (normalise) {
    f1=f1-fp; f2=f2-fp; f3=f3-fp; fp=fp-fp
  }
  
  prange=pmax(0,range(c(fp,f1,f2,f3)));
  if (normalise) yx=expression(paste(Delta,"(TDR)")) else yx="TDR"
  plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
  
  
  
  pse(xf,fp,sep,col="black",lty=2,lwd=mwd) # BH
  pse(xf,f1,se1,col="blue",lty=2,lwd=mwd) # cFDRn (adjusted)
  pse(xf,f2,se2,col="black",lty=1,lwd=mwd) # cFDR (unadjusted)
  pse(xf,f3,se3,col="red",lty=2,lwd=mwd) # cFDR (old method)
  
  abline(h=0)
  
  legend("bottomright",bg="white",
    c("P-val.",
      expression(widehat("cFDR")^n),
      expression(widehat("cFDR")),
      "Orig."),
    col=c("black","blue","black","red"),
    lty=c(2,2,1,2))
  
  
  
  # exceedances
  if (ia==1) {
    print(paste("TDR for cFDRn exceeds cFDR: ",exceeds(xx,tx1,tx2))) # does adjustment lead to improved TDR?
    print(paste("TDR using new FDR control method exceeds TDR using old method: ",exceeds(xx,tx2,tx3))) 
  }
  
  
  if (save_pdf) dev.off()
  
}






###########################################################################
## Power using PDF methods                                               ##
###########################################################################
## Figure 4, supp fig 10


for (ia in 1:2) { # use fdr=0.1 (ia=1) or fdr=0.01 (ia==2)
  for (id in 1:2) { # parametric assumptions satisfied, ie normal alternative (id==1) or not satisfied (id==2)
    
    atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
    if (ia==1) {
      alpha_cut=0.1
      attach(sim_gen_hfdr) 
    } else {
      alpha_cut=0.01
      attach(sim_gen_lfdr)
    }
    
    if (save_pdf) pdf(paste0(output_dir,"power3_alpha",ia,"_dist",id,".pdf"),width=5,height=5)
    
    
    
    tp=tdr_p
    t1=tdr_cf1_fdr4_adj1_dist1 # cFDR ECDF
    t4=tdr_cf4_fdr3_adj0_dist1 # cfdr KDE
    t5=tdr_cf5_fdr3_adj0_dist1 # cfdr param
    t6=tdr_cf6_fdr4_adj0_dist1 # cfdr oracle
    t7=tdr_cf7_fdr4_adj0_dist1 # cFDR oracle
    x=n1p+n1pq
    
    if (id==1) xdist=1 else xdist=2:3
    w=which(x>0 & pmin(tp,t1,t4,t5,t6,t7) > -1  & (dist1 %in% xdist) ) 
    txp=tp[w]; 
    tx1=t1[w]; tx4=t4[w]; tx5=t5[w]; tx6=t6[w]; tx7=t7[w]; 
    xx=x[w]
    
    
    sdx=sdx_general # gaussian smoothing
    mwd=lwd_general # line width
    
    sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)
    
    
    x0=seq(min(xx),max(xx),length.out=100)
    wt=outer(x0,xx,function(x,y) dnorm(x-y,sd=sdx))
    wt=wt/rowSums(wt)
    
    xf= (wt %*% xx);
    fp= (wt %*% txp); sep=sqrt(vw(txp,wt))
    f1= (wt %*% tx1); se1=sqrt(vw(tx1,wt))
    f4= (wt %*% tx4); se4=sqrt(vw(tx4,wt))
    f5= (wt %*% tx5); se5=sqrt(vw(tx5,wt))
    f6= (wt %*% tx6); se6=sqrt(vw(tx6,wt))
    f7= (wt %*% tx7); se7=sqrt(vw(tx7,wt))
    
    if (normalise) {
      f1=f1-fp; f4=f4-fp; f5=f5-fp;  f6=f6-fp;  f7=f7-fp; fp=fp-fp
    }
    
    
    
    prange=pmax(0,range(c(fp,f1,f4,f5,f6,f7)));
    if (normalise) yx=expression(paste(Delta,"(power)")) else yx="Power"
    plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
    
    
    
    pse(xf,fp,sep,col="black",lty=2,lwd=mwd) # BH
    pse(xf,f1,se1,col="black",lty=1,lwd=mwd) # cFDR ECDF 
    pse(xf,f4,se4,col="red",lty=3,lwd=mwd) # cfdr KDE
    pse(xf,f5,se5,col="blue",lty=3,lwd=mwd) # cfdr param
    pse(xf,f6,se6,col="gray",lty=3,lwd=mwd) # cfdr oracle
    pse(xf,f7,se7,col="gray",lty=1,lwd=mwd) # cFDR oracle
    
    abline(h=0)
    
    legend("bottomright",bg="white",
      c("P-val",
        expression(widehat("cFDR")^n),
        "CDF oracle",
        "PDF, param.",
        "PDF, KDE",
        "PDF oracle"),
      col=c("black","black","gray","blue","red","gray"),
      lty=c(2,1,1,3,3,3))
    
    # exceedances
    if (ia==1 & id==1) {
      print(paste("TDR using ECDF exceeds TDR using param. cfdr where param. assumptions hold: ",exceeds(xx,tx1,tx5))) # does adjustment lead to improved TDR?
      print(paste("TDR using ECDF exceeds TDR using KDE. cfdr where param. assumptions hold: ",exceeds(xx,tx1,tx4))) # does adjustment lead to improved TDR?
    }
    if (ia==1 & id==2) {
      print(paste("TDR using ECDF exceeds TDR using param. cfdr where param. assumptions do not hold: ",exceeds(xx,tx1,tx5))) # does adjustment lead to improved TDR?
      print(paste("TDR using ECDF exceeds TDR using KDE. cfdr where param. assumptions do not hold: ",exceeds(xx,tx1,tx4))) # does adjustment lead to improved TDR?
    }
    
    
    
    if (save_pdf) dev.off()
    
  }
}


###########################################################################
## Power using CDF methods                                               ##
###########################################################################
## Figure 5, supp fig 9

for (id in 1:2) {
  for (ia in 1:2) {
    
    atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
    if (ia==1) {
      alpha_cut=0.1
      attach(sim_gen_hfdr) 
    } else {
      alpha_cut=0.01
      attach(sim_gen_lfdr)
    }
    
    if (id==1) dists=1 else dists=2:3
    
    if (save_pdf) pdf(paste0(output_dir,"power2_dist",dists,"_alpha",ia,".pdf"),width=5,height=5)
    
    tp=tdr_p
    t1=tdr_cf1_fdr4_adj1_dist1 # cFDR ECDF
    t2=tdr_cf2_fdr3_adj1_dist1 # cFDR KDE
    t3=tdr_cf3_fdr3_adj1_dist1 # cFDR param
    t6=tdr_cf7_fdr4_adj0_dist1 # cFDR oracle
    x=n1p+n1pq
    
    w=which(pmin(tp,t1,t2,t3,t6) > -1 & (dist1 %in% dists)) 
    txp=tp[w]; 
    tx1=t1[w]; tx2=t2[w]; tx3=t3[w]; tx6=t6[w];
    xx=x[w]
    
    
    sdx=sdx_general # gaussian smoothing
    mwd=lwd_general # line width
    
    
    sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)
    
    x0=seq(min(xx),max(xx),length.out=100)
    wt=outer(x0,xx,function(x,y) dnorm(x-y,sd=sdx))
    wt=wt/rowSums(wt)
    
    xf= (wt %*% xx);
    fp= (wt %*% txp); sep=sqrt(vw(txp,wt))
    f1= (wt %*% tx1); se1=sqrt(vw(tx1,wt))
    f2= (wt %*% tx2); se2=sqrt(vw(tx2,wt))
    f3= (wt %*% tx3); se3=sqrt(vw(tx3,wt))
    f6= (wt %*% tx6); se6=sqrt(vw(tx6,wt))
    
    if (normalise) {
      f1=f1-fp; f2=f2-fp; f3=f3-fp; f6=f6-fp;  fp=fp-fp
    }
    
    prange=pmax(0,range(c(fp,f1,f2,f3,f6)));
    if (normalise) yx=expression(paste(Delta,"(TDR)")) else yx="TDR"
    plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
    
    
    pse(xf,fp,sep,col="black",lty=2,lwd=mwd,env=F) # BH
    pse(xf,f1,se1,col="black",lty=1,lwd=mwd,env=F) # cFDR ECDF 
    pse(xf,f2,se2,col="red",lty=1,lwd=mwd,env=F) # cFDR KDE
    pse(xf,f3,se3,col="blue",lty=1,lwd=mwd,env=F) # cFDR param
    pse(xf,f6,se6,col="gray",lty=1,lwd=mwd,env=F) # cfdr oracle
    
    abline(h=0)
    
    legend("bottomright",bg="white",
      c("P-val",
        expression(widehat("cFDR")^n),
        "CDF param.",
        "CDF KDE",
        "CDF Oracle"),
      col=c("black","black","blue","red","gray"),
      lty=c(2,1,1,1,1))
    
    # exceedances
    if (ia==1 & id==1) {
      print(paste("TDR using ECDF exceeds TDR using param. cFDR where param. assumptions hold: ",exceeds(xx,tx1,tx3))) # does adjustment lead to improved TDR?
      print(paste("TDR using ECDF exceeds TDR using KDE. cFDR where param. assumptions hold: ",exceeds(xx,tx1,tx2))) # does adjustment lead to improved TDR?
    }
    if (ia==1 & id==2) {
      print(paste("TDR using ECDF exceeds TDR using param. cFDR where param. assumptions do not hold: ",exceeds(xx,tx1,tx3))) # does adjustment lead to improved TDR?
      print(paste("TDR using ECDF exceeds TDR using KDE. cFDR where param. assumptions do not hold: ",exceeds(xx,tx1,tx2))) # does adjustment lead to improved TDR?
    }
    
    
    if (save_pdf) dev.off()
    
  }
}




###########################################################################
## Covariance between observations - FDR                                 ##
###########################################################################
## Figure 6, supp fig 11

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
attach(sim_cov) 

for (cor_level in 1:2) {
  xrho=c(min(rho),max(rho))[cor_level]
  
  suffix=c("lowcov","highcov")[cor_level]
  if (save_pdf) pdf(paste0(output_dir,"fdr_cov_",cor_level,".pdf"),width=5,height=5)
  
  wn=1:dim(sim_gen_hfdr)[1] # no correlation; use sim_gen_hfdr and restrict to continuous simulations
  wb=which(cor_mode==1 & rho==xrho) # block structure
  we=which(cor_mode==2 & rho==xrho) # equicorrelation
  
  fp=fdp_p
  f3=fdp_cf1_fdr3_adj1_dist1 # leave-out-block
  f4=fdp_cf1_fdr4_adj1_dist1 # leave-one-out
  xx=n1p+n1pq
  
  
  x0=seq(min(xx),max(xx),length.out=100)
  
  sdx=sdx_general # gaussian smoothing
  mwd=lwd_general # line width
  
  
  sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)
  
  
  ### no correlation (use sim_gen_hfdr)
  w=wn
  xw=(sim_gen_hfdr$n1p+sim_gen_hfdr$n1pq)[w]; 
  fpw=sim_gen_hfdr$fdp_p[w]; 
  f3w=sim_gen_hfdr$fdp_cf1_fdr3_adj1_dist1[w];
  f4w=sim_gen_hfdr$fdp_cf1_fdr4_adj1_dist1[w];
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfn=(wt %*% xw)
  ypn=(wt %*% fpw); spn=sqrt(vw(fpw,wt))
  y3n=(wt %*% f3w); s3n=sqrt(vw(f3w,wt))
  y4n=(wt %*% f4w); s4n=sqrt(vw(f4w,wt))
  
  # Values at x=0
  w0=1:dim(sim_gen_hfdr0)[1]
  ypn0=sim_gen_hfdr0$fdp_p[w0]
  y3n0=sim_gen_hfdr0$fdp_cf1_fdr3_adj1_dist1[w0];
  y4n0=sim_gen_hfdr0$fdp_cf1_fdr4_adj1_dist1[w0];
  
  
  
  
  ### block structure
  w=wb
  xw=xx[w]; fpw=fp[w]; f3w=f3[w]; f4w=f4[w]
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfb=(wt %*% xw)
  ypb=(wt %*% fpw); spb=sqrt(vw(fpw,wt))
  y3b=(wt %*% f3w); s3b=sqrt(vw(f3w,wt))
  y4b=(wt %*% f4w); s4b=sqrt(vw(f4w,wt))
  
  # Values at x=0
  w0=which(sim_cov0$cor_mode==1 & sim_cov0$rho==xrho) # block structure
  ypb0=sim_cov0$fdp_p[w0]
  y3b0=sim_cov0$fdp_cf1_fdr3_adj1_dist1[w0];
  y4b0=sim_cov0$fdp_cf1_fdr4_adj1_dist1[w0];
  
  
  
  ### equicorrelation
  w=we
  xw=xx[w]; fpw=fp[w]; f3w=f3[w]; f4w=f4[w]
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfe=(wt %*% xw)
  ype=(wt %*% fpw); spe=sqrt(vw(fpw,wt))
  y3e=(wt %*% f3w); s3e=sqrt(vw(f3w,wt))
  y4e=(wt %*% f4w); s4e=sqrt(vw(f4w,wt))
  
  # Values at x=0
  w0=which(sim_cov0$cor_mode==2 & sim_cov0$rho==xrho) # block structure
  ype0=sim_cov0$fdp_p[w0]
  y3e0=sim_cov0$fdp_cf1_fdr3_adj1_dist1[w0];
  y4e0=sim_cov0$fdp_cf1_fdr4_adj1_dist1[w0];
  
  
  
  ## short function to plot mean+/-2SE
  line0=function(offset, y,conf=0.95,...) {
    m=mean(y,na.rm=T); s=sd(y,na.rm=T)/sqrt(length(which(!is.na(y)))); z=-qnorm((1-conf)/2)
    lines(rep(offset,2),m + c(-1,1)*z*s,...)
  }
  
  
  prange=c(0.05,0.15);
  x0=0.01*max(xx)
  
  yx="FDR"
  plot(0,type="n",xlim=range(c(-5*x0,xx)),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
  
  pse(xfn,ypn,spn,col="black",lty=1,lwd=mwd,env=F) # BH, no correlation
  pse(xfb,ypb,spb,col="black",lty=2,lwd=mwd,env=F) # BH, block correlation
  pse(xfe,ype,spe,col="black",lty=3,lwd=mwd,env=F) # BH, equicorrelation
  pse(xfn,y3n,s3n,col="red",lty=1,lwd=mwd,env=F) # cFDR method 3, no correlation
  pse(xfb,y3b,s3b,col="red",lty=2,lwd=mwd,env=F) # cFDR method 3, block correlation
  pse(xfe,y3e,s3e,col="red",lty=3,lwd=mwd,env=F) # cFDR method 3, equicorrelation
  pse(xfn,y4n,s4n,col="blue",lty=1,lwd=mwd,env=F) # cFDR method 4, no correlation
  pse(xfb,y4b,s4b,col="blue",lty=2,lwd=mwd,env=F) # cFDR method 4, block correlation
  pse(xfe,y4e,s4e,col="blue",lty=3,lwd=mwd,env=F) # cFDR method 4, equi correlation
  
  line0(-4*x0,ypn0,col="black",lty=1,lwd=mwd) # BH, no correlation
  line0(-3*x0,ypb0,col="black",lty=2,lwd=mwd) # BH, block correlation
  line0(-2*x0,ype0,col="black",lty=3,lwd=mwd) # BH, equicorrelation
  line0(-1*x0,y3n0,col="red",lty=1,lwd=mwd) # cFDR method 3, no correlation
  line0( 0*x0,y3b0,col="red",lty=2,lwd=mwd) # cFDR method 3, block correlation
  line0( 1*x0,y3e0,col="red",lty=3,lwd=mwd) # cFDR method 3, equicorrelation
  line0( 2*x0,y4n0,col="blue",lty=1,lwd=mwd) # cFDR method 4, no correlation
  line0( 3*x0,y4b0,col="blue",lty=2,lwd=mwd) # cFDR method 4, block correlation
  line0( 4*x0,y4e0,col="blue",lty=3,lwd=mwd) # cFDR method 4, equi correlation
  
  abline(h=0.1)
  
  legend("topright",bg="white",
    c("P-val",
      "LOB",
      "LOO",
      "No cor.",
      "Block cor.",
      "Equicor."),
    col=c("black","red","blue","black","black","black"),
    lty=c(NA,NA,NA,1,2,3),pch=c(16,16,16,NA,NA,NA))
  
  
  if (save_pdf) dev.off()
  
}





###########################################################################
## Covariance between observations - TDR                                 ##
###########################################################################
## Figure 7, supp fig 12

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
attach(sim_cov) 

for (cor_level in 1:2) {
  xrho=c(min(rho),max(rho))[cor_level]
  
  suffix=c("lowcov","highcov")[cor_level]
  if (save_pdf) pdf(paste0(output_dir,"tdr_cov_",cor_level,".pdf"),width=5,height=5)
  
  wn=1:dim(sim_gen_hfdr)[1] # no correlation; use sim_gen_hfdr and restrict to continuous simulations
  wb=which(cor_mode==1 & rho==xrho) # block structure
  we=which(cor_mode==2 & rho==xrho) # equicorrelation
  
  fp=tdr_p
  f3=tdr_cf1_fdr3_adj1_dist1 # leave-out-block
  f4=tdr_cf1_fdr4_adj1_dist1 # leave-one-out
  xx=n1p+n1pq
  
  
  x0=seq(min(xx),max(xx),length.out=100)
  sdx=sdx_general # gaussian smoothing
  mwd=lwd_general # line width
  
  
  sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)
  
  
  ### no correlation (use sim_gen_hfdr)
  w=wn
  xw=(sim_gen_hfdr$n1p+sim_gen_hfdr$n1pq)[w]; 
  fpw=sim_gen_hfdr$tdr_p[w]; 
  f3w=sim_gen_hfdr$tdr_cf1_fdr3_adj1_dist1[w];
  f4w=sim_gen_hfdr$tdr_cf1_fdr4_adj1_dist1[w];
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfn=(wt %*% xw)
  ypn=(wt %*% fpw); spn=sqrt(vw(fpw,wt))
  y3n=(wt %*% f3w); s3n=sqrt(vw(f3w,wt))
  y4n=(wt %*% f4w); s4n=sqrt(vw(f4w,wt))
  
  print(paste("TDR using block-leave-out exceeds TDR using",
    "leave_one_out under independence: ",exceeds(xw,f3w,f4w))) 
  
  
  
  ### block structure
  w=wb
  xw=xx[w]; fpw=fp[w]; f3w=f3[w]; f4w=f4[w]
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfb=(wt %*% xw)
  ypb=(wt %*% fpw); spb=sqrt(vw(fpw,wt))
  y3b=(wt %*% f3w); s3b=sqrt(vw(f3w,wt))
  y4b=(wt %*% f4w); s4b=sqrt(vw(f4w,wt))
  
  print(paste("TDR using block-leave-out exceeds TDR using",
    "leave_one_out under block correlation ",
    "when rho=",xrho,": ",exceeds(xw,f3w,f4w))) 
  
  print(paste("Paired Wilcoxon rank-sum test for TDR using block-leave-out against TDR using",
    "leave_one_out under block correlation ",
    "when rho=",xrho,": p-value: ",wilcox.test(f3w,f4w,paired=T)$p.value,", ",
    "TDR using block-leave out larger in ", 100*length(which(f3w>f4w))/length(f3w),"% of cases, ",
    "TDR using block-leave out larger in ", 100*length(which(f4w>f3w))/length(f3w),"% of cases")) 
  
  
  
  
  ### equicorrelation
  w=we
  xw=xx[w]; fpw=fp[w]; f3w=f3[w]; f4w=f4[w]
  w2=which(is.finite(xw+fpw+f3w+f4w))
  xw=xw[w2]; fpw=fpw[w2]; f3w=f3w[w2]; f4w=f4w[w2]; 
  
  wt=outer(x0,xw,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xfe=(wt %*% xw)
  ype=(wt %*% fpw); spe=sqrt(vw(fpw,wt))
  y3e=(wt %*% f3w); s3e=sqrt(vw(f3w,wt))
  y4e=(wt %*% f4w); s4e=sqrt(vw(f4w,wt))
  
  print(paste("TDR using block-leave-out exceeds TDR using",
    "leave_one_out under equicorrelation ",
    "when rho=",xrho,": ",exceeds(xw,f3w,f4w))) 
  
  
  
  prange=range(c(ypn,ypb,ype,y3n,y3b,y3e,y4n,y4b,y4e));
  
  yx="TDR"
  plot(0,type="n",xlim=range(c(0,xx)),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")
  
  pse(xfn,ypn,spn,col="black",lty=1,lwd=mwd,env=F) # BH, no correlation
  pse(xfb,ypb,spb,col="black",lty=2,lwd=mwd,env=F) # BH, block correlation
  pse(xfe,ype,spe,col="black",lty=3,lwd=mwd,env=F) # BH, equicorrelation
  pse(xfn,y3n,s3n,col="red",lty=1,lwd=mwd,env=F) # cFDR method 3, no correlation
  pse(xfb,y3b,s3b,col="red",lty=2,lwd=mwd,env=F) # cFDR method 3, block correlation
  pse(xfe,y3e,s3e,col="red",lty=3,lwd=mwd,env=F) # cFDR method 3, equicorrelation
  pse(xfn,y4n,s4n,col="blue",lty=1,lwd=mwd,env=F) # cFDR method 4, no correlation
  pse(xfb,y4b,s4b,col="blue",lty=2,lwd=mwd,env=F) # cFDR method 4, block correlation
  pse(xfe,y4e,s4e,col="blue",lty=3,lwd=mwd,env=F) # cFDR method 4, equi correlation
  
  abline(h=0.1)
  
  legend("topleft",bg="white",
    c("P-val",
      "LOB",
      "LOO",
      "No cor.",
      "Block cor.",
      "Equicor."),
    col=c("black","red","blue","black","black","black"),
    lty=c(NA,NA,NA,1,2,3),pch=c(16,16,16,NA,NA,NA))
  
  if (save_pdf) dev.off()
  
}




###########################################################################
## Better to adjust empirically than by parametrisation                  ##
###########################################################################
## Supplementary material, figure 2

w1a=which(sim_gen_hfdr$dist1==1) # Non-parametric adjustment, parametric assumptions hold
w1b=which(sim_gen_hfdr$dist1!=1) # Non-parametric adjustment, parametric assumptions do not hold
w2a=which(sim_paradj$dist1==1) # Parametric adjustment, parametric assumptions hold
w2b=which(sim_paradj$dist1!=1) # Parametric adjustment, parametric assumptions do not hold

t1a=sim_gen_hfdr$tdr_cf3_fdr3_adj1_dist1[w1a]; x1a=(sim_gen_hfdr$n1p+sim_gen_hfdr$n1pq)[w1a]
t1b=sim_gen_hfdr$tdr_cf3_fdr3_adj1_dist1[w1b]; x1b=(sim_gen_hfdr$n1p+sim_gen_hfdr$n1pq)[w1b]
t2a=sim_paradj$tdr_cf3_fdr3_adj1_dist1[w2a]; x2a=(sim_paradj$n1p+sim_paradj$n1pq)[w2a]
t2b=sim_paradj$tdr_cf3_fdr3_adj1_dist1[w2b]; x2b=(sim_paradj$n1p+sim_paradj$n1pq)[w2b]

c1a=which(is.finite(t1a+x1a)); t1a=t1a[c1a]; x1a=x1a[c1a]
c1b=which(is.finite(t1b+x1b)); t1b=t1b[c1b]; x1b=x1b[c1b]
c2a=which(is.finite(t2a+x2a)); t2a=t2a[c2a]; x2a=x2a[c2a]
c2b=which(is.finite(t2b+x2b)); t2b=t2b[c2b]; x2b=x2b[c2b]

xx=c(x1a,x1b,x2a,x2b)

sdx=sdx_general # gaussian smoothing
mwd=lwd_general # line width

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

x0=seq(min(xx),max(xx),length.out=100)
wt1a=outer(x0,x1a,function(x,y) dnorm(x-y,sd=sdx)); wt1a=wt1a/rowSums(wt1a)
wt1b=outer(x0,x1b,function(x,y) dnorm(x-y,sd=sdx)); wt1b=wt1b/rowSums(wt1b)
wt2a=outer(x0,x2a,function(x,y) dnorm(x-y,sd=sdx)); wt2a=wt2a/rowSums(wt2a)
wt2b=outer(x0,x2b,function(x,y) dnorm(x-y,sd=sdx)); wt2b=wt2b/rowSums(wt2b)

f1a= (wt1a %*% t1a); s1a=sqrt(vw(t1a,wt1a)); xf1a=(wt1a %*% x1a)
f1b= (wt1b %*% t1b); s1b=sqrt(vw(t1b,wt1b)); xf1b=(wt1b %*% x1b)
f2a= (wt2a %*% t2a); s2a=sqrt(vw(t2a,wt2a)); xf2a=(wt2a %*% x2a)
f2b= (wt2b %*% t2b); s2b=sqrt(vw(t2b,wt2b)); xf2b=(wt2b %*% x2b)

prange=c(0,0.4)


if (save_pdf) pdf(paste0(output_dir,"adjustment_par_vs_emp.pdf"),width=5,height=5)

plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),
  ylab="TDR",bty="n",xaxs="i",yaxs="i")


pse(xf1a,f1a,s1a,col="black",lty=1,lwd=mwd) # Non-parametric adjustment, parametric assumptions hold
pse(xf1b,f1b,s1b,col="red",lty=1,lwd=mwd) # Non-parametric adjustment, parametric assumptions do not hold
pse(xf2a,f2a,s2a,col="black",lty=2,lwd=mwd) # Parametric adjustment, parametric assumptions hold
pse(xf2b,f2b,s2b,col="red",lty=2,lwd=mwd) # Parametric adjustment, parametric assumptions do not hold

legend("bottomright",bg="white",
  c("Emp. adj., par. hold",
    "Emp. adj., par. not hold",
    "Par. adj., par. hold",
    "Par. adj., par. not hold"),
  col=c("black","red","black","red"),
  lty=c(1,1,2,2))


if (save_pdf) dev.off()







###########################################################################
## How many shared non-null hypotheses do we need for a cFDR advantage?  ##
###########################################################################
## Figure 7, supp fig 12


for (ia in 1:2) {
  
  atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
  
  if (ia==1) attach(sim_gen_hfdr) else attach(sim_gen_lfdr)
  
  tp=tdr_p
  t1=tdr_cf1_fdr4_adj1_dist1 # cFDRn 
  t2=tdr_cf7_fdr4_adj0_dist1 # cFDR oracle
  x=n1pq/(n1p+n1pq)
  
  w=which(pmin(tp,t1,t2) > -1)
  txp=tp[w]; 
  tx1=t1[w];
  tx2=t2[w];
  xx=x[w]
  
  tx1=tx1-txp
  tx2=tx2-txp
  txp=0*txp
  
  sdx=sdx_general # gaussian smoothing
  mwd=lwd_general # line width
  
  sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)
  
  x0=seq(min(xx),max(xx),length.out=100)
  wt=outer(x0,xx,function(x,y) dnorm(x-y,sd=sdx))
  wt=wt/rowSums(wt)
  
  xf= (wt %*% xx);
  fp= (wt %*% txp); sp=sqrt(vw(txp,wt))
  f1= (wt %*% tx1); s1=sqrt(vw(tx1,wt))
  f2= (wt %*% tx2); s2=sqrt(vw(tx2,wt))
  
  
  # Values at 0
  n0=dim(sim_ind)[1]
  fpn=sim_ind$tdr_p
  f1n=sim_ind$tdr_cf1_fdr3_adj1_dist1 
  f2n=sim_ind$tdr_cf7_fdr4_adj0_dist1
  s1n=sd(f1n-fpn)/sqrt(n0)
  z0=-qnorm( (1-0.95)/2)
  c1n=mean(f1n-fpn) + c(-1,1)*z0*s1n
  
  if (save_pdf) pdf(paste0(output_dir,"min_percentage_alpha",c("_high","_low")[ia],".pdf"),width=5,height=5)
  
  prange=range(c(fp,f1,f2,c1n));
  yx=expression(paste(Delta,"(TDR)"))
  plot(0,type="n",xlim=c(-0.05,max(xx)),ylim=prange,
    xlab=expression(paste("n"[1]^{pq},"/(",n[1]^p, "+ n"[1]^{pq},")")),
    ylab=yx,bty="n")
  
  pse(xf,fp,0,col="black",lty=2,lwd=mwd) # BH
  pse(xf,f1,s1,col="black",lty=1,lwd=mwd) # cFDR ECDF
  pse(xf,f2,s2,col="red",lty=1,lwd=mwd) # cFDR oracle
  
  points(0,0,col="red",pch=16) # oracle performance at n1pq=0; identical to p-value
  lines(c(0,0),c1n,col="black",lty=2,lwd=mwd)
  points(0,mean(f1n-fpn),col="black",pch=16)
  
  legend("topleft",bg="white",
    c("P-val",
      expression(widehat("cFDR")^n),
      "CDF oracle"),
    col=c("black","black","red"),
    lty=c(2,1,1))
  
  if (save_pdf) dev.off()
  
}


###########################################################################
## Table of fixed-value simualation results                              ##
###########################################################################
## Table 2, supplementary table 1

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
attach(sim_fixed) 

pars_id=paste(sim_fixed$nhyp,sim_fixed$n1p,sim_fixed$n1q,sim_fixed$n1pq, 
   sim_fixed$sp,sim_fixed$sq,sim_fixed$dist1,cor_mode,rho,sep="_")
par_opt=unique(pars_id)


tab_out=c()
for (i in 1:length(par_opt)) {
  sub=which(pars_id==par_opt[i])
  pvec=as.numeric(unlist(strsplit(par_opt[i],"_")))
  rvec=c(mci(fdp_p[sub]),mci(tdr_p[sub]), # p-value
    mci(fdp_cf1_fdr3_adj1_dist1[sub]), mci(tdr_cf1_fdr3_adj1_dist1[sub]), # cFDR (ECDF, adjusted)
    mci(fdp_cf6_fdr4_adj0_dist1[sub]), mci(tdr_cf6_fdr4_adj0_dist1[sub])) # cfdr oracle (best possible)
  tab_out=rbind(tab_out,c(pvec,rvec))
}

description=rbind(
  c("5000_100_100_100_2_2_2_0_0", "Reference"),
  c("5000_0_0_0_2_2_2_0_0", "No effects"),
  c("5000_100_100_100_1.5_1.5_1_0_0", "Weak effects"),
  c("5000_100_100_100_3_3_3_0_0", "Large variance in effect sizes"),
  c("10000_100_100_100_2_2_2_0_0", "Larger n"),
  c("1000_100_100_100_2_2_2_0_0", "Smaller n"),
  c("5000_150_150_0_2_2_2_0_0", "No non-null shared hypotheses"),
  c("5000_0_0_200_2_2_2_0_0", "All non-null hypotheses shared"),
  c("5000_2000_2000_0_2_2_2_0_0", "Negative information"),
  c("5000_100_100_100_2_2_2_1_0.05", "Block correlation"),
  c("5000_100_100_100_2_2_2_2_0.05", "Equicorrelation"))

# Row and column names for raw table
rownames(tab_out)=description[,2][match(par_opt,description[,1])]
tab_out=tab_out[description[,2],]
colnames(tab_out)=c("nhyp","n1p","n1q","n1pq","sp","sq","dist","cor_mode","rho",
  "fdr_p","fdr_p_ci1","fdr_p_ci2","tdr_p","tdr_p_ci1","tdr_p_ci2",
  "fdr_cfdr","fdr_cfdr_ci1","fdr_cfdr_ci2","tdr_cfdr","tdr_cfdr_ci1","tdr_cfdr_ci2",
  "fdr_orac","fdr_orac_ci1","fdr_orac_ci2","tdr_orac","tdr_orac_ci1","tdr_orac_ci2")

# Write raw table
write.table(tab_out,file=paste0(output_dir,"fixed_tab.txt"),row.names=F,col.names=T,quote=F)


# Dress table for main text
tab_maintext=tab_out[,c("tdr_p","fdr_cfdr","tdr_cfdr","tdr_orac")]
colnames(tab_maintext)=c("TDR(P)", "FDR(cFDR)", "TDR(cFDR)", "TDR(oracle) \\\\")
tab_maintext=signif(tab_maintext,digits=3)
tab_maintext=cbind(Description=rownames(tab_maintext),tab_maintext)
tab_maintext[,dim(tab_maintext)[2]]=paste0(tab_maintext[,dim(tab_maintext)[2]]," \\\\")

# Write table for main text in a way that can be imported straight into TeX, and save
write.table(tab_maintext, row.names=F,quote=F,sep=" & ",col.names=T)
write.table(tab_maintext, file=paste0(output_dir,"fixed_tab_maintext.txt"),
  row.names=F,quote=F,sep=" & ",col.names=T)


# Dress tables for supplement
tab_supp1=tab_out[,c("nhyp", "n1p", "n1q", "n1pq", "sp", "sq", "dist", "rho")]
colnames(tab_supp1)=c("$n$","$n_1^p$", "$n_1^q$", "$n_1^{pq}$", "$s_p$", "$s_q$", "$d$","$\\rho$ \\\\")
tab_supp1=cbind(Description=rownames(tab_supp1),tab_supp1)
tab_supp1[,dim(tab_supp1)[2]]=paste0(tab_supp1[,dim(tab_supp1)[2]]," \\\\")

# Write table for supplement in a way that can be imported straight into TeX
write.table(tab_supp1,row.names=F,quote=F,sep=" & ",col.names=T)
write.table(tab_supp1, file=paste0(output_dir,"fixed_tab_supp1.txt"),
  row.names=F,quote=F,sep=" & ",col.names=T)


tab_supp2=tab_out[,c("fdr_p", "fdr_p_ci1", "fdr_p_ci2", "tdr_p", "tdr_p_ci1", "tdr_p_ci2")]
tab_supp2=signif(tab_supp2,digits=3)
f1=paste0(tab_supp2[,1]," (",tab_supp2[,2],",",tab_supp2[,3],")")
f2=paste0(tab_supp2[,4]," (",tab_supp2[,5],",",tab_supp2[,6],") \\\\")
tab_supp2=cbind(rownames(tab_supp2),f1,f2)
colnames(tab_supp2)=c("Description","FDR(P)","TDR(P) \\\\")
write.table(tab_supp2,quote=F,row.names=F,col.names=T,sep=" & ")
write.table(tab_supp2, file=paste0(output_dir,"fixed_tab_supp2.txt"),
  quote=F,row.names=F,col.names=T,sep=" & ")


tab_supp3=tab_out[,c("fdr_cfdr", "fdr_cfdr_ci1", "fdr_cfdr_ci2", "tdr_cfdr", "tdr_cfdr_ci1", "tdr_cfdr_ci2")]
tab_supp3=signif(tab_supp3,digits=3)
f1=paste0(tab_supp3[,1]," (",tab_supp3[,2],",",tab_supp3[,3],")")
f2=paste0(tab_supp3[,4]," (",tab_supp3[,5],",",tab_supp3[,6],") \\\\")
tab_supp3=cbind(rownames(tab_supp3),f1,f2)
colnames(tab_supp3)=c("Description","FDR(cFDR)","TDR(cFDR) \\\\")
write.table(tab_supp3,quote=F,row.names=F,col.names=T,sep=" & ")
write.table(tab_supp3, file=paste0(output_dir,"fixed_tab_supp3.txt"),
  quote=F,row.names=F,col.names=T,sep=" & ")


tab_supp4=tab_out[,c("fdr_orac", "fdr_orac_ci1", "fdr_orac_ci2", "tdr_orac", "tdr_orac_ci1", "tdr_orac_ci2")]
tab_supp4=signif(tab_supp4,digits=3)
f1=paste0(tab_supp4[,1]," (",tab_supp4[,2],",",tab_supp4[,3],")")
f2=paste0(tab_supp4[,4]," (",tab_supp4[,5],",",tab_supp4[,6],") \\\\")
tab_supp4=cbind(rownames(tab_supp4),f1,f2)
colnames(tab_supp4)=c("Description","FDR(oracle)","TDR(oracle) \\\\")
write.table(tab_supp4,quote=F,row.names=F,col.names=T,sep=" & ")
write.table(tab_supp4, file=paste0(output_dir,"fixed_tab_supp4.txt"),
  quote=F,row.names=F,col.names=T,sep=" & ")




### end of 'discovery' simulations



#############################################################################
## Generate and draw simulations for iterated cFDR                         ##
#############################################################################
## Fig 8

iterated_data="./data/iterated_cfdr_data.RData"

if (!file.exists(iterated_data)) {
  
  set.seed(1)
  
  N=1000 # Number of variables
  pi0=0.9 # Proportion of null hypotheses
  iter=500 # Total number of iterations to perform
  
  f1s=3 # Z-score distribution under HP1/HQ1 is normal with this standard deviation
  
  vmat=matrix(0,N,1+iter) # values of p/v with each iteration
  hmat=vmat # Hypothesis indicators
  pmat=vmat # P-values/Q-values
  
  n0=round(N*pi0); n1=N-n0
  P=c(2*pnorm(-abs(rnorm(n1,sd=f1s))),runif(n0))
  HP=c(rep(1,n1),rep(0,n0))
  pmat[,1]=P
  vmat[,1]=P
  hmat[,1]=HP
  
  for (i in 1:iter) {
    
    # If an even iteration, random distribution of HQ1; if odd, 15 times more likely to occur when HP=1 than when HP=0
    if ((i %% 2)==0) phen_overlap=1 else phen_overlap=15
    
    h1i=sample(N,n1,prob=c(rep(phen_overlap,n1),rep(1,n0)))
    Q=runif(N); Q1=2*pnorm(-abs(rnorm(n1,sd=f1s)))
    Q[h1i]=Q1
    HQ=rep(0,N); HQ[h1i]=1
    vlx=vl(vmat[,i],Q,mode=1,indices=1:N,nv=500,nt=1000,adj=F)
    V=il(vlx,sigma_null=f1s,pi0_null=1- (length(which(h1i>n1))/n0))
    ww=which(V>0.99); V[ww]=vmat[ww,i]
    pmat[,i+1]=Q
    vmat[,i+1]=V
    hmat[,i+1]=HQ
  }
  
  save(N,pi0,iter,f1s,pmat,vmat,hmat,n0,n1,file=iterated_data)
  
} else load(iterated_data)


## Analyse average amount of overlap for odd and even iterations
ovx=rep(0,iter); 
for (i in 1:iter) ovx[i]=sum(hmat[,1]*hmat[,i])

# Average overlap in associations for odd iterations
mean(ovx[2*(1:(iter/2)) - 1])
## around 10% of associations shared: shared at random

# Average overlap in associations for even iterations
mean(ovx[2*(1:(iter/2))])
## about half of association shared



############### First plot - used in former version of manuscript ##################
# Draw plot

if (save_pdf) pdf("./outputs/mix_condition.pdf",width=10,height=5)
par(mfrow=c(1,2))
i1=round(iter/3); i2=2*i1
h1p=1:n1; h0p=(n1+1):N

cc=colorRampPalette(c("red","purple","blue"))(iter)

plot(sort(vmat[h1p,1]),type="l",
  xlab="rank",ylab=expression(paste("p or v")),
  main=expression(paste("H"^1)[P]),ylim=c(0,1))
for (i in 2:iter) lines(sort(vmat[h1p,i]),col=cc[i])
lines(sort(vmat[h1p,1]),lwd=3);
legend("topleft",c("p","","","v","",""),
  lty=c(1,NA,1,1,1,NA),lwd=c(3,NA,1,1,1,NA),col=c("black",NA,"red","purple","blue",NA),
  y.intersp=0.5,bg="white")

plot(sort(vmat[h0p,1]),type="l",
  xlab="rank",ylab=expression(paste("p or v")),
  main=expression(paste("H"^0)[P]),ylim=c(0,1))
for (i in 2:iter) lines(sort(vmat[h0p,i]),col=cc[i])
lines(sort(vmat[h0p,1]),lwd=3);
legend("topleft",c("p","","","v","",""),
  lty=c(1,NA,1,1,1,NA),lwd=c(3,NA,1,1,1,NA),col=c("black",NA,"red","purple","blue",NA),
  y.intersp=0.5,bg="white")

if (save_pdf) dev.off()



######################################################
## Second plot

pdf("outputs/fig_iterated.pdf",width=8,height=6)

## NOT NICE CODE

plot(0,type="n",xlim=c(-100,2750),ylim=c(-1,12),ann=F,xaxt="n",yaxt="n",bty="n")

imax=500
ht=0.3; lx=1;
w=which(hmat[,1]>0.5)
segments(5*w,rep(10,length(w)),5*w,rep(10,length(w))+ht,lwd=lx,col="blue")
lines(c(0,500),c(10,10),lwd=lx); lines(c(550,1050),c(10,10),lwd=lx)
text(-100,10+ht/2,"p",col="blue")

for (i in 1:5) {
  w=which(hmat[,i+1]>0.5); w1=w[which(w<101)]; w2=w[which(w>101)]
  segments(5*w1,rep(10-i,length(w)),5*w1,rep(10-i,length(w))+ht,lwd=lx,col="red")
  segments(550+ (w2-100)*5/9,rep(10-i,length(w)),550+ (w2-100)*5/9,rep(10-i,length(w))+ht,lwd=lx,col="red")
  
  lines(c(0,500),c(10,10)-i,lwd=lx); lines(c(550,1050),c(10,10)-i,lwd=lx)
  text(-100,10-i+ht/2,TeX(paste0("q_{",i,"}")))
}
text(250,11,TeX("H_1^P (100)"),cex=); text(800,11,TeX("H_0^P (900)"))
text(500,4,"..."); text(500,2,"...")
text(500,3,"..."); text(500,1,"...")

i=imax
w=which(hmat[,i]>0.5); w1=w[which(w<101)]; w2=w[which(w>101)]
segments(5*w1,rep(0,length(w)),5*w1,rep(0,length(w))+ht,lwd=lx,col="red")
segments(550+ (w2-100)*5/9,rep(0,length(w)),550+ (w2-100)*5/9,rep(0,length(w))+ht,lwd=lx,col="red")

lines(c(0,500),c(0,0),lwd=lx); lines(c(550,1050),c(0,0),lwd=lx)
text(-100,ht/2,TeX(paste0("q^{",i,"}")))


xpos1=1350; xpos2=1950; xpos3=2650; xsc=2500/5; ysc=12/5; 

# Labels
text(xpos1 + xsc/2,11,TeX("H_1^P"),cex=); 
text(xpos2 + xsc/2,11,TeX("H_0^P"))
text(xpos3,adj=0.5,11,TeX("# < 5 x 10^{-5}"))


# Top row of plots
ypos=8

polygon(xpos1 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,1])),col="blue",lwd=lx)
lines(xpos1 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2,lwd=lx)
text(xpos1-100,ypos + ysc/2,"P",srt=90)

polygon(xpos2 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,1])),col="blue",lwd=lx)
lines(xpos2 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2,lwd=lx)

text(xpos1-30,ypos+0.2,"0",cex=0.7)
text(xpos1-30,ypos+ysc-0.2,"1",cex=0.7)

nx=length(which(vmat[,1]< (0.05/1000)))
text(xpos3,ypos+(ysc/2),paste0(nx,"/100"),col="blue")


# Second row of plots
ypos=5.5; ix=10

polygon(xpos1 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,1])),col="blue",lwd=lx)
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,ix])),col="red",lwd=lx)
lines(xpos1 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)
text(xpos1-100,ypos + ysc/2,c("P   ","  /  ",expression("       v"[10])),col=c("blue","black","red"),srt=90)


polygon(xpos2 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,1])),col="blue",lwd=lx)
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,ix])),col="red",lwd=lx)
lines(xpos2 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)

text(xpos1-30,ypos+0.2,"0",cex=0.7)
text(xpos1-30,ypos+ysc-0.2,"1",cex=0.7)

nx=length(which(vmat[,ix]< (0.05/1000)))
text(xpos3,ypos+(ysc/2),paste0(nx,"/100"),col="red")



# Third row of plots
ypos=3; ix=50

polygon(xpos1 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,1])),col="blue",lwd=lx)
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,ix])),col="red",lwd=lx)
lines(xpos1 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)
text(xpos1-100,ypos + ysc/2,c("P   ","  /  ",expression("       v"[50])),col=c("blue","black","red"),srt=90)


polygon(xpos2 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,1])),col="blue",lwd=lx)
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,ix])),col="red",lwd=lx)
lines(xpos2 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)

text(xpos1-30,ypos+0.2,"0",cex=0.7)
text(xpos1-30,ypos+ysc-0.2,"1",cex=0.7)

nx=length(which(vmat[,ix]< (0.05/1000)))
text(xpos3,ypos+(ysc/2),paste0(nx,"/100"),col="red")


# Fourth row of plots
ypos=-0.5; ix=500

polygon(xpos1 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,1])),col="blue",lwd=lx)
lines(xpos1 + xsc*((1:100)/101),ypos + ysc*(sort(vmat[1:100,ix])),col="red",lwd=lx)
lines(xpos1 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)
text(xpos1-100,ypos + ysc/2,c("P   ","  /  ",expression("        v"[500])),col=c("blue","black","red"),srt=90)

polygon(xpos2 + xsc*c(0,1,1,0),ypos + ysc*c(0,0,1,1))
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,1])),col="blue",lwd=lx)
lines(xpos2 + xsc*((1:900)/901),ypos + ysc*(sort(vmat[101:1000,ix])),col="red",lwd=lx)
lines(xpos2 + xsc*c(0,1),ypos + ysc*c(0,1),lty=2)

text(xpos1,ypos-0.3,"0",cex=0.7); 
text(xpos1+xsc,ypos-0.3,"1",cex=0.7)
text(xpos2,ypos-0.3,"0",cex=0.7); 
text(xpos2+xsc,ypos-0.3,"1",cex=0.7)

text(xpos1-30,ypos+0.2,"0",cex=0.7)
text(xpos1-30,ypos+ysc-0.2,"1",cex=0.7)

nx=length(which(vmat[,ix]< (0.05/1000)))
text(xpos3,ypos+(ysc/2),paste0(nx,"/100"),col="red")


dev.off()





#############################################################################
## Supplementary figure showing theorem on convergence of L-curve          ##
##  intersections with a line                                              ##
#############################################################################
## Figure 9 (appendix)


if (save_pdf) pdf("./outputs/convergence.pdf",width=5,height=5)

set.seed(1)
q0=1 # line of constant q
n=200 # number of points total
delta=0.1
delta2=0.05
p_eps=0.2 #p_epsilon

Fpq=function(p,q) sqrt(p)*q # CDF of P,Q
Fq=function(q) q # CDF of Q
Cpq=function(p,q) p*Fq(q)/Fpq(p,q)

P=runif(n)^2 # random values of P
px=seq(0,1,length.out=200) # values of p to evaluate function

cfdr=px/ ecdf(P)(px) # q0 is set to 1 in this case, so q is ignored
cfdrt=rev(cummin(rev(cfdr)))

plot(0,xlim=c(0,1),ylim=c(0.3,1),
  xlab="p",ylab=expression(paste("C(p,q"[0],")")),
  xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,at=c(p_eps,1),labels=c(expression("p"[epsilon]),1))
axis(2,at=c(Cpq(p_eps,q0)+delta,1),
  labels=c(expression(paste("C(p"[epsilon],",q"[0],") + ",delta)),1))

w1=which(px < p_eps); w2=which(px>= p_eps)
lines(px,Cpq(px,q0),col="red")
lines(px,Cpq(px,q0)-delta2,col="red",lty=2)
lines(px,Cpq(px,q0)+delta2,col="red",lty=2)
lines(px[w1],cfdrt[w1],col="black",lty=2)
lines(px[w2],cfdrt[w2],col="black",lty=1)

abline(v=p_eps,lty=2)
abline(h=Cpq(p_eps,q0)+delta,lty=2)
lines(c(0,p_eps),c(Cpq(p_eps,q0),Cpq(p_eps,q0)),lty=2)
arrows(p_eps/3,Cpq(p_eps,q0),p_eps/3,Cpq(p_eps,q0)+delta,lty=1,code=3,length=0.05)
text(p_eps/3+0.02,Cpq(p_eps,q0)+delta/2,expression(delta))

legend("bottomright",
  c(expression(paste("C(p,q"[0],")")),
    expression(paste("C(p,q"[0],")±",delta[2])),
    expression(paste(widehat(cFDRt),"(p,q"[0],")"))),
  col=c("red","red","black"),lty=c(1,2,1))

if (save_pdf) dev.off()





#############################################################################
## Supplementary figure for one-point-added section                        ##
#############################################################################
## Figure 10 (appendix)

n=20
q0=0.5
p1=0.3; q1=0.3

set.seed(1)
P=runif(n)^2; Q=runif(n)
xx=seq(0,1,length.out=1000)

m=length(which(Q <= q0))
cf=xx; cf2=xx
for (i in 1:length(xx)) cf[i]=xx[i]*m/length(which(P <= xx[i] & Q <= q0))
cfd=rev(cummin(rev(cf)))

for (i in 1:length(xx)) cf2[i]=xx[i]*(m+1)/length(which(c(P,p1) <= xx[i] & c(Q,q1) <= q0))
cfd2=rev(cummin(rev(cf2)))


if (save_pdf) pdf("./outputs/onepoint_demo.pdf",width=5,height=5)
plot(0,xlim=c(0,1),ylim=c(0,1.2),type="n",xlab="p",ylab="cFDR")
lines(xx,cf,lty=2)
lines(xx,cfd,lty=1)

lines(xx,cf2,col="red",lty=2)
lines(xx,cfd2,col="red",lty=1)
abline(v=q1,col="gray")

htest=0.79
points(htest*5/m,htest)
points(htest*8/(m+1),htest,col="red")

abline(h=htest,col="gray")

legend("bottomright",c(TeX("\\widehat{cFDR}"),TeX("\\widehat{cFDRt}"),
  TeX("\\widehat{cFDR} incl (p^*,q^*)"),TeX("\\widehat{cFDRt}  incl (p^*,q^*)")),
  lty=c(2,1,2,1),col=c("black","black","red","red"))
text(q1,1.1,"q*")

if (save_pdf) dev.off()



# Absolute dfference in cFDRt is less than the difference in cFDR:
if (save_pdf) pdf("./outputs/onepoint_diff.pdf",width=5,height=5)
plot(xx,cf2-cf,type="l",xlab="p",ylab="Change")
lines(xx,cfd2-cfd,col="red")
abline(h=0)
legend("topright",c(TeX("\\widehat{cFDR}"),TeX("\\widehat{cFDRt}")),col=c("black","red"),lty=1)
if (save_pdf) dev.off()



#############################################################################
## Supplementary figure showing similarity of PDF and CDF based L-curves   ##
#############################################################################
## Figure 11 (appendix)


if (save_pdf) pdf("./outputs/lconvergence.pdf",width=6,height=6)

par(mar=c(2,2,1,1))
par(mfrow=c(2,2))


for (i in 1:4) {
  if (i==1) par=c(0.7,0.1,0.1,2,3,1.5,1.5)
  if (i==2) par=c(0.7,0.1,0.1,2,2,3,3)
  if (i==3) par=c(0.99,0.0005,0.0003,2,3,4,3)
  if (i==4) par=c(0.7,0.1,0.05,3,2,2,2)
  #par=c(0.7,0.05,0.2,2,3,3,3)
  
  pi0=par[1];
  pi1=par[2]; 
  pi2=par[3]
  qx1=par[4]; 
  qy1=par[5]
  qx2=par[6];
  qy2=par[7];
  
  Fxy=function(x,y) 
    4*(pi0*pnorm(-abs(x))*pnorm(-abs(y)) + 
        pi1*pnorm(-abs(x),sd=qx1)*pnorm(-abs(y)) +
        pi2*pnorm(-abs(x))*pnorm(-abs(y),sd=qy1) + 
        (1-pi0-pi1-pi2)*pnorm(-abs(x),sd=qx2)*pnorm(-abs(y),sd=qy2))
  fxy=function(x,y) 
    pi0*dnorm(-abs(x))*dnorm(-abs(y)) + 
    pi1*dnorm(-abs(x),sd=qx1)*dnorm(-abs(y)) +
    pi2*dnorm(-abs(x))*dnorm(-abs(y),sd=qy1) + 
    (1-pi0-pi1-pi2)*dnorm(-abs(x),sd=qx2)*dnorm(-abs(y),sd=qy2)
  
  F0xy=function(x,y) 
    (4/(pi0+pi2))*(
      pi0*pnorm(-abs(x))*pnorm(-abs(y)) + 
        pi2*pnorm(-abs(x))*pnorm(-abs(y),sd=qy1))
  f0xy=function(x,y) 
    (1/(pi0+pi2))*(
      pi0*dnorm(-abs(x))*dnorm(-abs(y)) + 
        pi2*dnorm(-abs(x))*dnorm(-abs(y),sd=qy1))
  
  F1xy=function(x,y) Fxy(x,y)-F0xy(x,y)
  f1xy=function(x,y) fxy(x,y)-f0xy(x,y)
  
  Cxy=function(x,y) F0xy(x,y)/Fxy(x,y)
  cxy=function(x,y) f0xy(x,y)/fxy(x,y)
  
  #par(mfrow=c(2,1))
  xmax=10; ymax=15
  xx=seq(0,xmax,length.out=300); yy=seq(0,ymax,length.out=300)
  Fz=outer(xx,yy,Cxy); fz=outer(xx,yy,cxy)
  #image(xx,yy,(-log(Fz))^(1/3),col=rainbow(20))
  #image(xx,yy,(-log(fz))^(1/3),col=rainbow(20))
  
  
  atx=seq(2,xmax,length.out=6)
  atf=cxy(atx,3); atF=Cxy(atx,3)
  contour(xx,yy,Fz,levels=atF,drawlabels=FALSE)
  contour(xx,yy,fz,add=TRUE,col="blue",levels=atf,drawlabels=FALSE)
}

if (save_pdf) dev.off()




#############################################################################
## Supplementary figure showing failure of cFDR<alpha control at alpha     ##
#############################################################################
## Figure 12 (appendix)


# This figure does not require any data

if (save_pdf) pdf("./outputs/ml.pdf",width=5,height=5)

set.seed(1)

plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="P",ylab="Q",xaxs="i",yaxs="i")

q=seq(0,1,length.out=500); k=1/10; cf=0.5; p=cf*k/ (q - cf*q + k)
lines(p,q,col="red")

polygon(c(0,0,p,0,0),c(0,0,q,1,1),col="lightgray",border="red")
polygon(c(0,p[200],p[200],0),c(q[200],q[200],0,0),col="darkgray",border="black")

points(abs(rnorm(50))/50,abs(rnorm(50))/50,pch=16,cex=0.5)
points(p[200],q[200])

text(p[200]/2,q[200]/2,"M")
text(0.05,0.7,"L")


text(p[200]+0.03,q[200]+0.03,"p,q")

if (save_pdf) dev.off()





#############################################################################
## Figure showing chaoticity of cFDR1                                      ##
#############################################################################
## Supplementary material, figure 1

# This figure does not require any data

if (save_pdf) pdf("./outputs/chaoticity.pdf",width=5,height=5)

x1=0.001; y1=0.0025;
x2=0.0025; y2=0.001;
m=max(c(x1,y1,x2,y2))*1.5
m0=min(c(x1,y1,x2,y2))*0.75
m2=m0/2

plot(0,type="n",xlim=c(0,m),ylim=c(0,m),bty="n",xaxs="i",yaxs="i",xlab="P",ylab="Q",xaxt="n",yaxt="n")
axis(1,at=c(0,x1,x2,1),labels=c("0",expression(paste(p[1], " (min P)")),expression(p[2]),""))
axis(2,at=c(0,x1,x2,1),labels=c("0",expression(paste(q[2]," (min Q)")),expression(q[1]),""))


s=seq(0,m0,length.out=50); rs=50:1
z=outer(s,s,function(x,y) sqrt(x^2 + y^2))
cx=function(gl) c(gray(1 - gl*(100:1)/101),rep(gray(1),100))
image(x1+s,y1+s,z,col=cx(0.4),add=T)
image(rev(x1-s),rev(y1-s),z[rs,rs],col=cx(0.9),add=T)
image(rev(x1-s),y1+s,z[rs,],col=cx(1),add=T)
image(x1+s,rev(y1-s),z[,rs],col=cx(0.9),add=T)

image(x2+s,y2+s,z,col=cx(0.2),add=T)
image(rev(x2-s),rev(y2-s),z[rs,rs],col=cx(0.2),add=T)
image(rev(x2-s),y2+s,z[rs,],col=cx(0.6),add=T)
image(x2+s,rev(y2-s),z[,rs],col=cx(0.2),add=T)



points(c(x1,x2),c(y1,y2),pch=16)

np=200; Z=cbind(rnorm(np,mean=max(x1,x2),sd=m0*0.4),rnorm(np,mean=max(y1,y2),sd=m0*0.4)); 
#Z=Z[which(Z[,1]> min(x1,x2)+m2 & Z[,1]< max(x1,x2)+m2 & Z[,2]>min(y1,y2)+ m2 & Z[,2]<max(y1,y2)- m2),]
ccol=c(gray(0.8),gray(0.2))[1+ (Z[,2] < max(y1,y2))]
points(Z,pch=16,col=ccol,cex=0.5)

arrows(m*0.9,y1,m*0.9,y2,code=3); text(m*0.95,(y1+y2)/2,expression(N[Q]))
#arrows(x1,m*0.9,x2,m*0.9,code=3); text((x1+x2)/2,m*0.95,expression(N[P]))

segments(c(x1,x2,x1,x2),c(y1,y2,y1,y2),c(x1,x2,0,0),c(0,0,y1,y2),lty=3)

points(x1+ c(-1,-1,1,1)*m2/2,y1+c(-1,1,-1,1)*m2/2,col="red",pch=16)
points(x2+ c(-1,-1,1,1)*m2/2,y2+c(-1,1,-1,1)*m2/2,col="red",pch=16)

sc=1.2
text(x1+sc*m2/2,y1+sc*m2/2,expression(paste("p'(N"[Q],"+1)/2")),col="red",adj=c(0,0))
text(x1-sc*m2/2,y1+sc*m2/2,expression(paste("p'(N"[Q],"+1)")),col="red",adj=c(1,0))
text(x1+sc*m2/2,y1-sc*m2/2,expression(paste("p'N"[Q])),col="red",adj=c(0,1))
text(x1-sc*m2/2,y1-sc*m2/2,expression(paste("p'N"[Q])),col="red",adj=c(1,1))


text(x2+sc*m2/2,y2+sc*m2/2,"p'",col="red",adj=c(0,0))
text(x2-sc*m2/2,y2+sc*m2/2,"2p'",col="red",adj=c(1,0))
text(x2+sc*m2/2,y2-sc*m2/2,"p'",col="red",adj=c(0,1))
text(x2-sc*m2/2,y2-sc*m2/2,"p'",col="red",adj=c(1,1))


if (save_pdf) dev.off()





#############################################################################
## Supplementary figure showing convergence of PDF and CDF based L-regions ##
#############################################################################
## Figure 3,4,5 in supplementary material

set.seed(1)

convergence_file="./data/convergence_data.RData"

tp=c(0.1,0.1) # We will draw L-curves through this point
pars=c(0.7,0.1,0.15,1.5,3,1.5,3)
n_opts=c(1000,10000,100000)

if (!file.exists(convergence_file)) {
  
  pi0=pars[1]; pi1=pars[2]; pi2=pars[3]; pi3=1-pi0-pi1-pi2; s1=pars[4]; s2=pars[5]; t1=pars[6]; t2=pars[7]
  
  
  
  # True distributions
  pi1q=pi2/(pi0+pi2); pi0q=pi0/(pi0+pi2)
  
  dbase=dnorm
  pdfx0=function(zp,zq) dnorm(zp)*(pi0q*dnorm(zq) + pi1q*dbase(zq/t1)/t1)
  pdfx=function(zp,zq)  pi0*dnorm(zp)*dnorm(zq) +
    pi1*dbase(zp/s1)*dnorm(zq)/s1 +
    pi2*dnorm(zp)*dbase(zq/t1)/t1 +
    pi3*dbase(zp/s2)*dbase(zq/t2)/(s2*t2)
  # have to call this pdfx to avoid confusion with pdf function!
  
  pbase=pnorm
  cdf0=function(zp,zq) pnorm(-zp)*(pi0q*pnorm(-zq) + pi1q*pbase(-zq/t1))
  cdf=function(zp,zq)  pi0*pnorm(-zp)*pnorm(-zq) +
    pi1*pbase(-zp/s1)*pnorm(-zq) +
    pi2*pnorm(-zp)*pbase(-zq/t1) +
    pi3*pbase(-zp/s2)*pbase(-zq/t2)
  
  
  for (i in 1:length(n_opts)) {
    
    N=n_opts[i]
    n0=round(N*pi0); n1=round(N*pi1); n2=round(N*pi2); n3=round(N*pi3)
    ZP=abs(c(rnorm(n0,sd=1),rnorm(n1,sd=s1),rnorm(n2,sd=1),rnorm(n3,sd=s2))) # absolute Z scores (simulated)
    ZQ=abs(c(rnorm(n0,sd=1),rnorm(n1,sd=1),rnorm(n2,sd=t1),rnorm(n3,sd=t2)))
    
    P=2*pnorm(-abs(ZP)); Q=2*pnorm(-abs(ZQ))
    
    fit_pars=fit.4g(cbind(P,Q))$pars
    
    v0=vlo(c(tp[1],P),c(tp[2],Q),cdf0,cdf,indices=1)
    v1=vl(c(tp[1],P),c(tp[2],Q),adj=T,mode=1,indices=1)
    v2=vlx(c(tp[1],P),c(tp[2],Q),adj=T,pars=fit_pars,indices=1)
    v3=vly(c(tp[1],P),c(tp[2],Q),adj=T,mode=2,indices=1)
    
    v0l=vlo(c(tp[1],P),c(tp[2],Q),pdfx0,pdfx,indices=1)
    v2l=vlxl(c(tp[1],P),c(tp[2],Q),pars=fit_pars,indices=1)
    v3l=vlyl(c(tp[1],P),c(tp[2],Q),mode=2,indices=1)
    
    assign(paste0("v0_",i),v0)
    assign(paste0("v1_",i),v1)
    assign(paste0("v2_",i),v2)
    assign(paste0("v3_",i),v3)
    
    assign(paste0("v0l_",i),v0l)
    assign(paste0("v2l_",i),v2l)
    assign(paste0("v3l_",i),v3l)
    
    assign(paste0("P",i),P)
    assign(paste0("Q",i),Q)
    assign(paste0("fit_pars",i),fit_pars)
    
  }
  
  save(pars,n_opts,pdfx0,pdfx,cdf0,cdf,N,v0_1,v1_1,v2_1,v3_1,v0l_1,v2l_1,v3l_1,P1,Q1,fit_pars1,
    v0_2,v1_2,v2_2,v3_2,v0l_2,v2l_2,v3l_2,P2,Q2,fit_pars2,
    v0_3,v1_3,v2_3,v3_3,v0l_3,v2l_3,v3l_3,P3,Q3,fit_pars3,file=convergence_file)
  
} else load(convergence_file)


for (i in 1:3) {
  
  v0=get(paste0("v0_",i))
  v1=get(paste0("v1_",i))
  v2=get(paste0("v2_",i))
  v3=get(paste0("v3_",i))
  
  v0l=get(paste0("v0l_",i))
  v2l=get(paste0("v2l_",i))
  v3l=get(paste0("v3l_",i))
  
  P=get(paste0("P",i))
  Q=get(paste0("Q",i))
  fit_pars=get(paste0("fit_pars",i))
  
  
  if (save_pdf) pdf(paste0("./outputs/lcomp_",n_opts[i],"_cdf.pdf"),width=4,height=4)
  plot(0,type="n",xlab="P",ylab="Q",xaxs="i",yaxs="i",xlim=c(0,0.3),ylim=c(0,1))
  points(tp[1],tp[2],pch=16)
  lines(v0$x[1,],v0$y,col="red")
  lines(v1$x[1,],v1$y,col="black")
  lines(v2$x[1,],v2$y,col="black",lty=2)
  lines(v3$x[1,],v3$y,col="black",lty=3)
  legend("topright",c(expression(widehat(cFDR)^n),expression(widehat(cFDR)^p),expression(widehat(cFDR)^k),"Oracle"),
    lty=c(1,2,3,1),col=c("black","black","black","red"))
  if (save_pdf) dev.off()
  
  
  if (save_pdf) pdf(paste0("./outputs/lcomp_",n_opts[i],"_pdf.pdf"),width=4,height=4)
  plot(0,type="n",xlab="P",ylab="Q",xaxs="i",yaxs="i",xlim=c(0,0.3),ylim=c(0,1))
  points(tp[1],tp[2],pch=16)
  lines(v0l$x[1,],v0$y,col="red")
  lines(v2l$x[1,],v2$y,col="black",lty=2)
  lines(v3l$x[1,],v3$y,col="black",lty=3)
  legend("topright",c(expression(widehat(cfdr)^p),expression(widehat(cfdr)^k),"Oracle"),
    lty=c(2,3,1),col=c("black","black","red"))
  if (save_pdf) dev.off()
  
  
  print(i)
}

