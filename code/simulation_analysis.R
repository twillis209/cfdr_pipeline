###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Analyse results of simulations and draw plots                         ##
##                                                                       ##
## James Liley and Chris Wallace                                         ##
##                                                                       ##
###########################################################################
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

source("../cfdr/R/functions.R")
#library(cfdr)
library(mnormt)
library(mgcv)
library(pbivnorm)
library(MASS)
library(fields)
library(matrixStats)
library(latex2exp)

# Directories

# This variable gives the location of the folder cfdr_control. It should 
#  contain a document README, and four folders  named ./code, ./data, 
#  ./simulations, and ./outputs. Here it is assumed to be the current working 
#  directory
cfdr_dir="./"

# This variable gives the location to save outputs (plots, tables) to.
output_dir=paste0(cfdr_dir,"outputs/")



# Switches

# set to F to draw plots rather than saving them as PDFs
save_pdf=T 




###########################################################################
## Read simulation data                                                  ##
###########################################################################

# Read tables. They DO have headers but we regenerate names.
sx0=read.table(paste0(cfdr_dir,"./data/simmat0.1.txt"),fill=TRUE,header=T)
sx1=read.table(paste0(cfdr_dir,"./data/simmat0.01.txt"),fill=TRUE,header=T)

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

nsv=c("seed","alpha","N","dist1","dist2","n1p","n1q","n1pq","sp","sq",
paste0("fit_",c("pi0","pi1","pi2","tau1","tau2","s1","s2","conv")),"pi0_null","sigma_null")
nfdp=gsub("hit","fdp",vars)
ntdr=gsub("hit","tdr",vars)

nx=c(nsv, nfdp,ntdr)
 

colnames(sx0)=nx; colnames(sx1)=nx



###########################################################################
## FDR control using CDF methods                                         ##
###########################################################################

for (xdist in 1:2) {
for (ia in 1:2) {

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
if (ia==1) {
 alpha=0.1
 attach(sx0) 
} else {
 alpha=0.01
 attach(sx1)
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

sdx=0.15 # gaussian smoothing
mwd=2 # line width

xf=unique(xx)

fp=rep(0,length(xf));
f1=fp; f2=fp; f3=fp; f4=fp;

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

for (i in 1:length(xf)) {
   wts=dnorm(xx-xf[i],sd=sdx); wts=wts/sum(wts)
   fp[i]=weighted.mean(txp,wts,na.rm=T)
   f1[i]=weighted.mean(tx1,wts,na.rm=T)
   f2[i]=weighted.mean(tx2,wts,na.rm=T)
   f3[i]=weighted.mean(tx3,wts,na.rm=T)
   f4[i]=weighted.mean(tx4,wts,na.rm=T)
}

prange=c(0,alpha*1.5) #max(c(fp,f1,f2,f3,f4)));
yx="FDR"
plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")

lines(sort(xf),fp[order(xf)],col="black",lty=1,lwd=mwd) # BH
lines(sort(xf),f1[order(xf)],col="red",lty=1,lwd=mwd) # cFDR original
lines(sort(xf),f2[order(xf)],col="blue",lty=1,lwd=mwd) # cFDR naive
lines(sort(xf),f3[order(xf)],col="black",lty=2,lwd=mwd) # cFDR leave-out-block
lines(sort(xf),f4[order(xf)],col="black",lty=3,lwd=mwd) # cFDR leave-one-out

abline(h=alpha,col="gray")

legend("bottomright",bg="white",
     c("P-val.",
     	"Orig.",
       "Naive",
       "LOO",
       "LOB"),
     col=c("black","red","blue","black","black"),
     lty=c(1,1,1,2,3))

if (save_pdf) dev.off()

}
}


###########################################################################
## Power of ECDF estimator                                               ##
###########################################################################

for (ia in 1:2) {

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
if (ia==1) {
 alpha=0.1
 attach(sx0) 
} else {
 alpha=0.01
 attach(sx1)
}

if (save_pdf) pdf(paste0(output_dir,"power1_alpha",ia,".pdf"),width=5,height=5)

normalise=F # set to F to show absolute power rather than relative to fp
tp=tdr_p
t1=tdr_cf1_fdr3_adj1_dist1 # cFDRn 
t2=tdr_cf1_fdr3_adj0_dist1 # cFDR
t3=tdr_cf1_fdr1_adj1_dist1 # cFDR, old method
x=n1p+n1pq

w=which(x>0 &
		    pmin(tp,t1,t2,t3)>-1 & #txp-tx5 < 0.02 &
		    is.finite(rowSums(cbind(tp,t1,t2,t3))) &
        !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions
        #)
txp=tp[w]; 
tx1=t1[w]; tx2=t2[w]; tx3=t3[w]; 
xx=x[w]

# fix an issue with tx5
#tx5[which(fit_s1[w]<fit_s2[w])]=NA


ci=0.9 # confidence envelope width
sdx=0.20 # gaussian smoothing
mwd=2 # line width

xf=unique(xx)

fp=rep(0,length(xf));
f1=fp; f2=fp; f3=fp; 

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

for (i in 1:length(xf)) {
   wts=dnorm(xx-xf[i],sd=sdx); wts=wts/sum(wts)
   fp[i]=weighted.mean(txp,wts,na.rm=T)
   f1[i]=weighted.mean(tx1,wts,na.rm=T)
   f2[i]=weighted.mean(tx2,wts,na.rm=T)
   f3[i]=weighted.mean(tx3,wts,na.rm=T)
 }

if (normalise) {
  f1=f1-fp; f2=f2-fp; f3=f3-fp; fp=fp-fp
}

prange=pmax(0,range(c(fp,f1,f2,f3)));
if (normalise) yx=expression(paste(Delta,"(TDR)")) else yx="TDR"
plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")


lines(sort(xf),fp[order(xf)],col="black",lty=1,lwd=mwd) # BH
lines(sort(xf),f1[order(xf)],col="red",lty=2,lwd=mwd) # cFDR ECDF 
lines(sort(xf),f2[order(xf)],col="black",lty=2,lwd=mwd) # cFDR KDE
lines(sort(xf),f3[order(xf)],col="blue",lty=2,lwd=mwd) # cFDR param

abline(h=0)

legend("bottomright",bg="white",
     c("P-values",
     	expression(widehat("cFDR")^n),
       expression(widehat("cFDR")),
       "Old method"),
     col=c("black","red","black","blue"),
     lty=c(1,2,2,2))


if (save_pdf) dev.off()

}






###########################################################################
## FDR control using PDF methods                                         ##
###########################################################################

for (ia in 1:2) {

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
if (ia==1) {
 alpha=0.1
 attach(sx0) 
} else {
 alpha=0.01
 attach(sx1)
}

if (save_pdf) pdf(paste0(output_dir,"power3_alpha",ia,".pdf"),width=5,height=5)



normalise=F # set to F to show absolute power rather than relative to fp
tp=tdr_p
t1=tdr_cf1_fdr4_adj1_dist1 # cFDR ECDF
t4=tdr_cf4_fdr3_adj0_dist1 # cfdr KDE
t5=tdr_cf5_fdr3_adj0_dist1 # cfdr param
t6=tdr_cf6_fdr4_adj0_dist1 # cfdr oracle
t7=tdr_cf7_fdr4_adj0_dist1 # cFDR oracle
x=n1p+n1pq

w=which(x>0 & # pmin(tx1,tx2,tx3,tx4,tx5,tx6,tx7) >= 0.1 & # comment this line in to only include simulations for which the fitted distribution was close to correct.
		    pmin(tp,t1,t4,t5,t6,t7)>-1 & #txp-tx5 < 0.02 &
		    is.finite(rowSums(cbind(tp,t1,t4,t5,t6,t7))) &
        (dist1 == 1) & !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions
        #)
txp=tp[w]; 
tx1=t1[w]; tx4=t4[w]; tx5=t5[w]; tx6=t6[w]; tx7=t7[w]; 
xx=x[w]


ci=0.9 # confidence envelope width
sdx=0.15 # gaussian smoothing
mwd=2 # line width

xf=unique(xx)

fp=rep(0,length(xf));
f1=fp; f4=fp; f5=fp; f6=fp; f7=fp

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

for (i in 1:length(xf)) {
   wts=dnorm(xx-xf[i],sd=sdx); wts=wts/sum(wts)
   fp[i]=weighted.mean(txp,wts,na.rm=T)
   f1[i]=weighted.mean(tx1,wts,na.rm=T)
   f4[i]=weighted.mean(tx4,wts,na.rm=T)
   f5[i]=weighted.mean(tx5,wts,na.rm=T)
   f6[i]=weighted.mean(tx6,wts,na.rm=T)
   f7[i]=weighted.mean(tx7,wts,na.rm=T)
}

if (normalise) {
  f1=f1-fp; f4=f4-fp; f5=f5-fp;  f6=f6-fp;  f7=f7-fp; fp=fp-fp
}

prange=pmax(0,range(c(fp,f1,f4,f5,f6,f7)));
if (normalise) yx=expression(paste(Delta,"(power)")) else yx="Power"
plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")


lines(sort(xf),fp[order(xf)],col="black",lty=2,lwd=mwd) # BH
lines(sort(xf),f1[order(xf)],col="black",lty=1,lwd=mwd) # cFDR ECDF 
lines(sort(xf),f4[order(xf)],col="gray",lty=3,lwd=mwd) # cfdr KDE
lines(sort(xf),f5[order(xf)],col="blue",lty=3,lwd=mwd) # cfdr param
lines(sort(xf),f6[order(xf)],col="red",lty=3,lwd=mwd) # cfdr oracle
lines(sort(xf),f7[order(xf)],col="red",lty=1,lwd=mwd) # cFDR oracle

abline(h=0)

legend("bottomright",bg="white",
     c("P-value",
     	 expression("cFDR"^n),
       "CDF oracle",
       "PDF, param.",
       "PDF, KDE",
       "PDF oracle"),
     col=c("black","black","red","blue","gray","red"),
     lty=c(2,1,1,3,3,3))


if (save_pdf) dev.off()

}



###########################################################################
## FDR control using CDF methods                                         ##
###########################################################################

for (id in 1:2) {
for (ia in 1:2) {

atx=intersect(search(),ls()); if (length(atx)>0) for (i in 1:length(atx)) detach(atx[i],character=TRUE)
if (ia==1) {
 alpha=0.1
 attach(sx0) 
} else {
 alpha=0.01
 attach(sx1)
}

if (id==1) dists=1 else dists=2:3

if (save_pdf) pdf(paste0(output_dir,"power2_dist",dists,"_alpha",ia,".pdf"),width=5,height=5)

normalise=F # set to F to show absolute power rather than relative to fp
tp=tdr_p
t1=tdr_cf1_fdr4_adj1_dist1 # cFDR ECDF
t2=tdr_cf2_fdr3_adj1_dist1 # cFDR KDE
t3=tdr_cf3_fdr3_adj1_dist1 # cFDR param
t6=tdr_cf6_fdr4_adj0_dist1 # cfdr oracle
x=n1p+n1pq

w=which(x>0 & # pmin(tx1,tx2,tx3,tx4,tx5,tx6,tx7) >= 0.1 & # comment this line in to only include simulations for which the fitted distribution was close to correct.
		    pmin(tp,t1,t2,t3,t6)>-1 & #txp-tx5 < 0.02 &
		    is.finite(rowSums(cbind(tp,t1,t2,t3,t6))) &
        (dist1 %in% dists) & !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions
        #)
txp=tp[w]; 
tx1=t1[w]; tx2=t2[w]; tx3=t3[w]; tx6=t6[w];
xx=x[w]


ci=0.9 # confidence envelope width
sdx=0.15 # gaussian smoothing
mwd=2 # line width

xf=unique(xx)

fp=rep(0,length(xf));
f1=fp; f2=fp; f3=fp;  f6=fp;

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

for (i in 1:length(xf)) {
   wts=dnorm(xx-xf[i],sd=sdx); wts=wts/sum(wts)
   fp[i]=weighted.mean(txp,wts,na.rm=T)
   f1[i]=weighted.mean(tx1,wts,na.rm=T)
   f2[i]=weighted.mean(tx2,wts,na.rm=T)
   f3[i]=weighted.mean(tx3,wts,na.rm=T)
   f6[i]=weighted.mean(tx6,wts,na.rm=T)
}

if (normalise) {
  f1=f1-fp; f2=f2-fp; f3=f3-fp; f6=f6-fp;  fp=fp-fp
}

prange=pmax(0,range(c(fp,f1,f2,f3,f6)));
if (normalise) yx=expression(paste(Delta,"(TDR)")) else yx="TDR"
plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n",xaxs="i",yaxs="i")


lines(sort(xf),fp[order(xf)],col="black",lty=2,lwd=mwd) # BH
lines(sort(xf),f1[order(xf)],col="black",lty=1,lwd=mwd) # cFDR ECDF 
lines(sort(xf),f2[order(xf)],col="gray",lty=1,lwd=mwd) # cFDR KDE
lines(sort(xf),f3[order(xf)],col="blue",lty=1,lwd=mwd) # cFDR param
lines(sort(xf),f6[order(xf)],col="red",lty=1,lwd=mwd) # cfdr oracle

abline(h=0)

legend("bottomright",bg="white",
     c("P-value",
     	 "ECDF",
       "Param.",
       "KDE",
       "Oracle"),
     col=c("black","black","blue","gray","red"),
     lty=c(2,1,1,1,1))


if (save_pdf) dev.off()

}
}





#############################################################################
## Generate and draw simulations for iterated cFDR                         ##
#############################################################################

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





#############################################################################
## Supplementary figure showing theorem on convergence of L-curve          ##
##  intersections with a line                                              ##
#############################################################################

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

if (save_pdf) pdf("lconvergence.pdf",width=6,height=6)

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
sx1=par[4]; 
sy1=par[5]
sx2=par[6];
sy2=par[7];

Fxy=function(x,y) 
  4*(pi0*pnorm(-abs(x))*pnorm(-abs(y)) + 
  pi1*pnorm(-abs(x),sd=sx1)*pnorm(-abs(y)) +
  pi2*pnorm(-abs(x))*pnorm(-abs(y),sd=sy1) + 
  (1-pi0-pi1-pi2)*pnorm(-abs(x),sd=sx2)*pnorm(-abs(y),sd=sy2))
fxy=function(x,y) 
  pi0*dnorm(-abs(x))*dnorm(-abs(y)) + 
  pi1*dnorm(-abs(x),sd=sx1)*dnorm(-abs(y)) +
  pi2*dnorm(-abs(x))*dnorm(-abs(y),sd=sy1) + 
  (1-pi0-pi1-pi2)*dnorm(-abs(x),sd=sx2)*dnorm(-abs(y),sd=sy2)

F0xy=function(x,y) 
  (4/(pi0+pi2))*(
    pi0*pnorm(-abs(x))*pnorm(-abs(y)) + 
    pi2*pnorm(-abs(x))*pnorm(-abs(y),sd=sy1))
f0xy=function(x,y) 
  (1/(pi0+pi2))*(
    pi0*dnorm(-abs(x))*dnorm(-abs(y)) + 
    pi2*dnorm(-abs(x))*dnorm(-abs(y),sd=sy1))

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



