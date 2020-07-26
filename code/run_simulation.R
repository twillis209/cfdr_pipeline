###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Script to run a single simulation                                     ##
##                                                                       ##
## James Liley and Chris Wallace, 2020                                   ##
## Correspondence: JL, james.liley@igmm.ed.ac.uk                         ##
##                                                                       ##
###########################################################################
##
## This script does the following:
##  1. Simulates a set of P and Q, given a range of parameters
##  2. For potential 'hits', evaluates either CDF-based cFDR, or PDF-based
##   'local' cFDR in a range of ways (details below)
##  3. Given these cFDR values, determines 'hits' using a range of 
##   putatively FDR-controlling methods
##  4. Evaluates the FDP (false discovery proportion) and TDP (true 
##   discovery proportion) for each set of hits.
##  5. Save results to a file
##
## This script can be invoked with or without arguments. Command arguments
##  are interpreted as follows:
##   1. Random seed. Set to 0 to choose based on time.
##   2. 0-1 string indicating which cFDR types to compute; 0 for 'do not 
##    compute', 1 for 'compute'. String is length 10. The p-value alone is
##    trivially always computed. Digits stand for: 
##     i)    Empirical CDF-based cFDR without estimation of Pr(HP0|Q)
##     ii)   Empirical CDF-based cFDRn with estimation of Pr(HP0|Q)
##     iii)  KDE-estimated CDF based cFDR without estimation of Pr(HP0|Q)
##     iv)   KDE-estimated CDF based cFDR with estimation of Pr(HP0|Q)
##     v)    Parametrically-estimated CDF based cFDR without estimation of Pr(HP0|Q)
##     vi)   Parametrically-estimated CDF based cFDR with estimation of Pr(HP0|Q)
##     vii)  KDE-estimated PDF based (local) cfdr
##     viii) Parametrically-estimated PDF based (local) cfdr
##     ix)    Oracle PDF-based (local) cfdr
##     x)   Oracle CDF-based cFDR
##   3. 0-1 string indicating which FDR-controlling methods to use.
##     i)    Original FDR-controlling method
##     ii)   New naive FDR controlling method (no leave-out)
##     iii)  New FDR controlling method, leave-out-block
##     iv)   New FDR controlling method, leave-one-out
##   4. Indicator variable for whether to use fixed correct distribution 
##    of Q|HP0 or estimate empirically. 0: true; 1: empirical; 2: both
##   5. Correlation indicator between p_i,p_j and between q_i,q_j. 
##    0: independent; 1:correlated when i,j are in same fold; 
##    2: equicorrelated.
##   6. Value of correlation for 5, ignored if argument 5 is 0.
##   7. Value at which to control cFDR (alpha)
##   8. nhyp: number of hypotheses total
##   9. n1p:  number of non-null hypotheses for only P
##   10. n1q: number of non-null associations for only Q
##   11. n1pq: number of non-null hypotheses for both P and Q
##   12. sp: scale of effect size distribution for associations with P
##   13. sq: scale of effect size distribution for associations with Q
##   14. distx: distribution governing z-scores for non-null hypotheses: 1: Gaussian; 2: t-distribution (df=3); 3: Cauchy distribution
## Arguments 8-12 can be omitted, in which case they are chosen randomly.
##  If only one parameter is given, it is interpreted as a random seed.
##  If no command arguments are given, this is equivalent to command args:
##   [random] 1111111111 1111 2 0 0 0.1
##  that is, choose all parameters randomly (either uniformly with prob
##  0.9, or from given fixed values with prob 0.1) and evaluate every 
##  cFDR type and every possible FDR-controlling method. This can be slow.
###


###########################################################################
## Directory for results ##################################################
###########################################################################

# Assumes current working directory is folder containing file ./data, 
#  ./code, ./output and ./simulations. 
out_dir="./simulations/" 
 
###########################################################################
## Packages, scripts, options #############################################
###########################################################################

source("code/functions.R") ## or library(cfdr)
library(mnormt)
library(mvtnorm)
library(MASS)
library(fields)
library(matrixStats)
options(digits.secs=8)



###########################################################################
## Interpret command arguments ############################################
###########################################################################

args=commandArgs(trailingOnly=TRUE)
print(args)
if (!(length(args) %in% c(0,1,7,14))) stop("Incorrect number of command arguments: should be 0, 1, 7 or 14. See script header")

if (length(args)==0) {
 seed=as.numeric(substr(Sys.time(),21,27))
 cft=rep(1,10)
 fdrt=rep(1,4)
 qnt=2
 cort=0
 rho=0
 alpha=0.1
 par_rand=TRUE
} else {
if (length(args)==1) {
 seed=as.numeric(args[1])
 if (seed==0) seed=as.numeric(substr(Sys.time(),21,27))
 cft=rep(1,10)
 fdrt=rep(1,4)
 qnt=2
 cort=0
 rho=0
 alpha=0.1
 par_rand=TRUE
} else {
 seed=as.numeric(args[1])
 if (seed==0) seed=as.numeric(substr(Sys.time(),21,27))
 cft=as.numeric(unlist(strsplit(args[2],"")))
 fdrt=as.numeric(unlist(strsplit(args[3],"")))
 qnt=as.numeric(args[4])
 cort=as.numeric(args[5])
 rho=as.numeric(args[6])
 alpha=as.numeric(args[7])
 if (length(args)==7) par_rand=TRUE else {
  par_rand=FALSE
  nhyp=as.numeric(args[8])
  n1p=as.numeric(args[9])
  n1q=as.numeric(args[10])
  n1pq=as.numeric(args[11])
  sp=as.numeric(args[12])
  sq=as.numeric(args[13])
  distx=as.numeric(args[14])
 }
}
}


###########################################################################
## Simulation parameters ##################################################
###########################################################################

# random seed
set.seed(seed) 


if (par_rand) { # choose parameters randomly

 distx=sample(3,1)  # 1 for normal, 2 for t (3df), 3 for Cauchy

 throwaway=runif(1); # This parameter does not do anything, but it used to, so is kept to keep old simulations reproducible.
   
   nhyp=round(10^runif(1,3,4)) 
  
   n1p=sample(200,1) # number of variables associated ONLY with P (principal)
   n1q=sample(200,1) # number of variables associated ONLY with Q (conditional)
   n1pq=sample(200,1) # number of variables associated with both
  
   sp=runif(1,1.5,3) # scale of effect size distribution for associations with A
   sq=runif(1,1.5,3) # scale of effect size distribution for associations with B
}



###########################################################################
## Simulation data ########################################################
###########################################################################

if (cort==0) { # Select zp and zq independently
 if (distx==1) {
   zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(nhyp-n1p-n1pq,sd=1))
   zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(nhyp-n1p-n1pq-n1q,sd=1))
 }
 if (distx==2) {
   zp=c(rt(n1pq,df=3)*sp,rt(n1p,df=3)*sp,rnorm(nhyp-n1p-n1pq,sd=1))
   zq=c(rt(n1pq,df=3)*sq,rnorm(n1p,sd=1),rt(n1q,df=3)*sq,rnorm(nhyp-n1p-n1pq-n1q,sd=1))
 }
 if (distx==3) {
   zp=c(rcauchy(n1pq,sc=sp),rcauchy(n1p,sc=sp),rnorm(nhyp-n1p-n1pq,sd=1))
   zq=c(rcauchy(n1pq,sc=sq),rnorm(n1p,sd=1),rcauchy(n1q,sc=sq),rnorm(nhyp-n1p-n1pq-n1q,sd=1))
 }
}


# Fold assignments
nfold=3
fold=rep(1:nfold,1+floor(nhyp/nfold))[1:nhyp][order(runif(nhyp))]



if (cort==1) { # Select zp and zq dependently within-fold

 zp=rep(0,nhyp); zq=rep(0,nhyp) # initialise

 for (i in 1:nfold) {
  sub=which(fold==i)
  
  if (n1pq>0) indpq=1:n1pq else indpq=c()
  if (n1p>0) indp=(n1pq+1):(n1pq+n1p) else indp=c()
  if (n1q>0) indq=(n1pq+n1p+1):(n1pq+n1p+n1q) else indq=c()
  ind0=setdiff(1:nhyp,c(indpq,indp,indq))
  
  subpq=sub[which(sub %in% indpq)] # fold=i, associated with both P and Q
  subp=sub[which(sub %in% indp)] # fold=i, associated with P
  subq=sub[which(sub %in% indq)]
  sub0=sub[which(sub %in% ind0)]
  
  cpq=diag(length(subpq)) + rho*(1-diag(length(subpq))) # Correlation matrix for samples subpq
  cp=diag(length(subp)) + rho*(1-diag(length(subp))) # Correlation matrix for samples subp
  cq=diag(length(subq)) + rho*(1-diag(length(subq))) 
  c0=diag(length(sub0)) + rho*(1-diag(length(sub0))) 
  
  z0pi=rmnorm(1,varcov=c0) # values P, within-fold, for null hyps for both P and Q
  z0qi=rmnorm(1,varcov=c0) # values Q, within fold, ...
  
  if (dim(cq)[1]>0) zqpi=rmnorm(1,varcov=cq) else zqpi=c() # values P, within fold, for hyps non-null for only Q
  if (dim(cp)[1]>0) zpqi=rmnorm(1,varcov=cp) else zpqi=c() # values Q, within fold, for hyps non-null for only P
  
  if (distx==1) {
    if (dim(cp)[1]>0) zppi=sp*rmnorm(1,varcov=cp) else zppi=c() # values P for hyps non-null for only P
    if (dim(cq)[1]>0) zqqi=sq*rmnorm(1,varcov=cq) else zqqi=c() # values Q for hyps non-null for only Q
    if (dim(cpq)[1]>0) zxpi=sp*rmnorm(1,varcov=cpq) else zxpi=c() # values P for hyps non-null for both
    if (dim(cpq)[1]>0) zxqi=sq*rmnorm(1,varcov=cpq) else zxqi=c() # values Q for hyps non-null for both
  }
  if (distx==2) {
    if (dim(cp)[1]>0) zppi=sp*rmvt(1,sigma=cp,df=3)  else zppi=c()
    if (dim(cq)[1]>0) zqqi=sq*rmvt(1,sigma=cq,df=3)  else zqqi=c()
    if (dim(cpq)[1]>0) zxpi=sp*rmvt(1,sigma=cpq,df=3) else zxpi=c()
    if (dim(cpq)[1]>0) zxqi=sq*rmvt(1,sigma=cpq,df=3) else zxqi=c()
  }
  if (distx==3) {
    if (dim(cp)[1]>0) zppi=sp*rmvt(1,sigma=cp,df=1)  else zppi=c()
    if (dim(cq)[1]>0) zqqi=sq*rmvt(1,sigma=cq,df=1)  else zqqi=c()
    if (dim(cpq)[1]>0) zxpi=sp*rmvt(1,sigma=cpq,df=1) else zxpi=c()
    if (dim(cpq)[1]>0) zxqi=sq*rmvt(1,sigma=cpq,df=1) else zxqi=c()
  }
  
  zp[sub0]=z0pi; zp[subp]=zppi; zp[subq]=zqpi; zp[subpq]=zxpi
  zq[sub0]=z0qi; zq[subp]=zpqi; zq[subq]=zqqi; zq[subpq]=zxqi    
 }
}







if (cort==2) { # Select zp and zq dependently overall (conditioned on class)

 zp=rep(0,nhyp); zq=rep(0,nhyp) # initialise

 if (n1pq>0) subpq=1:n1pq else subpq=c()
  if (n1p>0) subp=(n1pq+1):(n1pq+n1p) else subp=c()
  if (n1q>0) subq=(n1pq+n1p+1):(n1pq+n1p+n1q) else subq=c()
  sub0=setdiff(1:nhyp,c(subpq,subp,subq))
    
  cpq=diag(length(subpq)) + rho*(1-diag(length(subpq))) # Correlation matrix for samples subpq
  cp=diag(length(subp)) + rho*(1-diag(length(subp))) # Correlation matrix for samples subp
  cq=diag(length(subq)) + rho*(1-diag(length(subq))) 
  c0=diag(length(sub0)) + rho*(1-diag(length(sub0))) 
  
  z0pi=rmnorm(1,varcov=c0) # values P, within-fold, for null hyps for both P and Q
  z0qi=rmnorm(1,varcov=c0) # values Q, within fold, ...
  
  if (dim(cq)[1]>0) zqpi=rmnorm(1,varcov=cq) else zqpi=c() # values P, within fold, for hyps non-null for only Q
  if (dim(cp)[1]>0) zpqi=rmnorm(1,varcov=cp) else zpqi=c() # values Q, within fold, for hyps non-null for only P
  
  if (distx==1) {
    if (dim(cp)[1]>0) zppi=sp*rmnorm(1,varcov=cp) else zppi=c() # values P for hyps non-null for only P
    if (dim(cq)[1]>0) zqqi=sq*rmnorm(1,varcov=cq) else zqqi=c() # values Q for hyps non-null for only Q
    if (dim(cpq)[1]>0) zxpi=sp*rmnorm(1,varcov=cpq) else zxpi=c() # values P for hyps non-null for both
    if (dim(cpq)[1]>0) zxqi=sq*rmnorm(1,varcov=cpq) else zxqi=c() # values Q for hyps non-null for both
  }
  if (distx==2) {
    if (dim(cp)[1]>0) zppi=sp*rmvt(1,sigma=cp,df=3)  else zppi=c()
    if (dim(cq)[1]>0) zqqi=sq*rmvt(1,sigma=cq,df=3)  else zqqi=c()
    if (dim(cpq)[1]>0) zxpi=sp*rmvt(1,sigma=cpq,df=3) else zxpi=c()
    if (dim(cpq)[1]>0) zxqi=sq*rmvt(1,sigma=cpq,df=3) else zxqi=c()
  }
  if (distx==3) {
    if (dim(cp)[1]>0) zppi=sp*rmvt(1,sigma=cp,df=1)  else zppi=c()
    if (dim(cq)[1]>0) zqqi=sq*rmvt(1,sigma=cq,df=1)  else zqqi=c()
    if (dim(cpq)[1]>0) zxpi=sp*rmvt(1,sigma=cpq,df=1) else zxpi=c()
    if (dim(cpq)[1]>0) zxqi=sq*rmvt(1,sigma=cpq,df=1) else zxqi=c()
  }
  
  zp[sub0]=z0pi; zp[subp]=zppi; zp[subq]=zqpi; zp[subpq]=zxpi
  zq[sub0]=z0qi; zq[subp]=zpqi; zq[subq]=zqqi; zq[subpq]=zxqi    

}




# P-values
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

mp=min(p[which(p>0)]); p[which(p==0)]=mp
mq=min(q[which(q>0)]); q[which(q==0)]=mq

zp=-qnorm(p/2); zq=-qnorm(q/2)

# H (hypothesis indicators)
if (n1pq+n1p > 0) h1a=1:(n1pq+n1p) else h1a=c()
h0a=setdiff(1:nhyp, h1a)





###########################################################################
## Properties #############################################################
###########################################################################


# Fit four-groups model
if (any(cft[c(5,6,8)]>0)) {
parst=c((nhyp-n1p-n1q)/nhyp,n1p/nhyp,n1q/nhyp,sq,sq,sp,sp) # true parameter values
fit.4gs=function(P,pars=parst,...) {
f1=fit.4g(P,pars=c(0.7,0.1,0.1,2,2,2,2),...); f2=fit.4g(P,pars=pars,...)
if (f1$lhood> f2$lhood) return(f1) else return(f2)
 }
mxit=1000
fit1=fit.4gs(cbind(zp,zq),maxit=mxit,pars=parst)
pars=fit1$pars
conv=dim(fit1$hist)[1] < mxit
} else {
pars=rep(NA, 7)
conv=NA
}



# Parmameters for underlying distribution of zq|HP0
if (qnt>0) {
ff=tryCatch(fit.2g(q[which(p>0.5)]),
  error=function(e) list(pars=c(0.5,1)), warning=function(w) list(pars=c(0.5,1)))

pi0_null = c(1-(n1q/(nhyp-n1p-n1pq)),ff$pars[1])
sigma_null = c(sq,ff$pars[2])
dist_null = c(distx,1)
} else {
pi0_null = rep(1-(n1q/(nhyp-n1p-n1pq)),2)
sigma_null = c(sq,sq)
dist_null = c(distx,distx)
}


# Potential rejections
if (n1p + n1q < nhyp/2) {
subx=which(log10(p)+log10(q) < -2)
cf=cfdr(p,q,sub=subx)
sub=which(cf<6*alpha)
if (length(sub)<30) sub=order(cf)[1:30]
} else {
sub=which(p<0.5)
}




###########################################################################
## True PDFs and CDFs #####################################################
###########################################################################



if (distx==1) {
   dbase=dnorm; pbase=pnorm
}
if (distx==2) {
   dbase=function(x) dt(x,df=3); pbase=function(x) pt(x,df=3)
}
if (distx==3) {
   dbase=dcauchy; pbase=pcauchy
}

pi1p=n1p/nhyp; pi1q=n1q/nhyp; pi1pq=n1pq/nhyp; pi0=1-pi1p-pi1q-pi1pq
pi1q0=n1q/(nhyp-n1p-n1pq); pi0q0=1-pi1q0

pdf0=function(zp,zq) dnorm(zp)*(pi0q0*dnorm(zq) + pi1q0*dbase(zq/sq)/sq)
pdf=function(zp,zq)  pi0*dnorm(zp)*dnorm(zq) +
                     pi1p*dbase(zp/sp)*dnorm(zq)/sp +
                     pi1q*dnorm(zp)*dbase(zq/sq)/sq +
                     pi1pq*dbase(zp/sp)*dbase(zq/sq)/(sp*sq)
  
cdf0=function(zp,zq) pnorm(-zp)*(pi0q0*pnorm(-zq) + pi1q0*pbase(-zq/sq))
cdf=function(zp,zq)  pi0*pnorm(-zp)*pnorm(-zq) +
                     pi1p*pbase(-zp/sp)*pnorm(-zq) +
                     pi1q*pnorm(-zp)*pbase(-zq/sq) +
                     pi1pq*pbase(-zp/sp)*pbase(-zq/sq)




###########################################################################
## Save procedure #########################################################
###########################################################################


vars=c("hit_p",as.vector(
  outer(as.vector(
    outer(as.vector(
      outer(as.vector(
        outer(0:1,1:4,function(x,y) paste0("_fdr",y,"_adj",x))),
        1:2,function(x,y) paste0(x,"_dist",y))),
      1:7,function(x,y) paste0("hit_cf",y,x))),
    c("",paste0("_fold",1:nfold)),function(x,y) paste0(x,y))))

save_vec=function() {
  
  fdp_vec=c()
  tdr_vec=c()
  
  for (i in 1:length(vars)) {
    if (exists(vars[i])) {
      hitx=get(vars[i])
      
      # false discovery proportions
      if (length(hitx)>0) {
        if (hitx[1]>0) {
          fdp=length(intersect(hitx,h0a))/length(hitx) 
        } else fdp=-1
      } else fdp=0 # false-discovery proportion by p-value
      fdp_vec=c(fdp_vec,fdp)
      
      # type-2 error rate
      tdr=length(intersect(hitx,h1a))/length(h1a)
      tdr_vec=c(tdr_vec,tdr)
    } else {
      fdp_vec=c(fdp_vec,NA)
      tdr_vec=c(tdr_vec,NA)
    }
  }
  
  state_vec=c(seed,alpha,nhyp,dist_null,n1p,n1q,n1pq,sp,sq,pars,conv,ff$pars,cort,rho)
  names(state_vec)=c("seed","alpha","N","dist1","dist2","n1p","n1q","n1pq","sp","sq",
    paste0("fit_",c("pi0","pi1","pi2","tau1","tau2","s1","s2","conv")),"pi0_null",
    "sigma_null","cor_mode","rho")
  names(fdp_vec)=gsub("hit","fdp",vars)
  names(tdr_vec)=gsub("hit","tdr",vars)
  
  out=c(state_vec,fdp_vec,tdr_vec)
  
  write(out,file=paste0(out_dir,"cfdrsim",seed,".txt"),ncol=length(out))
}

# Initialise
for (i in 1:length(vars)) {
    if (exists(vars[i])) rm(list=vars[i])
}

save_vec()


###########################################################################
## Rejections with p-value: Benjamini-Hochberg procedure ##################
###########################################################################

# Rejections with p-value: Benjamini-Hochberg procedure
rp=rank(p)
hit_p0=which(p< alpha*rank(p)/length(p))
if (length(hit_p0)> 0) hit_p=which(rp <= max(rp[hit_p0])) else hit_p=c()

save_vec()


###########################################################################
## Standard non-parametric cFDR ###########################################
###########################################################################

# Needed several times
if (cft[1]==1) vLf=vl(p,q,indices=sub,mode=0,adj=F); 
if (cft[2]==1) vLt=vl(p,q,indices=sub,mode=0,adj=T); 

if (cft[1]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
if (cft[2]==1) iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

### FDR control method 1
# No adjustment for Pr(H0|Q), using true null distribution
if (cft[1]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr1_adj0_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

# using estimated null distribution
if (cft[1]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr1_adj0_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}



# Adjustment for Pr(H0|Q)
if (cft[2]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr1_adj1_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

if (cft[2]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr1_adj1_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}




save_vec()


### FDR control, new naive method
# No adjustment for Pr(H0|Q)
if (cft[1]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[1]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

save_vec()



# Adjustment for Pr(H0|Q)
if (cft[2]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[2]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

save_vec()



### FDR control, new method (blocks)
# No adjustment for Pr(H0|Q)
if (cft[1]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr3_adj0_dist1 = ({
  hit=c()
  vv_all=rep(1,length(p)) 
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=F)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    vv=rep(1,length(p)); vv[sub2]=vvx[,1]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}

if (cft[1]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr3_adj0_dist2 = ({
  hit=c()
  vv_all=rep(1,length(p)) 
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=F)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    vv=rep(1,length(p)); vv[sub2]=vvx[,2]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}

save_vec()




# Adjustment for Pr(H0|Q)
if (cft[2]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr3_adj1_dist1 = ({
  vv_all=rep(1,length(p)) 
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=T)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)    
    vv=rep(1,length(p)); vv[sub2]=vvx[,1]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
save_vec()
}

if (cft[2]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr3_adj1_dist2 = ({
  vv_all=rep(1,length(p)) 
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=T)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)    
    vv=rep(1,length(p)); vv[sub2]=vvx[,2]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}
save_vec()





### FDR control, new method (leave-one-out)
# No adjustment for Pr(H0|Q)
if (cft[1]==1 & fdrt[4]==1) vLf=vl(p,q,indices=sub,adj=F,mode=1)
if (cft[1]==1 & fdrt[4]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)


if (cft[1]==1 & fdrt[4]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[1]==1 & fdrt[4]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}
save_vec()



# Adjustment for Pr(H0|Q)
if (cft[2]==1 & fdrt[4]==1) vLt=vl(p,q,indices=sub,adj=T,mode=1)
if (cft[2]==1 & fdrt[4]==1) iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

if (cft[2]==1 & fdrt[4]==1 & qnt %in% c(0,2)) {
hit_cf1_fdr4_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[1]==1 & fdrt[4]==1 & qnt %in% c(1,2)) {
hit_cf1_fdr4_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}
save_vec()




###########################################################################
## KDE-based cFDR #########################################################
###########################################################################

## general
if (cft[3]==1) vLf=vly(p,q,indices=sub,adj=F); 
if (cft[4]==1) vLt=vly(p,q,indices=sub,adj=T); 

if (cft[3]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
if (cft[4]==1) iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)



### FDR control method 1
# No adjustment for Pr(H0|Q)
if (cft[3]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr1_adj0_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

if (cft[3]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr1_adj0_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}


# Adjustment for Pr(H0|Q)
if (cft[4]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr1_adj1_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

if (cft[4]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr1_adj1_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

save_vec()




### FDR control method 2
# No adjustment for Pr(H0|Q)
if (cft[3]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[3]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}



# Adjustment for Pr(H0|Q)
if (cft[4]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[4]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

save_vec()


### FDR control method 3
# No adjustment for Pr(H0|Q)
if (cft[3]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr3_adj0_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=F)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,1])
    out=vv[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})
}

if (cft[3]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr3_adj0_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=F)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,2])
    out=vv[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}


# Adjustment for Pr(H0|Q)
if (cft[4]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf2_fdr3_adj1_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=T)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,1])
    out=vv[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf2_fdr3b_adj1_dist1",hit_all) 
  hit
})
}

if (cft[4]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf2_fdr3_adj1_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=T)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,2])
    out=vv[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}

save_vec()




###########################################################################
## Parametric cFDR, four-groups ###########################################
###########################################################################

## general
if (cft[5]==1) vLf=vlx(p,q,pars,indices=sub,adj=F); 
if (cft[6]==1) vLt=vlx(p,q,pars,indices=sub,adj=T); 

if (cft[5]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
if (cft[6]==1) iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)


### FDR control method 1
# No adjustment for Pr(H0|Q)
if (cft[5]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf3_fdr1_adj0_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

if (cft[5]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr1_adj0_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}
save_vec()



# Adjustment for Pr(H0|Q)
if (cft[6]==1 & fdrt[1]==1 & qnt %in% c(0,2)) {
hit_cf3_fdr1_adj1_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)])  # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

if (cft[6]==1 & fdrt[1]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr1_adj1_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
}

save_vec()



### FDR control method 2
# No adjustment for Pr(H0|Q)
if (cft[5]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf3_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[5]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

save_vec()


# Adjustment for Pr(H0|Q)
if (cft[6]==1 & fdrt[2]==1 & qnt %in% c(0,2)) {
hit_cf3_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[6]==1 & fdrt[2]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}
save_vec()




### FDR control method 3. 
# No adjustment for Pr(H0|Q)
if (cft[5]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
parsx=matrix(0,nfold,7)
hit_cf3_fdr3_adj0_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=F)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}

if (cft[5]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr3_adj0_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=F)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}


# Adjustment for Pr(H0|Q)
if (cft[6]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf3_fdr3_adj1_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=T)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})
}

if (cft[6]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf3_fdr3_adj1_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=T)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})
}



###########################################################################
## KDE-based local cfdr ###################################################
###########################################################################

if (cft[7]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
hit_cf4_fdr3_adj0_dist1=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vlyl(p,q,indices=inf,fold=inf,mode=2)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})
}

if (cft[7]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf4_fdr3_adj0_dist2=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vlyl(p,q,indices=inf,fold=inf,mode=2)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}
save_vec()


###########################################################################
## Parametric local cfdr, four-groups #####################################
###########################################################################


if (cft[8]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
parsx=matrix(0,nfold,7)
hit_cf5_fdr3_adj0_dist1=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars,sgm=c(1,1000))$pars #only need to do this once per fold
    vl0=vlxl(p,q,indices=inf,pars=parsx[i,])
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})
}

if (cft[8]==1 & fdrt[3]==1 & qnt %in% c(1,2)) {
hit_cf5_fdr3_adj0_dist2=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars,sgm=c(1,1000))$pars #only need to do this once per fold
    vl0=vlxl(p,q,indices=inf,pars=parsx[i,])
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})
}
save_vec()



###########################################################################
## Oracle local cfdr ######################################################
###########################################################################

### Leave-one-out method (4)
if (cft[9]==1) vLf=vlo(p,q,f0=pdf0,f=pdf,indices=sub)
if (cft[9]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

if (cft[9]==1 & qnt %in% c(0,2)) {
hit_cf6_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[9]==1 & qnt %in% c(1,2)) { 
hit_cf6_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}
save_vec()




###########################################################################
## Oracle CDF-based cFDR ##################################################
###########################################################################

### Leave-one-out method (4)
if (cft[10]==1) vLf=vlo(p,q,f0=cdf0,f=cdf,indices=sub)
if (cft[10]==1) iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

if (cft[10]==1 & qnt %in% c(0,2)) { 
hit_cf7_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}

if (cft[10]==1 & qnt %in% c(1,2)) { 
hit_cf7_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
}
save_vec()
