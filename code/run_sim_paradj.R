###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Subversion of run_simulation used for a specific test                 ##
##                                                                       ##
## James Liley and Chris Wallace, 2020                                   ##
## Correspondence: JL, james.liley@igmm.ed.ac.uk                         ##
##                                                                       ##
###########################################################################
#
# This script runs analogously to run_simulation.R, but only analysing 
#  parametric cFDR, and using a a parametric rather than ECDF based 
#  adjustment. It is not intended to be used other than to reproduce 
#  simulations


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
subx=which(log10(p)+log10(q) < -2)
cf=cfdr(p,q,sub=subx)
sub=which(cf<6*alpha)
if (length(sub)<30) sub=order(cf)[1:30]


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
## Parametric cFDR, four-groups ###########################################
###########################################################################

## general
if (cft[6]==1) vLt=vlx(p,q,pars,indices=sub,adj=2); 

if (cft[6]==1) iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)


### FDR control method 3. 
# Adjustment for Pr(H0|Q)
if (cft[6]==1 & fdrt[3]==1 & qnt %in% c(0,2)) {
parsx=matrix(0,nfold,7)
hit_cf3_fdr3_adj1_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=2)
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
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=2)
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

save_vec()