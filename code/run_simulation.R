###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Run simulation                                                        ##
##                                                                       ##
## James Liley and Chris Wallace                                         ##
##                                                                       ##
###########################################################################

###########################################################################
## Directory to write results to ##########################################
###########################################################################

# Assumes current working directory is folder containing file ./data, 
#  ./code, ./output and ./simulations. An example output with seed=1 is 
#  shown in ./simulations.
out_dir="./simulations/"

###########################################################################
## Packages and scripts ###################################################
###########################################################################

library(cfdr)
library(mnormt)
library(MASS)
library(fields)
library(matrixStats)


###########################################################################
## Simulation parameters ##################################################
###########################################################################

if (!exists("seed")) { #choose seed based on system time if it is not already set
  options(digits.secs=8)
  seed=as.numeric(substr(Sys.time(),21,27))
}
set.seed(seed) # random seed

distx=sample(3,1)  # 1 for normal, 2 for t (3df), 3 for Cauchy

par_cont=runif(1);
if (par_cont>0.9) { # select pars from a discrete distribution
  
  alpha= 0.1 # universal
  nsnp=sample(c(1000,10000),1) # number of variables
  
  n1p=sample(c(0,10,200),1) # number of variables associated ONLY with P (principal)
  n1q=sample(c(0,10,200),1) # number of variables associated ONLY with Q (conditional)
  n1pq=sample(c(0,10,200),1) # number of variables associated with both
  
  sp=sample(c(1.5,3),1) # scale of effect size distribution for associations with P
  sq=sample(c(1.5,3),1) # scale of effect size distribution for associations with Q
  
} else {
  
  alpha= 0.1 
  nsnp=round(10^runif(1,3,4)) 
  
  n1p=sample(200,1) # number of variables associated ONLY with P (principal)
  n1q=sample(200,1) # number of variables associated ONLY with Q (conditional)
  n1pq=sample(200,1) # number of variables associated with both
  
  sp=runif(1,1.5,3) # scale of effect size distribution for associations with A
  sq=runif(1,1.5,3) # scale of effect size distribution for associations with B
  
}



###########################################################################
## Simulation data ########################################################
###########################################################################

if (distx==1) {
  zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}
if (distx==2) {
  zp=c(rt(n1pq,df=3)*sp,rt(n1p,df=3)*sp,rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rt(n1pq,df=3)*sq,rnorm(n1p,sd=1),rt(n1q,df=3)*sq,rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}
if (distx==3) {
  zp=c(rcauchy(n1pq,sc=sp),rcauchy(n1p,sc=sp),rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rcauchy(n1pq,sc=sq),rnorm(n1p,sd=1),rcauchy(n1q,sc=sq),rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}

# P-values
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

mp=min(p[which(p>0)]); p[which(p==0)]=mp
mq=min(q[which(q>0)]); q[which(q==0)]=mq

zp=-qnorm(p/2); zq=-qnorm(q/2)

# H (hypothesis indicators)
if (n1pq+n1p > 0) h1a=1:(n1pq+n1p) else h1a=c()
h0a=setdiff(1:nsnp, h1a)

# Fit four-groups model
parst=c((nsnp-n1p-n1q)/nsnp,n1p/nsnp,n1q/nsnp,sq,sq,sp,sp) # true parameter values
fit.4gs=function(P,pars=parst,...) {
	f1=fit.4g(P,pars=c(0.7,0.1,0.1,2,2,2,2),...); f2=fit.4g(P,pars=pars,...)
	if (f1$lhood> f2$lhood) return(f1) else return(f2)
} # quick function to scan likelihood surface
mxit=1000
fit1=fit.4gs(cbind(zp,zq),maxit=mxit,pars=parst)
pars=fitx$pars
conv=dim(fitx$hist)[1] < mxit

# Folds for method 3
nfold=3
fold=rep(1:nfold,1+floor(nsnp/nfold))[1:nsnp][order(runif(nsnp))]

# Parmameters for underlying distribution of zq|HP0
ff=tryCatch(fit.2g(q[which(p>0.5)]),
  error=function(e) list(pars=c(0.5,1)), warning=function(w) list(pars=c(0.5,1)))

pi0_null = c(1-(n1q/(nsnp-n1p-n1pq)),ff$pars[1])
sigma_null = c(sq,ff$pars[2])
dist_null = c(distx,1)

# Potential rejections
subx=which(log10(p)+log10(q) < -2)
cf=cfdr(p,q,sub=subx)
sub=which(cf<6*alpha)
if (length(sub)<30) sub=order(cf)[1:30]


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

pi1p=n1p/nsnp; pi1q=n1q/nsnp; pi1pq=n1pq/nsnp; pi0=1-pi1p-pi1q-pi1pq
pi1q0=n1q/(nsnp-n1p-n1pq); pi0q0=1-pi1q0

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
  
  state_vec=c(seed,alpha,nsnp,dist_null,n1p,n1q,n1pq,sp,sq,pars,conv,ff$pars)
  names(state_vec)=c("seed","alpha","N","dist1","dist2","n1p","n1q","n1pq","sp","sq",
    paste0("fit_",c("pi0","pi1","pi2","tau1","tau2","s1","s2","conv")),"pi0_null","sigma_null")
  names(fdp_vec)=gsub("hit","fdp",vars)
  names(tdr_vec)=gsub("hit","tdr",vars)
  
  out=c(state_vec,fdp_vec,tdr_vec)
  
  write(out,file=paste0(out_dir,"cfdrsim",seed,".txt"),ncol=length(out))
}

# Initialise
for (i in 1:length(vars)) {
	if (exists(vars[i])) rm(list=vars[i])
}




###########################################################################
## Rejections with p-value: Benjamini-Hochberg procedure ##################
###########################################################################

# Rejections with p-value: Benjamini-Hochberg procedure
rp=rank(p)
hit_p0=which(p< alpha*rank(p)/length(p))
if (length(hit_p0)> 0) hit_p=which(rp <= max(rp[hit_p0])) else hit_p=c()



###########################################################################
## Standard non-parametric cFDR ###########################################
###########################################################################

# Needed several times
vLf=vl(p,q,indices=sub,mode=0,adj=F); 
vLt=vl(p,q,indices=sub,mode=0,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

### FDR control method 1
# No adjustment for Pr(H0|Q), using true null distribution
hit_cf1_fdr1_adj0_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

# using estimated null distribution
hit_cf1_fdr1_adj0_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf1_fdr1_adj1_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf1_fdr1_adj1_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


save_vec()


### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf1_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()



# Adjustment for Pr(H0|Q)
hit_cf1_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()




### FDR control method 3/3b
# No adjustment for Pr(H0|Q)
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


save_vec()

# Adjustment for Pr(H0|Q)
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
save_vec()





### FDR control method 4
# No adjustment for Pr(H0|Q)
vLf=vl(p,q,indices=sub,adj=F,mode=1)
iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf1_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()



# Adjustment for Pr(H0|Q)
vLt=vl(p,q,indices=sub,adj=T,mode=1)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf1_fdr4_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr4_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()




###########################################################################
## KDE-based cFDR #########################################################
###########################################################################

## general
vLf=vly(p,q,indices=sub,adj=F); 
vLt=vly(p,q,indices=sub,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)



### FDR control method 1
# No adjustment for Pr(H0|Q)
hit_cf2_fdr1_adj0_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf2_fdr1_adj0_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf2_fdr1_adj1_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf2_fdr1_adj1_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

save_vec()


### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf2_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf2_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf2_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf2_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()


### FDR control method 3
# No adjustment for Pr(H0|Q)
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


# Adjustment for Pr(H0|Q)
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

save_vec()




###########################################################################
## Parametric cFDR, four-groups ###########################################
###########################################################################

## general
vLf=vlx(p,q,pars,indices=sub,adj=F); 
vLt=vlx(p,q,pars,indices=sub,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)


### FDR control method 1
# No adjustment for Pr(H0|Q)
hit_cf3_fdr1_adj0_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf3_fdr1_adj0_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
save_vec()


# Adjustment for Pr(H0|Q)
hit_cf3_fdr1_adj1_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)])  # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf3_fdr1_adj1_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

save_vec()



### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf3_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf3_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()

# Adjustment for Pr(H0|Q)
hit_cf3_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf3_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()




### FDR control method 3. 
# No adjustment for Pr(H0|Q)
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

# Adjustment for Pr(H0|Q)
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


###########################################################################
## KDE-based local cfdr ###################################################
###########################################################################


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

save_vec()


###########################################################################
## Parametric local cfdr, four-groups #####################################
###########################################################################

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

save_vec()



###########################################################################
## Oracle local cfdr ######################################################
###########################################################################

### Leave-one-out method (4)
vLf=vlo(p,q,f0=pdf0,f=pdf,indices=sub)
iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf6_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf6_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()




###########################################################################
## Oracle CDF-based cFDR ##################################################
###########################################################################

### Leave-one-out method (4)
vLf=vlo(p,q,f0=cdf0,f=cdf,indices=sub)
iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf7_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf7_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()


