###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Reproducibility check                                                 ##
##                                                                       ##
## James Liley and Chris Wallace, 2020                                   ##
## Correspondence: JL, james.liley@igmm.ed.ac.uk                         ##
##                                                                       ##
###########################################################################
#
# This R script will read simulation results, and for each matrix, choose
#  a random row to recreate.
#
# This is not done for matrix sim_parametric_adjustment.txt as it involves
#  changing the script run_simulation.R
#  
# It will typically take about two hours to run. No outputs are generated 
#  except for results printed to terminal


###########################################################################
## Packages, directories and switches                                    ##
###########################################################################

# Packages

source("code/functions.R") # or library(cfdr)
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



###########################################################################
## Check package and R versions                                          ##
###########################################################################

rcheck=(R.Version()$major=="3" & R.Version()$minor=="3.3")

pkgs=c("mnormt", "mgcv", "pbivnorm", "MASS", "fields", "matrixStats", 
   "latex2exp", "maps", "spam", "grid", "nlme")
vsns=c("1.5.5", "1.8.17", "0.6.0", "7.3.45", "8.10", "0.51.0", "0.4.0", 
   "3.1.1", "1.4.0", "3.3.3", "3.1.131.1")

loaded_versions=unlist(sapply(pkgs,function(x) as.character(packageVersion(x))))

if (all(loaded_versions==vsns) & rcheck) {
 cat("R version and package versions correct.\n") 
} else {
 cat("R version or packages below do not match those used in generation. Simulations may not reproduce\n")
 if (!rcheck) cat(paste0("R version is ",R.Version()$major,".",R.Version()$minor," but should be 3.3.3\n"))
 w=which(!(loaded_versions==vsns))
 str=paste0(paste0(pkgs[w]," is version ",loaded_versions[w]," but should be version ",vsns[w],"\n"),collapse="")
 cat(str)
}



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


# Column names for tables
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
## Check each matrix                                                     ##
###########################################################################


## sim_gen_high_fdr.txt

i0=sample(nrow(sim_gen_hfdr),1)
X=sim_gen_hfdr[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1])
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_gen_high_fdr.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_gen_high_fdr.txt correctly"))






## sim_gen_high_fdr_null.txt

i0=sample(nrow(sim_gen_hfdr0),1)
X=sim_gen_hfdr0[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 1111111111 1111 2 0 0 0.1 ", 
 paste(X[c(3,6,7,8,9,10,4)],collapse=" "))
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_gen_high_fdr_null.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_gen_high_fdr_null.txt correctly"))




## sim_gen_low_fdr.txt

i0=sample(nrow(sim_gen_lfdr),1)
X=sim_gen_lfdr[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 1111111111 1111 2 0 0 0.01 ")
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_gen_low_fdr.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_gen_low_fdr.txt correctly"))




## sim_gen_low_fdr_null.txt

i0=sample(nrow(sim_gen_lfdr0),1)
X=sim_gen_lfdr0[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 1111111111 1111 2 0 0 0.01 ", 
 paste(X[c(3,6,7,8,9,10,4)],collapse=" "))
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_gen_low_fdr_null.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_gen_low_fdr_null.txt correctly"))



## sim_fixed.txt

i0=sample(nrow(sim_fixed),1)
X=sim_fixed[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 1111111111 1111 2 ",X[21]," ",X[22]," 0.1 ", 
 paste(X[c(3,6,7,8,9,10,4)],collapse=" "))
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_fixed.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_fixed.txt correctly"))



## sim_cov.txt

i0=sample(nrow(sim_cov),1)
X=sim_cov[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 0100000000 0011 2 ",X[21]," ",X[22]," 0.1 ")
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_cov.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_cov.txt correctly"))




## sim_cov_null.txt

i0=sample(nrow(sim_cov0),1)
X=sim_cov0[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 0100000000 0011 2 ",X[21]," ",X[22]," 0.1 ", 
 paste(X[c(3,6,7,8,9,10,4)],collapse=" "))
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_cov_null.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_cov_null.txt correctly"))



## sim_unrelated.txt

i0=sample(nrow(sim_ind),1)
X=sim_ind[i0,]
cmd=paste0("Rscript code/run_simulation.R ",X[1]," 0100000000 0011 2 ",X[21]," ",X[22]," ",X[2]," ", 
 paste(X[c(3,6,7,8,9,10,4)],collapse=" "))
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_unrelated.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_unrelated.txt correctly"))




## sim_parametric adjustment.txt
# Note that this does not use run_simulation.R

i0=sample(nrow(sim_paradj),1)
X=sim_paradj[i0,]
cmd=paste0("Rscript code/run_sim_paradj.R ",X[1]," 0000010000 0011 2 0 0 0.1 ")
system(cmd)
r0=read.table(paste0("simulations/cfdrsim",X[1],".txt"))
if (max(abs(r0-X),na.rm=T)< 1e-5) print(paste0("Reproduced row ",i0," of sim_parametric_adjustment.txt correctly")) else 
  print(paste0("Failed to reproduce row ",i0," of sim_parametric_adjustment.txt correctly"))


