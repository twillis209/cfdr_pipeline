###########################################################################
##                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
## Analysis of TWAS data in breast and ovarian cancer                    ##
## James Liley and Chris Wallace                                         ##
##                                                                       ##
###########################################################################

# This script will deterministically generate plots related to to the 
#  transcriptome-wide association study analysis used in the motivating 
#  example for the paper above. 

# The working directory should be one level above this 

###########################################################################
## Packages, scripts and functions ########################################
###########################################################################


library(mgcv)
library(plotrix)
library(latex2exp)
library(cfdr)


###########################################################################
## Switches ###############################################################
###########################################################################

# Draw PDFs: set to F to display plots directly
pdfx=T


###########################################################################
## Dataset ################################################################
###########################################################################

## Raw datasets ./data/TWAS/raw/BCAC.dat and ./data/TWAS/raw/OCAC.dat can be 
##  downloaded from http://twas-hub.org/.
## Initial processing of data is fairly slow so initial results are saved in
##  ./data/TWAS/twas_xy_summary.RData

twasfile="./data/TWAS/twas_summary.RData"

if (!file.exists(twasfile)) {

## Read data
rxbrca0=read.table(paste0(home_dir,"./data/TWAS/raw/BCAC.dat"),header=T)
l1=levels(rxbrca0$PANEL); l1x=l1[which(grepl("GTEx",l1))]; rxbrca=rxbrca0[which(rxbrca0$PANEL %in% l1x),]
cname1=paste0(rxbrca$PANEL,rxbrca$ID); rxbrca=rxbrca[match(unique(cname1),cname1),]; cname1=unique(cname1); rownames(rxbrca)=cname1

rxoca0=read.table(paste0(home_dir,"./data/TWAS/raw/OCAC.dat"),header=T)
l1=levels(rxoca0$PANEL); l1x=l1[which(grepl("GTEx",l1))]; rxoca=rxoca0[which(rxoca0$PANEL %in% l1x),]
cname1=paste0(rxoca$PANEL,rxoca$ID); rxoca=rxoca[match(unique(cname1),cname1),]; cname1=unique(cname1); rownames(rxoca)=cname1

# Extract and name p-values
p_brca=rxbrca$TWAS.P; names(p_brca)=paste0(rxbrca$PANEL,rxbrca$ID)
p_ov=rxov$TWAS.P; names(p_ov)=paste0(rxov$PANEL,rxov$ID)

# Find variables common to both datasets
inx=intersect(names(p_brca),names(p_oca))
p_brca=p_brca[inx]; p_oca=p_oca[inx]

# Remove variables for which p_brca is 0 or NA
w=which(is.finite(-qnorm(p_brca/2)+ -qnorm(p_oca/2)))
p_brca=p_brca[w]; p_oca=p_oca[w]

# Panel and ID
panel=as.character(rxbrca[names(p_brca),]$PANEL); 
id=as.character(rxbrca[names(p_brca),]$ID)

# Folds for CV: same gene is always in the same fold
fold=as.numeric(as.factor(id))
nfold=max(fold)


# Other data
chr=as.character(rxbrca[names(p_brca),]$chr)
pos=as.numeric(rxbrca[names(p_brca),]$pos)
gene=as.character(rxbrca[names(p_brca),]$ID)
tissue=as.character(rxbrca[names(p_brca),]$PANEL)


###########################################################################
## BRCA - determine L-regions #############################################
###########################################################################

# We determine L-regions corresponding to cFDR methods both with and without
#  our adjust for estimating the value of P(Q<q|HP0). In these plots, we 
#  only plot rejections resulting from adjusted cFDR.

sub_brca=which(z_oca^2 + z_brca^2 > 4^2 )

xbrca1=matrix(-1,length(sub_brca),2004); xbrca2=xbrca1

r_brca=rank(p_brca)

totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_brca,inx)
if (length(sub)>0) {
  v1=vl(p_brca,p_ov,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_brca,p_ov,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xbrca1[match(sub,sub_brca),]=v1$x; ybrca1=v1$y
  xbrca2[match(sub,sub_brca),]=v2$x; ybrca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}





###########################################################################
## OCA - determine L-regions ##############################################
###########################################################################


sub_oca=which(z_oca^2 + z_brca^2 > 4^2 )

xoca1=matrix(-1,length(sub_oca),2004); xoca2=xoca1

r_oca=rank(p_oca)

totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_oca,inx)
if (length(sub)>0) {
  v1=vl(p_oca,p_brca,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_oca,p_brca,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xoca1[match(sub,sub_oca),]=v1$x; yoca1=v1$y
  xoca2[match(sub,sub_oca),]=v2$x; yoca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}


###########################################################################
## Integrate over L-regions ###############################################
###########################################################################

# Parametrisation of Q|H0: BRCA|OCA
pars_brca=fit.2g(p_oca[which(p_brca> 0.5)])$pars
pi0_brca=pars_brca[1]
sigma_brca=pars_brca[2]

# Parametrisation of Q|H0: OCA|BRCA
pars_oca=fit.2g(p_brca[which(p_oca> 0.5)])$pars
pi0_oca=pars_oca[1]
sigma_oca=pars_oca[2]

# Integrate over L: BRCA|OCA
i_brca1=il(xbrca1,ybrca1,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # adjusted cFDR
i_brca2=il(xbrca2,ybrca2,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # unadjusted cFDR

# Integrate over L: OCA|BRCA
i_oca1=il(xoca1,yoca1,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # adjusted cFDR
i_oca2=il(xoca2,yoca2,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # unadjusted cFDR




###########################################################################
## Save everything ########################################################
###########################################################################


save(p_brca,p_oca,fold,sub_brca,sub_oca,
     chr,pos,gene,tissue,
     xbrca1,xbrca2,xoca1,xoca2,ybrca1,ybrca2,yoca1,yoca2,
     pi0_brca,sigma_brca,pi0_oca,sigma_oca,
     i_brca1,i_brca2,i_oca1,i_oca2,
     file=twasfile)

} else load(twasfile)





###########################################################################
## Load data ##############################################################
###########################################################################

# Databases

xv1=xbrca1; xv2=xbrca2; yv1=ybrca1; yv2=ybrca2; sub=sub_brca
p=p_brca; q=p_oca; v1=i_brca1; v2=i_brca2; pi0=pi0_brca; sigma=sigma_brca; 




###########################################################################
## General plot parameters ################################################
###########################################################################


k=length(sub); # number of values for which L-regions were computed
xsub=3:2002 # set of indices of finite-valued points of L-curves

# FDR control level
alpha=1e-6

# Benjamini-Hochberg method on p-values
r=rank(p)
hp=which(r <= max(r[which(p/r < alpha/length(p))]))

# Benjamini-Hochberg method on v-values
vbc1=rep(1,length(p)); vbc1[sub]=v1; rvbc1=rank(vbc1)
mvx1=max(rvbc1[which(vbc1/rvbc1 < alpha/length(vbc1))]); mvxl1=which(sub==match(mvx1,rvbc1))
h1=which(rvbc1 <= mvx1)

# Plot labels and parameters
zp=-qnorm(p/2); zq=-qnorm(q/2)
plot_name=paste0("twas_regions_brca")
rlim=3; osub=which(zp^2 + zq^2 > rlim^2)






###########################################################################
## B-H on all p-values rejection region ###################################
###########################################################################


if (pdfx) pdf(paste0("./outputs/",plot_name,"_bh.pdf"),width=5,height=5)

header=expression(paste("A. Rej. reg. from B-H on p"[BRCA]," only"))

xl=suppressWarnings(TeX("$z_{BRCA}$ ($=-\\Phi^{-1}(p_{BRCA}/2)$)"))
yl=suppressWarnings(TeX("$z_{OCA}$ ($=-\\Phi^{-1}(p_{OCA}/2)$)") )
zx=min(zp[hp])

plot(0,type="n",xlim=c(0,8),ylim=c(0,8),xaxs="i",yaxs="i",main=header,xlab=xl,ylab=yl)
polygon(c(zx,zx,100,100),c(0,100,100,0),col="lightblue",border=NA)

points(zp[osub],zq[osub],cex=0.3)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

points(zp[hp],zq[hp],col="blue",pch=5,cex=0.5)
abline(v=zx,lwd=3,lty=2,col="blue")

legend("bottomleft",c("Obs. Z scores", expression(paste("H"[0]," rejected")),
       "Rej. reg. border"),bg="white",
       col=c("black","blue","blue"),pch=c(1,5,NA),lty=c(NA,NA,1),lwd=c(NA,NA,3),
       pt.cex=c(0.3,0.75,NA))

if (pdfx) dev.off()





###########################################################################
## Independent filtering rejection region #################################
###########################################################################

if (pdfx) pdf(paste0("./outputs/",plot_name,"_if.pdf"),width=5,height=5)

qcut=4
vx=which(zq>qcut)
pvx=p[vx]; sqvx=pvx/rank(pvx); pcut=-qnorm(sort(pvx)[max(rank(pvx)[which(sqvx < alpha/length(pvx))])]/2)


header=expression(paste("B. Rej. reg: B-H on p"[BRCA],"|p"[OCA],"<10"^-4))

xl=suppressWarnings(TeX("$z_{BRCA}$"))
yl=suppressWarnings(TeX("$z_{OCA}$") )

plot(0,type="n",xlim=c(0,8),ylim=c(0,8),xaxs="i",yaxs="i",main=header,xlab=xl,ylab=yl)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")
polygon(c(pcut,pcut,100,100),c(qcut,100,100,qcut),col="lightblue",border=NA)

points(zp[osub],zq[osub],cex=0.3)

abline(h=qcut,col="red",lwd=3)
lines(c(pcut,pcut),c(qcut,100),col="blue",lty=2,lwd=3)
lines(c(pcut,100),c(qcut,qcut)+0.01,col="blue",lty=2,lwd=3)
abline(v=min(zp[hp]),lwd=3,lty=2,col="darkgray")

wp=which(zp>pcut & zq>qcut)
points(zp[wp],zq[wp],col="blue",pch=5,cex=0.5)

legend("bottomleft",c("Obs. Z scores", expression(paste("H"[0]," rejected")),
       "Rej. reg. border","Prev. rej. reg."),bg="white",
       col=c("black","blue","blue","darkgray"),pch=c(1,5,NA,NA),lty=c(NA,NA,2,2),lwd=c(NA,NA,3,3),
       pt.cex=c(0.3,0.75,0.75))

if (pdfx) dev.off()






###########################################################################
## cFDR-based rejection region ############################################
###########################################################################

if (pdfx) pdf(paste0("./outputs/",plot_name,".pdf"),width=5,height=5)

border_check=rep(1000,length(v1))
for (i in (mvx1-20):(mvx1+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h1) - in.out(cbind(xv1[order(v1)[i],],yv1),cbind(p,q))))
bmin=order(v1)[which.min(border_check)]
xb=-qnorm(xv1[bmin,xsub]/2); yb=-qnorm(yv1[xsub]/2)

header=expression(paste("C. Rej. reg: cFDR-based"))

xl=suppressWarnings(TeX("$z_{BRCA}$"))
yl=suppressWarnings(TeX("$z_{OCA}$") )

plot(0,type="n",xlim=c(0,8),ylim=c(0,8),xaxs="i",yaxs="i",main=header,xlab=xl,ylab=yl)
polygon(c(xb,xb[length(xb)],100,100),c(yb,100,100,0),col="lightblue",border=NA)

points(zp[osub],zq[osub],cex=0.3)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

points(zp[h1],zq[h1],col="blue",pch=5,cex=0.5)
lines(xb,yb,col="blue",lwd=3,lty=2)

lines(c(pcut,pcut),c(qcut,100),col="darkgray",lty=2,lwd=3)
lines(c(pcut,100),c(qcut,qcut)+0.01,col="darkgray",lty=2,lwd=3)
abline(v=min(zp[hp]),lwd=3,lty=2,col="darkgray")

border_check=rep(1000,length(v1))
for (i in (mvx1-20):(mvx1+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h1) - in.out(cbind(xv1[order(v1)[i],],yv1),cbind(p,q))))
bmin=order(v1)[which.min(border_check)]
lines(-qnorm(xv1[bmin,xsub]/2),-qnorm(yv1[xsub]/2),col="blue",lwd=3,lty=2)

legend("bottomleft",c("Obs. Z scores", expression(paste("H"[0]," rejected")),
       "Rej. reg. border","Prev. rej. regions"),bg="white",
       col=c("black","blue","blue","darkgray"),pch=c(1,5,NA,NA),lty=c(NA,NA,2,2),lwd=c(NA,NA,3,3),
       pt.cex=c(0.3,0.75,0.75))


if (pdfx) dev.off()



###########################################################################
## L-region on p-value scale ##############################################
###########################################################################

if (pdfx) pdf(paste0("./outputs/",plot_name,"_pval.pdf"),width=5,height=5)


xl=suppressWarnings(TeX("$p_{BRCA}$"))
yl=suppressWarnings(TeX("$p_{OCA}$") )

header="D. L-reg. , p-val scale"

plot(0,type="n",xlim=c(0,1e-7),ylim=c(0,1),xaxs="i",yaxs="i",main=header,xlab=xl,ylab=yl)

polygon(c(0,0,xv1[bmin,xsub],0),c(0,1,yv1[xsub],0),col="lightblue",border=NA)
points(p[osub],q[osub],cex=0.3)

points(p[h1],q[h1],col="blue",pch=5,cex=0.5)

lines(xv1[bmin,xsub],yv1[xsub],col="blue",lwd=3,lty=2)

legend("topright",c("Obs. P-vals", expression(paste("H"[0]," rejected")),
       "Rej. region (L-reg.)"),bg="white",
       col=c("black","blue","blue"),pch=c(1,5,NA),lty=c(NA,NA,2),lwd=c(NA,NA,3),
       pt.cex=c(0.3,0.75,0.75))

if (pdfx) dev.off()






###########################################################################
## OCA rejection region (cFDR) ############################################
###########################################################################


# Data
xv1=xoca1; xv2=xoca2; yv1=yoca1; yv2=yoca2; sub=sub_oca
p=p_oca; q=p_brca; v1=i_oca1; v2=i_oca2; pi0=pi0_oca; sigma=sigma_oca; 


k=length(sub); # number of values for which L-regions were computed
xsub=3:2002 # set of indices of finite-valued points of L-curves

# FDR control level
alpha=1e-6

# Benjamini-Hochberg method on p-values
r=rank(p)
hp=which(r <= max(r[which(p/r < alpha/length(p))]))

# Benjamini-Hochberg method on v-values, adjusted
vbc1=rep(1,length(p)); vbc1[sub]=v1; rvbc1=rank(vbc1)
mvx1=max(rvbc1[which(vbc1/rvbc1 < alpha/length(vbc1))]); mvxl1=which(sub==match(mvx1,rvbc1))
h1=which(rvbc1 <= mvx1)

# Plot labels and parameters
zp=-qnorm(p/2); zq=-qnorm(q/2)
plot_name="twas_regions_oca"
rlim=3; osub=which(zp^2 + zq^2 > rlim^2)



# save as pdf if switch is set
if (pdfx) pdf(paste0("./outputs/",plot_name,".pdf"),width=5,height=5)

border_check=rep(1000,length(v1))
for (i in (mvx1-20):(mvx1+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h1) - in.out(cbind(xv1[order(v1)[i],],yv1),cbind(p,q))))
bmin=order(v1)[which.min(border_check)]
xb=-qnorm(xv1[bmin,xsub]/2); yb=-qnorm(yv1[xsub]/2)

header=expression(paste("Rej. reg: cFDR-based, OCA|BRCA"))

xl=suppressWarnings(TeX("$z_{OCA}$") )
yl=suppressWarnings(TeX("$z_{BRCA}$"))

plot(0,type="n",xlim=c(0,8),ylim=c(0,8),xaxs="i",yaxs="i",main=header,xlab=xl,ylab=yl)
polygon(c(xb,xb[length(xb)],100,100),c(yb,100,100,0),col="lightblue",border=NA)

points(zp[osub],zq[osub],cex=0.3)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

points(zp[h1],zq[h1],col="blue",pch=5,cex=0.5)
lines(xb,yb,col="blue",lwd=3,lty=2)

border_check=rep(1000,length(v1))
for (i in (mvx1-20):(mvx1+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h1) - in.out(cbind(xv1[order(v1)[i],],yv1),cbind(p,q))))
bmin=order(v1)[which.min(border_check)]
lines(-qnorm(xv1[bmin,xsub]/2),-qnorm(yv1[xsub]/2),col="blue",lwd=3,lty=2)

legend("bottomleft",c("Obs. Z scores", expression(paste("H"[0]," rejected")),
       "Rej. reg. border"),bg="white",
       col=c("black","blue","blue"),pch=c(1,5,NA),lty=c(NA,NA,2),lwd=c(NA,NA,3),
       pt.cex=c(0.3,0.75,0.75))

if (pdfx) dev.off()

