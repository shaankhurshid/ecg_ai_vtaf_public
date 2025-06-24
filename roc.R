# Depends
library(data.table)
library(plyr)
library(timeROC)
library(survival)
library(APtools)
library(prodlim)
library(Cairo)
library(polspline)
library(nricens)
library(dcurves)
source(file='~functions/functions.R')

# Load data
analysis <- fread(file='analysis.csv')

################################# ROC curves
charge <- timeROC(T=analysis$time_to_af, delta=analysis$incd_af,
                  marker=analysis$CHARGE_AF,cause=1,times=c(1,1.999))
ecgai <- timeROC(T=analysis$time_to_af, delta=analysis$incd_af,
                 marker=analysis$ecg_ai_logit,cause=1,times=c(1,1.999))
chai <- timeROC(T=analysis$time_to_af, delta=analysis$incd_af,
                marker=analysis$ch_ai,cause=1,times=c(1,1.999))

pdf(file='roc.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(charge,1.999,add=T,col='#d95f028C',lwd=1.2)
par(new=TRUE)
plot(ecgai,1.999,add=T,col='#1b9e778C',lwd=1.2)
par(new=TRUE)
plot(chai,1.999,add=T,col='#7570b38C',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.3,0.2,legend=c('CH-AI (0.781)','ECG-AI (0.780)','CHARGE-AF (0.700)'),
       col=c('#7570b38C','#1b9e778C','#d95f028C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

### AUPRC
points_charge <- auprc(data=analysis,time='time_to_af',status='incd_af',marker='CHARGE_AF',eval.t=1.999,tolerance=2)
points_ai <- auprc(data=analysis,time='time_to_af',status='incd_af',marker='ecg_ai_logit',eval.t=1.999,tolerance=2)
points_chai <- auprc(data=analysis,time='time_to_af',status='incd_af',marker='ch_ai',eval.t=1.999,tolerance=2)

pr_no_skill <- APSurv(stime=analysis$time_to_af,status=analysis$incd_af,
                      marker=analysis$ecg_ai_logit,t0.list=c(1,1.999))$ap_summary[4]

# Plot
pdf(file='auprc.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))

plot(x=points_ai$sens,y=points_ai$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#1b9e778C',type='l')
par(new=TRUE)
plot(x=points_charge$sens,y=points_charge$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#d95f028C',type='l')
par(new=TRUE)
plot(x=points_chai$sens,y=points_chai$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#7570b38C',type='l')
par(new=TRUE)

axis(1,at=seq(0,1,0.2),cex.axis=1.6)
axis(2,at=seq(0,1,0.1),las=2,cex.axis=1.6)

mtext("Sensitivity/Recall",1,line=2.8,cex=1.8)
mtext("Precision/PPV",2,line=3.5,cex=1.8)

segments(0,pr_no_skill,1,pr_no_skill,lty=5)

legend(0.3,1,legend=c('CH-AI (0.131)','ECG-AI (0.129)','CHARGE-AF (0.0935)'),
       col=c('#7570b38C','#1b9e778C','#d95f028C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)

dev.off()