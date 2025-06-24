# Depends
library(data.table)
library(stringr)
library(plyr)
library(sktools)
library(lubridate)
library(APtools)

# Load baseline
deciles <- fread(file='output_decile_long.csv')

# Split
setorder(deciles,decile,-Intervention)

## PR
pdf(file='barplot.pdf',height=6,width=16,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(6,6,5,1),xpd=TRUE)

coords <- barplot(deciles$af*100,
                  space=c(0.05,0.05,rep(c(1,0.05),9)),border=NA,
                  col=rep(c('#1f78b48C','#33a02c8C'),10),ylim=c(0,20),
                  yaxt='n',xaxt='n')

avg_coords <- c()
for (i in seq(1,19,2)){
  avg_coords <- c(avg_coords,(coords[i]+coords[i+1])/2)
}

for (i in 0:9){
  mtext(paste0(i*10,'-',i*10+10),side=1,cex=3,line=1.8,at=avg_coords[i+1])
}

for (i in 1:20){
  arrows(coords[i],deciles$lci[i]*100,coords[i],deciles$uci[i]*100,
         col='black',code=3,lwd=1.4,angle=90,length=0.05)
}

axis(2,cex.axis=2.5,las=2,at=seq(0,20,5),pos=coords[1]-0.6)
mtext("New AF diagnosis rate (%)",side=2,cex=3,line=1)
mtext("CH-AI AF risk percentile",side=1,cex=3,line=5.5)

legend(x=min(coords),y=20,legend=c("Screening",
                                   "Control"),
       col=c('#1f78b48C','#33a02c8C'),lwd=5,bty='n',cex=3)

dev.off()