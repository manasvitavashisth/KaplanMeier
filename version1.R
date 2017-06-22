rm(list=ls())
setwd("~/Box Sync/Manasvita Vashisth/KM Plot/Melanoma")
df.input = read.table("Phenotype", header = TRUE,sep="\t")
df.map=read.table("Proteins.tsv",header=TRUE,sep="\t",na.strings=c("", "NA"))
df.map <- na.omit(df.map)
#censor=df.input[2]
name=df.input[1]
namepro=df.map[1]
#df.input.transpose = t(df.input) #transpose
x <- colnames(df.map)
#Trk <- data.frame(matrix(ncol = dim(df.map)[1]+2, nrow = 0))
#Trk=matrix(0,dim(df.map)[1],dim(df.map)[2]+1)
Trk=df.map
for(i in 1:dim(df.map)[1])
{
  for(j in 1:dim(df.input)[1])
  {
    if(name[j,1]==namepro[i,1])
       {
         Trk[i,dim(df.map)[2]+1]=df.input[j,11]/365
         Trk[i,dim(df.map)[2]+2]=df.input[j,2]
         break
       }
  }
}
Trk <- na.omit(Trk)

library(maxstat)
library(survival)

MLANAtest=maxstat.test(Surv(V67, V68) ~ MLANA,data=Trk, smethod="LogRank",pmethod="condMC")
splitMLANA <- rep(1, nrow(Trk))
Trk<- cbind(Trk, splitMLANA)
Trk$splitMLANA[Trk$MLANA <= MLANAtest$estimate] <- 0
plot(survfit(Surv(V67, V68) ~ splitMLANA, data=Trk),
     xlab = "Survival time in years",
     ylab="Probability",mark=3,col=c(1,2))
x=MLANAtest$estimate
lgn=c(paste("MLANA<",x),paste("MLANA>",x))
legend("topright",legend=lgn,text.col=c(1,2),lty = c(1,2))
ttl=paste("Log rank p-value=",MLANAtest$p.value)
title(ttl)
