library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(gplots)
library(timeROC)
library(ggsci)

tcga<-read.table("TCGA_exp_sig.txt",header = T,sep = "\t",quote = "",check.names = F)
GSE57303<-read.table("GSE57303os.txt",header = T,sep = " ",quote = "\"")
names(GSE57303)[1:2]<-c("OS","OS.time")
GSE57303$OS.time<-GSE57303$OS.time/12
GSE57303<-cbind(sample=rownames(GSE57303),GSE57303)
#write.table(GSE57303,"GSE57303os.txt",col.names = T,row.names = F,sep = "\t",quote = F)

GSE62254<-read.table("GSE62254os.txt",header = T,sep = " ",quote = "\"")
names(GSE62254)[1:2]<-c("OS","OS.time")
GSE62254$OS.time<-GSE62254$OS.time/12
GSE62254<-cbind(sample=rownames(GSE62254),GSE62254)
#write.table(GSE62254,"GSE62254os.txt",col.names = T,row.names = F,sep = "\t",quote = F)

mm<-list(TCGA=tcga,
         GSE57303=GSE57303,GSE62254=GSE62254)

mm <- lapply(mm,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
est_data <- mm$TCGA
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
#rm(mm)
rf_nodesize <- 5
seed <- 123
###StepCox[forward]+Lasso###
fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = "forward")
rid <- names(coef(fit))
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))
set.seed(seed)
library(ggsci)
mypal<-pal_npg("nrc")(10)
fit = cv.glmnet(x1, x2,
                nfold=10, 
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='response',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
###timeROC
rt<-rs$GSE62254
library(timeROC)
ROC.bili.marginal<-timeROC(T=rt$OS.time,
                           delta=rt$OS,marker=rt$RS,
                           cause=1,weighting="marginal",
                           times=quantile(rt$OS.time,probs=seq(0,1,0.1)),
                           iid=TRUE)
plotAUCcurve(ROC.bili.marginal,conf.int=TRUE,conf.band=TRUE,col="#00A087FF")
###C指数条形图
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}),
                 se=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[2])}))%>%
  rownames_to_column('ID')
library(gplots)
Cindex<-cc$Cindex
names(Cindex)<-cc$ID
write.table(cc,"Cindex.txt",col.names = T,row.names = F,sep = "\t",quote = F)

barplot2(Cindex,plot.ci=TRUE,ci.l=c(cc$Cindex-cc$se),ci.u=c(cc$Cindex+cc$se),
         col=mypal[1:3],plot.grid = TRUE,border = "black",ci.lwd=3,grid.lwd	=3,
         cex.axis=1.5,xpd=FALSE,main = "C-index")
library(compareC)
tcga_cli<-read.table("TCGA_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(tcga$sample,tcga_cli$sample)
tcga_cli<-tcga_cli[loc,]
GSE57303_cli<-read.table("GSE57303_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE57303$sample,GSE57303_cli$sample)
GSE57303_cli<-GSE57303_cli[loc,]
GSE62254_cli<-read.table("GSE62254_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE62254$sample,GSE62254_cli$sample)
GSE62254_cli<-GSE62254_cli[loc,]

###TCGA
rt<-tcga_cli
tcga_rs<-list(RS=rs$TCGA,Age=rt[,c(2:4)],Grade=rt[,c(2:3,5)],pM=rt[,c(2:3,6)],
              pN=rt[,c(2:3,7)],pT=rt[,c(2:3,8)],Gender=rt[,c(2:3,9)],Stage=rt[,c(2:3,10)])

cc <- data.frame(Cindex=sapply(tcga_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[1])}),
                 se=sapply(tcga_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[2])}))%>%
  rownames_to_column('ID')
write.table(cc,"TCGA_cli_Cindex.txt",col.names = T,row.names = F,sep = "\t",quote = F)
Cindex<-cc$Cindex
names(Cindex)<-cc$ID
barplot2(Cindex,plot.ci=TRUE,ci.l=c(cc$Cindex-cc$se),ci.u=c(cc$Cindex+cc$se),
         col=mypal,plot.grid = TRUE,border = "black",ci.lwd=3,grid.lwd	=3,
         cex.axis=1.5,xpd=FALSE,main = "TCGA")
tcga_compareC_p<-data.frame(Var=colnames(rt[,4:10]),pval=c(1:length(colnames(rt[,4:10]))))
for (i in colnames(rt[,4:10])) {
  p<-compareC(rt$OS.time, rt$OS, rs$TCGA$RS, rt[,i])$pval
  tcga_compareC_p[which(tcga_compareC_p$Var==i),2]<-p
}
write.table(tcga_compareC_p,"TCGA_cli_Cindex_p.txt",col.names = T,row.names = F,sep = "\t",quote = F)
###GSE57303
rt<-GSE57303_cli
GSE57303_rs<-list(RS=rs$GSE57303,Age=rt[,c(2:3,6)],Grade=rt[,c(2:4)],pM=rt[,c(2:3,9)],
              pN=rt[,c(2:3,8)],pT=rt[,c(2:3,7)],Gender=rt[,c(2:3,5)])

cc <- data.frame(Cindex=sapply(GSE57303_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[1])}),
                 se=sapply(GSE57303_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[2])}))%>%
  rownames_to_column('ID')
write.table(cc,"GSE57303_cli_Cindex.txt",col.names = T,row.names = F,sep = "\t",quote = F)
Cindex<-cc$Cindex
names(Cindex)<-cc$ID
barplot2(Cindex,plot.ci=TRUE,ci.l=c(cc$Cindex-cc$se),ci.u=c(cc$Cindex+cc$se),
         col=mypal,plot.grid = TRUE,border = "black",ci.lwd=3,grid.lwd	=3,
         cex.axis=1.5,xpd=FALSE,main = "GSE57303")
GSE57303_compareC_p<-data.frame(Var=colnames(rt[,4:9]),pval=c(1:length(colnames(rt[,4:9]))))
for (i in colnames(rt[,4:9])) {
  p<-compareC(rt$OS.time, rt$OS, rs$GSE57303$RS, rt[,i])$pval
  GSE57303_compareC_p[which(GSE57303_compareC_p$Var==i),2]<-p
}
write.table(GSE57303_compareC_p,"GSE57303_cli_Cindex_p.txt",col.names = T,row.names = F,sep = "\t",quote = F)
###GSE62254
rt<-GSE62254_cli

rt$Stage[rt$Stage==1]<-3
rt$Stage[rt$Stage==2]<-4
rt$T[rt$T==2]<-3
rt$N[rt$N==1]<-3
rt$N[rt$N==2]<-4

GSE62254_rs<-list(RS=rs$GSE62254,Age=rt[,c(2:3,7)],pM=rt[,c(2:3,10)],
                  pN=rt[,c(2:3,9)],pT=rt[,c(2:3,8)],Gender=rt[,c(2:3,6)],Stage=rt[,c(2:3,11)])

cc <- data.frame(Cindex=sapply(GSE62254_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[1])}),
                 se=sapply(GSE62254_rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~.,x))$concordance[2])}))%>%
  rownames_to_column('ID')
write.table(cc,"GSE62254_cli_Cindex.txt",col.names = T,row.names = F,sep = "\t",quote = F)
Cindex<-cc$Cindex
names(Cindex)<-cc$ID
barplot2(Cindex,plot.ci=TRUE,ci.l=c(cc$Cindex-cc$se),ci.u=c(cc$Cindex+cc$se),
         col=mypal,plot.grid = TRUE,border = "black",ci.lwd=3,grid.lwd	=3,
         cex.axis=1.5,xpd=FALSE,main = "GSE62254")
GSE62254_compareC_p<-data.frame(Var=colnames(rt[,6:11]),pval=c(1:length(colnames(rt[,6:11]))))
for (i in colnames(rt[,6:11])) {
  p<-compareC(rt$OS.time, rt$OS, rs$GSE62254$RS, rt[,i])$pval
  GSE62254_compareC_p[which(GSE62254_compareC_p$Var==i),2]<-p
}
write.table(GSE62254_compareC_p,"GSE62254_cli_Cindex_p.txt",col.names = T,row.names = F,sep = "\t",quote = F)

tcga_cli<-read.table("TCGA_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(tcga$sample,tcga_cli$sample)
tcga_cli<-tcga_cli[loc,]
tcga_cli$RS<-rs$TCGA$RS
GSE57303_cli<-read.table("GSE57303_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE57303$sample,GSE57303_cli$sample)
GSE57303_cli<-GSE57303_cli[loc,]
GSE57303_cli$RS<-rs$GSE57303$RS
GSE62254_cli<-read.table("GSE62254_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE62254$sample,GSE62254_cli$sample)
GSE62254_cli<-GSE62254_cli[loc,]
GSE62254_cli$RS<-rs$GSE62254$RS

library(ggplot2)
library(dplyr)
dat <-tcga_cli
dat$Age<-ifelse(dat$Age>65,">65","<=65")
dat$RS<-ifelse(dat$RS>median(dat$RS),"High","Low")
names(dat)[6:8]<-c("pM","pN","pT")
dat$Grade<-ifelse(dat$Grade=="G3","G3","G1+G2")
dat$pN<-ifelse(dat$pN=="N0","N0+N1",ifelse(dat$pN=="N1","N0+N1","N2+N3"))
dat$pT<-ifelse(dat$pT=="T1","T1+T2",ifelse(dat$pT=="T2","T1+T2","T3+T4"))
dat$Stage<-ifelse(dat$Stage=="stage i",'stage i+ii',
                  ifelse(dat$Stage=="stage ii","stage i+ii","stage iii+iv"))

fit <- coxph(Surv(OS.time,OS)~.,dat[,-1])
multiCoxSum=summary(fit)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab<-as.data.frame(outTab)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="TCGA_cli_multiCox.txt",sep="\t",row.names=F,quote=F)
library(forestplot)
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  data=as.matrix(rt)
  HR=data[,1:3]
  hr=sprintf("%.3f",HR[,"HR"])
  hrLow=sprintf("%.3f",HR[,"HR.95L"])
  hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  pVal=data[,"pvalue"]
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue",zero = "black")     
  tabletext <- list(c(NA, rownames(HR)),append("pvalue", pVal),append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")))  
  pdf(file=forestFile,width =10,height =6,onefile = FALSE)
  print(forestplot(tabletext, 
                   rbind(rep(NA, 3), HR),
                   col=clrs,
                   graphwidth=unit(50, "mm"),
                   xlog=T,
                   lwd.ci=4,
                   lwd.zero=2,
                   col.zero="black",
                   boxsize=0.3,
                   hrzl_lines=T,
                   xlab="Hazard ratio",
                   txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))))
  dev.off()
}

bioForest(coxFile="TCGA_aaa.txt",forestFile="TCGA_cli_multiCox.pdf",forestCol="#B24745FF")

dat <-GSE57303_cli
dat$Age<-ifelse(dat$Age>65,">65","<=65")
dat$RS<-ifelse(dat$RS>median(dat$RS),"High","Low")
names(dat)[7:9]<-c("pT","pN","pM")
dat$Grade<-ifelse(dat$Grade=="I","G1+G2",ifelse(dat$Grade=="II","G1+G2","G3+G4"))
dat$pN<-ifelse(dat$pN=="N0","N0+N1",ifelse(dat$pN=="N1","N0+N1","N2+N3"))
dat$pT<-ifelse(dat$pT=="T1","T1+T2",ifelse(dat$pT=="T2","T1+T2","T3+T4"))

fit <- coxph(Surv(OS.time,OS)~.,dat[,-1])
multiCoxSum=summary(fit)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab<-as.data.frame(outTab)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="GSE57303_cli_multiCox.txt",sep="\t",row.names=F,quote=F)

bioForest(coxFile="GSE57303_aaa.txt",forestFile="GSE57303_cli_multiCox.pdf",forestCol="#B24745FF")

dat <-GSE62254_cli
dat$Age<-ifelse(dat$Age>65,">65","<=65")
dat$RS<-ifelse(dat$RS>median(dat$RS),"High","Low")
names(dat)[c(6,8:10)]<-c("Gender","pT","pN","pM")
dat$pN<-ifelse(dat$pN=="N0","N0+N1",ifelse(dat$pN=="N1","N0+N1","N2+N3"))
dat$pT<-ifelse(dat$pT=="T1","T1+T2",ifelse(dat$pT=="T2","T1+T2","T3+T4"))
dat$Stage<-ifelse(dat$Stage=="I",'stage i+ii',
                  ifelse(dat$Stage=="II","stage i+ii","stage iii+iv"))

fit <- coxph(Surv(OS.time,OS)~.,dat[,-c(1,4,5)])
multiCoxSum=summary(fit)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab<-as.data.frame(outTab)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="GSE62254_cli_multiCox.txt",sep="\t",row.names=F,quote=F)

bioForest(coxFile="GSE62254_aaa.txt",forestFile="GSE62254_cli_multiCox.pdf",forestCol="#B24745FF")
###GSE62254中的DFS多因素
fit <- coxph(Surv(DFE.time,DFS)~.,dat[,-c(1:3)])
multiCoxSum=summary(fit)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab<-as.data.frame(outTab)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="GSE62254_cli_multiCox_DFS.txt",sep="\t",row.names=F,quote=F)

bioForest(coxFile="GSE62254_DFS_aaa.txt",forestFile="GSE62254_cli_multiCox_DFS.pdf",forestCol="#B24745FF")
library(survminer)
fit<-survfit(Surv(DFE.time,DFS)~RS,data =dat)
pdf("GSE62254_DFS_surv.pdf",width = 8,height = 6,onefile = T)
ggsurvplot(fit,data = dat, pval = TRUE,risk.table = T, 
           risk.table.col="strata",
           risk.table.fontsize=5,
           palette = c("#DC0000B2","#00A087B2"),
           legend.title="",
           pval.size=8,
           pval.method = T,
           surv.median.line = "hv",
           legend=c(0.7,0.9),
           legend.labs=c("High risk","low risk"),
           ggtheme=theme_survminer(font.x = c(18, "bold", "red"),
                                   font.y = c(18, "bold", "black"),
                                   #font.main = c(16, "bold", "black"),
                                   #font.submain = c(15, "bold", "black"),
                                   #font.caption = c(14, "bold", "blue"),
                                   font.tickslab = c(12, "bold", "black"),
                                   font.legend = c(20, "bold", "black")),
           tables.theme = theme_survminer(font.x = c(18, "bold", "red"),
                                          font.y = c(18, "bold", "black"),
                                          font.main = c(16, "bold", "black"),
                                          text=element_text(face = "bold"),
                                          #font.submain = c(15, "bold", "black"),
                                          #font.caption = c(14, "bold", "blue"),
                                          font.tickslab = c(12, "bold", "black"),
                                          font.legend = c(18, "bold", "black")),
           conf.int = F)+xlab("Times(years)")
dev.off()
write.table(tcga_cli,"TCGA_cli.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(GSE57303_cli,"GSE57303_cli.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(GSE62254_cli,"GSE62254_cli.txt",col.names = T,row.names = F,sep = "\t",quote = F)

tr<-read.table("TCGA_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(tcga$sample,tr$sample)
tr<-tr[loc,]
tr$RS<-rs$TCGA$RS
ROC.a <- timeROC(T=tr$OS.time, 
                 delta=tr$OS, marker=tr$RS,
                 other_markers=as.matrix(tr[,c("Age","Gender")]),
                 cause=1,
                 weighting="marginal",
                 times=quantile(tr$OS.time,probs=seq(0,1,0.1)),
                 iid=T)

tr<-read.table("GSE57303_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE57303$sample,tr$sample)
tr<-tr[loc,]
tr$RS<-rs$GSE57303$RS
ROC.b <- timeROC(T= tr$OS.time, delta= tr$OS, marker=tr$RS,
                 other_markers=as.matrix(tr[,c("M")]),
                 cause=1,weighting="marginal",
                 times=quantile(tr$OS.time,probs=seq(0,1,0.1)),
                 iid=TRUE)

tr<-read.table("GSE62254_cli_number.txt",header = T,sep = "\t",quote = "",check.names = F)
loc<-match(GSE62254$sample,tr$sample)
tr<-tr[loc,]
tr$RS<-rs$GSE62254$RS

ROC.c <- timeROC(T=tr$OS.time, delta= tr$OS, marker=tr$RS,
                 other_markers=as.matrix(tr[,c("Age","T","N","M")]),
                 cause=1,weighting="marginal",
                 times=quantile(tr$OS.time,probs=seq(0,1,0.1)),
                 iid=TRUE)
pdf("timeROC.pdf", 6, 5)
plotAUCcurve(ROC.a, conf.int=FALSE, col="#E64B35FF")
plotAUCcurve(ROC.b, conf.int=FALSE, col="#4DBBD5FF", add=TRUE)
plotAUCcurve(ROC.c, conf.int=FALSE, col="#00A087FF", add=TRUE)
title(main = "Time-dependent AUC")
# 图例设置
legend("topright", c("TCGA-adjust(Age+Gender)",
                     "GSE57303-adjust(M)","GSE62254-adjust(Age+T+N+M)"),
       col=c("#E64B35FF","#4DBBD5FF","#00A087FF"),
       bty='n', lty=1, lwd=2, cex=1)
dev.off()
###timeCindex拟合所有预后显著的变量
###TCGA
dat<-tcga_cli
cox1 <- coxph(Surv(OS.time,OS)~RS,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~Age,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Gender,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~RS+Age + Gender,data = dat,x=TRUE,y=TRUE)

eval.time <- seq(1,floor(max(dat$OS.time)),0.5) 
obj <- list("cox1"=cox1,
            "cox2"=cox2,
            "cox3"=cox3,
            "cox4"=cox4)

timeC <- pec::cindex(object = obj,
                     #formula=form,
                     formula=Surv(OS.time,OS)~.,
                     data=dat,
                     eval.times=eval.time, 
                     splitMethod = "BootCv")
timeC.mat <- do.call(cbind,timeC$BootCvCindex) 
ymin <- min(timeC.mat) 

mycol <- mypal

pdf("TCGA_timeCindex.pdf",width = 6,height = 5.5)
par(bty="l", 
    mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25) 

for (i in 1:ncol(timeC.mat)) { 
  if(i == 1){ 
    plot(eval.time,timeC.mat[,i],main = "TCGA Time-dependent C index",
         type="l",
         col = mycol[i],
         lwd = 2,
         ylim = c(ymin,1),xlim = range(dat$OS.time),
         xaxt = "n",
         xlab="Time (Years)",ylab = "Concordance index")
    axis(side = 1,
         at = seq(0,max(eval.time),1),
         labels = seq(0,max(eval.time),1))
  } else { 
    lines(eval.time,timeC.mat[,i],
          col = mycol[i],
          lwd = 2)
  }
}
if(ymin < 0.5) {abline(h = 0.5,lty = 4,col = "grey50",lwd = 2)}

legend("topright",
       legend = c("RS","Age","Gender","RS+Age+Gender"),
       col = mycol,
       lty = 1,
       lwd = 2,
       y.intersp = 1, x.intersp = 0.5, 
       bty = "o") 
invisible(dev.off()) 
###GSE57303
dat<-GSE57303_cli
cox1 <- coxph(Surv(OS.time,OS)~RS,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~M,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~RS+M,data = dat,x=TRUE,y=TRUE) 

eval.time <- seq(1,floor(max(dat$OS.time)),0.5) 
obj <- list("cox1"=cox1,
            "cox2"=cox2,
            "cox3"=cox3)

timeC <- pec::cindex(object = obj,
                     #formula=form,
                     formula=Surv(OS.time,OS)~.,
                     data=dat,
                     eval.times=eval.time, 
                     splitMethod = "BootCv") 
timeC.mat <- do.call(cbind,timeC$BootCvCindex) 
ymin <- min(timeC.mat) 

mycol <- mypal

pdf("GSE57303_timeCindex.pdf",width = 6,height = 5.5)
par(bty="l", 
    mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25) 

for (i in 1:ncol(timeC.mat)) { 
  if(i == 1){ 
    plot(eval.time,timeC.mat[,i],main = "GSE57303 Time-dependent C index",
         type="l",
         col = mycol[i],
         lwd = 2,
         ylim = c(ymin,1),xlim = range(dat$OS.time),
         xaxt = "n",
         xlab="Time (Years)",ylab = "Concordance index")
    axis(side = 1,
         at = seq(0,max(eval.time),1),
         labels = seq(0,max(eval.time),1))
  } else { 
    lines(eval.time,timeC.mat[,i],
          col = mycol[i],
          lwd = 2)
  }
}
if(ymin < 0.5) {abline(h = 0.5,lty = 4,col = "grey50",lwd = 2)} 

legend("topright",  
       legend = c("RS","M","RS+M"),
       col = mycol,
       lty = 1,
       lwd = 2,
       y.intersp = 1, x.intersp = 0.5, 
       bty = "o") 
invisible(dev.off()) 
###GSE62254
dat<-GSE62254_cli
dat$Stage[dat$Stage=="I"]<-"III"
dat$Stage[dat$Stage=="II"]<-"IV"
dat$T[dat$T=="T2"]<-"T3"
dat$N[dat$N=="N0"]<-"N2"
dat$N[dat$N=="N1"]<-"N3"

cox1 <- coxph(Surv(OS.time,OS)~RS,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~Age,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~T,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~N,data = dat,x=TRUE,y=TRUE) 
cox5 <- coxph(Surv(OS.time,OS)~M,data = dat,x=TRUE,y=TRUE) 
cox6 <- coxph(Surv(OS.time,OS)~RS+Age+T+N+M,data = dat,x=TRUE,y=TRUE) 

eval.time <- seq(1,floor(max(dat$OS.time)),0.5) 
obj <- list("cox1"=cox1,
            "cox2"=cox2,
            "cox3"=cox3,
            "cox4"=cox4,
            "cox5"=cox5,
            "cox6"=cox6)

timeC <- pec::cindex(object = obj,
                     #formula=form,
                     formula=Surv(OS.time,OS)~.,
                     data=dat,
                     eval.times=eval.time, 
                     splitMethod = "BootCv") 
timeC.mat <- do.call(cbind,timeC$BootCvCindex) 
ymin <- min(timeC.mat) 

mycol <- mypal

pdf("GSE62254_timeCindex.pdf",width = 6,height = 5.5)
par(bty="l", 
    mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25) 

for (i in 1:ncol(timeC.mat)) { 
  if(i == 1){ 
    plot(eval.time,timeC.mat[,i],main = "GSE62254 Time-dependent C index",
         type="l",
         col = mycol[i],
         lwd = 2,
         ylim = c(ymin,1),xlim = range(dat$OS.time),
         xaxt = "n",
         xlab="Time (Years)",ylab = "Concordance index")
    axis(side = 1,
         at = seq(0,max(eval.time),1),
         labels = seq(0,max(eval.time),1))
  } else { 
    lines(eval.time,timeC.mat[,i],
          col = mycol[i],
          lwd = 2)
  }
}
if(ymin < 0.5) {abline(h = 0.5,lty = 4,col = "grey50",lwd = 2)} 

legend("topright", 
       legend = c("RS","Age","T","N","M","RS+Age+T+N+M"),
       col = mycol,
       lty = 1,
       lwd = 2,
       y.intersp = 1, x.intersp = 0.5, 
       bty = "o") 
invisible(dev.off()) 


