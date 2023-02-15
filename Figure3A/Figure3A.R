library(survival)
library(forestplot)
library(stringr)
est_dd<-read.table("est_dd.txt",header = T,sep = "\t",quote = "",check.names = F)
fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = "forward")
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
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
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
  pdf(file=forestFile,width =8,height =12,onefile = FALSE)
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

bioForest(coxFile="multiCox.txt",forestFile="multiCox.pdf",forestCol="#B24745FF")
