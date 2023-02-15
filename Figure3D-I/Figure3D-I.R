library(survivalROC)
rt<-read.table("rt.txt",header = T,sep = "\t",quote = "",check.names = F)
rt<-rs$TCGA
nobs <- NROW(rt)
roc1=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =1, span = 0.25*nobs^(-0.20))
roc2=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =3, span = 0.25*nobs^(-0.20))
roc3=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =5, span = 0.25*nobs^(-0.20))


plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="#E64B35B2",
     xlab="False positive rate", ylab="True positive rate",
     main=paste("TCGA ROC"),smooth=T,
     lwd = 4, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

abline(0,1)
lines(roc2$FP,roc2$TP,type = "l",xlim=c(0,1), ylim=c(0,1),col="#4DBBD5B2",lwd = 4)
lines(roc3$FP,roc3$TP,type = "l",col="#00A087B2",xlim=c(0,1), ylim=c(0,1),lwd = 4)
roc1=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =1, method = "KM")
roc2=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =3, method = "KM")
roc3=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$RS, 
                 predict.time =5, method = "KM")

legend(0.5,0.4,c(paste("1year:AUC=",round(roc1$AUC,3)),
                 paste("3year:AUC =",round(roc2$AUC,3)),
                 paste("5year:AUC =",round(roc3$AUC,3))
),
x.intersp = 1,y.intersp = 0.8,
lty = 1,lwd = 2,col = c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#7A142C","#5D90BA","#431A3D"),
bty="n",seg.len=1,cex=1.5)
###生存分析
library(survminer)
rt$risk<-ifelse(rt$RS>median(rt$RS),"high","low")
fit<-survfit(Surv(OS.time,OS)~risk,data =rt)
pdf("TCGA_surv.pdf",width = 8,height = 6,onefile = T)
ggsurvplot(fit,data = rt, pval = TRUE,risk.table = T, #线的类型
           risk.table.col="strata",
           risk.table.fontsize=5,
           palette = c("#DC0000B2","#00A087B2"),#线的颜色
           legend.title="",
           pval.size=8,
           pval.method = T,
           surv.median.line = "hv",
           legend=c(0.7,0.9),
           legend.labs=c("High risk","low risk"),
           #title    = "Survival curves", #主标题
           #subtitle = "Based on Kaplan-Meier estimates", #副标题
           #caption  = "created with survminer",#图底部添加说明
           ggtheme=theme_survminer(font.x = c(18, "bold", "red"),
                                   font.y = c(18, "bold", "black"),
                                   #font.main = c(16, "bold", "black"),
                                   #font.submain = c(15, "bold", "black"),
                                   #font.caption = c(14, "bold", "blue"),
                                   font.tickslab = c(12, "bold", "black"),#坐标轴
                                   font.legend = c(20, "bold", "black")),#图例
           tables.theme = theme_survminer(font.x = c(18, "bold", "red"),
                                          font.y = c(18, "bold", "black"),
                                          font.main = c(16, "bold", "black"),
                                          text=element_text(face = "bold"),
                                          #font.submain = c(15, "bold", "black"),
                                          #font.caption = c(14, "bold", "blue"),
                                          font.tickslab = c(12, "bold", "black"),#坐标轴
                                          font.legend = c(18, "bold", "black")),#图例
           conf.int = F)+xlab("Times(years)")
dev.off()
