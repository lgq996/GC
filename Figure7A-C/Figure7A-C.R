data<-data.table::fread("unicox.txt",header = T,sep = "\t",quote = "",check.names = F)
data<-as.data.frame(data)
rownames(data)<-data$sample
data<-data[,-c(1:3)]
data<-as.data.frame(t(data))
cli<-read.table("TCGA_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
gene<-c("SNHG5",
        "LINC01270",
        "CHKB-AS1",
        "NUTM2A-AS1",
        "MIR181A2HG",
        "CCNT2-AS1",
        "DLG3-AS1",
        "LINC01134",
        "NIFK-AS1",
        "RP11-443B7.1",
        "LSAMP-AS1",
        "HMGN3-AS1",
        "LPP-AS2",
        "RP11-710C12.1",
        "RP11-155O18.6",
        "CASC15",
        "RP11-449P15.2",
        "FLJ16779"
)
library(IOBR)
library(tibble)
library(tidyr)
library(stringr)
library(reshape2)
####28
ss<-read.table("ssGSEA_28cell.txt",header = T,sep = "\t",quote = "",check.names = F)
ss<-as.data.frame(t(ss))
ss$sample<-rownames(ss)
ss<-merge(cli[,c(1,11)],ss,by="sample")
rownames(ss)<-ss$sample
ss$sample<-NULL

boxplot<-ss
boxplot$RS<-ifelse(boxplot$RS>median(boxplot$RS),"High","Low")
melt <- reshape2:: melt(boxplot,
                        id.vars = c("RS"),
                        variable.name ="Signature",
                        value.name = "Signature_score")

box<- ggplot(melt,
             aes(x=Signature, y=Signature_score, 
                 fill = RS 
             )) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", 
               outlier.size = 0.65) +
  coord_flip()+
  scale_fill_manual(values= c("#00A087B2","#DC0000B2")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size =15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold")
  ) +
  xlab("")+ylab("Enrichment Score")+
  theme(legend.position = "top")

box + stat_compare_means(label = "p.signif",size=5)
ggsave("boxplot_28.pdf",width = 10,height = 8)
library(tidyverse)
library(corrplot)
library(RColorBrewer)
boxplot<-ss
corr <- cor(boxplot, method = "pearson")
p.corr<- cor.mtest(boxplot, conf.level = .95) 
pdf("corrplot_28.pdf",width = 8,height = 8)
corrplot(corr,title = "",
         p.mat = p.corr$p,
         sig.level = .05,
         insig = "pch",
         pch.cex = 1.5,
         pch.col = "red",
         #plotCI ="circle",
         #lowCI.mat =res1$lowCI,
         #uppCI.mat = res1$uppCI,
         col = brewer.pal(10,"RdYlBu"),
         #bg = "gold2",
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="original", 
         addrect = 3, 
         #mar = c(4,0,4,0), 
         rect.col = "#00A087B2", rect.lwd = 3, cl.pos = "r", 
         tl.offset=0.4,
         cl.offset = 0.4,
         cl.ratio = 0.2,
         tl.col = "black", tl.cex =1, cl.cex =1, tl.srt=60)
dev.off()
tcga_gsva<-ss
tcga_gsva<-tibble::column_to_rownames(tcga_gsva,"ID")
hub<-gene
tcga_expr<-data[which(rownames(data)%in%hub),]
genelist <- hub
gene1 <- genelist
immuscore <- function(gene1){
  y <- as.numeric(tcga_expr[gene1,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene1,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
data1 <- do.call(rbind,lapply(genelist,immuscore))
head(data1)
write.csv(data1, "correlation_28.csv", quote = F, row.names = F)
data1$pstar <- ifelse(data1$p.value < 0.05,
                      ifelse(data1$p.value < 0.01,"**","*"),
                      "")
data1$pstar[1:20]
library(ggplot2)
#data1$immune_cells<-stringr::str_replace_all(data1$immune_cells,"_CIBERSORT","")
ggplot(data1, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#00AF66FF",mid = "white",high = "#CC0C00FF")+
  geom_text(aes(label=pstar),col ="black",size = 10)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 15,face = "bold"),
        axis.text.y = element_text(size =20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("gene_cell_cor_28.pdf",width = 12,height = 8)








