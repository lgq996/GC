data<-read.table("unicox.txt",header = T,sep = "\t",quote = "",check.names = F)
rownames(data)<-data$sample
data<-as.data.frame(t(data[,-c(1:3)]))
cli<-read.table("TCGA_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
cli$Risk=ifelse(cli$RS>median(cli$RS),"High","Low")
library(limma)
group_list <- cli$Risk
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(data)
contrast.matrix <- makeContrasts(paste0(unique(group_list), collapse = '-'),
                                 levels = design)
fit <- lmFit(data, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)
DEG <- nrDEG  
DEG$symbol <- rownames(DEG)
sig<-DEG[which(DEG$P.Value < 0.05&abs(DEG$logFC)>1), ]
write.table(sig, 'DEG.txt',
            sep = '\t', quote = FALSE, col.names = T, row.names = FALSE)

tcga_gsva<-data[sig$symbol,]
tcga_gsva<-as.data.frame(t(tcga_gsva))
hub<-read.table("18.txt",header = F)
hub<-hub$V1
tcga_expr<-data[hub,]
genelist <- hub
gene1 <- genelist
immuscore <- function(gene1){
  y <- as.numeric(tcga_expr[gene1,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="pearson")
    data.frame(gene=gene1,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
data1 <- do.call(rbind,lapply(genelist,immuscore))
head(data1)
write.csv(data1, "correlation_86_18.csv", quote = F, row.names = F)
data1$pstar <- ifelse(data1$p.value < 0.05,
                      ifelse(data1$p.value < 0.01,"**","*"),
                      "")
data1$pstar[1:20]
library(ggplot2)
ggplot(data1, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#00AF66FF",mid = "white",high = "#CC0C00FF")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12,face = "bold"),
        axis.text.y = element_text(size =15,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("gene_cell_cor_86_18.pdf",width = 24,height = 8)

DEG$State = ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) >1, 
                   ifelse(DEG$logFC> 0.5 ,'Up','Down'),
                   'None')
library(ggplot2)
p <- ggplot(data = DEG, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=State)) +
  scale_color_manual(values=c("#00468BB2", "grey","#ED0000B2"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size =21,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 21,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 21,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 21,face = "bold"),
        plot.title = element_text(size = 21,face = "bold",hjust = 0.5),
        legend.title = element_text(size = 21,face =  "bold"),
        legend.text = element_text(size = 14,face = "bold"))

p
for_label <-DEG[which(DEG$symbol%in%sig$symbol),]
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    label.size = 1,
    data = for_label,
    color="black"
  )
ggsave("Volcano.pdf",width = 8,height = 6)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
cli<-rbind(cli[which(cli$Risk=="High"),],cli[which(cli$Risk=="Low"),])
loc<-match(cli$sample,colnames(data))
data<-data[,loc]
col_ha<-HeatmapAnnotation(which = "col",'RS'=rep(c("High","Low"),
                                                          times=c(148,148)),
                          annotation_name_gp = gpar(fontsize = 15,fontface = "bold"),
                          annotation_name_side = "left",
                          OSstat=ifelse(cli$OS==0,"Alive","Dead"),
                          OStime=cli$OS.time,
                          Gender=cli$Gender,
                          Stage=cli$Stage,
                          Age=cli$Age,
                          Grade=cli$Grade,
                          pM=cli$pathologic_M,
                          pN=cli$pathologic_N,
                          pT=cli$pathologic_T,
                          col = list(RS=c("High"="#CC0C00FF","Low"="#5C88DAFF")),
                          annotation_legend_param=list(RS=list(title="RS",
                                                                    title_position="topleft",
                                                                    title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                    labels_rot =0,
                                                                    legend_height=unit(1,"cm"),
                                                                    legend_width =unit(5,"mm"),
                                                                    labels_gp = gpar(fontsize = 15,
                                                                                     fontface = "bold")),
                                                       OSstat=list(title="OSstat",
                                                                   title_position="topleft",
                                                                   title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                   labels=c("Alive","Dead"),
                                                                   labels_rot =0,
                                                                   legend_height=unit(1,"cm"),
                                                                   legend_width =unit(5,"mm"),
                                                                   labels_gp = gpar(fontsize = 15,
                                                                                    fontface = "bold")),
                                                       OStime=list(title="OStime(years)",
                                                                   title_position="topleft",
                                                                   title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                   labels_rot =0,
                                                                   legend_height=unit(3,"cm"),
                                                                   legend_width =unit(5,"mm"),
                                                                   labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       Age=list(title="Age",
                                                                title_position="topleft",
                                                                title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                                labels_rot =0,
                                                                legend_height=unit(3,"cm"),
                                                                legend_width =unit(5,"mm"),
                                                                labels_gp = gpar(fontsize = 9,fontface = "bold")),
                                                       Grade=list(title="Grade",
                                                                title_position="topleft",
                                                                title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                labels_rot =0,
                                                                legend_height=unit(1.5,"cm"),
                                                                legend_width =unit(5,"mm"),
                                                                labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       pM=list(title="pM",
                                                                title_position="topleft",
                                                                title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                labels_rot =0,
                                                                legend_height=unit(1,"cm"),
                                                                legend_width =unit(5,"mm"),
                                                                labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       pN=list(title="pN",
                                                                title_position="topleft",
                                                                title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                labels_rot =0,
                                                                legend_height=unit(2,"cm"),
                                                                legend_width =unit(5,"mm"),
                                                                labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       pT=list(title="pT",
                                                                   title_position="topleft",
                                                                   title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                   labels_rot =0,
                                                                   legend_height=unit(2,"cm"),
                                                                   legend_width =unit(5,"mm"),
                                                                   labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       Gender=list(title="Gender",
                                                                  title_position="topleft",
                                                                  title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                  labels_rot =0,
                                                                  legend_height=unit(1,"cm"),
                                                                  legend_width =unit(5,"mm"),
                                                                  labels_gp = gpar(fontsize = 15,fontface = "bold")),
                                                       Stage=list(title="Stage",
                                                                title_position="topleft",
                                                                title_gp = gpar(fontsize = 21, fontface = "bold"),
                                                                labels_rot =0,
                                                                legend_height=unit(2,"cm"),
                                                                legend_width =unit(5,"mm"),
                                                                labels_gp = gpar(fontsize = 15,fontface = "bold"))
                          )
)
exp<-data[which(rownames(data)%in%sig$symbol),]
scale_exp <- apply(exp, 1, scale)
rownames(scale_exp) <- colnames(exp)
scale_exp <- t(scale_exp)
pdf("ComplexHeatmap.pdf",width = 8,height = 8)
heatmap<-Heatmap(scale_exp,name = " ",
                 heatmap_legend_param = list(title="",title_position="topleft",labels_rot =0,
                                             legend_height=unit(8,"cm"),
                                             legend_width =unit(5,"mm"),
                                             labels_gp = gpar(fontsize = 21,fontface = "bold")),
                 show_column_names = F,
                 show_row_names = T,
                 col = colorRamp2(c(-2,0,4),c("#00468BB2","white", "#ED0000B2")),
                 column_title ="",
                 column_split = c(rep("High",148),rep("Low",148)),
                 column_gap = unit(1, "mm"),
                 #row_title ="Intersect Gene",
                 column_title_side = "top",
                 row_title_side = "left",
                 row_title_rot = 90, 
                 column_title_gp = gpar(fontsize = 21, fontface = "bold",col = "black"), 
                 row_names_gp = gpar(fontsize = 3, fontface = "bold",col = "black"),
                 #row_title_gp = gpar(fontsize = 15, fontface = "bold",col = "black"),
                 cluster_columns =F,
                 cluster_rows = T,
                 column_order=c(colnames(scale_exp)),
                 show_row_dend = T,
                 top_annotation = col_ha
)
heatmap + rowAnnotation(link = anno_mark(at = which(rownames(exp)%in%sig$symbol), 
                                         labels = sig$symbol, labels_gp = gpar(fontsize = 8, fontface = "bold",col = "#ED0000B2")))

dev.off()
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05         
qvalueFilter=0.25        
rt<-sig
names(rt)[7]<-"SYMBOL"
m<-bitr(geneID = rt$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
rt<-merge(rt,m,by="SYMBOL")
gene=rt$ENTREZID
geneFC=rt$logFC
names(geneFC)=gene

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05,
            qvalueCutoff =1,  readable =T,ont = "BP")
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO-MF.txt",sep="\t",quote=F,row.names = F)
library(ggplot2)
library(tidyverse)

pdf(file="barplot_GOBP.pdf",width =10,height = 8)
ggplot(GO[1:5,],aes(Count,Description,color=pvalue))+
  geom_point(aes(size=Count))+
  scale_size_area(max_size =10)+
  scale_colour_gradient(low="#ED0000B2",high="#0099B4B2",
                        guide = guide_colorbar(reverse = TRUE))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+
  labs(x="GeneRatio")+
  theme_bw()+
  theme(axis.text.y = element_text(size=18,face = "bold"),
        axis.text.x = element_text(size=21,face = "bold"),
        axis.title = element_text(size = 21,face = "bold"),
        legend.title = element_text(size=21,face = "bold"),
        legend.text = element_text(size=15,face = "bold"))
dev.off()
R.utils::setOption("clusterProfiler.download.method",'auto') 

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx,layout = "star",colorEdge =F,color_category = "#0099B4B2",color_gene = "#ED0000B2",
         cex_category=2,cex_gene=1,cex_label_category=1,cex_label_gene=1)
ggsave("network_KEGG.pdf",width = 10,height = 8)

KEGG=as.data.frame(kkx)
KEGG=KEGG[(KEGG$pvalue<pvalueFilter),]
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

pdf(file="barplot_KEGG.pdf",width = 10,height =8)
ggplot(KEGG,aes(Count,Description,color=pvalue))+
  geom_point(aes(size=Count))+
  scale_size_area(max_size =10)+
  scale_colour_gradient(low="#ED0000B2",high="#0099B4B2",
                        guide = guide_colorbar(reverse = TRUE))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+
  labs(x="GeneRatio")+
  theme_bw()+
  theme(axis.text.y = element_text(size=18,face = "bold"),
        axis.text.x = element_text(size=21,face = "bold"),
        axis.title = element_text(size = 21,face = "bold"),
        legend.title = element_text(size=21,face = "bold"),
        legend.text = element_text(size=15,face = "bold"))
dev.off()

