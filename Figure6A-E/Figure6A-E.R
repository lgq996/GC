data<-read.table("unicox.txt",header = T,sep = "\t",quote = "",check.names = F)
rownames(data)<-data$sample
data<-as.data.frame(t(data[,-c(1:3)]))
cli<-read.table("TCGA_cli.txt",header = T,sep = "\t",quote = "",check.names = F)
cli$Risk=ifelse(cli$RS>median(cli$RS),"High","Low")
library(IOBR)
library(ggsci)
mypal<-pal_npg("nrc")(10)

sig<-calculate_sig_score(pdata= NULL,
                             eset= data,
                             signature= signature_collection,
                             method= "ssgsea",
                             mini_gene_count = 2)
sig$Index<-NULL
names(sig)[1]<-"sample"
wt<-merge(cli[,c(1,11)],sig,by="sample")
wt$RS<-ifelse(wt$RS>median(wt$RS),"High","Low")

out<-data.frame(sig=colnames(wt)[3:ncol(wt)],p=rep(1,length(3:ncol(wt))))
for (i in colnames(wt)[3:ncol(wt)]) {
  p=wilcox.test(get(i) ~ RS, data = wt, var.equal = TRUE)$p.value
  out[which(out$sig==i),2]<-p
}
out<-out[out$p<0.05,]
names(signature_tme)
names(signature_metabolism)
names(signature_tumor)
boxplot<-wt[,c("RS",intersect(out$sig,names(signature_tumor)))]
melt <- reshape2:: melt(boxplot,
                        id.vars = c("RS"),
                        variable.name ="Signature",
                        value.name = "Signature_score")

box<- ggplot(melt,
             aes(x=Signature, y=Signature_score, 
                 fill = RS 
             )) + ɫ
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier???ú?ɫ
               outlier.size = 0.65) +
  #coord_flip()+
  scale_fill_manual(values= c("#E64B35FF","#4DBBD5FF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold")
  ) +
  xlab("Signature tumor")+ylab("Enrichment Score")+
  theme(legend.position = "top")

box + stat_compare_means(label = "p.signif",size=5)
ggsave("boxplot_tumor.pdf",width = 6,height = 8)
library(cBioPortalData)
cbio <- cBioPortal()
studies = getStudies(cbio)
head(studies$studyId)
id = "stad_tcga"
clinical = clinicalData(cbio, id)
colnames(clinical)

df = na.omit(clinical[,c("patientId","DFS_MONTHS","DFS_STATUS")])
names(df)<-c("sample","DFS.time","DFS")
df$sample<-paste(df$sample,"-01A",sep="")
df$DFS.time<-as.numeric(df$DFS.time)
df$DFS.time<-df$DFS.time/12
df$DFS<-ifelse(df$DFS=="0:DiseaseFree",0,1)
write.table(df,"TCGA_DFS_cli.txt",col.names = T,row.names = F,sep = "\t",quote = F)

df = na.omit(clinical[,c("patientId","TMB_NONSYNONYMOUS")])
names(df)<-c("sample","TMB")
df$sample<-paste(df$sample,"-01A",sep="")
df<-merge(df,cli[,c(1:3,11)],by="sample",all=F)
cor.test(df$TMB,df$RS,method = "pearson")
boxplot<-df[,c(5,2)]
boxplot$TMB<-as.numeric(boxplot$TMB)
boxplot$TMB<-log2(boxplot$TMB+1)
library(ggpubr)
p<-ggplot(boxplot,aes(x=Risk, y=TMB)) + 
  geom_violin(aes(fill = Risk),trim = F,alpha=0.3,width=0.8)+
  geom_boxplot(aes(fill=Risk),notch = T,width=0.5) +
  geom_jitter(aes(fill = Risk),position = position_jitter(0.1),shape=21, size = 4,alpha=0.9)+
  scale_fill_manual(values= c("#E64B35FF","#4DBBD5FF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 90, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 15,face =  "bold"),
        legend.text = element_text(size = 15,face = "bold")
  ) +
  xlab("")+ylab("log2(TMB+1)")+guides(fill=F)+
  theme(legend.position = "top")

p + stat_compare_means(size=7,label.y =max(boxplot[,2]))
ggsave(paste0("boxplot_",i,".pdf"),width =8,height = 8)

library(maftools)
tmp<-read.table("TCGA-STAD.mutect2_snv.tsv",header = T,sep = "\t",quote = "",check.names = F)
tmp<-tmp[which(tmp$Sample_ID%in%cli$sample),]

tmp<-tmp[which(tmp$Sample_ID%in%cli[which(cli$Risk=="High"),"sample"]),]
tmp<-tmp[which(tmp$Sample_ID%in%cli[which(cli$Risk=="Low"),"sample"]),]

colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )
tcga.brca = read.maf(maf = tmp,
                     vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = names(tail(sort(table(tmp$Variant_Classification ))))

oncoplot(maf = tcga.brca, top =20,colors = col) 
tmp<-read.table("TCGA-STAD.mutect2_snv.tsv",header = T,sep = "\t",quote = "",check.names = F)
tmp<-tmp[which(tmp$Sample_ID%in%cli$sample),]

tmp<-tmp[,c(1,2,9)]
tmp$effect<-ifelse(tmp$effect=="stop_gained",1,
                   ifelse(tmp$effect=="3_prime_UTR_variant",2,
                          ifelse(tmp$effect=="frameshift_variant",3,
                                 ifelse(tmp$effect=="synonymous_variant",4,
                                        ifelse(tmp$effect=="missense_variant",
                                               5,6)))))
tmp$effect<-1
melt <- reshape2::acast(tmp,Sample_ID~gene,max)
melt<-as.data.frame(melt)
melt[melt=="-Inf"]<-0
melt$sample<-rownames(melt)
melt<-merge(cli[,c(1,11)],melt,by="sample")
tcga_gsva<-melt[,3:ncol(melt)]
tcga_expr<-melt[,1:2]
gene1 <- c("RS")
immuscore <- function(gene1){
  y <- as.numeric(tcga_expr[,gene1])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene1,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
data1 <- do.call(rbind,lapply(gene1,immuscore))
head(data1)
write.csv(data1, "correlation_RS_genemut.csv", quote = F, row.names = F)
