INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

mONCTSG_Gene<-read.csv(file.path(INPUT_PATH, "mONCTSG_Gene.csv")) 
mONCTSG_MutFree<-read.csv(file.path(INPUT_PATH, "mONCTSG_MutFree.csv")) 
mONCTSG_Pathway<-read.csv(file.path(INPUT_PATH, "mONCTSG_Pathway.csv")) 

if (file.exists(file.path(OUTPUT_PATH, "mONCTSG_Gene_Result.txt"))) {file.remove(file.path(OUTPUT_PATH, "mONCTSG_Gene_Result.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "mONCTSG_MutFree_Result.txt"))) {file.remove(file.path(OUTPUT_PATH, "mONCTSG_MutFree_Result.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.txt"))) {file.remove(file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "mOT_gene_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "mOT_gene_q.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "mOT_mutfree_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "mOT_mutfree_q.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "mOT_pathway_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "mOT_pathway_q.csv"))}

for(i in 5:ncol(mONCTSG_Gene)) {
  mwut<-wilcox.test(mONCTSG_Gene[,i]~mONCTSG_Gene[,4])
  mwup<-mwut$p.value
  ttest<-t.test(mONCTSG_Gene[,i]~mONCTSG_Gene[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "mONCTSG_Gene_Result.txt"), append=TRUE)
}

for(i in 5:ncol(mONCTSG_MutFree)) {
  mwut<-wilcox.test(mONCTSG_MutFree[,i]~mONCTSG_MutFree[,4])
  mwup<-mwut$p.value
  ttest<-t.test(mONCTSG_MutFree[,i]~mONCTSG_MutFree[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "mONCTSG_MutFree_Result.txt"), append=TRUE)
}

for(i in 5:ncol(mONCTSG_Pathway)) {
  mwut<-wilcox.test(mONCTSG_Pathway[,i]~mONCTSG_Pathway[,4])
  mwup<-mwut$p.value
  ttest<-t.test(mONCTSG_Pathway[,i]~mONCTSG_Pathway[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.txt"), append=TRUE)
}

mOT_p<-read.csv(file.path(INPUT_PATH, "mONCTSG_p.csv"))

mOT_gene_q<-p.adjust(mOT_p$gene, method="BH")
mOT_mutfree_q<-p.adjust(mOT_p$mutfree, method="BH")
mOT_pathway_q<-p.adjust(mOT_p$pathway, method="BH")
write.csv(mOT_gene_q,file.path(OUTPUT_PATH, "mOT_gene_q.csv"))
write.csv(mOT_mutfree_q,file.path(OUTPUT_PATH, "mOT_mutfree_q.csv"))
write.csv(mOT_pathway_q,file.path(OUTPUT_PATH, "mOT_pathway_q.csv"))

mOT_gene2<-read.csv(file.path(INPUT_PATH, "mONCTSG_Gene_Result2.csv"))
mOT_pathway<-read.csv(file.path(INPUT_PATH, "mONCTSG_Pathway_Result.csv"))

ggplot(mOT_gene2,aes(x=mean_2DS,y=mean_3DO))+geom_point(shape=21, alpha=0.5, aes(fill=oncTSG,colour=oncTSG,size=log10P))+scale_color_manual(values=c("#931751","#005492"))+scale_fill_manual(values=c("#931751","#005492"))+scale_alpha_manual(values=c(0.6,0.3,0.3))+scale_size(range=c(4,16))+theme(aspect.ratio=1)+xlim(-2.2,0.4)+ylim(-2.2,0.4)

ggplot(mOT_pathway,aes(x=mean_2DS,y=mean_3DO))+geom_point(alpha=0.6, aes(shape=significance, colour=oncTSG,fill=oncTSG,size=log10P))+scale_shape_manual(values=c(1,21))+scale_color_manual(values=c("#931751","#005492"))+scale_fill_manual(values=c("#931751","#005492"))+scale_size(range=c(3,30))+theme(aspect.ratio=1)+xlim(-0.6,0.2)+ylim(-0.6,0.2)

