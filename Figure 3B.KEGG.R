library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       
qvalueFilter=0.05       


colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("C:\\Users")          
rt=read.table("diffgenes.txt", header=T, sep="\t", check.names=F)     


genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        


kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$id[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)


showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}


pdf(file="barplot.pdf", width=10, height=15)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()


pdf(file="bubble.pdf", width = 10, height = 15)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()