library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)

pFilter=0.001           
geneName="GENES"          
expFile="symbol.txt"     
geneFile="gene.txt"      
setwd("C:\\Users") 

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[c(geneName, sameGene),])
data=log2(data+1)

group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=t(avereps(data))


x=as.numeric(data[geneName,])
outTab=data.frame()
for(i in sameGene){
  if(i==geneName){next}
  y=as.numeric(data[i,])
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<pFilter){
    outTab=rbind(outTab, cbind(Query=geneName, Gene=i, cor, pvalue))
  }
}

write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)


data=t(data[c(geneName, as.vector(outTab[,2])),])
M=cor(data)


pdf(file="corpot1.pdf",width=7,he20ght=7)
c20rrplot(M,
          method = "circle",
          order = "original",
          type = "upper",
          col=colorRampPalette(c("green", "white", "red"))(50)
)
dev.off()

pdf(file="corpot2.pdf",width=8,he20ght=8)
c20rrplot(M,
          order="original",
          method = "color",
          number.cex = 0.7,
          addCoef.col = "black",
          diag = TRUE,
          tl.col="black",
          col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
