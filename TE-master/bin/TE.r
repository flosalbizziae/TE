#cite: Jain, A, Tuteja, G. (2018) TissueEnrich: Tissue-specific gene enrichment analysis. Bioinformatics, bty890. 10.1093/bioinformatics/bty890
args<-commandArgs(T) #从命令行传参

type<-args[1]
#传入参数：
#1. 选用的工作种类
#2. 基因列表
#3. 参数
#4. 生成的特定标签

if(type=="type=1"){
  #teEnrichment
  genelist<-args[2]
  par<-args[3]
  inp1="TRUE"
  inp2=1
  inputGenes=unlist(strsplit(genelist,","))#genelist即参数args[2]格式为“g1,g2,g3,...,gn"
  library("GSEABase")
  gs<-GeneSet(geneIds=inputGenes,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  library("TissueEnrich")
  parameters=substr(par,5,nchar(par))
  if(nchar(parameters)==0){
    output<-teEnrichment(gs)
  }else{
    tmplist=strsplit(parameters,",")
    for( i in 1:length(tmplist[[1]])){
      p=strsplit(tmplist[[1]][i],"=")
      if(p[[1]][1]=="BH"){
        inp1=as.logical(p[[1]][2])
      }
      if(p[[1]][1]=="db"){
        inp2=as.numeric(p[[1]][2])
      }
    }
    output=teEnrichment(gs,multiHypoCorrection=inp1,rnaSeqDataset=inp2)
  }
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]),colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  write.csv(enrichmentOutput,paste(args[4],"_out.csv",sep=''))
  
  #the barplot of the -log10(pvalue) of the tissue
  library(ggplot2)
  #png("barplot.png")
  p<-ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
    geom_bar(stat = 'identity')+
    labs(x='', y = '-LOG10(P-Adjusted)')+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
    ggsave(p,filename=paste(args[4],"_barplot.png",sep=''))
  #dev.off()
  
  #the heatmap of the expression profiles of tissue-specific genes using
  library(tidyr)
  #png("heatmap.png")
  if(length(which(enrichmentOutput[,1]==max(enrichmentOutput[,1])))==1){
  seExp<-output[[2]][[row.names(enrichmentOutput[which(enrichmentOutput[,1]==max(enrichmentOutput[,1])),])]]
  exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
  exp$Gene<-row.names(exp)
  exp<-exp %>% gather(key = "Tissue", value = "expression",1:(ncol(exp)-1))
  
  p<-ggplot(exp, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                             colour = "white") + scale_fill_gradient(low = "white",
                                                                                     high = "steelblue")+
    labs(x='', y = '')+
    theme_bw()+
    guides(fill = guide_legend(title = "Log2(TPM)"))+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
	ggsave(p,filename=paste(args[4],"_heatmap.png",sep=''))
	}
  #dev.off()
}else if(type=="type=2"){
  #teGeneRetrieval
  #teEnrichmentCustom
#传入参数：
#1. 选用的工作种类
#2. 基因列表文件
#3. 组织特异表达列表
#4. 参数
#5. 生成的特定标签
  genelist<-args[2]
  tg<-args[3]
  par<-args[4]
  inp1=5
  inp2=1
  inp3=7
  inp4=1
  in5=as.logical("TRUE")
  inputGenes=unlist(strsplit(genelist,","))#genelist即参数args[2]格式为“g1,g2,g3,...,gn"
  library("GSEABase")
  gs<-GeneSet(geneIds=inputGenes,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  expressionData<-read.table(tg,header=TRUE,row.names=1,sep='\t')
  library("SummarizedExperiment")
  se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
  library("TissueEnrich")
  parameters=substr(par,5,nchar(par))
  if(nchar(parameters)==0){
    output1<-teGeneRetrieval(se)
    output2<-teEnrichmentCustom(gs,output1)
  }else{
    tmplist=strsplit(parameters,",")
    for( i in 1:length(tmplist[[1]])){
      p=strsplit(tmplist[[1]][i],"=")
      if(p[[1]][1]=="fc"){
        inp1=as.numeric(p[[1]][2])
      }
      if(p[[1]][1]=="minexp"){
        inp2=as.numeric(p[[1]][2])
      }
      if(p[[1]][1]=="maxnt"){
        inp3=as.numeric(p[[1]][2])
      }
      if(p[[1]][1]=="tg"){
        inp4=as.numeric(p[[1]][2])
      }
      if(p[[1]][1]=="BH"){
        inp5=as.logical(p[[1]][2])
      }
    }
    output1<-teGeneRetrieval(se,foldChangeThreshold=inp1,expressedGeneThreshold=inp2,maxNumberOfTissues=inp3)
    output2<-teEnrichmentCustom(gs,output1,tissueSpecificGeneType=inp4,multiHypoCorrection=inp5)
  }
  write.csv(assay(output1),paste(args[5],"_tissue-specific-genes.csv",sep=''))
  
  enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),
                                        row.names = rowData(output2[[1]])[,1]),
                             colData(output2[[1]])[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  write.csv(enrichmentOutput,paste(args[5],"_out.csv",sep=''))
 
  library(ggplot2)
  #png("barplot.png")
  p<-ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
    geom_bar(stat = 'identity')+
    labs(x='', y = '-LOG10(P-Adjusted)')+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
    ggsave(p,filename=paste(args[5],"_barplot.png",sep=''))
  #dev.off()
  
  #the heatmap of the expression profiles of tissue-specific genes using
  library(tidyr)
  #png("heatmap.png")
  if(length(which(enrichmentOutput[,1]==max(enrichmentOutput[,1])))==1){
  seExp<-output2[[2]][[row.names(enrichmentOutput[which(enrichmentOutput[,1]==max(enrichmentOutput[,1])),])]]
  exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
  exp$Gene<-row.names(exp)
  exp<-exp %>% gather(key = "Tissue", value = "expression",1:(ncol(exp)-1))
  
  p<-ggplot(exp, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                             colour = "white") + scale_fill_gradient(low = "white",
                                                                                     high = "steelblue")+
    labs(x='', y = '')+
    theme_bw()+
    guides(fill = guide_legend(title = "Log2(TPM)"))+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
	ggsave(p,filename=paste(args[5],"_heatmap.png",sep=''))
	}
  #dev.off()
  
}else{
  cat("The type is not recognizable!")
}

