library(Seurat)
library(tidyverse)

PatientList <- grep("/filtered_feature_bc_matrix$",list.dirs(),value=T)
scRNAlist<-list()

for(i in 1:length(PatientList)){
          counts<-Read10X(data.dir=PatientList[i])
          samplename=strsplit(PatientList[i],'/')[[1]][3]
          scRNAlist[[i]]<-CreateSeuratObject(counts,project=samplename,min.features=200,min.cells=3)
          scRNAlist[[i]]@meta.data$samplename<-samplename
          cat(samplename,dim(counts),dim(scRNAlist[[i]]),sep='\t',end='\n')
          }

for(i in 1:length(scRNAlist)){
           sc<-scRNAlist[[i]]
           sc[["mt_percent"]]<-PercentageFeatureSet(sc,pattern="^MT-")
           sc[["rb_percent"]]<-PercentageFeatureSet(sc,pattern="^RP[SL]")
           HB_genes <-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
           HB_m<-match(HB_genes,row.names(sc@assays$RNA))
           HB_genes<-rownames(sc@assays$RNA)[HB_m]
           HB_genes<-HB_genes[!is.na(HB_genes)]
           sc[["hb_percent"]]<-PercentageFeatureSet(sc,features=HB_genes)
           count_q=quantile(sc$nCount_RNA, 0.98)
           sc=subset(sc,subset = nFeature_RNA>500&nFeature_RNA<8000&mt_percent<15&rb_percent<40&hb_percent<3&nCount_RNA<count_q)
           print(dim(sc))
           scRNAlist[[i]]<-sc
}

scdata=merge(x=scRNAlist[[1]],y=scRNAlist[-1])

counts <- GetAssayData(scdata, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"))),]
filter<-grep("^MT-|^RP[SL]", rownames(counts))
counts1<-rownames(counts)[-filter]
counts2 <- counts[which(rownames(counts) %in% counts1),]
scdata <- subset(scdata, features = rownames(counts2))
