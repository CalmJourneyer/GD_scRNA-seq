library(Seurat)
library(tidyverse)

datalist<-list()

label<-list('GD-B-10'='1','GD-B-11'='2','GD-B-13'='4','GD-B-14'='5','GD-B-15'='6','GD-B-17'='7','GD-B-19'='8','GD-B-7'='11','GD-B-9'='12','p02n'='13','p04n'='14','p08n'='15','p11n'='16','p12n'='17','p15n'='18','p16n'='19')
for(i in names(label)){
vdj<-read.csv(sprintf('/data3/home/congjia/scRNAseq/newest_version/TCR/%s/outs/filtered_contig_annotations.csv',i))
vdj[['rename_barcode']]<-paste0(vdj[['barcode']],'_',label[i])
vdj[['filename']]<-i
datalist[[i]]<-vdj
}
vdj <- do.call(rbind, datalist)

vdj <- vdj %>% filter(high_confidence =="true" & 
                  chain %in% c("TRA","TRB") &
                  productive =="true")

vdj_a <- vdj %>% filter(chain =="TRA") %>% arrange(desc(umis), desc(reads)) %>% distinct(rename_barcode,.keep_all = TRUE)
vdj_b <- vdj %>% filter(chain =="TRB") %>% arrange(desc(umis), desc(reads)) %>% distinct(rename_barcode,.keep_all = TRUE)

final_vdj = full_join(x = vdj_a, y=vdj_b, by = c("rename_barcode"), suffix = c(".TRA",".TRB"))

subdata_T<-read_rds('/data3/home/congjia/scRNAseq/newest_version/result/T.rds.gz')
subdata_T@meta.data <- subdata_T@meta.data %>% 
  mutate(rename_barcode = rownames(subdata_T@meta.data)) %>% 
  inner_join(final_vdj, by = "rename_barcode")

subdata_T@meta.data$Clone_AA = paste(subdata_T@meta.data$cdr3.TRA, subdata_T@meta.data$cdr3.TRB, sep="_")

subT= subset(subdata_T@meta.data,productive.TRA == "true" & productive.TRB == "true" & subgroup!='γδT')
tmp = subT%>% group_by(Clone_AA) %>%
  summarize(Clone_NUM = n()) %>%
  mutate(Clone_ID = paste0("Clone_",rownames(.)))
subT = merge.data.frame(subT, tmp)
names<-c('CD4 naïve','Th1-like','Th17','Tfh','Tph','Treg')
subT$stype <- ifelse(subT$subgroup%in%names,"CD4","CD8")

filter_subT<-subT %>% select(rename_barcode,Clone_AA,Clone_ID,samplename,stype,subgroup,condition)
write.table(filter_subT,file='./result/clonaltype.txt',sep='\t',quote=FALSE,row.names = FALSE)
