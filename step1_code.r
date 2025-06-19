setwd("E:\\SYSU\\SLH\\GSE156928_RAW")
file_list=read.table("file_list.txt")

files <- list.files(path = "./", 
                    pattern = "quant.genes.sf", 
                    full.names = TRUE, 
                    recursive = TRUE)
matrix_A=NULL
for(i in 1:length(files)){
  data=read.table(files[i],sep="\t",header=T)
  matrix_A=cbind(matrix_A,as.integer(data$NumReads))
}
colnames(matrix_A)=file_list$V1
rownames(matrix_A)=data$Name


gene_I_S_list=read.delim("../IDtoSym.GRCh38-2024-A.txt",header=F,sep="\t")

ENSG_ID=sub("\\.\\d+$", "",data$Name)
loc=match(ENSG_ID,gene_I_S_list$V1)
ENSG_ID1=ENSG_ID
for(i in 1:length(ENSG_ID1)){
  if(!is.na(loc[i])){
  ENSG_ID1[i]=gene_I_S_list[loc[i],2]
    }
  
}

rownames(matrix_A)=ENSG_ID1
normalize_exp=RNAseq.Normalize(matrix_A,log2=T,method = "DESeq2")

