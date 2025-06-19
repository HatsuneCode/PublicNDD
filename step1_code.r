## tar -xvf GSE156928_RAW.tar
## for i in *.tar.gz; do tar -xvf $i; done

## 1. read tables
fs = list.files(pattern = 'quant.genes.sf', full.names = T, recursive = T)
expr = do.call(cbind, lapply(fs, function(f) {
  exp = read.table(f, sep = '\t', header = T)[, c('Name', 'NumReads')]
  setNames(data.frame(exp$NumReads, row.names = exp$Name), unlist(strsplit(f, '/'))[2] )
}))
rownames(expr) = sub('\\..*', '', rownames(expr))

## 2. id map
ID  = read.table('IDtoSym.GRCh38-2024-A.txt', sep = '\t')
Sym = ID$V2[match(rownames(expr), ID$V1)] 
rownames(expr) = ifelse(is.na(Sym), rownames(expr), paste0(rownames(expr), '_', Sym))
expr = RNAseq.checkDupRow(expr)
expr = round(expr)
write.table(cbind(Gene = rownames(expr), expr), 'Raw.counts.txt', sep = '\t', row.names = F)

## 3. normalize
expr = RNAseq.Normalize(expr)
write.table(cbind(Gene = rownames(expr), expr), 'Normalized.counts.txt', sep = '\t', row.names = F)

###PCA
library("ggrepel")
library(ggplot2)
colnames(expr)=sapply(colnames(expr), function(x) {paste0(strsplit(x, "_")[[1]][1],"_",strsplit(x, "_")[[1]][3])})
expr.f = expr[order(apply(expr, 1, mad), decreasing = T)[1:3e3], ]
pca = PCA(expr.f)
pca$group = ifelse(grepl("PD",colnames(expr)),"PD","N")
pca$group=factor(pca$group,levels=c("PD","N"))
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 4, max.overlaps = 50) +
  labs(title = "PCA of Top3k Variable Genes",
       x = 'PC1', y = 'PC2', color = 'Group') +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    axis.title = element_text(size = 16), axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave('PCA.Top3k.png', width = 9, height = 7, dpi = 300)

# 差异分析
raw_expr=read.table("Raw.counts.txt",header=T,sep="\t",row.names = 1)
colnames(raw_expr)=sapply(colnames(raw_expr), function(x) {paste0(strsplit(x, "_")[[1]][1],"_",strsplit(x, "_")[[1]][3])})
pos=colnames(raw_expr)[grep("PD",colnames(raw_expr))]
neg = colnames(raw_expr)[-grep("PD",colnames(raw_expr))]
DEG = RNAseq.DESeq2(raw_expr, pos, neg)
write.table(DEG, 'DEG.DvsCtrl.txt', sep = '\t', row.names = F)
