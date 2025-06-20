## tar -xvf GSE156928_RAW.tar
## for i in *.tar.gz; do tar -xvf $i; done

## 1. Read tables
fs = list.files(pattern = 'quant.genes.sf', full.names = T, recursive = T)
expr = do.call(cbind, lapply(fs, function(f) {
  exp = read.table(f, sep = '\t', header = T)[, c('Name', 'NumReads')]
  n   = paste(unlist(strsplit(unlist(strsplit(f, '/'))[2], '_'))[2:3], collapse = '_')
  setNames(data.frame(exp$NumReads, row.names = exp$Name), n)
}))
rownames(expr) = sub('\\..*', '', rownames(expr))

## 2. ID map
ID  = read.table('IDtoSym.GRCh38-2024-A.txt', sep = '\t')
Sym = ID$V2[match(rownames(expr), ID$V1)] 
rownames(expr) = ifelse(is.na(Sym), rownames(expr), paste0(rownames(expr), '_', Sym))
expr = RNAseq.checkDupRow(expr)
expr = round(expr)
write.table(cbind(Gene = rownames(expr), expr), 'Raw.counts.txt', sep = '\t', row.names = F)

## 3. Normalize
expr = RNAseq.Normalize(expr)
write.table(cbind(Gene = rownames(expr), expr), 'Normalized.counts.txt', sep = '\t', row.names = F)

## 4. PCA
library(ggrepel)
library(ggplot2)
expr = read.table('Normalized.counts.txt', sep = '\t', header = T, row.names = 1, check.names = F)
expr.f = expr[order(apply(expr, 1, mad), decreasing = T)[1:3e3], ]
pca    = PCA(expr.f)
groups = c('PD', 'Ctrl')
pca$group = ifelse(grepl('PD', colnames(expr)), groups[1], groups[2])
pca$group = factor(pca$group, groups)
p = ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = group), size = 4) +
  geom_text_repel(aes(label = sample), size = 4) +
  labs(title = 'PCA of Top3k Variable Genes', x = 'PC1', y = 'PC2', color = 'Group') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), size = 18, face = 'bold',
        legend.text = element_text(size = 14));p
ggsave('PCA.Top3k.png', p, w = 9, h = 7)

## 5. DEG
raw = read.table('Raw.counts.txt', sep = '\t', header = T, row.names = 1, check.names = F)
pos = grep('PD', colnames(raw), value = T)
neg = setdiff(colnames(raw), pos)
DEG = RNAseq.DESeq2(raw, pos, neg, name = 'GSE156928.PD')
write.table(DEG, 'DEG.PDvsCtrl.txt', sep = '\t', row.names = F)

## 6. GSEA
#### command
expr = read.table('Normalized.counts.txt', sep = '\t', header = T, row.names = 1, check.names = F)
expr.p = expr[,  grepl('PD', colnames(expr))]  
expr.n = expr[, !grepl('PD', colnames(expr))]
GSEA.prepare(expr.p, expr.n, 'PD', 'Ctrl', 
             cmd.out = '/mnt/e/Project/NDD-Database/hsa/GSE156928')
#### results
res  = list.files('GSEA', 'gsea_report_for.*.tsv', full.names = T, recursive = T)
GESA = do.call(rbind, lapply(res, function(f) {
  gsea = read.table(f, sep = '\t', header = T)
  gsea = data.frame(Terms = gsea$NAME, NES = gsea$NES, Pvalue = gsea$NOM.p.val)
  gsea$Sig = ifelse(gsea$Pvalue < .05, 'Yes', 'No')
  gsea
}))
write.table(GSEA, 'GSEA.results.txt', sep = '\t', row.names = F)
