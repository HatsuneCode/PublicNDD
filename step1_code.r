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
