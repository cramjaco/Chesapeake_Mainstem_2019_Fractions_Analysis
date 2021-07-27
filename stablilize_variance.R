counts <- read_tsv('ASVs_counts.tsv') %>% rename(ASV = "X1")
taxa0 <- read_tsv("ASVs_taxonomy.tsv") %>% rename(ASV = "X1")
library(DESeq2)

counts_mtx <- counts %>% column_to_rownames("ASV") %>% as.matrix()

deseq_pre <- DESeqDataSetFromMatrix(counts_mtx, sample, ~ID)
dds = deseq_pre[rowSums(counts(deseq_pre)) > 5,]

dds_esf <- estimateSizeFactors(dds, type = "poscounts")
dds01 <- DESeq(dds_esf)
dds_res <- results(dds01)

dds_counts <- counts(dds01, normalized = TRUE)

## why aren't the flagged samples removed? Dear god.