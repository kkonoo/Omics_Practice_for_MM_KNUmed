if (!require("tidyverse")) install.packages("tidyverse")
if (!require("edgeR")) BiocManager::install("edgeR")
if (!require("DESeq2")) BiocManager::install("DESeq2")
# 라이브러리 로드
library(tidyverse)
library(edgeR)
library(DESeq2)


set.seed(123)
setwd("T:/homes/eha/교육/[학부] 4-1 분자의학 선택연구/2026/")


####
# 1. creating example data ####
## 1-1. fake 1000 genes, 6 samples ####
n_genes <- 1000
n_samples <- 6
sample_names <- c(paste0("Ctrl_", 1:3), paste0("Treat_", 1:3))
sample_names
gene_names <- paste0("Gene_", 1:n_genes)
head(gene_names)

## 1-2. creating Raw Count Matrix (Negative Binomial distribution) ####
rnorm(10,mean=0,sd=1) # making mock data
x = rnorm(1000,mean=0,sd=1)
mean(x); sd(x)
rbinom(); rpois(); rbeta()

raw_counts <- matrix(rnbinom(n_genes * n_samples, # total 1000*6 values
                             mu = 100, size = 5), # mean: 100, dispersion: 5
                     nrow = n_genes # making 1000 genes
                     )
head(raw_counts)

### varying Library Size(Sequencing depth) ####
raw_counts[, 4:6] <- raw_counts[, 4:6] * 3
head(raw_counts)

colnames(raw_counts) <- sample_names
rownames(raw_counts) <- gene_names
head(raw_counts)

## 1-3. Gene Length in bp (for RPKM, TPM) ####
gene_lengths <- round(runif(n_genes, min = 200, max = 5000))
head(gene_lengths)
length(gene_lengths)


####
# 2. CPM (counts per million) ####
## library size: Total Counts per sample 
lib_sizes <- colSums(raw_counts)
lib_sizes

## calculating CPM
prop_counts <- sweep(raw_counts, 
                     MARGIN = 2, # 1: by row(gene), 2: by column(sample)
                     STATS = lib_sizes, 
                     FUN = "/"
                     )
head(raw_counts)
head(prop_counts)
lib_sizes
79/99904
152/99904

cpm_counts <- prop_counts * 1e6
head(cpm_counts)


####
# 3. RPKM (reads per killobase million) ####
## gene length (kb)
kb_lengths <- gene_lengths / 1000
head(kb_lengths)

## calculating RPKM
## 3-1. correcting for gene length (RPK) ####
rpk <- sweep(raw_counts, 
             MARGIN = 1, # 1: by row(gene)
             STATS = kb_lengths, 
             FUN = "/"
)
head(rpk)
## 3-2. correcting for Library Size (million-based) (RPKM) ####
rpkm_counts <- sweep(rpk, 
                     MARGIN = 2, # 1: by column(sample)
                     STATS = lib_sizes, 
                     FUN = "/"
                     )*1e6
head(rpkm_counts)


####
# 4. TPM (transcripts per million) ####
## 4-1. RPK (Reads Per Kilobase) ####
rpk <- sweep(raw_counts, 1, kb_lengths, "/")

## 4-2. per-sample RPK ####
rpk_sums <- colSums(rpk)
head(rpk_sums)

## 4-3. calculating TPM (correcting for per-sample RPK) ####
tpm_counts <- sweep(rpk, 
                    MARGIN = 2, 
                    STATS = rpk_sums, # instead of library size
                    FUN = "/")*1e6

## check: TPM sum = 1M
colSums(tpm_counts)


####
# 5. TMM from edgeR ####
## 5-1. make edgeR object (DGEList) ####
dge <- DGEList(counts = raw_counts)
dge
str(dge)

## 5-2. make Normalization Factor from TMM ####
dge2 <- calcNormFactors(dge, method = "TMM")
str(dge2)
dge2@.Data[[2]] %>% head
dge@.Data[[2]] %>% head

dge2 <- calcNormFactors(dge, method = "TMM", 
                        refColumn=3,
                        logratioTrim=0.4,
                        sumTrim=0
                          )
dge2@.Data[[2]] %>% head

## 5-3. normalizing counts (CPM, but using TMM factor) ####
tmm_counts <- cpm(dge2, normalized.lib.sizes = TRUE)
head(tmm_counts)

cpm_counts2 <- cpm(dge2, normalized.lib.sizes = FALSE)
head(cpm_counts2)
head(cpm_counts)


####
# 6. RLE from DESeq2 ####
## 6-1. metadata for DESeq2 object ####
colData <- data.frame(condition = factor(rep(c("Ctrl", "Treat"), each = 3)))
colData
rownames(colData) <- sample_names
colData
head(raw_counts)

## 6-2. make DESeqDataSet (dds) ####
dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = colData, 
                              design = ~ condition)
dds
str(dds)

## 6-3. estimating Size Factor ####
dds2 <- estimateSizeFactors(dds)
dds2 # creating sizeFactor in colData
dds2@colData

## 6-4. normalized counts ####
deseq_counts <- counts(dds2, normalized = TRUE)
head(deseq_counts)
raw_counts2 <- counts(dds2, normalized = FALSE)
head(raw_counts2)


# 7. visualization ####
## 7-1. reshaping counts (Long format for ggplot) ####
transform_to_long <- function(count_matrix, method_name) {
  as.data.frame(count_matrix) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Count") %>%
    mutate(Method = method_name,
           # Log2 transformation (pseudocount: +1 for 0)
           LogCount = log2(Count + 1))
}

## applying the function ####
transform_to_long(raw_counts, "1. Raw") %>% head
transform_to_long(cpm_counts, "2. CPM") %>% head

df_all <- rbind(
  transform_to_long(raw_counts, "1. Raw"),
  transform_to_long(cpm_counts, "2. CPM"),
  transform_to_long(rpkm_counts, "3. RPKM"),
  transform_to_long(tpm_counts, "4. TPM"),
  transform_to_long(tmm_counts, "5. TMM (edgeR)"),
  transform_to_long(deseq_counts, "6. DESeq2")
)
head(df_all)
tail(df_all)

## 7-2. visualization (Boxplot) ####
df_all %>%
ggplot(aes(x = Sample, y = LogCount, fill = Sample)) +
  # box plot
  geom_boxplot() + 
  # plot by method
  facet_wrap(~Method, scales = "free_y") + 
  # theme
  theme_bw() + 
  # x-axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # plot title
  labs(title = "Normalization Method Comparison",
       y = "Log2(Counts + 1)",
       x = "")

# TMM / DESeq2: less highly expressed genes 



