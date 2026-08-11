library(tidyverse)    # dplyr, ggplot2 등 포함
library(DESeq2)

setwd("T:/homes/eha/교육/[학부] 4-1 분자의학 선택연구/2026/")

# 왜 EDA를 하는가? ####
# - 데이터 구조 확인
# - 이상 샘플(outlier) 탐지
# - Batch effect 존재 여부 확인
# - 기대한 대로 그룹이 분리되는지 시각적 확인


# 1. read data (count table and metadata) ####
countData <- readRDS('data/0_raw_counts_final.rds')
dim(countData)
head(countData[,1:5])

colData <- readRDS("data/0_meta_data_final.rds")
head(colData)
table(colData$tissue_type)

# check sample ids of the count table and metadata
all(colnames(countData) == colData$barcode) # should be TRUE

# set the Factor level (condition)
# default: alphabetical ... ref = Control/Normal
colData$condition <- factor(colData$tissue_type, levels = c("Normal", "Tumor")) 
colData$condition
# case / control ... treatment / control ... KO / control



# 2. DESeq2 object (dds) ####
# countData, colData, design
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ patient + condition) # DEG by condition
colData$patient = gsub('-','_',colData$patient)
dds
saveRDS(dds, file='data/1_1_dds.rds')

# 3. Pre-filtering (de-noising) ####
dim(counts(dds))

head(counts(dds)[,1:4])
head(counts(dds)[,1:4]>=1)
head(rowSums(counts(dds)[,1:4]>=1))
head(rowSums(counts(dds)>=1))
rowSums(counts(dds)) %>% head

## 3-1. lowly expressed genes (for speed-up & statistical power) ####
### option 1: simple (Total Count < 10?) ####
keep <- rowSums(counts(dds)) >= 10
head(keep)
cat(paste("   Removed", sum(!keep), "low count genes.\n"))
cat(paste("   Remaining genes:", sum(keep), "\n"))
dds_filtered <- dds[keep, ]

### option 2: considring group size (recommended) ####
# Total Count < 10 in the smallest group?
# ex: Control(3), Case(10) -> non-zero counts in at least 3 samples?
main_condition = 'condition'
smallestGroupSize <- min(table(colData[[main_condition]]))
smallestGroupSize*0.4
keep <- rowSums(counts(dds) >= 1) >= (smallestGroupSize*0.5)
cat(paste("   Removed", sum(!keep), "low count genes.\n"))
cat(paste("   Remaining genes:", sum(keep), "\n"))
dds_filtered <- dds[keep, ]

## 3-2. confounding genes ####
all_genes <- rownames(dds_filtered)
head(all_genes)
# Mitochodiral: starting with "^MT-" 
mt_genes <- grep("^MT-", all_genes, value = TRUE)
mt_genes
# Ribosomal: starting with "^RPS|^RPL"
ribo_genes <- grep("^RPS|^RPL", all_genes, value = TRUE)
ribo_genes
# Immunoglobulin: starting with "^IG[HKL]" : "IGH", "IGK", "IGL"
ig_genes <- grep("^IG[HKL]", all_genes, value = TRUE)
ig_genes
# filtering gene list
contam_list <- c(mt_genes, ribo_genes, ig_genes)
length(contam_list)
# !: 'NOT'
keep_indices <- !(all_genes %in% contam_list)
dds_filtered <- dds_filtered[keep_indices, ]
dds_filtered

saveRDS(dds_filtered, file='data/1_2_dds_filtered.rds')


# 4. Variance Stabilizing Transformation (vst) ####
## 4-1. vst (normalization and transformation) ####
vsd <- vst(dds_filtered, blind=FALSE) 
# blind = FALSE: using design (recommended), blind=TRUE: using only technical noises
saveRDS(vsd, file='data/1_3_vsd.rds')

# check the result visually
vsd_data = assay(vsd)
head(vsd_data[,1:4]) # similar to log-scale
df_vst <- data.frame(mean = rowMeans(vsd_data),
                     sd = apply(vsd_data, 1, sd)
                     )
head(df_vst)

ggplot(df_vst, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.1) +
  geom_smooth(color = "red") + # trend line
  theme_minimal() +
  labs(title = "Variance Stabilization Check (VST)",
       x = "Mean Expression", y = "Standard Deviation")

# no transformation
count_data = assay(dds_filtered)
df_count <- data.frame(mean = rowMeans(count_data),
                       sd = apply(count_data, 1, sd)
)
ggplot(df_count, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.1) +
  geom_smooth(color = "red") + #  trend line
  theme_minimal() +
  labs(title = "Variance Stabilization Check (no transf)",
       x = "Mean Expression", y = "Standard Deviation")


## 4-2. rlog (regularized log) ####
rld <- rlog(dds_filtered, blind=FALSE)



# 5. PCA ####
# ✅ 해석 포인트:
# 1) Tumor와 Normal이 잘 분리되는가? → 생물학적 신호가 있다!
# 2) 성별이나 나이 같은 변수에 의해 분리되는가? → Batch effect?
# 3) 이상하게 동떨어진 샘플은 없는가? → Outlier 확인

vsd_data = assay(vsd)

## 5-1. PCA ####
pca_res <- prcomp(t(vsd_data))

pca_res$x %>% head

# PC1 and PC2
pca_coords <- as.data.frame(pca_res$x)[,1:2]
pca_coords %>% head

ggplot(pca_coords, aes(x=PC1, y=PC2)) + 
  geom_point(size=3) + ggtitle("PCA Result") + theme_bw()


## 5-2. using DESeq2 function (using top 100) ####
plotPCA(vsd, intgroup = "condition") +
  theme_classic() +
  labs(title = "PCA Plot")

plotPCA(vsd, ntop=nrow(vsd), intgroup = "condition") +
  theme_classic() +
  labs(title = "PCA Plot")

# by other covariates
plotPCA(vsd, intgroup = "race") +
  theme_classic() +
  labs(title = "PCA Plot")

plotPCA(vsd, intgroup = "gender") +
  theme_classic() +
  labs(title = "PCA Plot")

plotPCA(vsd, intgroup = "age_at_index") +
  theme_classic() +
  labs(title = "PCA Plot")


# better visulation by ggplot2
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - VST data") +
  theme_bw()
