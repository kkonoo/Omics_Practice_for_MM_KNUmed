if (!require("apeglm")) BiocManager::install("apeglm")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("pheatmap")) install.packages("pheatmap")


library(tidyverse)    # dplyr, ggplot2 등 포함
library(DESeq2)
library(apeglm)
library(ggrepel) # for labeling
library(pheatmap)

setwd("T:/homes/eha/교육/[학부] 4-1 분자의학 선택연구/2026/")


# 0. load DESeq data ####
dds <- readRDS('data/2_1_dds_DESeq2.rds')
dds@design
dds$condition

# 1. Extracting Results ####
res <- results(dds, contrast=c("condition", "Tumor", "Normal"))
head(res)
# 결과 구조 이해
# | 열            | 의미                                   |
# | baseMean      | 전체 샘플의 평균 발현량                |
# | log2FoldChange| log2(Tumor / Normal) = LFC             |
# | lfcSE         | LFC의 표준오차                         |
# | stat          | Wald test statistic                    |
# | pvalue        | 원래 p-value                           |
# | padj          | FDR 보정된 p-value (adjusted p-value)  |
summary(res)

# another option
dds@colData$condition
dds@colData$condition = relevel(dds@colData$condition, ref='Normal')
dds@colData$condition
dds@colData$condition = relevel(dds@colData$condition, ref='Tumor')

## 1-2. Independent Filtering ####
# x: Mean Normalized Counts, y: number of the significant
plot(metadata(res)$filterNumRej, 
     type="b", ylab="Number of rejections", xlab="Quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta, lty=3)
title("Independent Filtering Threshold Optimization")

res %>% as_tibble() %>% arrange(padj) %>% tail
res %>% as_tibble() %>% filter(is.na(padj)) %>% dim


## 1-3. check results ####
goi = c('TOP2A','SFTPA1','SFTPB','SFTPC','AGER')
res[goi,]

dds = readRDS('data/2_2_dds_DESeq2_simple.rds')
res2 <- results(dds, contrast = c("condition", "Tumor", "Normal"))
summary(res2)
res2[goi,]

head(res[order(res$padj), ]) # order by signifiance


# 2. Visualization ####
res_df <- res %>% as.data.frame() %>% mutate(gene_id=rownames(.)) %>% arrange(padj)
head(res_df)

## 2-1. MA plot ####
# X축: Mean Expression (평균 발현량)
# Y축: Log2 Fold Change
# → 전체 데이터의 분포와 bias 확인
plotMA(res, ylim=c(-5,5), main="Before Shrinkage")

### LFC Shrinkage to control low count genes ####
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Tumor_vs_Normal", type="apeglm")
summary(resLFC)
saveRDS(resLFC, file='data/3_dds_DESeq2_shr.rds')

res_df <- as.data.frame(resLFC) %>% 
  rownames_to_column("gene_id") %>%
  mutate(is_significant = padj < 0.05 & abs(log2FoldChange) > 1)

# check the shrinkage effect
par(mfrow=c(1,2))
plotMA(res, main="Before Shrinkage", ylim=c(-5,5))
plotMA(resLFC, main="After Shrinkage (Recommended)", ylim=c(-5,5))
par(mfrow=c(1,1))


## 2-2. volcano ####
p_thr=0.05
lfc_thr = 1
vol_data <- res_df %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > lfc_thr & padj < p_thr ~ "UP",
    log2FoldChange < -lfc_thr & padj < p_thr ~ "DOWN",
    TRUE ~ "NO"
  ))
head(vol_data)

ggplot(vol_data, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("DOWN"="steelblue", "NO"="grey", "UP"="indianred")) +
  geom_vline(xintercept=c(-lfc_thr, lfc_thr), linetype="dashed") +
  geom_hline(yintercept=-log10(p_thr), linetype="dashed") +
  # labeling the top 10
  geom_text_repel(data=head(arrange(vol_data, padj), 10), 
                  aes(label=gene_id), color="black") +
  theme_minimal() +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-value")

## using lfc (recommended)
vol_data2 <- resLFC %>% as.data.frame() %>% mutate(gene_id=rownames(.))  %>% arrange(padj) %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > lfc_thr & padj < p_thr ~ "UP",
    log2FoldChange < -lfc_thr & padj < p_thr ~ "DOWN",
    TRUE ~ "NO"
  ))
head(vol_data2)

ggplot(vol_data2, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("DOWN"="steelblue", "NO"="grey", "UP"="indianred")) +
  geom_vline(xintercept=c(-lfc_thr, lfc_thr), linetype="dashed") +
  geom_hline(yintercept=-log10(p_thr), linetype="dashed") +
  # labeling the top 10
  geom_text_repel(data=head(arrange(vol_data2, padj), 10), 
                  aes(label=gene_id), color="black") +
  theme_minimal() +
  labs(title="Volcano Plot (shr)", x="Log2 Fold Change", y="-Log10 Adjusted P-value")


## 2-3. heatmap ####
### 1) vst matrix ####
vsd <- readRDS('data/1_3_vsd.rds')

### 2) top 20 ####
top_genes <- res_df %>% 
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20) %>% 
  .$gene_id
head(top_genes)

### 3) vst matrix ####
mat <- assay(vsd)[top_genes, ]
head(mat)

## 4) annotations ####
df <- as.data.frame(colData(dds)[,c("gender","condition")])
head(df)

# Z-score scaling
pheatmap(mat, 
         scale="row",                    # z-score
         show_rownames=TRUE,             # label gene names
         show_colnames=FALSE,            # label sample names
         annotation_col=df,              # sample annotation  
         cluster_cols = TRUE,            # sample clustering
         main = "Top 20 DEGs Heatmap"    # plot title
         )


# FYI. FDR ####
## 1) example data (5 genes) ####
pvalues <- c(0.04, 0.001, 0.05, 0.0001, 0.2)
genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
df <- data.frame(Gene=genes, P_raw=pvalues)
df

## 2) sorting p-values ####
df <- df[order(df$P_raw), ] 
df

## 3) ranking * BH cutoff (m=5, alpha=0.05) ####
m <- 5
alpha <- 0.05
df$Rank <- 1:m
df$BH_Threshold <- (df$Rank / m) * alpha
df

## 4) pass or not ####
df$Pass <- df$P_raw <= df$BH_Threshold
df

## 5) padj ####
df$padj_calc <- df$P_raw * (m / df$Rank)
df
df$padj_calc <- pmin(df$padj_calc, 1)
df

# Cumulative Minimum from bottom
df$padj_final <- rev(cummin(rev(df$padj_calc)))
df

df$R_function <- p.adjust(df$P_raw, method="BH")
df #%>% select(Gene, P_raw, Rank, BH_Threshold, padj_final)
