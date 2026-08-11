if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("ReactomePA"))     BiocManager::install("ReactomePA")
#if (!require("msigdbr")) BiocManager::install("msigdbr")
if (!require("enrichplot"))     BiocManager::install("enrichplot")
# if you have any issues, please run R as administrator

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)  # 시각화 전용 패키지
# library(org.Hs.eg.db) # Gene ID 변환이 필요한 경우 사용

setwd("T:/homes/eha/교육/[학부] 4-1 분자의학 선택연구/2026/")



# 1. load results (with entrez_id) ####
res_df <- readRDS('data/4_1_res_DESeq2_with_entrez_and_rank.rds')
head(res_df)

## 1-1. input: ranked gene list ####
cor(res_df$rank_lfc, res_df$rank_pi)
cor(res_df$rank_lfc, res_df$rank_wald)
cor(res_df$rank_pi, res_df$rank_wald)

### (optional) removing duplicates
dim(res_df)
res_df <- res_df %>% distinct(entrez, .keep_all = TRUE) %>% filter(!is.na(entrez))
dim(res_df)

## 1-2. ranked gene list ####
geneList <- res_df$rank_wald
names(geneList) <- res_df$entrez 
geneList <- sort(geneList, decreasing = TRUE)

head(geneList) # top up-regulated
tail(geneList) # top down-regulated

names(geneList) <- res_df$gene 
geneList <- sort(geneList, decreasing = TRUE)

head(geneList) # top up-regulated
tail(geneList) # top down-regulated


# 2. GSEA ####
## Gene ontology (GO) ####
gsea_res <- gseGO(geneList = geneList,
                  OrgDb        = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont          = "BP",            
                  minGSSize    = 10,              # excluding too small pathways
                  maxGSSize    = 500,             # excluding too large pathways
                  pvalueCutoff = 1,               # pvalue cutoff to show terms
                  by = "fgsea" )                  # fgsea: pvalues by tail approximate (fast), DOSE: permutation
head(gsea_res)
tail(gsea_res)
saveRDS(gsea_res,'data/4_2_gsea_GO.rds')

## REACTOME ####
gsea_res <- gsePathway(geneList, 
                pvalueCutoff = 1)
head(gsea_res)
saveRDS(gsea_res,'data/4_2_gsea_RE.rds')


# 3. visualization ####
gsea_res <- readRDS('data/4_2_gsea_GO.rds')
## 3-1. Dot Plot (summary) ####
# x: GeneRatio, dot size: Gene Count, color: p.adjust
dotplot(gsea_res, showCategory = 5, split = ".sign") + 
  ggtitle("GO BP Enrichment")

dotplot(gsea_res, showCategory = 5, split = ".sign") + 
  facet_grid(.~.sign) + # Up/Down Pathway 
  ggtitle("GO BP Enrichment (split)")


## 3-2. GSEA Enrichment Plot (for the specific pathway) ####
p1 <- gseaplot2(gsea_res, geneSetID = 1, title = gsea_res$Description[1])
print(p1)

### multiple (Top 3)
p2 <- gseaplot2(gsea_res, geneSetID = 1:3)
print(p2)

# 3-3. Ridge Plot (expression level), requring ggridges
ridgeplot(gsea_res, showCategory = 5) + labs(x = "Enrichment Distribution")

# 3-4. Network plot
cnetplot(gsea_res, foldChange=geneList)


# 4. Core Enrichment Genes ####
gsea_res_df <- as.data.frame(gsea_res)
head(gsea_res_df)

# save the result
write.csv(gsea_res_df, "data/4_2_gsea_GO_results.csv")

# show core genes for a pathway
core_genes <- strsplit(gsea_res_df$core_enrichment[1], "/")[[1]]
print(core_genes)


# 5. 자유 탐색 아이디어: 관심 유전자 탐색 ####
# 특정 유전자가 종양에서 어떻게 변하는지 확인
my_gene <- "EGFR"  # 원하는 유전자로 바꿔보세요
my_gene <- "ARHGAP44"  # 원하는 유전자로 바꿔보세요

# 발현량 시각화
plotCounts(dds, gene = my_gene, intgroup = "condition", returnData = FALSE)

# 더 예쁘게
plot_data <- plotCounts(dds, gene = my_gene, intgroup = "condition", returnData = TRUE)
plot_data$count %>% head
plot_data$count %>% max
counts(dds,normalized=T)[my_gene,] %>% max
log2(counts(dds,normalized=T)[my_gene,] %>% max)

ggplot(plot_data, aes(x = condition, y = log2(count+1), fill = condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = paste0(my_gene, " Expression"),
       y = "log2(Normalized Counts)") +
  theme_minimal()


