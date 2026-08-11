########################################################################
# 분자의학 선택실습: 오믹스 데이터 분석 기술
# Day 2 - DEG 분석, 시각화, Pathway 분석
# 담당교수: 하은지 (eha@knu.ac.kr)
########################################################################

# Day 1 결과 불러오기
# load("data/day1_results.RData")

library(tidyverse)
library(DESeq2)

# =====================================================================
# Session 1: DEG 분석 이론 + DESeq2 실습 (10:00 ~ 11:00)
# =====================================================================

## 1-1. DEG 분석이란? ####
# DEG = Differentially Expressed Genes (차등 발현 유전자)
# "Tumor에서 Normal 대비 유의하게 발현이 변한 유전자 찾기"
# 
# 핵심 개념: Variance Decomposition (분산 분해)
# Total Variance = Condition + Batch + Residual
#                 (신호)     (노이즈) (랜덤 오차)
#
# 우리의 목표: Condition의 효과를 정확하게 추정하기

## 1-2. DESeq2의 통계 모델 ####
# (1) RNA-seq 데이터의 특성
#   - Count data: 이산적 정수 (0, 1, 2, ...)
#   - Overdispersion: 분산 >> 평균
#   → 일반 t-test / linear regression으로는 부적절!
#
# (2) Negative Binomial (음이항) 분포 사용
#   Var(μ) = μ + α·μ²
#   - μ: 평균 발현량
#   - α: dispersion (생물학적 변동)
#
# (3) Generalized Linear Model (GLM)
#   log(μ) = β₀ + β₁·condition + β₂·batch + ...
#   - β₁이 바로 log2 Fold Change!

## 1-3. Design Formula 설정 ####
# Simple model: ~ condition
# Paired model: ~ patient + condition (환자별 변동을 보정)
#               → 더 정밀한 분석!

# ▶ 우리 데이터: paired tumor-normal → patient 변수 포함
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData   = metadata,
  design    = ~ patient + condition  # Paired model
)
dds$condition <- relevel(dds$condition, ref = "normal")  # 기준: normal

## 1-4. DESeq2 실행! ####
# 이 한 줄이 모든 것을 해줍니다:
# ① size factor 계산 (normalization)
# ② dispersion 추정
# ③ GLM fitting
# ④ Wald test (통계 검정)

dds <- DESeq(dds)

# 시간이 좀 걸립니다... (약 1-3분)

## 1-5. 결과 추출 ####
res <- results(dds)
res

# 결과 구조 이해
# | 열            | 의미                                    |
# |---------------|----------------------------------------|
# | baseMean      | 전체 샘플의 평균 발현량                    |
# | log2FoldChange| log2(Tumor / Normal) = LFC             |
# | lfcSE         | LFC의 표준오차                           |
# | stat          | Wald test statistic                     |
# | pvalue        | 원래 p-value                            |
# | padj          | FDR 보정된 p-value (adjusted p-value)   |

summary(res)

## 1-6. Multiple Testing Correction ####
# 문제: 20,000개 유전자를 검정 → p < 0.05로 해도 1,000개가 우연히 유의!
# 해결: FDR (False Discovery Rate) 보정
#   - Bonferroni: α/n ... 너무 엄격
#   - BH (Benjamini-Hochberg): FDR 조절 → padj 사용!
#   - FDR = FP / (FP + TP)

# 유의한 DEG 수
sum(res$padj < 0.05, na.rm = TRUE)
sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)

## 1-7. 특정 유전자 확인 ####
# 폐 관련 Gold Standard 유전자
genes_of_interest <- c("TOP2A", "SFTPA1", "SFTPB", "SFTPC", "AGER")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df %>%
  filter(gene %in% genes_of_interest) %>%
  select(gene, baseMean, log2FoldChange, padj) %>%
  arrange(padj)

# TOP2A: 세포 분열 관련 → Tumor에서 UP
# SFTPA1/B/C: 폐 surfactant → Tumor에서 DOWN
# AGER: 폐 상피세포 마커 → Tumor에서 DOWN

# ✅ 여기까지 되셨으면 손 들어주세요!


# =====================================================================
# Session 2: 결과 시각화 (11:00 ~ 12:00)
# =====================================================================

## 2-1. MA Plot ####
# X축: Mean Expression (평균 발현량)
# Y축: Log2 Fold Change
# → 전체 데이터의 분포와 bias 확인

plotMA(res, main = "MA Plot: Tumor vs. Normal", ylim = c(-5, 5))

# LFC shrinkage: 발현량 낮은 유전자의 LFC 과대추정 보정
res_shrunk <- lfcShrink(dds, coef = "condition_tumor_vs_normal", type = "apeglm")
plotMA(res_shrunk, main = "MA Plot (shrunken LFC)", ylim = c(-5, 5))

## 2-2. Volcano Plot ####
# X축: Log2 Fold Change
# Y축: -Log10(adjusted p-value)
# → 어떤 유전자가 크게 변하면서 통계적으로도 유의한지 한눈에!

res_df <- res_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"  # Not Significant
    )
  )

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "Log2 Fold Change", y = "-Log10(adjusted p-value)",
       title = "Volcano Plot: Tumor vs. Normal",
       color = "Significance") +
  theme_minimal() +
  xlim(c(-10, 10))

# Gold standard 유전자 표시
library(ggrepel)

highlight_df <- res_df %>%
  filter(gene %in% genes_of_interest)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey70")) +
  geom_point(data = highlight_df, color = "purple", size = 3) +
  geom_text_repel(data = highlight_df, aes(label = gene),
                  color = "purple", fontface = "bold", size = 4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "Log2 Fold Change", y = "-Log10(adjusted p-value)",
       title = "Volcano Plot with Key Lung Genes",
       color = "Significance") +
  theme_minimal() +
  xlim(c(-10, 10))

## 2-3. Heatmap ####
# 상위 DEG들의 발현 패턴을 heatmap으로 시각화

library(pheatmap)

# 상위 50개 DEG 선택 (padj 기준)
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

# rlog 변환된 데이터에서 추출
rld <- rlog(dds, blind = FALSE)
mat <- assay(rld)[top_genes, ]

# Z-score 변환 (행 기준: 유전자별 평균 0, 표준편차 1)
mat_z <- t(scale(t(mat)))

# annotation 준비
ann_col <- data.frame(
  Condition = metadata$condition,
  row.names = colnames(dds)
)

pheatmap(mat_z,
         annotation_col = ann_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         main = "Top 50 DEGs: Tumor vs. Normal",
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100))


# =====================================================================
# 점심시간 (12:00 ~ 13:30)
# =====================================================================


# =====================================================================
# Session 3: Over-Representation Analysis (ORA) (13:30 ~ 14:30)
# =====================================================================

## 3-1. Functional Enrichment 개념 ####
# DEG 리스트를 얻었다 → "이 유전자들이 어떤 기능/경로에 관여하는가?"
#
# 주요 데이터베이스 (사전):
# - Gene Ontology (GO)
#   · BP (Biological Process): 무슨 일을 하는가? (e.g., cell division)
#   · MF (Molecular Function): 분자적 활성 (e.g., kinase activity)
#   · CC (Cellular Component): 어디에 있는가? (e.g., nucleus)
# - KEGG Pathway: 대사/신호전달 경로
# - Reactome: peer-reviewed, manually curated pathways

## 3-2. ORA (Over-Representation Analysis) 원리 ####
# "내 DEG 리스트에 특정 pathway의 유전자가 기대 이상으로 많은가?"
# → Fisher's exact test (2x2 분할표)
#
#              | In Pathway | Not in Pathway |
# | DEG        |     a      |       b        |
# | Background |     c      |       d        |
#
# ⚠️ ORA의 한계:
# - 유전자를 DEG / non-DEG 이분법으로 나눔 → threshold 의존적!
# - padj = 0.049 vs. 0.051 → 하나는 포함, 하나는 제외

## 3-3. clusterProfiler 실습 ####
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db"))   BiocManager::install("org.Hs.eg.db")
if (!require("enrichplot"))     BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# (1) DEG 리스트 준비
deg_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

length(deg_genes)  # 몇 개의 DEG?

# (2) Gene ID 변환 (Gene Symbol → Entrez ID)
# GO/KEGG는 Entrez ID를 사용합니다
gene_ids <- bitr(deg_genes,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
head(gene_ids)

# Background genes (전체 발현 유전자)
bg_genes <- bitr(rownames(filtered_counts),
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

# (3) GO Enrichment (Biological Process)
ego <- enrichGO(
  gene          = gene_ids$ENTREZID,
  universe      = bg_genes$ENTREZID,  # background!
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",               # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE                # Entrez → Symbol 변환
)

head(ego, 10)

# (4) GO 결과 시각화
# Dot plot
dotplot(ego, showCategory = 15, title = "GO Biological Process: DEGs")

# Bar plot
barplot(ego, showCategory = 15)

# Gene-Concept Network
cnetplot(ego, showCategory = 5, categorySize = "pvalue")

## 3-4. KEGG Pathway 분석 ####
ekk <- enrichKEGG(
  gene         = gene_ids$ENTREZID,
  universe     = bg_genes$ENTREZID,
  organism     = "hsa",          # human
  pvalueCutoff = 0.05
)

head(ekk)
dotplot(ekk, showCategory = 15, title = "KEGG Pathways: DEGs")


# =====================================================================
# Session 4: GSEA 개념 소개 + 결과 해석 (14:30 ~ 15:30)
# =====================================================================

## 4-1. GSEA (Gene Set Enrichment Analysis) 개념 ####
# ORA의 한계를 극복한 방법!
# 
# 핵심 차이:
# ORA: DEG/non-DEG 이분법 → threshold 의존적
# GSEA: 모든 유전자를 rank → threshold-free!
#
# 원리:
# 1) 모든 유전자를 LFC 또는 통계량으로 정렬 (rank)
# 2) 특정 gene set이 리스트의 위쪽(up) / 아래쪽(down)에 몰려있는가?
# 3) Enrichment Score (ES) 계산
# 4) Permutation test로 유의성 검정
#
# 장점:
# - 모든 유전자를 사용 (정보 손실 없음)
# - padj = 0.06인 유전자도 기여할 수 있음
# - 방향성 (up/down) 동시에 파악

## 4-2. GSEA 실습 ####

# (1) Ranked gene list 만들기
# 모든 유전자의 통계량으로 rank

# Gene Symbol → Entrez ID 변환 (전체 유전자)
all_genes_df <- res_df %>%
  filter(!is.na(padj)) %>%
  inner_join(bg_genes, by = c("gene" = "SYMBOL"))

# Ranking: Wald statistic 사용 (방향성 + 유의성 반영)
gene_list <- all_genes_df$stat
names(gene_list) <- all_genes_df$ENTREZID

# 중복 Entrez ID 제거 (최대 절대값 유지)
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)  # 내림차순 정렬!

head(gene_list, 5)   # 가장 UP
tail(gene_list, 5)   # 가장 DOWN

# (2) GSEA 실행 (GO BP)
gsea_res <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)

head(gsea_res, 10)

# (3) GSEA 결과 시각화
# GSEA Enrichment Plot (classic mountain plot)
gseaplot2(gsea_res, geneSetID = 1, title = gsea_res$Description[1])

# 상위 3개 pathway
gseaplot2(gsea_res, geneSetID = 1:3)

# Dot plot (양방향)
dotplot(gsea_res, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA: GO BP (Up vs. Down)")

# Ridge plot
ridgeplot(gsea_res, showCategory = 15) +
  labs(x = "Enrichment Distribution")

## 4-3. 결과 해석: 임상적 의미 ####
# ✅ 해석 연습 (토론)
#
# Q1. GO BP에서 가장 유의한 pathway는 무엇인가?
#     → 폐암의 생물학적 특성과 일치하는가?
#
# Q2. KEGG에서 "Pathways in cancer"가 나타나는가?
#     → 구체적으로 어떤 sub-pathway가 관련되는가?
#
# Q3. Up-regulated DEGs vs. Down-regulated DEGs
#     → 각각 어떤 기능에 enriched 되어 있는가?
#     → (예: cell cycle ↑, immune response ↓?)
#
# Q4. 이 분석 결과가 실제 환자 치료에 어떤 시사점을 줄 수 있는가?
#     → 잠재적 drug target? biomarker?


# =====================================================================
# Session 5: 자유 탐색 + 마무리 (15:30 ~ 17:00)
# =====================================================================

## 5-1. 자유 탐색 아이디어 ####

# 아이디어 A: 관심 유전자 탐색
# 특정 유전자가 종양에서 어떻게 변하는지 확인
my_gene <- "EGFR"  # 원하는 유전자로 바꿔보세요

# 발현량 시각화
plotCounts(dds, gene = my_gene, intgroup = "condition", returnData = FALSE)

# 더 예쁘게
plot_data <- plotCounts(dds, gene = my_gene, intgroup = "condition", returnData = TRUE)
ggplot(plot_data, aes(x = condition, y = count, fill = condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_log10() +
  labs(title = paste0(my_gene, " Expression"),
       y = "Normalized Count (log10)") +
  theme_minimal()

# 아이디어 B: Threshold 바꿔보기
# padj < 0.01 & |LFC| > 2 로 엄격하게 하면?
deg_strict <- res_df %>%
  filter(padj < 0.01, abs(log2FoldChange) > 2) %>%
  pull(gene)
length(deg_strict)

# 아이디어 C: 다른 데이터셋 찾아보기
# GEO에서 관심 있는 질환 데이터셋을 검색해보세요!
# https://www.ncbi.nlm.nih.gov/geo/
# 검색 예시: "Colorectal Cancer AND RNA-seq"

## 5-2. 전체 워크플로우 요약 ####
#
# ┌─────────────────┐
# │  GEO/TCGA에서    │
# │  데이터 탐색     │  ← Day 1 Session 2
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  Metadata 정리   │
# │  + Count 매칭    │  ← Day 1 Session 3
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  Filtering      │
# │  Normalization  │  ← Day 1 Session 4
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  EDA: PCA, QC   │  ← Day 1 Session 5
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  DEG Analysis   │
# │  (DESeq2)       │  ← Day 2 Session 1
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  시각화          │
# │  Volcano/Heatmap │ ← Day 2 Session 2
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  Enrichment     │
# │  ORA + GSEA     │  ← Day 2 Session 3-4
# └────────┬────────┘
#          ▼
# ┌─────────────────┐
# │  임상적 해석     │
# │  + 보고서 작성   │
# └─────────────────┘

## 5-3. 추가 학습 자원 ####
# - DESeq2 Vignette: https://bioconductor.org/packages/DESeq2/
# - clusterProfiler Book: https://yulab-smu.top/biomedical-knowledge-mining-book/
# - 통계 개념: StatQuest (YouTube)
# - R 기초: R for Data Science (https://r4ds.had.co.nz/)

########################################################################
# END OF DAY 2
########################################################################
