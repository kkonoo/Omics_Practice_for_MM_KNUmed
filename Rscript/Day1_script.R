########################################################################
# 분자의학 선택실습: 오믹스 데이터 분석 기술
# Day 1 - 공공 오믹스 데이터 탐색 및 전처리
# 담당교수: 하은지 (eha@knu.ac.kr)
########################################################################

# =====================================================================
# Session 1: R/RStudio 환경 세팅 + R 기초 (10:00 ~ 11:00)
# =====================================================================

## 1-1. R과 RStudio 소개 ####
# R: 통계 분석 및 시각화를 위한 프로그래밍 언어
# RStudio: R을 편리하게 사용하기 위한 IDE (통합 개발 환경)
# 화면 구성: Source(코드 작성) | Console(실행) | Environment | Files/Plots

## 1-2. 기본 연산 ####
1 + 1
3 * 5
10 / 3
2^10  # 2의 10제곱

## 1-3. 변수 할당 ####
x <- 10        # 변수에 값 할당 (단축키: Alt + -)
y = 20         # 이것도 가능하지만, <- 를 권장
x + y

my_name <- "홍길동"
my_name

## 1-4. 데이터 타입 ####

# (1) Vector: 같은 타입의 값들의 모음
ages <- c(25, 30, 28, 35, 27)
names <- c("김", "이", "박", "최", "정")
ages
names

# Vector 연산
mean(ages)     # 평균
length(ages)   # 길이
ages[1]        # 첫 번째 원소
ages[c(1,3)]   # 1번째, 3번째

# (2) Data frame: 표 형태의 데이터 (엑셀 시트와 비슷)
patients <- data.frame(
  pt_id   = c("P001", "P002", "P003", "P004", "P005"),
  gender  = c("M", "F", "M", "F", "M"),
  height  = c(175, 162, 178, 158, 180),
  weight  = c(85, 52, 90, 60, 75),
  disease = c("Diabetes", "None", "Hypertension", "Diabetes", "None")
)
patients

# Data frame 접근
patients$height           # 열 선택: $ 사용
patients[1, ]             # 1번째 행
patients[, "disease"]     # "disease" 열
dim(patients)             # 행, 열 수
head(patients, 3)         # 처음 3행

## 1-5. 패키지 설치 및 로드 ####
# 패키지 = 다른 사람이 만든 함수 모음 (앱스토어의 앱과 비슷)
# install.packages("패키지이름")  → 한 번만 하면 됨
# library(패키지이름)             → 매번 R 시작할 때마다 해야 함

# 오늘 사용할 패키지 설치 (처음 한 번만!)
if (!require("tidyverse"))   install.packages("tidyverse")
if (!require("GEOquery"))    BiocManager::install("GEOquery")
if (!require("DESeq2"))      BiocManager::install("DESeq2")

# BiocManager가 없으면 먼저 설치
if (!require("BiocManager")) install.packages("BiocManager")

# 패키지 로드
library(tidyverse)    # dplyr, ggplot2 등 포함
library(GEOquery)
library(DESeq2)

## 1-6. Pipe 연산자 %>% (chain) ####
# "그리고 나서" 라고 읽으면 됩니다
# 단축키: Ctrl + Shift + M

# 기존 방식
head(patients, 3)

# Pipe 방식 (같은 결과)
patients %>% head(3)

# 여러 작업을 연결할 때 강력함
patients %>%
  filter(gender == "M") %>%   # 남성만 골라서
  select(pt_id, height) %>%   # id와 키만 선택해서
  arrange(desc(height))       # 키 내림차순 정렬

# ✅ 여기까지 되셨으면 손 들어주세요!


# =====================================================================
# Session 2: GEO 탐색 + 데이터 다운로드 (11:00 ~ 12:00)
# =====================================================================

## 2-1. 공개 데이터베이스 소개 ####
# GEO (Gene Expression Omnibus): NCBI에서 운영하는 유전체 데이터 저장소
# - https://www.ncbi.nlm.nih.gov/geo/
# - 논문에서 "Data availability" 섹션에 GEO 번호 (GSExxxxxx)
#
# 주요 용어:
# - GSE: Series (연구/논문 단위)
# - GSM: Sample (개별 샘플)
# - GPL: Platform (분석 장비)
#
# 좋은 데이터셋 기준:
# ✓ sample size: 충분한 N (최소 30개 이상)
# ✓ experimental design: control vs. case (정상 vs. 종양)
# ✓ metadata availability: 임상정보 (나이, 성별, 병기 등)

## 2-2. 오늘의 데이터: GSE159857 ####
# 제목: Lung adenocarcinoma (폐선암)
# 디자인: Paired tumor vs. normal (같은 환자의 종양-정상 쌍)
# 샘플 수: 57쌍 (tumor 57 + normal 57 = 114)
# 플랫폼: RNA-seq (bulk)

# GEO에서 직접 검색해봅시다!
# 1) https://www.ncbi.nlm.nih.gov/geo/ 접속
# 2) "GSE159857" 검색
# 3) 어떤 정보가 있는지 살펴보기

## 2-3. GEOquery로 데이터 가져오기 ####
gse_id <- "GSE159857"

# (1) 메타데이터 가져오기
gset_data <- getGEO(gse_id, GSEMatrix = TRUE)
gset_data <- gset_data[[1]]  # 리스트에서 첫 번째 원소 추출

# (2) 메타데이터 (임상 정보) 확인
meta <- pData(gset_data)
dim(meta)                    # 몇 개의 샘플, 몇 개의 변수?
colnames(meta)               # 어떤 변수가 있는지?

# 중요한 변수들 살펴보기
meta %>%
  select(title, geo_accession,
         `age_at_index:ch1`,
         `gender:ch1`,
         `tissue_type:ch1`) %>%
  head(10)

## 2-4. count matrix 다운로드 ####
# GEO의 supplementary files에서 count matrix 다운로드
# 실습에서는 미리 다운받은 데이터를 사용합니다

# 방법 1: getGEOSuppFiles() - 자동 다운로드
# getGEOSuppFiles(gse_id)

# 방법 2: 수동 다운로드 후 읽기 (권장)
# GEO 페이지 하단 Supplementary Files에서 count matrix 다운로드

# 실습용 데이터 읽기 (교수님이 미리 준비한 파일)
# raw_counts <- read.table("data/GSE159857_raw_counts.tsv.gz",
#                          header = TRUE, sep = "\t", row.names = 1)

# ▶ 시간 절약을 위해, 준비된 RData 파일을 로드합니다
# load("data/day1_data.RData")
# 포함: raw_counts, metadata

# [교수님 참고] 아래 코드로 미리 데이터 준비해두세요:
# 실제 수업에서는 위의 load() 한 줄로 진행합니다

## 2-5. 데이터 구조 파악 ####
# (실제 데이터가 로드된 후 실행)
# dim(raw_counts)             # 행(유전자) x 열(샘플)
# head(raw_counts[, 1:5])     # 첫 5개 샘플의 처음 몇 유전자
# rownames(raw_counts)[1:10]  # 유전자 이름 확인

# ✅ 점심시간 전 체크포인트:
# Q1. 이 데이터셋에는 몇 개의 유전자가 있나요?
# Q2. 샘플은 총 몇 개인가요?
# Q3. Tumor와 Normal은 각각 몇 개인가요?


# =====================================================================
# 점심시간 (12:00 ~ 13:30)
# =====================================================================


# =====================================================================
# Session 3: 데이터 구조 파악 + Metadata 매칭 (13:30 ~ 14:30)
# =====================================================================

## 3-1. Metadata 정리하기 ####
# GEO metadata는 보통 지저분 → 깔끔하게 정리 필요

# 필요한 변수만 추출하여 깔끔한 metadata 만들기
metadata <- meta %>%
  select(
    geo_accession,
    title,
    age       = `age_at_index:ch1`,
    gender    = `gender:ch1`,
    condition = `tissue_type:ch1`
  ) %>%
  mutate(
    age       = as.numeric(age),
    condition = factor(condition, levels = c("normal", "tumor")),
    # 환자 ID 추출 (paired design이므로)
    patient   = str_extract(title, "LUAD-\\d+")
  )

head(metadata, 10)
table(metadata$condition)   # normal vs. tumor 수
table(metadata$gender)      # 성별 분포

## 3-2. count matrix와 metadata 매칭 (중요!) ####
# ⚠️ 가장 흔한 실수: count matrix의 열 순서와 metadata의 행 순서가 불일치

# 순서 확인
all(colnames(raw_counts) == metadata$geo_accession)  # TRUE여야 함!

# 만약 FALSE라면 재정렬
# metadata <- metadata[match(colnames(raw_counts), metadata$geo_accession), ]

## 3-3. 기본 탐색 ####
# (1) Library size (총 read 수) 확인
lib_sizes <- colSums(raw_counts)
summary(lib_sizes)

# Library size 시각화
lib_df <- data.frame(
  sample    = colnames(raw_counts),
  lib_size  = lib_sizes,
  condition = metadata$condition
)

ggplot(lib_df, aes(x = reorder(sample, lib_size), y = lib_size / 1e6,
                   fill = condition)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Sample", y = "Library Size (Millions)", 
       title = "Library Size per Sample") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5))

# (2) 발현 분포 확인
# 특정 유전자의 count 분포
hist(as.numeric(raw_counts["TP53", ]), 
     main = "TP53 Expression", xlab = "Raw Count", breaks = 20)


# =====================================================================
# Session 4: 전처리 + Normalization (14:30 ~ 15:30)
# =====================================================================

## 4-1. RNA-seq 데이터 분석 워크플로우 ####
# ① Read trimming      (FASTQ → trimmed FASTQ)   ← 우리는 여기부터 안 함
# ② Alignment          (trimmed FASTQ → BAM)
# ③ Quantification     (BAM → Count Matrix)
# ④ Normalization      (Count Matrix → Normalized)  ← 여기부터 시작!
# ⑤ Differential Analysis (DEG 분석)
# ⑥ Functional Profiling  (Pathway 분석)

## 4-2. Gene Filtering ####
# 낮은 발현 유전자 제거 → 통계적 검정력 향상

# 현재 유전자 수
nrow(raw_counts)

# 기준: 전체 샘플의 50% 이상에서 count > 0인 유전자만 유지
keep <- rowSums(raw_counts > 0) >= ncol(raw_counts) * 0.5
filtered_counts <- raw_counts[keep, ]
nrow(filtered_counts)
cat(sprintf("Filtered: %d → %d genes\n", nrow(raw_counts), nrow(filtered_counts)))

## 4-3. Normalization 개념 ####
# 왜 normalization이 필요한가?
# - 샘플마다 sequencing depth(library size)가 다름
# - 이 차이를 보정해야 샘플 간 비교가 가능

# 주요 normalization 방법:
# | 방법        | 보정           | 비교     | 비고              |
# |-------------|----------------|----------|-------------------|
# | Raw Count   | 없음           | 불가     | DESeq2/edgeR 입력 |
# | RPKM/FPKM   | depth + length | 샘플 내  | 사용하지 마세요   |
# | TPM         | depth + length | 샘플 간  | 발현량 비교용     |
# | TMM / RLE   | depth + 조성   | 샘플 간  | DEG 분석 표준!    |

## 4-4. DESeq2로 Normalization 실습 ####

# (1) DESeq2 입력 객체 만들기 (DESeqDataSet)
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData   = metadata,
  design    = ~ condition  # 어떤 조건을 비교할 건지
)
dds

# (2) Size Factor 계산 (RLE normalization)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Size Factor 시각화
sf_df <- data.frame(
  sample = colnames(dds),
  sf     = sizeFactors(dds),
  condition = metadata$condition
)
ggplot(sf_df, aes(x = reorder(sample, sf), y = sf, fill = condition)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Sample", y = "Size Factor", title = "DESeq2 Size Factors") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5))

# (3) Normalized count 추출
normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts[, 1:5])

# Raw vs. Normalized 비교
par(mfrow = c(1, 2))
boxplot(log2(counts(dds, normalized = FALSE) + 1)[, 1:20],
        main = "Raw Counts (log2)", las = 2, cex.axis = 0.6)
boxplot(log2(counts(dds, normalized = TRUE) + 1)[, 1:20],
        main = "Normalized Counts (log2)", las = 2, cex.axis = 0.6)
par(mfrow = c(1, 1))


# =====================================================================
# Session 5: EDA - PCA Plot + QC (15:30 ~ 17:00)
# =====================================================================

## 5-1. 왜 EDA를 하는가? ####
# - 데이터 구조 확인
# - 이상 샘플(outlier) 탐지
# - Batch effect 존재 여부 확인
# - 기대한 대로 그룹이 분리되는지 시각적 확인

## 5-2. rlog 변환 ####
# DESeq2는 log2 변환에 rlog 또는 vst를 권장
# rlog: 작은 데이터셋에 적합 (느리지만 정확)
# vst: 큰 데이터셋에 적합 (빠름)

rld <- rlog(dds, blind = TRUE)
# blind = TRUE: design 정보를 무시하고 변환 (탐색적 분석용)

## 5-3. PCA Plot ####
# PCA (Principal Component Analysis): 고차원 데이터를 2차원으로 축소

# DESeq2 내장 함수 사용
plotPCA(rld, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA: Tumor vs. Normal")

# 성별로도 확인
plotPCA(rld, intgroup = "gender") +
  theme_minimal() +
  ggtitle("PCA: by Gender")

## 5-4. PCA 해석 ####
# ✅ 해석 포인트:
# 1) Tumor와 Normal이 잘 분리되는가? → 생물학적 신호가 있다!
# 2) 성별이나 나이 같은 변수에 의해 분리되는가? → Batch effect?
# 3) 이상하게 동떨어진 샘플은 없는가? → Outlier 확인

## 5-5. Sample-to-Sample Distance ####
# Heatmap으로 샘플 간 유사도 확인

library(pheatmap)
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)

# annotation 추가
ann_col <- data.frame(
  Condition = metadata$condition,
  Gender    = metadata$gender,
  row.names = colnames(dds)
)

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = ann_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample-to-Sample Distance")

## 5-6. Day 1 마무리 ####
# 
# 오늘 배운 것:
# ✓ R 기본 문법 (변수, 벡터, 데이터프레임, pipe)
# ✓ GEO에서 데이터 탐색 및 다운로드
# ✓ Metadata 정리 및 count matrix 매칭
# ✓ Gene filtering & Normalization (DESeq2)
# ✓ EDA: PCA plot, sample distance heatmap
#
# 내일 할 것:
# → DEG 분석 (DESeq2)
# → 결과 시각화 (Volcano plot, Heatmap)
# → Functional Enrichment (GO/KEGG pathway)
# → 임상적 해석

# DESeq2 객체 저장 (내일 이어서 사용)
# save(dds, metadata, filtered_counts, file = "data/day1_results.RData")

########################################################################
# END OF DAY 1
########################################################################
