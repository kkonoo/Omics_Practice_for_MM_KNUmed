library(tidyverse)    # dplyr, ggplot2 등 포함
library(DESeq2)

setwd("T:/homes/eha/교육/[학부] 4-1 분자의학 선택연구/2026/")

# 1. read data (dds) ####
dds <- readRDS('data/1_2_dds_filtered.rds')
dds@colData %>% head
dds@design


# 2. DESeq2 running ####
# Dispersion estimation -> GLM Fitting -> Wald Test
# Simple model: ~ condition
# Paired model: ~ patient + condition (환자별 변동을 보정)
#               → 더 정밀한 분석!

dds <- DESeq(dds)
# ① size factor 계산 (normalization)
# ② dispersion 추정
# ③ GLM fitting
# ④ Wald test (통계 검정)
saveRDS(dds, file='data/2_1_dds_DESeq2.rds')


dds@design = as.formula('~ condition')
dds_smpl <- DESeq(dds)
saveRDS(dds_smpl, file='data/2_2_dds_DESeq2_simple.rds')


# 3. check Dispersion Plot ####
dds <- readRDS('data/2_1_dds_DESeq2.rds')
plotDispEsts(dds)

dds_smpl <- readRDS(file='data/2_2_dds_DESeq2_simple.rds')
plotDispEsts(dds_smpl)

