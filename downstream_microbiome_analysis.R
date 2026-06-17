################################################################################
## Downstream phyloseq import and analyses
################################################################################

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(picante)
  library(ggpubr)
  library(ANCOMBC)
  library(mia)
  library(microViz)
  library(TreeSummarizedExperiment)
  library(ggord)
  library(microbiome)
  library(caret)
  library(randomForest)
  library(pROC)
})


## ---------------------------------------------------------------------------
## Alpha diversity
## ---------------------------------------------------------------------------

evenness <- read_qza("data/core-metrics-results/evenness_vector.qza")$data

adiv <- data.frame(
  Evenness  = evenness,
  Observed  = phyloseq::estimate_richness(physeq, measures = "Observed")[, 1],
  Shannon   = phyloseq::estimate_richness(physeq, measures = "Shannon")[, 1],
  Faith_PD  = picante::pd(
    samp = data.frame(t(data.frame(phyloseq::otu_table(physeq)))),
    tree = phyloseq::phy_tree(physeq)
  )[, 1],
  Group_A = phyloseq::sample_data(physeq)$Group_A,
  Group_B = phyloseq::sample_data(physeq)$Group_B
)

wilcox.test(Observed ~ Group_A, data = subset(adiv, Group_A %in% c("HC", "RA")))
wilcox.test(Faith_PD ~ Group_A, data = subset(adiv, Group_A %in% c("HC", "RA")))

## ---------------------------------------------------------------------------
## Beta diversity: unweighted UniFrac PCoA and PERMANOVA
## ---------------------------------------------------------------------------

uw_pcoa <- read_qza("data/core-metrics-results/unweighted_unifrac_pcoa_results.qza")
uw_dist <- read_qza("data/core-metrics-results/unweighted_unifrac_distance_matrix.qza")

metadata <- data.frame(phyloseq::sample_data(physeq))
metadata$SampleID <- rownames(metadata)

adonis2(
  uw_dist$data ~ Group_A,
  data = metadata,
  permutations = 99999
)

## Optional plotting code retained as rough outline only
pcoa_df <- uw_pcoa$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  dplyr::left_join(metadata, by = "SampleID")

ggplot(pcoa_df, aes(PC1, PC2, color = Group_A)) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Group_A), linetype = 2) +
  theme_bw()

## ---------------------------------------------------------------------------
## ANCOM-BC2: taxonomic differential abundance
## ---------------------------------------------------------------------------

## HC vs RA, genus level
ancombc.genus.HCRA <- ancombc2(
  data = physeq.rmOA,
  fix_formula = "Group_A",
  assay_name = "counts",
  tax_level = "Genus",
  pseudo = 0,
  pseudo_sens = TRUE,
  p_adj_method = "BH",
  prv_cut = 0.2,
  lib_cut = 0,
  group = "Group_A",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  verbose = TRUE
)

res.genus.HCRA <- ancombc.genus.HCRA$res

res.genus.HCRA.sig <- res.genus.HCRA %>%
  dplyr::filter(diff_Group_ARA == 1) %>%
  dplyr::arrange(desc(lfc_Group_ARA))

write.csv(
  res.genus.HCRA.sig,
  "results/ancombc_genus_HC_vs_RA.csv",
  row.names = FALSE
)

## SPRA vs SNRA, genus level
ancombc.genus.RA <- ancombc2(
  data = physeq.RA,
  fix_formula = "Group_B",
  assay_name = "counts",
  tax_level = "Genus",
  pseudo = 0,
  pseudo_sens = TRUE,
  p_adj_method = "BH",
  prv_cut = 0.2,
  lib_cut = 0,
  group = "Group_B",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  verbose = TRUE
)

res.genus.RA <- ancombc.genus.RA$res

write.csv(
  res.genus.RA,
  "results/ancombc_genus_SPRA_vs_SNRA.csv",
  row.names = FALSE
)


## ---------------------------------------------------------------------------
## HC vs RA, RDA
## ---------------------------------------------------------------------------

rdadata <- makeTreeSummarizedExperimentFromPhyloseq(physeq)
rdadata <- mia:::transformAssay(rdadata, method = "relabundance")

rdadata <- mia:::runRDA(
  rdadata,
  assay.type = "relabundance",
  formula = assay ~ Periodontal_condition + Group_A + Sex,
  distance = "bray",
  na.action = na.exclude
)

rda_info <- attr(reducedDim(rdadata, "RDA"), "significance")
rda_info$permanova
rda_info$homogeneity

write.csv(
  rda_info$permanova,
  "results/dbRDA_HC_RA_permanova.csv"
)

## ---------------------------------------------------------------------------
## Representative HC-vs-RA random forest/RFE block
## ---------------------------------------------------------------------------

rf_data <- read.csv("results/RF_metadata_genus.csv", check.names = FALSE)
rf_data <- rf_data %>%
  dplyr::filter(Group_A %in% c("HC", "RA"))

rf_data$Group_A <- as.factor(rf_data$Group_A)

## In the original analysis, the input feature range was the filtered genus table.
## Column positions should be checked against the private metadata table.
feature_cols <- setdiff(
  colnames(rf_data),
  c("#SampleID", "Group_A", "Group_B")
)

feature_cols <- feature_cols[sapply(rf_data[, feature_cols, drop = FALSE], is.numeric)]

ctrl <- rfeControl(
  functions = rfFuncs,
  method = "repeatedcv",
  number = 10,
  repeats = 25,
  verbose = FALSE
)

subsets <- seq_along(feature_cols)

## One representative split.
## Full repeated objects and fitted models are intentionally not deposited here.
randomsamples <- rf_data$Group_A
rss <- split(seq_along(randomsamples), randomsamples)

idx <- sort(as.numeric(unlist(
  sapply(rss, function(x) sample(x, length(x) * 0.8))
)))

train_df <- rf_data[idx, ]
test_df  <- rf_data[-idx, ]

x_train <- train_df[, feature_cols, drop = FALSE]
y_train <- train_df$Group_A

set.seed(121)

rf_rfe_fit <- rfe(
  x = x_train,
  y = y_train,
  sizes = subsets,
  rfeControl = ctrl
)

selected_features <- predictors(rf_rfe_fit)

write.csv(
  data.frame(selected_features = selected_features),
  "results/rf_selected_features_representative.csv",
  row.names = FALSE
)

## Training ROC
train_prob <- as.data.frame(rf_rfe_fit$fit$votes)
train_prob$observed <- y_train

roc_train <- pROC::roc(
  response = ifelse(train_prob$observed == "RA", "RA", "HC"),
  predictor = as.numeric(train_prob$RA),
  direction = "<",
  levels = c("HC", "RA"),
  percent = TRUE
)

pROC::auc(roc_train)

## Test ROC
test_prob <- as.data.frame(
  predict(
    rf_rfe_fit$fit,
    test_df[, selected_features, drop = FALSE],
    type = "prob"
  )
)

test_prob$observed <- test_df$Group_A

roc_test <- pROC::roc(
  response = ifelse(test_prob$observed == "RA", "RA", "HC"),
  predictor = as.numeric(test_prob$RA),
  direction = "<",
  levels = c("HC", "RA"),
  percent = TRUE
)

pROC::auc(roc_test)

write.csv(
  data.frame(
    set = c("training", "test"),
    auc = c(as.numeric(pROC::auc(roc_train)), as.numeric(pROC::auc(roc_test)))
  ),
  "results/rf_auc_representative.csv",
  row.names = FALSE
)


## ---------------------------------------------------------------------------
## Import PICRUSt2 MetaCyc pathway abundance table
## ---------------------------------------------------------------------------

metacyc <- readr::read_tsv(
  gzfile("data/picrust2_pipeline/pathways_out/path_abun_descrip.tsv.gz")
)

metacyc_describ <- metacyc[, 1:2]
metacyc_otu <- metacyc[, -c(1:2)]
metacyc_otu <- as.matrix(metacyc_otu)

rownames(metacyc_otu) <- metacyc_describ$pathway
metacyc_otu2 <- otu_table(metacyc_otu, taxa_are_rows = TRUE)

metacyc_describ <- as.matrix(metacyc_describ)
rownames(metacyc_describ) <- rownames(metacyc_otu)

metacyc_physeq <- phyloseq(
  metacyc_otu2,
  tax_table(metacyc_describ),
  sample_data(physeq)
)

## ---------------------------------------------------------------------------
## Pathway-level beta diversity
## ---------------------------------------------------------------------------

metacyc_otu_tab <- t(data.frame(
  phyloseq::otu_table(metacyc_physeq),
  check.names = FALSE
))

metacyc_bc_dist <- vegan::vegdist(metacyc_otu_tab, method = "bray")

adonis2(
  metacyc_bc_dist ~ sample_data(metacyc_physeq)$Group_A,
  permutations = 999
)

## ---------------------------------------------------------------------------
## Pathway-level ANCOM-BC
## ---------------------------------------------------------------------------

ancombc.metacyc <- ancombc2(
  data = metacyc_physeq,
  fix_formula = "Group_A",
  assay_name = "counts",
  pseudo = 0,
  pseudo_sens = TRUE,
  p_adj_method = "BH",
  prv_cut = 0.2,
  lib_cut = 0,
  group = "Group_A",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  verbose = TRUE
)

res.metacyc <- ancombc.metacyc$res

res.metacyc.sig <- res.metacyc %>%
  dplyr::filter(diff_Group_ARA == 1) %>%
  dplyr::arrange(desc(lfc_Group_ARA)) %>%
  dplyr::mutate(Group = ifelse(lfc_Group_ARA > 0, "RA", "HC"))

write.csv(
  res.metacyc.sig,
  "results/ancombc_metacyc_HC_vs_RA.csv",
  row.names = FALSE
)

sessionInfo()