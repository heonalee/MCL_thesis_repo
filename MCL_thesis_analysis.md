MCL_thesis_analysis
================
Heona
2026-02-06

- [MCL 1-100 Data Analysis](#mcl-1-100-data-analysis)
- [RNA data](#rna-data)
  - [Import the data & merge](#import-the-data--merge)
  - [Merge with only shared genes](#merge-with-only-shared-genes)
  - [PCA](#pca)
  - [Data check](#data-check)
  - [Batch correction with ComBat](#batch-correction-with-combat)
  - [Import the data](#import-the-data)
  - [Quick raw data check](#quick-raw-data-check)
  - [Sample loading normalization](#sample-loading-normalization)
  - [Principal Component analysis
    (PPCA)](#principal-component-analysis-ppca)
  - [Internal reference scaling](#internal-reference-scaling)
  - [Average replicates and
    log2-transform](#average-replicates-and-log2-transform)
  - [HarmonizR](#harmonizr)
  - [Mapping and unifying colnames](#mapping-and-unifying-colnames)
  - [DreamAI](#dreamai)
  - [Cola clustering](#cola-clustering)

# MCL 1-100 Data Analysis

# RNA data

## Import the data & merge

I have RNA1 and RNA2 from 2 different sources and have to unite them.

``` r
## RNA data extraction

rna_1 <- readRDS("/Users/heona/git-repos/MCL_thesis_repo/data/raw_data/20250303_MCL_RNAseq1.rds")
rna_2 <- readRDS("/Users/heona/git-repos/MCL_thesis_repo/data/raw_data/20250303_MCL_RNAseq2.rds") 

print(paste("Missing values in rna_1:", sum(is.na(rna_1))))
```

    ## [1] "Missing values in rna_1: 0"

``` r
print(paste("Missing values in rna_2:", sum(is.na(rna_2))))
```

    ## [1] "Missing values in rna_2: 0"

``` r
rna_1_data <- rna_1 %>%
  as.data.frame() %>%
  dplyr::select(where(~ !any(is.na(.))))

rna_2_data <- rna_2 %>%
  as.data.frame() %>%
  dplyr::select(where(~ !any(is.na(.))))
```

First, I need to unify the column names.

dplyr::select(any_of(colnames(rna_data))) searches for overlapping
column names, we still have different column names in each matrix - I
have to change them to have uniform colnames.

``` r
mapping <- read_excel("/Users/heona/git-repos/MCL_thesis_repo/data/raw_data/20240830_MCL_IDs_Heona_V2_reduced.xlsx") 

#create a name vector, which replaces the RNAseq_ID names with the MS2_ID names
name_map <- setNames (mapping$MS2_ID, mapping$RNAseq_ID)

# rename columns in rna_2_data
colnames(rna_2_data) <- ifelse(colnames(rna_2_data) %in% names(name_map),
                               name_map[colnames(rna_2_data)], 
                               colnames(rna_2_data))

#create a name vector, which replaces the Cegat_id names with the MS2_ID names
name_map_2 <- setNames (mapping$MS2_ID, mapping$Cegat_id)

# only keep the columns that could be mapped, omit the missing ones that are still named "S..."
rna_2_data <- rna_2_data[, !grepl("^MCL", colnames(rna_2_data))]

# rename columns in rna_1_data
colnames(rna_1_data) <- ifelse(colnames(rna_1_data) %in% names(name_map_2),
                               name_map_2[colnames(rna_1_data)], 
                               colnames(rna_1_data))

# only keep the columns that could be mapped, omit the missing ones that are still named "S..."
rna_1_data <- rna_1_data[, !grepl("^S", colnames(rna_1_data))]
cat("dimension of rna_1_data:", dim(rna_1_data), "\n")
```

    ## dimension of rna_1_data: 14828 25

``` r
#hhow many gene names are overlapping
cat("shared rownames:", sum(rownames(rna_1) %in% rownames(rna_2)), "\n")
```

    ## shared rownames: 9995

``` r
cat("rows only in rna_1:", sum(!(rownames(rna_1) %in% rownames(rna_2))), "\n")
```

    ## rows only in rna_1: 4833

``` r
cat("rows only in rna_2:",sum(!(rownames(rna_2) %in% rownames(rna_1))), "\n")
```

    ## rows only in rna_2: 1851

``` r
#to make a full join of rna 1 and 2, move rownames to a new column called "gene" - bc full_join can only work with columns, not rownames
rna_1_data$gene <- rownames(rna_1_data)
rna_2_data$gene <- rownames(rna_2_data)

# now full_join by gene
rna_data <- full_join(rna_1_data, rna_2_data, by = "gene")
rownames(rna_data) <- rna_data$gene # set the column gene back to rownames
rna_data$gene <- NULL # delete the column gene
```

rna_1_data and rna_2_data are from different labs and measured different
genes. There are 9995 overlapping genes and about 5000 genes unique to
rna_1 and 2000 unique to rna_2.

## Merge with only shared genes

For further analysis, we need to proceed with only the data on the genes
that were shared between both datasets.

``` r
shared_genes <- intersect(rownames(rna_1_data), rownames(rna_2_data))

rna_data_shared <- cbind(rna_1_data[shared_genes, ], rna_2_data[shared_genes, ]) %>%
  select(-gene)
```

## PCA

``` r
sum(is.na(rna_data_shared))
```

    ## [1] 0

``` r
# 1. Filter out genes with too much missing data
threshold <- 0.9
rna_pca_df <- rna_data_shared %>%
  filter(rowMeans(is.na(.)) < threshold)

# Standard PCA (much faster)
rna_pca <- prcomp(t(as.matrix(rna_pca_df)), center = TRUE, scale. = TRUE)

# Extract scores
pca_out <- as.data.frame(rna_pca$x[, 1:2]) %>%
  rownames_to_column("Sample_ID")

# For variance explained, use:
var_explained <- summary(rna_pca)$importance[2, 1:2] * 100

# 5. for ex. "928_04" is plex_sample
pca_out <- pca_out %>%
  separate(Sample_ID, into = c("Plex", "Sample"), sep = "_", remove = FALSE)

# 6. Plot
ggplot(pca_out, aes(PC1, PC2, color = Plex)) +
  geom_point(size = 3) +
  theme_classic() +
  xlab(paste0("PC1 (", round(var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  labs(title = "PCA of shared RNA-seq Samples")
```

![](MCL_thesis_analysis_files/figure-gfm/check%20for%20rna%20batch%20effect-1.png)<!-- -->

The PCA shows a strong batch effect between the two datasets, and also
within the plexes. For BayesdeBulk, which we will perform with this rna
data, it is crucial to not have any batch effects!

## Data check

``` r
# Count how many values are below 0 and get the range of negative values
cat("=== rna_data_shared ===\n")
```

    ## === rna_data_shared ===

``` r
cat("Total sum:", sum(rna_data_shared), "\n")
```

    ## Total sum: 3322889

``` r
cat("Data range:", range(rna_data_shared, na.rm = TRUE), "\n")
```

    ## Data range: -3.525104 14.11628

``` r
cat("Negative values:", sum(rna_data_shared < 0, na.rm = TRUE), "\n")
```

    ## Negative values: 29161

``` r
cat("Range of negative values:", range(rna_data_shared[rna_data_shared < 0], na.rm = TRUE), "\n\n")
```

    ## Range of negative values: -3.525104 -0.0008044756

``` r
cat("=== rna_1_data ===\n")
```

    ## === rna_1_data ===

``` r
cat("Total sum:", sum(rna_1_data[ , sapply(rna_1_data, is.numeric)], na.rm = TRUE), "\n")
```

    ## Total sum: 1769492

``` r
cat("Data range:", range(rna_1_data[ , sapply(rna_1_data, is.numeric)], na.rm = TRUE), "\n")
```

    ## Data range: -2.811264 16.17664

``` r
cat("Negative values:", sum(rna_1_data < 0, na.rm = TRUE), "\n")
```

    ## Negative values: 6503

``` r
cat("Range of negative values:", range(rna_1_data[rna_1_data < 0], na.rm = TRUE), "\n\n")
```

    ## Range of negative values: -0.002074767 -2.81126414

``` r
cat("=== rna_2_data ===\n")
```

    ## === rna_2_data ===

``` r
cat("Total sum:", sum(rna_2_data[ , sapply(rna_2_data, is.numeric)], na.rm = TRUE), "\n")
```

    ## Total sum: 2319244

``` r
cat("Data range:", range(rna_2_data[ , sapply(rna_2_data, is.numeric)], na.rm = TRUE), "\n")
```

    ## Data range: -3.525104 18.29288

``` r
cat("Negative values:", sum(rna_2_data < 0, na.rm = TRUE), "\n")
```

    ## Negative values: 35482

``` r
cat("Range of negative values:", range(rna_2_data[rna_2_data < 0], na.rm = TRUE), "\n")
```

    ## Range of negative values: -0.0008044756 -3.52510448

Expression values ranged from approximately -2,8 to 16,2 in samples from
lab 1 and from -3,5 to 18,3 in samples from lab 2.

### Normalization pattern

Since I have negative values in my RNA datasets that hinder batch
correction methods to function, I want to check if my RNA data is
normalized and that is the cause of the negative values.

``` r
library(tidyverse)
library(cowplot)

rna_data_shared %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  separate(sample, into = c("plex", "number"), sep = "_", remove = FALSE) %>%
  ggplot(aes(x = sample, y = expression, fill = plex)) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, hjust = 1)) +
  labs(x = "", y = "log2 (expression)") +
  scale_fill_brewer(palette = "Set3")
```

![](MCL_thesis_analysis_files/figure-gfm/check%20normalization%20in%20rna-1.png)<!-- -->

``` r
# Density plots
rna_data_shared %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  separate(sample, into = c("plex", "number"), sep = "_", remove = FALSE) %>%
  ggplot(aes(x = expression, color = plex)) +
  geom_density() +
  theme_classic() +
  labs(x = "log2 (expression)", color = "Plex") +
  scale_color_brewer(palette = "Set3")
```

![](MCL_thesis_analysis_files/figure-gfm/check%20normalization%20in%20rna-2.png)<!-- -->

``` r
rna_data_shared %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  separate(sample, into = c("plex", "number"), sep = "_", remove = FALSE) %>%
  ggplot(aes(x = expression, color = plex, fill = plex, group = sample)) +
  geom_density(fill = NA) +
  theme_classic() +
  labs(x = "log2 (expression)", color = "Plex", fill = "Plex") +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3")
```

![](MCL_thesis_analysis_files/figure-gfm/check%20normalization%20in%20rna-3.png)<!-- -->

``` r
# Check column medians
rna_data_shared %>%
  summarise(across(everything(), ~median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "median") %>%
  summary()
```

    ##     sample              median     
    ##  Length:72          Min.   :4.288  
    ##  Class :character   1st Qu.:4.420  
    ##  Mode  :character   Median :4.473  
    ##                     Mean   :4.767  
    ##                     3rd Qu.:5.345  
    ##                     Max.   :5.495

The boxplot visualization shows aligned medians within samples of the
same lab origin, confirming normalization within each batch while
revealing differences in normalization between the two laboratories

## Batch correction with ComBat

with lab origin as batch variable

The samples starting with 9 are from rna_1, the samples starting with 7
are from rna_2.

``` r
batch_vector <- ifelse(grepl("^9", colnames(rna_data_shared)), "rna_1", "rna_2")
modcombat <- model.matrix(~1, data = data.frame(batch = batch_vector))

rna_combat <- ComBat(dat = as.matrix(rna_data_shared),
                        batch = batch_vector,
                        mod = modcombat,
                        par.prior = TRUE)
```

    ## Found2batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

#### PCA after Batch correction

``` r
# Standard PCA
rna_pca_corrected <- prcomp(t(as.matrix(rna_combat)), scale. = TRUE)

# Create data frame for plotting
pca_out_corrected <- as.data.frame(rna_pca_corrected$x[, 1:2]) %>%
  rownames_to_column("Sample_ID") %>%
  separate(Sample_ID, into = c("Plex", "Sample"), sep = "_", remove = FALSE)

# Plot
ggplot(pca_out_corrected, aes(PC1, PC2, color = Plex)) +
  geom_point(size = 3) +
  theme_classic() +
  xlab(sprintf("PC1 (%.2f%%)", summary(rna_pca_corrected)$importance[2,1] * 100)) +
  ylab(sprintf("PC2 (%.2f%%)", summary(rna_pca_corrected)$importance[2,2] * 100)) +
  labs(title = "PCA of rna_combat (After ComBat)")
```

![](MCL_thesis_analysis_files/figure-gfm/standard%20PCA%20of%20rna%20after%20combat-1.png)<!-- -->

color by lab origin:

``` r
# Standard PCA
pca_rna_combat <- prcomp(t(as.matrix(rna_combat)), scale. = TRUE)

# Create data frame for plotting
pca_df <- data.frame(
  PC1 = pca_rna_combat$x[,1],
  PC2 = pca_rna_combat$x[,2],
  Sample_ID = rownames(pca_rna_combat$x)
) %>%
  mutate(Lab = ifelse(substr(Sample_ID, 1, 1) == "9", "RNA_1", "RNA_2"))

# Plot
ggplot(pca_df, aes(PC1, PC2, color = Lab)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("RNA_1" = "#E41A1C", "RNA_2" = "#377EB8")) +
  theme_classic() +
  labs(
    x = sprintf("PC1 (%.1f%%)", summary(pca_rna_combat)$importance[2,1] * 100),
    y = sprintf("PC2 (%.1f%%)", summary(pca_rna_combat)$importance[2,2] * 100),
    title = "PCA of rna_combat (After ComBat) - colored by Lab"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/standard%20PCA%20of%20rna%20after%20combat%20colored%20by%20lab%20origin-1.png)<!-- -->

It looks like the batch correction was successful, the dominant
separation by lab has been eliminated.

\#Proteomic data

## Import the data

``` r
getwd()
```

    ## [1] "/Users/heona/git-repos/MCL_thesis_repo"

``` r
setwd("/Users/heona/git-repos/MCL_thesis_repo")
mcl_original <- read.table("/Users/heona/git-repos/MCL_thesis_repo/data/raw_data/MCL_1-100_Cologne_Essen.txt", header = TRUE, sep = "\t")%>%
  janitor::clean_names()
```

## Quick raw data check

``` r
dim(mcl_original)
```

    ## [1] 6327  253

The dataset contains 6327 rows and 253 columns. I select only the
relevant columns and rename them.

``` r
mcl_unfiltered <- mcl_original %>%
  dplyr::select(uniprot_i_ds, gene_names, contains("reporter_intensity_corrected")) %>%
  dplyr::rename("uniprot_id" = "uniprot_i_ds" ) %>%
  dplyr::rename("gene_id" = "gene_names") %>%
  dplyr::rename_with(.cols = -c("uniprot_id", "gene_id"), .fn = ~ stringr::str_remove(string = ., "reporter_intensity_corrected_"))

dim(mcl_unfiltered)
```

    ## [1] 6327  244

The transformed dataset mcl_unfiltered contains 6327 rows and 244
columns.

### Missing values

Check for number of missing values - R likes to work with NA as missing
values, but in our dataset missing values are NaN.

``` r
print(paste("Total number of missing values: ", sum(is.na(mcl_unfiltered))))
```

    ## [1] "Total number of missing values:  564006"

``` r
print(paste("Proportion of missing values:", 
            round(sum(is.na(mcl_unfiltered)) / length(as.matrix(mcl_unfiltered)) * 100, 2), "%"))
```

    ## [1] "Proportion of missing values: 36.53 %"

### Initial visualization of NAs

``` r
# Process both replicates together
missing_per_sample <- bind_rows(
  mcl_unfiltered %>%
    select(contains("_r1_")) %>%
    summarise(across(everything(), ~sum(is.na(.)) / n() * 100)) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "pct_missing") %>%
    mutate(replicate = "r1"),
  
  mcl_unfiltered %>%
    select(contains("_r2_")) %>%
    summarise(across(everything(), ~sum(is.na(.)) / n() * 100)) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "pct_missing") %>%
    mutate(replicate = "r2")
) %>%
  mutate(plex = str_extract(sample, "\\d{3}$"))

# Faceted plot
ggplot(missing_per_sample, aes(x = reorder(sample, as.numeric(plex)), y = pct_missing, fill = plex)) +
  geom_bar(stat = "identity") +
  facet_wrap(~replicate, scales = "free_x", ncol = 1) +
  labs(x = "Samples", 
       y = "% of proteins missing", 
       fill = "Plex",
       title = "Missing values per sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4))
```

![](MCL_thesis_analysis_files/figure-gfm/visualize%20NAs%20first%20for%20r1%20and%20r2-1.png)<!-- -->

Samples across plexes show varying % of missing values, while within the
plex usually have the same % of missingness. The spikes in missingness
are the “empty” plexes, those that were labeled but not filled with
samples.The mass spec does not measure perfectly, and slight overlap is
usual, explaining the measurements of some signals even in the empty
samples. I exclude replicate 1 and 2 of the TMT-labels that were not
used in the plexes. Also, I exclude the plexes 949 ad 783 because they
did not yield enough data (\> 50% missing values)

``` r
mcl <- mcl_unfiltered %>%
  select(-c("8_r1_775", "9_r1_775", "10_r1_775", # from Plex 5 from MCL1-54
                          "1_r1_930", "2_r1_930", "3_r1_930", "4_r1_930", "5_r1_930", # from Plex 3 from MCL56-100
                          "8_r2_775", "9_r2_775", "10_r2_775", # r2
                          "1_r2_930", "2_r2_930", "3_r2_930", "4_r2_930", "5_r2_930", # r2
                          )) %>%
  select(-contains("949")) %>%
  select(-contains("783"))
```

I convert NaN to NA and remove rows with only NA, then check the data
again.

### Data on protein identified

``` r
#convert NaN to NA
mcl[mcl == "NaN"] <- NA_integer_

#convert Zero values to NA
mcl[,-c(1:2)][mcl[,-c(1:2)] == "0"] <- NA_integer_

#remove rows with only NA
mcl <- mcl %>%
  filter(rowSums(is.na(mcl[,-c(1:2)])) != ncol(mcl[,-c(1:2)]))
```

Data check:

``` r
dim(mcl)
```

    ## [1] 6321  184

``` r
print(paste("Total number of proteins identified:", nrow(mcl)))
```

    ## [1] "Total number of proteins identified: 6321"

``` r
na_per_sample <- colSums(is.na(mcl[,-c(1:2)]))

print(paste("Median number of proteins identified per sample:", nrow(mcl) - median(na_per_sample)))
```

    ## [1] "Median number of proteins identified per sample: 4624"

``` r
print(paste("Minimum number of proteins identified per sample:", nrow(mcl) - max(na_per_sample)))
```

    ## [1] "Minimum number of proteins identified per sample: 3759"

``` r
print(paste("Maximum number of proteins identified per sample:", nrow(mcl) - min(na_per_sample)))
```

    ## [1] "Maximum number of proteins identified per sample: 5243"

``` r
print(paste("Total number of missing values:", sum(is.na(mcl))))
```

    ## [1] "Total number of missing values: 319128"

``` r
print(paste("Proportion of missing values:", 
            round(sum(is.na(mcl)) / length(as.matrix(mcl)) * 100, 2), "%"))
```

    ## [1] "Proportion of missing values: 27.44 %"

``` r
print(paste("Negative values:", sum(mcl[sapply(mcl, is.numeric)] < 0, na.rm = TRUE)))
```

    ## [1] "Negative values: 0"

Proteomics typically has 20-40% missing values, which is much higher
than RNA-seq. This is normal and expected due to the stochastic nature
of mass spectrometry.

### Duplicates

check for duplicated features and samples:

``` r
print(paste("Number of duplicated features: ", sum(duplicated(mcl$uniprot_id))))
```

    ## [1] "Number of duplicated features:  0"

``` r
print(paste("Number of duplicated samples: ", sum(duplicated(colnames(mcl[,-c(1:2)])))))
```

    ## [1] "Number of duplicated samples:  0"

### Raw data distribution of mcl

inspect how the data is distributed without any correction

``` r
order_vec <- colnames(mcl[,-c(1:2)])

# both plots without legend
p1 <- mcl %>%
  dplyr::select(-gene_id) %>%
  pivot_longer(!uniprot_id, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r1_", Prot_id)) %>%
  ggplot(aes(factor(Prot_id, levels = order_vec), log2(intensity), fill = plex)) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 1") +
  scale_fill_brewer(palette = "Set3")

p2 <- mcl %>%
  dplyr::select(-gene_id) %>%
  pivot_longer(!uniprot_id, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r2_", Prot_id)) %>%
  ggplot(aes(factor(Prot_id, levels = order_vec), log2(intensity), fill = plex)) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 2") +
  scale_fill_brewer(palette = "Set3")

# legend from one plot
legend <- get_legend(
  p1 + theme(legend.position = "right")
)
```

    ## Warning: Removed 157266 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

``` r
# combine plots and legend
plot_grid(
  plot_grid(p1, p2, nrow = 2),
  legend,
  rel_widths = c(1, 0.1)
)
```

    ## Warning: Removed 157266 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 161862 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](MCL_thesis_analysis_files/figure-gfm/raw%20data%20distribution-1.png)<!-- -->
This boxplot visualization of raw protein intensities on a log₂-scale
reveals consistent intensity distributions between replicates while
showing variations both within and across plexes indicating the need for
normalization.

### Replicate correlation

Each plex was measured in two replicates - they should be similar - so
we investigate the correlation between both replicates.

``` r
mcl_rep1 <- mcl %>%
  dplyr::select(contains("r1"))
mcl_rep2 <- mcl %>%
  dplyr::select(contains("r2"))

#set the color-vector 
color_fun_corr = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

cor(mcl_rep1, mcl_rep2, method = "spearman", use = "pairwise.complete.obs") %>%
  Heatmap(as.matrix(.),
                   column_title = "Replicate 1",
                   row_title = "Replicate 2",
                   col = color_fun_corr,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   name = "Spearman's R") %>%
  draw()
```

![](MCL_thesis_analysis_files/figure-gfm/Heatmap%20replicate%20correlation-1.png)<!-- -->
This shows an acceptable correlation, though optimal would be a higher
contrast between the diagonal and off-diagonal correlations. Still, we
can combine the replicate values to an average value.

### Combine replicates & reshape dataframe from wide to long

I combine the replicates to a dataframe with the average of r1 and r2.

``` r
mcl_r1 <- mcl %>%
  dplyr::select(uniprot_id, contains("r1")) %>%
    tidyr::pivot_longer(!uniprot_id,
                        names_to = "samples_r1",
                        values_to = "count_R1")

mcl_r2 <- mcl %>%
    dplyr::select(uniprot_id, contains("r2")) %>%
    tidyr::pivot_longer(!uniprot_id,
                        names_to = "samples_r2",
                        values_to = "count_R2")

# mearge mcl_r1 and mcl_r2 and add an average column
mcl_long <- cbind(mcl_r1, mcl_r2[, -1]) # [, -1] to leave out the second uniprot_id column

mcl_avg <- mcl_long %>%
  rowwise() %>%
  mutate(avg_count = mean(c(count_R1, count_R2), na.rm = TRUE)) 

#reformat the NaN in the new column to NA
mcl_avg[mcl_avg == "NaN"] <- NA_integer_
```

### Visualize NAs

``` r
comb <- mcl_avg %>%
  group_by(samples_r1) %>%
  summarise(na_rate = sum(is.na(avg_count))) %>%
  mutate(samples_r2 = samples_r1) %>% 
  separate(samples_r2, c("remove", "Plex"), sep = "r1_") %>% 
  mutate(na_perc = na_rate / nrow(mcl)) 

ggplot(comb, aes(
    x = forcats::fct_reorder(samples_r1, na_perc),
    y = na_perc,
    fill = Plex
  )) +
  geom_col() +
  theme_cowplot() +
  labs(x = "samples", y = "[%] of proteins missing") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y =  element_line(colour = "black"),
    axis.line.x =  element_line(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.title = element_blank(),
    axis.text.x = element_text(
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "none"
  ) 
```

![](MCL_thesis_analysis_files/figure-gfm/visualize%20NAs%20of%20mcl_avg-1.png)<!-- -->

### Reference channel correlation

The reference channel is the 11th sample, so I only select the 11th
sample of each plex as ref_cor.

Then, cor() calculates the pairwise Spearman correlations for the
selected columns starting with 11, considering only complete
observations (pairwise.complete.obs). The resulting correlation matrix
is converted to a data frame using as.data.frame(). The mutate()
function adds new column “type” with the value “ref” to distinguish it
as the reference-channel correlation. pivot_longer() converts the
correlation matrix from wide format to long format, resulting in a
tibble with two columns: “Prot_id” and “count” (the correlation values)
filter(count != 1) removes self-correlation values (correlation of a
variable with itself), as these values are always 1 in correlation
matrices.

``` r
mcl_wide <- mcl_avg %>%
  dplyr::select(c("uniprot_id", "samples_r1", "avg_count")) %>% #mcl_wide has the average count, we only take the name of samples_r1 (in case of confusion)
  pivot_wider(names_from = "samples_r1", values_from = avg_count)

#reference-channel correlation 
ref_cor11 <- mcl_wide %>%
  dplyr::select(starts_with("11_")) %>% 
  cor(., method = "spearman", use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  mutate(type = "ref") %>%
  pivot_longer(!c(type), names_to = "samples", values_to = "count") %>%
  filter(count != 1)

head(ref_cor11, n = 5)
```

    ## # A tibble: 5 × 3
    ##   type  samples   count
    ##   <chr> <chr>     <dbl>
    ## 1 ref   11_r1_757 0.808
    ## 2 ref   11_r1_764 0.807
    ## 3 ref   11_r1_772 0.798
    ## 4 ref   11_r1_775 0.760
    ## 5 ref   11_r1_920 0.808

``` r
ref_cor <- ref_cor11[duplicated(ref_cor11$count),] # remove duplicates that appear twice bc of pivot_longer

head(ref_cor, n = 5)
```

    ## # A tibble: 5 × 3
    ##   type  samples   count
    ##   <chr> <chr>     <dbl>
    ## 1 ref   11_r1_753 0.808
    ## 2 ref   11_r1_753 0.807
    ## 3 ref   11_r1_757 0.839
    ## 4 ref   11_r1_753 0.798
    ## 5 ref   11_r1_757 0.821

``` r
ggplot(ref_cor, aes(count)) +
  geom_density(fill = "#2a9d8f") +
  cowplot::theme_cowplot() +
  labs(x = "Spearman's R" ) 
```

![](MCL_thesis_analysis_files/figure-gfm/ref%20channel%20correlation-1.png)<!-- -->
The reference channel correlation ensures that measurements are accurate
and not biased by technical irregularities. Spearman’s R showed an
overall positive correlation for reference channels and a peak around
0.8, when ideally it would be at \> 0.9. This sub-optimal correlation
suggested the presence of technical variation between plexes, requiring
batch effect correction.

## Sample loading normalization

To correct for differences in total protein amount loaded across samples
and plexes, sample loading normalization is applied. A “target” is
defined as the median total intensity across all samples (median of all
column sums). For each sample, a normalization factor is calculated by
dividing the target by the sample’s total intensity. Each protein
intensity value is then multiplied by its corresponding normalization
factor, scaling all samples to a common total intensity level while
preserving the relative differences between individual proteins. To
verify successful normalization, protein intensities are again
visualized by boxplots on a log2-scale.

``` r
library(dplyr)

# create a vector of plexes
plex_vec <- unique(
  stringr::str_remove(
    colnames(mcl[,-c(1,2)]), 
    "._r1_|._r2_|.._r1_|.._r2_")
  ) 

# Store each plex into a list embedment
plex_list <- list() 

for(i in plex_vec) {
  plex_list[[i]] <- mcl %>% 
    dplyr::select(contains(i))
}

#create the target scaling factor 
colsum_vec <- c() 

for(i in plex_vec) {
  df <- plex_list[[i]]
  colsum_vec <- c(colsum_vec, colSums(df, na.rm = TRUE))
}

target <- median(colsum_vec, na.rm = TRUE)

#correct per plex with a helper function
run_sl_correct <- function(Y) { 
  norm_facs <- target / colSums(Y, na.rm = TRUE) # calculate normalization factor for each column in matrix Y
  output_sl <- sweep(Y, 2, norm_facs, FUN = "*") # multiplies norm_facs with each column
  return(output_sl) # output_sl is the normalized matrix
}

data_sl_repl <- plex_list%>%
  purrr::map(run_sl_correct) %>% 
  bind_cols()
```

visualize effect

``` r
# Plots without legends
p1 <- data_sl_repl %>%
  tibble::rownames_to_column("n") %>%
  pivot_longer(!n, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r1_", Prot_id)) %>%
  ggplot(aes(
    factor(Prot_id, levels = order_vec), log2(intensity), fill = plex
  )) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(
    size = 6,
    angle = 90,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 1") +
  scale_fill_brewer(palette = "Set3")

p2 <- data_sl_repl %>%
  tibble::rownames_to_column("n") %>%
  pivot_longer(!n, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r2_", Prot_id)) %>%
  ggplot(aes(
    factor(Prot_id, levels = order_vec), log2(intensity), fill = plex
  )) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(
    size = 6,
    angle = 90,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 2") +
  scale_fill_brewer(palette = "Set3")

# Extract legend from one plot
legend <- get_legend(
  p1 + theme(legend.position = "right")
)
```

    ## Warning: Removed 157266 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

``` r
# combine plots and legend
plot_grid(
  plot_grid(p1, p2, nrow = 2),
  legend,
  rel_widths = c(1, 0.1)
)
```

    ## Warning: Removed 157266 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 161862 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](MCL_thesis_analysis_files/figure-gfm/plot%20with%20one%20legend%20for%20thesis-1.png)<!-- -->
After sample loading normalization protein intensities are effectively
balanced across all samples. The boxplots show aligned medians and
comparable intensity distributions across both replicates and plexes,
indicating successful correction of technical variability in protein
loading.

### Correlation between random sets of channels

I am comparing the ones with TMT-label 129C since they are present in
every plex.

``` r
library(psych)

sl_data <- data_sl_repl %>%
  dplyr::select(contains("r1")) %>%
  dplyr::select(contains(plex_vec[1:9])) %>%
  dplyr::select(contains("7_"))

pairs.panels(log2(sl_data), lm = TRUE, main = "Random channel 129C over plexes")
```

![](MCL_thesis_analysis_files/figure-gfm/corr%20between%20random%20sets%20after%20normalization-1.png)<!-- -->

## Principal Component analysis (PPCA)

Let’s further have a look into the general structure in the dataset by
PCA. Since we are working with a dataset that contains missing data, we
probabilistic PCA (pPCA) that accounts for missing data
<http://www.cs.columbia.edu/~blei/seminar/2020-representation/readings/TippingBishop1999.pdf>.

Highly correlated samples cluster together in a 2D graph The axes are
ranked in order of importance - differences along x axis (PC1 =
principle component 1) are more important than differences along the 2nd
principal component on y-axis (PC2)

The scales: Large absolute values (e.g., -60 or 20) indicate that those
samples differ significantly from the average sample along that
principal component. Smaller values (closer to 0) suggest that a sample
is near the dataset’s mean in that dimension.

``` r
library(pcaMethods)
library(scrime)

threshold <- 0.5

ppca_df <- data_sl_repl %>% 
  filter(rowMeans(is.na(.)) < threshold) %>%
  #feature-wise scaling and centering 
  rowScales() %>%
  as.data.frame()

data_ppca <- pcaMethods::pca(t(as.matrix(ppca_df)), method = "ppca", nPcs = 2, seed = 123)

ppca_out <- as.data.frame(scores(data_ppca)) %>%
  rownames_to_column("Prot_id") %>%
  separate(Prot_id, c("Number", "Plex"), sep = "_r1_|_r2_")
ggplot(ppca_out, aes(PC1, PC2, col = Plex)) +
  geom_point(size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  xlab(paste("PC1 (", round(data_ppca@R2[1] * 100, digits = 2), "%)")) +
  ylab(paste("PC2 (", round(data_ppca@R2[2] * 100, digits = 2), "%)")) 
```

![](MCL_thesis_analysis_files/figure-gfm/initial%20pca-1.png)<!-- -->

This is the PCA of the mcl data that was corrected with the
normalization factor. The PCA reveals prominent clustering of samples
according to their plex membership, confirming that batch effects were a
major source of variance in the dataset.

## Internal reference scaling

The reference channel (channel 11, containing an identical pooled sample
in each plex) serves as the basis for reference for protein-wise scaling
within each plex. For each protein, the geometric mean of intensities
across all reference channels is defined as the reference average. A
scaling factor is calculated by dividing the ref-erence average by the
observed value in the reference channel of a plex. All values of that
protein within the plex were then multiplied with the corresponding
plex- and protein-specific scaling factor.

``` r
#make a dataframe of the reference channels per plex 
irs_factors <- data_sl_repl %>% 
  dplyr::select(contains("11_"))

#calculate the geometric mean per sample as global reference for scaling
irs_factors$geomean <- apply(irs_factors, 1, function(x) exp(mean(log(x), na.rm = TRUE)))

#pull out sample names to identify  !!! Adjust removal position to geomeam column
irs_factor_vec = as.vector(names(irs_factors[,-c((length(plex_vec)*2)+1)]))
irs_factor_vec
```

    ##  [1] "11_r1_753" "11_r2_753" "11_r1_757" "11_r2_757" "11_r1_764" "11_r2_764"
    ##  [7] "11_r1_772" "11_r2_772" "11_r1_775" "11_r2_775" "11_r1_920" "11_r2_920"
    ## [13] "11_r1_928" "11_r2_928" "11_r1_930" "11_r2_930" "11_r1_935" "11_r2_935"

``` r
#create the protein and reference channel wise scaling factor
for(i in irs_factor_vec) {
  factor <-  irs_factors$geomean / (irs_factors %>% dplyr::select(i))
  irs_factors[,paste0("fac_",i)] <- factor
}
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(i)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(i))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
irs_mult <- colnames(irs_factors %>% dplyr::select(contains("fac"))) %>% as.data.frame

irs_mult <- cbind(irs_mult, irs_mult) 

colnames(irs_mult) <- c("factor", "exp")

irs_mult<- irs_mult %>%
  separate(exp, c("discard", "plex"), sep = "\\_11_") %>%
  dplyr::select(-discard) %>%
  filter(plex !=  "r1_753") # i just took the first plex to start with

irs_mult_vec = as.vector(irs_mult$plex)

irs_fac_filt <- irs_factors %>% dplyr::select(contains("fac"))

all_irs <- (data_sl_repl %>% dplyr::select(contains( "r1_753"))) * unlist(irs_fac_filt %>% dplyr::select(contains( "r1_753"))) # starting with the first plex, the scaling is applied on all samples

for(i in irs_mult_vec) {
  all_irs <- cbind(all_irs, (data_sl_repl %>% dplyr::select(contains(i))) * unlist(irs_fac_filt %>% dplyr::select(contains(i))))
}
```

visualize effect

``` r
# Create plots without legends
p1 <- all_irs %>%
  tibble::rownames_to_column("n") %>%
  pivot_longer(!n, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r1_", Prot_id)) %>%
  ggplot(aes(
    factor(Prot_id, levels = order_vec), log2(intensity), fill = plex
  )) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(
    size = 6,
    angle = 90,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 1") +
  scale_fill_brewer(palette = "Set3")

p2 <- all_irs %>%
  tibble::rownames_to_column("n") %>%
  pivot_longer(!n, names_to = "Prot_id", values_to = "intensity") %>%
  mutate(splitter = Prot_id) %>%
  separate(splitter, c("remove", "plex"), sep = "_r1_|_r2_") %>%
  filter(grepl("_r2_", Prot_id)) %>%
  ggplot(aes(
    factor(Prot_id, levels = order_vec), log2(intensity), fill = plex
  )) +
  geom_boxplot() +
  theme_cowplot() +
  theme(axis.text.x = element_text(
    size = 6,
    angle = 90,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "none") +
  labs(x = "", y = "log2(intensity)", title = "Replicate 2") +
  scale_fill_brewer(palette = "Set3")

# Extract legend from one plot
legend <- get_legend(
  p1 + theme(legend.position = "right")
)
```

    ## Warning: Removed 157378 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

``` r
# Combine plots and legend
plot_grid(
  plot_grid(p1, p2, nrow = 2),
  legend,
  rel_widths = c(1, 0.15)
)
```

    ## Warning: Removed 157378 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 161910 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](MCL_thesis_analysis_files/figure-gfm/boxplots%20after%20internal%20reference%20scaling-1.png)<!-- -->

### Correlation between random set of channels after IRS

``` r
sl_test_data <- all_irs %>%
  dplyr::select(contains("r1")) %>%
  dplyr::select(contains(plex_vec[1:9])) %>%
  dplyr::select(contains("7_"))

pairs.panels(log2(sl_test_data), lm = TRUE, main = "Random channel 129C over plexes")
```

![](MCL_thesis_analysis_files/figure-gfm/correlation%20between%20random%20set%20of%20channels%202-1.png)<!-- -->

### PPCA after IRS

``` r
ppca_df <- all_irs %>% 
  filter(rowMeans(is.na(.)) < threshold) %>%
  #feature-wise scaling and centering 
  rowScales() %>%
  as.data.frame()
data_ppca <- pcaMethods::pca(t(as.matrix(ppca_df)), method = "ppca", nPcs = 2, seed = 123)
ppca_out <- as.data.frame(scores(data_ppca)) %>%
  rownames_to_column("Prot_id") %>%
  separate(Prot_id, c("Number", "Plex"), sep = "_r1_|_r2_")
ggplot(ppca_out, aes(PC1, PC2, col = Plex)) +
  geom_point(size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  xlab(paste("PC1 (", round(data_ppca@R2[1] * 100, digits = 2), "%)")) +
  ylab(paste("PC2 (", round(data_ppca@R2[2] * 100, digits = 2), "%)"))
```

![](MCL_thesis_analysis_files/figure-gfm/check%20ppca-1.png)<!-- -->

After internal reference scaling, data consistency was enhanced,
improving reliability of downstream analyses. Correlation analysis of a
random channel (channel 7) across plexes showed increased correlations
(ρ ≈ 0.91–0.98) compared to pre-IRS values (ρ ≈ 0.63–0.80), confirming
successful batch correction. However, PPCA reveals that samples still
clustered by plex, indicating residual batch effects requiring further
correction.

## Average replicates and log2-transform

I now average replicate values to a single value per protein per sample
to create mcl_proteome_final.

``` r
mcl_norm <- cbind(mcl[,c(1)], all_irs) %>%
  dplyr::select(-contains("11_"))

colnames(mcl_norm)[1] <- "uniprot_id"

# create the datasets for r1 and r2 first separately
mcl_long_norm_r1 <- mcl_norm %>%
    dplyr::select(uniprot_id, contains("r1")) %>%
    tidyr::pivot_longer(!uniprot_id,
                        names_to = "Prot_id1",
                        values_to = "count_R1")

mcl_long_norm_r2 <- mcl_norm %>%
    dplyr::select(uniprot_id, contains("r2")) %>%
    tidyr::pivot_longer(!uniprot_id,
                        names_to = "Prot_id2",
                        values_to = "count_R2") 

head(mcl_long_norm_r1, n = 3)
```

    ## # A tibble: 3 × 3
    ##   uniprot_id               Prot_id1 count_R1
    ##   <chr>                    <chr>       <dbl>
    ## 1 P60709;E7EVS6;A0A6Q8PFE4 1_r1_753 2608005.
    ## 2 P60709;E7EVS6;A0A6Q8PFE4 2_r1_753 2986740.
    ## 3 P60709;E7EVS6;A0A6Q8PFE4 3_r1_753 3179421.

``` r
head(mcl_long_norm_r2, n = 3)
```

    ## # A tibble: 3 × 3
    ##   uniprot_id               Prot_id2 count_R2
    ##   <chr>                    <chr>       <dbl>
    ## 1 P60709;E7EVS6;A0A6Q8PFE4 1_r2_753 2398527.
    ## 2 P60709;E7EVS6;A0A6Q8PFE4 2_r2_753 2961353.
    ## 3 P60709;E7EVS6;A0A6Q8PFE4 3_r2_753 3175631.

``` r
# now merge and leave out uniprot_id (first col) from r2 dataset

mcl_long_norm <- cbind(mcl_long_norm_r1, mcl_long_norm_r2[, -1]) 

mcl_long_norm <- mcl_long_norm %>%
  rowwise() %>%
  mutate(avg_count = mean(c(count_R1, count_R2), na.rm = TRUE))

#reformat the NA 
mcl_long_norm$avg_count[mcl_long_norm$avg_count == "NaN"] <- NA_integer_

mcl_wide_norm <- mcl_long_norm %>% 
  dplyr::select(uniprot_id, Prot_id1, avg_count) %>%
  pivot_wider(names_from = "Prot_id1", values_from = "avg_count")
```

Assign the gene_ids to the samples and log2-transform

``` r
m_comb <- as.matrix(mcl_wide_norm[,-1])

rownames(m_comb) <- mcl_wide_norm$uniprot_id

# t() function swaps rows and cols

t_comb <- t(m_comb) %>% 
  as.data.frame() %>%
  rownames_to_column("Pre_id") %>%
  separate(Pre_id, c("number", "plex"), sep = "_r1_") %>%
  mutate(plex = paste0("P", as.character(plex))) %>%
  mutate(plex = paste0(as.character(plex), as.character(number))) %>%
  dplyr::select(-number) %>% 
  column_to_rownames("plex")

mcl_proteome <-
  t(t_comb) %>% 
  as.data.frame() %>% 
  rownames_to_column("uniprot_id")

uniprot_gene_mcl <- mcl %>% dplyr::select(uniprot_id, gene_id) 

mcl_proteome_final <-  left_join(uniprot_gene_mcl, mcl_proteome, by  = "uniprot_id") 

mcl_proteome_final[,3:ncol(mcl_proteome_final)]<- sapply(mcl_proteome_final[,3:ncol(mcl_proteome_final)], as.numeric)

mcl_proteome_final[,3:ncol(mcl_proteome_final)] <- log2(mcl_proteome_final[,3:ncol(mcl_proteome_final)]) # log-transformed!
```

### PPCA before HarmonizR

``` r
ppca_df <- mcl_proteome_final[,-c(1,2)] %>% 
  filter(rowMeans(is.na(.)) < threshold) %>%
  #feature-wise scaling and centering 
  rowScales() %>%
  as.data.frame()

data_ppca <- pcaMethods::pca(t(as.matrix(ppca_df)), method = "ppca", nPcs = 2, seed = 123)

ppca_out <- as.data.frame(scores(data_ppca)) %>%
  rownames_to_column("Prot_id") %>%
  mutate(Plex = str_remove(str_sub(Prot_id, 1, 4), "P"))
# sep = 4 means separate at the 4th string character

ggplot(ppca_out, aes(PC1, PC2, col = Plex)) +
  geom_point(size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  xlab(paste("PC1 (", round(data_ppca@R2[1] * 100, digits = 2), "%)")) +
  ylab(paste("PC2 (", round(data_ppca@R2[2] * 100, digits = 2), "%)"))
```

![](MCL_thesis_analysis_files/figure-gfm/mcl_proteome_final%20ppca-1.png)<!-- -->
missingness distribution

``` r
comb_df <- comb %>%
  dplyr::select(Plex, na_perc) %>%
  group_by(Plex) %>%
  summarise(na = mean(na_perc)) %>%
  mutate(Plex = paste0("P", Plex, sep = ""))

ppca_out <- as.data.frame(scores(data_ppca)) %>%
  rownames_to_column("Prot_id") %>%
  separate(Prot_id, c("Plex", "Number"), sep = 4) %>%
  left_join(comb_df, 
            by = "Plex")

ggplot(ppca_out, aes(PC1, PC2, col = na)) +
  geom_point(size = 3) +
  scale_color_viridis_c()+
  theme_classic() +
  xlab(paste("PC1 (", round(data_ppca@R2[1] * 100, digits = 2), "%)")) +
  ylab(paste("PC2 (", round(data_ppca@R2[2] * 100, digits = 2), "%)")) +
  labs(color = "% missing")
```

![](MCL_thesis_analysis_files/figure-gfm/final%20data%20NA%20distribution%201-1.png)<!-- -->

## HarmonizR

with Plex = Batch, and with specifying parameters

HarmonizR integrates limma’s removeBatchEffect() function employing
linear regres-sion and the batch effect correction method ComBat which
uses an empirical Bayes framework. Unlike most harmonization tools that
require complete datasets, HarmonizR can handle missing values by matrix
dissection, without imputing or reducing data, which makes it suitable
for proteomic data.

I use the normalized log-transformed data mcl_proteome_final for
HarmonizR.

``` r
print(paste("Total number of missing values: ", sum(is.na(mcl_proteome_final))))
```

    ## [1] "Total number of missing values:  114545"

``` r
prod(dim(mcl_proteome_final))
```

    ## [1] 530964

``` r
print(sum(is.na(mcl_proteome_final)) / prod(dim(mcl_proteome_final)))
```

    ## [1] 0.2157303

We have 21,57% missing values before applying harmonizr.

``` r
# 1: Create the batch vector information
# Extract sample names (excluding uniprot_id and gene_name columns)
# make batch_df with cols: samplename, numbered ID and batch
batch_df <- tibble(samplename = colnames(mcl_proteome_final[,-c(1,2)])) %>%
  rowid_to_column("ID") %>%
  relocate(ID, .after = samplename) %>%
  mutate(
    plex_number = sub("P(.{3}).*", "\\1", samplename),
    batch = case_when(
      plex_number == "753" ~ 1,
      plex_number == "757" ~ 2,
      plex_number == "764" ~ 3,
      plex_number == "772" ~ 4,
      plex_number == "775" ~ 5,
      plex_number == "920" ~ 6,
      plex_number == "928" ~ 7,
      plex_number == "930" ~ 8,
      plex_number == "935" ~ 9,
      TRUE ~ NA_real_
    )
  ) %>%
  select(samplename, ID, batch)

# Check batch assignments
table(batch_df$batch)
```

    ## 
    ##  1  2  3  4  5  6  7  8  9 
    ## 10 10 10 10  7 10 10  5 10

``` r
# 2: Prepare dataframe for harmonization
# Remove gene_name column and set uniprot_id as rownames
harmonize_df <- mcl_proteome_final[,-2] %>%
  column_to_rownames("uniprot_id")

# 3: Run HarmonizR with explicit parameters
mcl_harmonized <- harmonizR(
  data_as_input = harmonize_df,
  description_as_input = batch_df,
  algorithm = "ComBat",
  ComBat_mode = 3,  # parametric adjustment (1 = non-parametric if needed)
  plot = "samplemeans"  # generates diagnostic plots
)
```

    ## Initializing HarmonizR...

    ## Reading the files...

    ## Preparing...

    ## Splitting the data using ComBat adjustment...

    ## Rebuilding...

    ## Writing file...

    ## Saving plot to pdf...

    ## Visualizing samplemeans...

![](MCL_thesis_analysis_files/figure-gfm/harmonizr%20specifying%20the%20correction%20parameters-1.png)<!-- -->

    ## Termination.

``` r
# 4: Convert back to dataframe with protein IDs
mcl_proteome_harmonized <- mcl_harmonized %>%
  as.data.frame() %>%
  rownames_to_column("uniprot_id") %>%
  left_join(mcl_proteome_final[, c("uniprot_id", "gene_id")], by = "uniprot_id") %>%
  relocate(gene_id, .after = uniprot_id)

# store files in the processed_data folder

file.rename("cured_data.tsv", "data/processed_data/mcl_harmonized.tsv")
```

    ## [1] TRUE

``` r
file.rename("cured_data.pdf", "data/processed_data/mcl_harmonized_diagnostics.pdf")
```

    ## [1] TRUE

\###PPCA after HarmonizR

``` r
ppca_df <- mcl_proteome_harmonized[,-c(1,2)] %>% 
  filter(rowMeans(is.na(.)) < threshold) %>%
  #feature-wise scaling and centering 
  rowScales() %>%
  as.data.frame()

data_ppca <- pcaMethods::pca(t(as.matrix(ppca_df)), method = "ppca", nPcs = 2, seed = 123)

ppca_out <- as.data.frame(scores(data_ppca)) %>%
  rownames_to_column("Prot_id") %>%
  separate(Prot_id, into = c("Plex", "Number"), sep = 4)
# sep = 4 means separate at the 4th string character

ggplot(ppca_out, aes(PC1, PC2, col = Plex)) +
  geom_point(size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  xlab(paste("PC1 (", round(data_ppca@R2[1] * 100, digits = 2), "%)")) +
  ylab(paste("PC2 (", round(data_ppca@R2[2] * 100, digits = 2), "%)"))
```

![](MCL_thesis_analysis_files/figure-gfm/ppca%20mcl_proteome_harmonized-1.png)<!-- -->

HarmonizR successfully corrected the batch effects. We proceed with
mcl_proteome_harmonized for our missing value imputation.

## Mapping and unifying colnames

``` r
# for protein data, create a name vector, which replaces the MS_Label names with the MS2_ID names
name_map_ms <- setNames (mapping$MS2_ID, mapping$MS_Label)

proteome_data_log <- mcl_proteome_harmonized

# rename columns in mcl_proteome_cured_data to "753_01" to fit the same format -> prot_data
colnames(proteome_data_log) <- ifelse(colnames(proteome_data_log) %in% names(name_map_ms),
                               name_map_ms[colnames(proteome_data_log)], 
                               colnames(proteome_data_log))
```

proteome_data_log now has the same colnames as the rna_data_shared

We still have multiple genes in some gene_id rownames, that are
seperated by “;” So we select one of the genes (the first) and drop the
rest. If the first of the row is already taken prior, the next is
selected as unique identifier.

### Create unique gene_id identifiers

``` r
# A vector to keep track of the genes that have already been used
used_genes <- c()

# A vector to store the new, unique gene IDs
unique_gene_id <- sapply(proteome_data_log$gene_id, function(genes_str) {
  
  # Split the string of genes by the semicolon
  candidate_genes <- str_split(genes_str, ";")[[1]]
  
  # Find the first gene in the list that has not been used yet
  for (gene in candidate_genes) {
    if (!gene %in% used_genes) {
      used_genes <<- c(used_genes, gene) # If found, add it to our list of used genes
      return(gene) # Return this gene as the chosen one for this row
    }
  }
  # If no unique gene is found for this row (all candidates were already used),
  # return NA as a placeholder.
  return(NA_character_)
}) 

# Add the new column of unique gene IDs to dataframe
proteome_data_log$unique_gene_id <- unique_gene_id

prot_data_log <- proteome_data_log %>%
  select(uniprot_id, gene_id, unique_gene_id, everything())

sum(is.na(unique_gene_id))
```

    ## [1] 165

165 rows cannot be matched with a unique_gene_id, and thus have to be
removed. We only work with those that have both gene_id and uniprot_id
because we focus on the features that are gene/protein-code.

``` r
prot_data_final_log <- prot_data_log %>%
  select(-uniprot_id, -gene_id) %>%
  filter(!is.na(unique_gene_id)) 

rownames(prot_data_final_log) <- prot_data_final_log$unique_gene_id
prot_data_final_log <- prot_data_final_log[,-1]
```

## DreamAI

In order to apply BayesDeBulk deconvolution analysis with RNA and
protein data, we need our datasets to be free of missing values. To
impute the missing values I use DreamAI.

### Missing values

``` r
na_percentage <- rowMeans(is.na(prot_data_final_log)) * 100 #
rows_above_50_na <- which(na_percentage > 50)

row_lowest_na <- which.min(na_percentage) 
min_na_percentage <- na_percentage[row_lowest_na] 

row_highest_na <- which.max(na_percentage) 
max_na_percentage <- na_percentage[row_highest_na] 

# Proportion of rows with >50% missing values
proportion_above_50_na <- sum(na_percentage > 50) / nrow(mcl_proteome_harmonized) *100

# Output results
list(
  num_rows_above_50_na = sum(na_percentage > 50),
  row_lowest_na = row_lowest_na,
  min_na_percentage = min_na_percentage,
  row_highest_na = row_highest_na,
  max_na_percentage = max_na_percentage,
  proportion_above_50_na = proportion_above_50_na
)
```

    ## $num_rows_above_50_na
    ## [1] 1273
    ## 
    ## $row_lowest_na
    ## ACTB 
    ##    1 
    ## 
    ## $min_na_percentage
    ## ACTB 
    ##    0 
    ## 
    ## $row_highest_na
    ## AAMP 
    ## 5449 
    ## 
    ## $max_na_percentage
    ## AAMP 
    ##  100 
    ## 
    ## $proportion_above_50_na
    ## [1] 20.14241

### Filtering for DreamAI

As a preparation, we need to filter the rows that have NA \> 50%, We
only take the columns that are also found in RNA data. For some, RNA
data was missing and we have to exclude these.

``` r
sum(is.na(prot_data_final_log))
```

    ## [1] 109320

``` r
# for DreamAI limit the threshold of missing values to 0.5
threshold <- 0.5

prot_data_final_dai_log <- prot_data_final_log %>%
  filter(rowMeans(is.na(.[,])) < threshold) %>%
  select(, any_of(colnames(rna_combat))) # only the columns which are found in rna_data. bc for some, rna data was missing and we have to exclude these

sum(is.na(prot_data_final_dai_log))
```

    ## [1] 29610

``` r
dim(prot_data_final_log)
```

    ## [1] 6155   82

``` r
dim(prot_data_final_dai_log)
```

    ## [1] 4882   72

``` r
dim(rna_combat)
```

    ## [1] 9995   72

The dimensions are now matched between RNA and protein data.

### Missing values after filtering

Before I run DreamAI, I want to see, how many missing values there are.

``` r
vis_miss(prot_data_final_dai_log[, order(colnames(prot_data_final_dai_log))]) +
  theme(
    axis.text.x = element_text(size = 4, angle = 90, vjust = 1, hjust = 1),
    plot.margin = ggplot2::margin(t = 50, r = 5, b = 5, l = 5)
  )
```

![](MCL_thesis_analysis_files/figure-gfm/missing%20values%20right%20before%20dreamai%20-%20log%20trasnformed%20data-1.png)<!-- -->

The graph shows how many missing values are in each sample now, it
varies from 3% to 22%, on average 8,4%. You can derive from this graph,
that the missing values are not missing at random, they follow a
pattern.

### Run DreamAI()

on prot_data_final_dai_log

Install DreamAI:

``` r
require("cluster")
require("survival")
require("randomForest")
require("missForest")
require("glmnet")
require("Rcpp")
```

    ## Loading required package: Rcpp

    ## Warning: package 'Rcpp' was built under R version 4.4.3

``` r
require("foreach")
require("itertools")
require("iterators")
require("Matrix")
require("devtools")
```

    ## Loading required package: devtools

    ## Loading required package: usethis

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("impute") #original code: BiocManager::install("impute", version = "3.8"), caused problems because my R version is newer (4.4) and is only compatible with BiocManager 3.19. If i need the 3.8 version i need to install R 3.5.x
```

    ## Bioconductor version 3.19 (BiocManager 1.30.27), R 4.4.1 (2024-06-14)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'impute'

    ## Old packages: 'cluster', 'doRNG', 'lubridate', 'skmeans', 'viridisLite'

``` r
require("impute")
```

    ## Loading required package: impute

``` r
BiocManager::version()
```

Run DreamAI

``` r
library(DreamAI)

imputed_prot_data_log <- DreamAI(prot_data_final_dai_log, k = 10, maxiter_MF = 10, ntree = 100,
  maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
  gamma_ADMIN = NA, gamma = 50, CV = FALSE,
  fillmethod = "row_mean", maxiter_RegImpute = 10,
  conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
  "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),
  out = c("Ensemble"))
```

    ## 
    ##  6 methods specified, ensemble imputation will be generated with those algorithms:
    ##  KNN, MissForest, ADMIN, Birnn, SpectroFM, RegImpute

    ## [1] "Method 1 complete"
    ## [1] "Method 2 complete"
    ## [1] "Method 3 complete"
    ## [1] "Method 4 complete"
    ## [1] "Method 5 complete"
    ## [1] "Method 6 complete"

``` r
#imputed_prot_data_log$Ensemble has the new data

#check NAs
sum(is.na(imputed_prot_data_log$Ensemble))
```

    ## [1] 0

``` r
#make it a dataframe (and later maybe a matrix)
prot_data_dreamai_log <- imputed_prot_data_log$Ensemble %>%
  as.data.frame()
```

prot_data_dreamai_log is the imputed dataframe to further work with.

``` r
order_asc = sort(colnames(prot_data_dreamai_log))

prot_data_dreamai_log <- prot_data_dreamai_log[, order_asc]
```

### PCA after DreamAI

of prot_data_dreamai_log

as a last check to not have batch effects standard PCA (bc we have a
complete dataset after DreamAI)

``` r
# Prepare data with scaling
pca_df_dreamai_log <- prot_data_dreamai_log %>%
  as.data.frame() %>%
  rowScales() %>%
  as.data.frame()

# Run standard PCA using prcomp (base R)
data_pca_dreamai_log <- prcomp(t(as.matrix(pca_df_dreamai_log)), 
                                center = FALSE,  # already centered via rowScales
                                scale. = FALSE)  # already scaled via rowScales

# Extract scores and separate sample IDs
pca_out_dreamai <- as.data.frame(data_pca_dreamai_log$x[, 1:2]) %>%
  rownames_to_column("Sample_id") %>%
  separate(Sample_id, into = c("Plex", "Number"), sep = "_")

# Calculate variance explained
var_explained <- summary(data_pca_dreamai_log)$importance[2, 1:2] * 100

# Plot
ggplot(pca_out_dreamai, aes(PC1, PC2, col = Plex)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  xlab(paste0("PC1 (", round(var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  ggtitle("PCA of prot_data_dreamai_log")
```

![](MCL_thesis_analysis_files/figure-gfm/standard%20PCA%20of%20prot_data_dreamai_log-1.png)<!-- -->
No batch effect visible - which means we can take prot_data_dreamai_log
for BayesdeBulk.

\##BayesDeBulk

<https://github.com/WangLab-MSSM/BayesDeBulk.git>

using prot_data_dreamai_log and rna_combat

BayesDeBulk expects MATRICES with: - Genes as ROWS - Samples as COLUMNS

``` r
# Check your data structure
class(rna_combat)     
```

    ## [1] "matrix" "array"

``` r
class(prot_data_dreamai_log)      
```

    ## [1] "data.frame"

``` r
dim(rna_combat)        
```

    ## [1] 9995   72

``` r
dim(prot_data_dreamai_log)        
```

    ## [1] 4882   72

``` r
head(rownames(prot_data_dreamai_log)) # Should be gene names
```

    ## [1] "ACTB"      "HIST2H2BE" "HIST1H4A"  "VIM"       "HBB"       "MYH9"

``` r
head(colnames(prot_data_dreamai_log)) # Should be sample names
```

    ## [1] "753_01" "753_02" "753_03" "753_04" "753_05" "753_06"

``` r
head(rownames(rna_combat))   # Should be gene names
```

    ## [1] "NOC2L"    "KLHL17"   "AGRN"     "C1orf159" "SDF4"     "UBE2J2"

``` r
head(colnames(rna_combat))   # Should be sample names
```

    ## [1] "928_04" "928_09" "935_01" "928_07" "935_02" "930_09"

run BayesDeBulk

``` r
out <-BayesDeBulk(n.iter=5000, burn.in=2500,
                  Y=list(rna_combat,prot_data_dreamai_log), 
                  markers = LM22_markers(list(rna_combat,prot_data_dreamai_log)))

bayesdb_out_log <- out$cell.fraction %>%
  as.data.frame()

write.csv(bayesdb_out_log, 
          file = "data/processed_data/bayesdb_out_log.csv", 
          row.names = TRUE) 
```

#### Cell composition data check

``` r
library(tidyverse)

# 1. Check structure
dim(bayesdb_out_log)
```

    ## [1] 72 22

``` r
colnames(bayesdb_out_log)
```

    ##  [1] "B.cells.naive"                "B.cells.memory"              
    ##  [3] "Plasma.cells"                 "T.cells.CD8"                 
    ##  [5] "T.cells.CD4.naive"            "T.cells.CD4.memory.resting"  
    ##  [7] "T.cells.CD4.memory.activated" "T.cells.follicular.helper"   
    ##  [9] "T.cells.regulatory..Tregs."   "T.cells.gamma.delta"         
    ## [11] "NK.cells.resting"             "NK.cells.activated"          
    ## [13] "Monocytes"                    "Macrophages.M0"              
    ## [15] "Macrophages.M1"               "Macrophages.M2"              
    ## [17] "Dendritic.cells.resting"      "Dendritic.cells.activated"   
    ## [19] "Mast.cells.resting"           "Mast.cells.activated"        
    ## [21] "Eosinophils"                  "Neutrophils"

``` r
# 2. Summary statistics per cell type
celltype_summary <- bayesdb_out_log %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "proportion") %>%
  group_by(cell_type) %>%
  summarise(
    median = median(proportion, na.rm = TRUE),
    mean = mean(proportion, na.rm = TRUE),
    sd = sd(proportion, na.rm = TRUE),
    min = min(proportion, na.rm = TRUE),
    max = max(proportion, na.rm = TRUE),
    cv = sd / mean * 100
  ) %>%
  arrange(desc(median))

print(celltype_summary, n = Inf)
```

    ## # A tibble: 22 × 7
    ##    cell_type                    median   mean     sd     min    max    cv
    ##    <chr>                         <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>
    ##  1 B.cells.naive                0.0533 0.234  0.289  0.00193 0.851  124. 
    ##  2 Neutrophils                  0.0431 0.161  0.263  0.00407 0.914  164. 
    ##  3 Macrophages.M0               0.0282 0.0603 0.152  0.00363 0.872  252. 
    ##  4 Dendritic.cells.activated    0.0259 0.0362 0.0606 0.00195 0.436  167. 
    ##  5 T.cells.regulatory..Tregs.   0.0255 0.0422 0.105  0.00179 0.847  248. 
    ##  6 T.cells.gamma.delta          0.0242 0.0521 0.133  0.00182 0.813  256. 
    ##  7 NK.cells.activated           0.0220 0.0485 0.116  0.00267 0.798  239. 
    ##  8 Mast.cells.resting           0.0164 0.0232 0.0172 0.00192 0.0603  74.1
    ##  9 T.cells.CD4.naive            0.0163 0.0366 0.103  0.00181 0.884  281. 
    ## 10 T.cells.follicular.helper    0.0153 0.0243 0.0190 0.00185 0.0750  78.2
    ## 11 B.cells.memory               0.0153 0.0258 0.0265 0.00191 0.192  103. 
    ## 12 Macrophages.M1               0.0153 0.0237 0.0181 0.00195 0.0592  76.1
    ## 13 T.cells.CD4.memory.activated 0.0149 0.0237 0.0182 0.00181 0.0578  76.7
    ## 14 T.cells.CD8                  0.0143 0.0243 0.0190 0.00174 0.0578  78.2
    ## 15 Eosinophils                  0.0143 0.0231 0.0174 0.00194 0.0667  75.5
    ## 16 Plasma.cells                 0.0142 0.0226 0.0170 0.00193 0.0499  75.4
    ## 17 Monocytes                    0.0142 0.0228 0.0169 0.00203 0.0460  74.1
    ## 18 Dendritic.cells.resting      0.0142 0.0228 0.0171 0.00199 0.0474  75.1
    ## 19 Mast.cells.activated         0.0141 0.0226 0.0168 0.00198 0.0477  74.3
    ## 20 Macrophages.M2               0.0138 0.0228 0.0172 0.00197 0.0476  75.5
    ## 21 NK.cells.resting             0.0135 0.0236 0.0182 0.00180 0.0494  77.1
    ## 22 T.cells.CD4.memory.resting   0.0135 0.0243 0.0189 0.00184 0.0559  77.8

``` r
# 3. Boxplot of cell type proportions
bayesdb_out_log %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "proportion") %>%
  ggplot(aes(x = reorder(cell_type, proportion, FUN = median), y = proportion)) +
  geom_boxplot() +
  coord_flip() +
  theme_classic() +
  labs(x = "Cell type", y = "Estimated proportion",
       title = "Cell type proportions across samples")
```

![](MCL_thesis_analysis_files/figure-gfm/bayesdb_out_log%20data%20check-1.png)<!-- -->

``` r
# 4. Heatmap
library(ComplexHeatmap)

Heatmap(
  t(as.matrix(bayesdb_out_log)),  # transpose so cell types are rows
  name = "Proportion",
  row_title = "Cell types",
  column_title = "Samples",
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
```

![](MCL_thesis_analysis_files/figure-gfm/bayesdb_out_log%20data%20check-2.png)<!-- -->

``` r
# 5. Variability plot (coefficient of variation)
celltype_summary %>%
  ggplot(aes(x = reorder(cell_type, cv), y = cv)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  labs(x = "Cell type", y = "Coefficient of variation (%)",
       title = "Variability of cell type proportions")
```

![](MCL_thesis_analysis_files/figure-gfm/bayesdb_out_log%20data%20check-3.png)<!-- -->

``` r
# 6. Quick summary for your thesis
cat("Number of cell types:", ncol(bayesdb_out_log), "\n")
```

    ## Number of cell types: 22

``` r
cat("Number of samples:", nrow(bayesdb_out_log), "\n\n")
```

    ## Number of samples: 72

``` r
cat("Highest median abundance:\n")
```

    ## Highest median abundance:

``` r
print(head(celltype_summary, 3))
```

    ## # A tibble: 3 × 7
    ##   cell_type      median   mean    sd     min   max    cv
    ##   <chr>           <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl>
    ## 1 B.cells.naive  0.0533 0.234  0.289 0.00193 0.851  124.
    ## 2 Neutrophils    0.0431 0.161  0.263 0.00407 0.914  164.
    ## 3 Macrophages.M0 0.0282 0.0603 0.152 0.00363 0.872  252.

``` r
cat("\nLowest median abundance:\n")
```

    ## 
    ## Lowest median abundance:

``` r
print(tail(celltype_summary, 3))
```

    ## # A tibble: 3 × 7
    ##   cell_type                  median   mean     sd     min    max    cv
    ##   <chr>                       <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>
    ## 1 Macrophages.M2             0.0138 0.0228 0.0172 0.00197 0.0476  75.5
    ## 2 NK.cells.resting           0.0135 0.0236 0.0182 0.00180 0.0494  77.1
    ## 3 T.cells.CD4.memory.resting 0.0135 0.0243 0.0189 0.00184 0.0559  77.8

``` r
cat("\nMost variable cell types (highest CV):\n")
```

    ## 
    ## Most variable cell types (highest CV):

``` r
print(celltype_summary %>% arrange(desc(cv)) %>% head(3))
```

    ## # A tibble: 3 × 7
    ##   cell_type           median   mean    sd     min   max    cv
    ##   <chr>                <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl>
    ## 1 T.cells.CD4.naive   0.0163 0.0366 0.103 0.00181 0.884  281.
    ## 2 T.cells.gamma.delta 0.0242 0.0521 0.133 0.00182 0.813  256.
    ## 3 Macrophages.M0      0.0282 0.0603 0.152 0.00363 0.872  252.

#### UMAP for plexes

``` r
umap_bdb_log <- umap(bayesdb_out_log)

umap_df <- data.frame(
  UMAP1 = umap_bdb_log$layout[,1],
  UMAP2 = umap_bdb_log$layout[,2],
  Sample = rownames(bayesdb_out_log)
) %>%
  separate(Sample, into = c("Plex", "Number"), sep = "_", remove = FALSE)

ggplot(umap_df, aes(UMAP1, UMAP2, color = Plex)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "UMAP of bayesdb_out_log colored by Plex"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/umap%20for%20plexes-1.png)<!-- -->

#### UMAP for single plexes

``` r
umap_bdb_log <- umap(bayesdb_out_log, n_neighbors = 15)

# Create data frame for plotting
umap_df <- data.frame(
  UMAP1 = umap_bdb_log$layout[,1],
  UMAP2 = umap_bdb_log$layout[,2],
  Sample = rownames(bayesdb_out_log)
) %>%
  separate(Sample, into = c("Plex", "Number"), sep = "_", remove = FALSE) %>%
  mutate(Highlight = ifelse(Plex %in% c("930", "935"), Plex, "Other"))

# Plot
ggplot(umap_df, aes(UMAP1, UMAP2, color = Highlight)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("930" = "#E41A1C", "935" = "#377EB8", "Other" = "grey80")) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "UMAP of bayesdb_out_log colored by Plex 930 and 935"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/umap%20single%20plexes-1.png)<!-- -->

UMAP of 928:

``` r
umap_bdb_log <- umap(bayesdb_out_log)

umap_df <- data.frame(
  UMAP1 = umap_bdb_log$layout[,1],
  UMAP2 = umap_bdb_log$layout[,2],
  Sample = rownames(bayesdb_out_log)
) %>%
  separate(Sample, into = c("Plex", "Number"), sep = "_", remove = FALSE) %>%
  mutate(Highlight = ifelse(Plex %in% c("928"), Plex, "Other"))

ggplot(umap_df, aes(UMAP1, UMAP2, color = Highlight)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("928" = "#E41A1C", "Other" = "grey80")) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "UMAP of bayesdb_out_log colored by Plex 928"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/UMAP%20of%20928-1.png)<!-- -->
UMAP of 775

``` r
umap_bdb_log <- umap(bayesdb_out_log)

umap_df <- data.frame(
  UMAP1 = umap_bdb_log$layout[,1],
  UMAP2 = umap_bdb_log$layout[,2],
  Sample = rownames(bayesdb_out_log)
) %>%
  separate(Sample, into = c("Plex", "Number"), sep = "_", remove = FALSE) %>%
  mutate(Highlight = ifelse(Plex %in% c("775"), Plex, "Other"))

ggplot(umap_df, aes(UMAP1, UMAP2, color = Highlight)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("775" = "#E41A1C", "Other" = "grey80")) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "UMAP of bayesdb_out_log colored by Plex 775"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/umap%20of%20775-1.png)<!-- -->

There are no big batch effects visible.

#### UMAP lab origin

``` r
umap_df <- data.frame(
  UMAP1 = umap_bdb_log$layout[,1],
  UMAP2 = umap_bdb_log$layout[,2],
  Sample = rownames(bayesdb_out_log)
) %>%
  mutate(Lab = ifelse(substr(Sample, 1, 1) == "9", "RNA_1", "RNA_2"))

ggplot(umap_df, aes(UMAP1, UMAP2, color = Lab)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("RNA_1" = "#E41A1C", "RNA_2" = "#377EB8")) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "UMAP of bayesdb_out_log - colored by RNA source lab"
  )
```

![](MCL_thesis_analysis_files/figure-gfm/umap%20lab%20origin-1.png)<!-- -->

#### Scaled UMAPs

``` r
# Scaled UMAP
set.seed(123)

bayesdb_log_scaled <- scale(bayesdb_out_log)
umap_scaled <- umap(bayesdb_log_scaled)

# Base data frame
umap_base <- data.frame(
  UMAP1 = umap_scaled$layout[,1],
  UMAP2 = umap_scaled$layout[,2],
  Sample = rownames(bayesdb_log_scaled)
) %>%
  separate(Sample, into = c("Plex", "Number"), sep = "_", remove = FALSE) %>%
  mutate(
    Lab = ifelse(substr(Sample, 1, 1) == "9", "RNA_1", "RNA_2"),
    Highlight_930_935 = case_when(
      Plex == "930" ~ "930",
      Plex == "935" ~ "935",
      TRUE ~ "Other"
    ),
    Highlight_928 = ifelse(Plex == "928", "928", "Other"),
    Highlight_775 = ifelse(Plex == "775", "775", "Other")
  )

# 1. Plex 930 and 935
p1 <- ggplot(umap_base, aes(UMAP1, UMAP2, color = Highlight_930_935)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("930" = "#E41A1C", "935" = "#377EB8", "Other" = "grey80"), name = "Plex") +
  theme_classic() +
  labs(title = "Plex 930 & 935")

# 2. Plex 928
p2 <- ggplot(umap_base, aes(UMAP1, UMAP2, color = Highlight_928)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("928" = "#E41A1C", "Other" = "grey80"), name = "Plex") +
  theme_classic() +
  labs(title = "Plex 928")

# 3. Plex 775
p3 <- ggplot(umap_base, aes(UMAP1, UMAP2, color = Highlight_775)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("775" = "#E41A1C", "Other" = "grey80"), name = "Plex") +
  theme_classic() +
  labs(title = "Plex 775")

# 4. All plexes
p4 <- ggplot(umap_base, aes(UMAP1, UMAP2, color = Plex)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  theme_classic() +
  labs(title = "All Plexes")

# 5. RNA source lab
p5 <- ggplot(umap_base, aes(UMAP1, UMAP2, color = Lab)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("RNA_1" = "#E41A1C", "RNA_2" = "#377EB8")) +
  theme_classic() +
  labs(title = "RNA Source Lab")

# Combine all
library(patchwork)
```

    ## 
    ## Attaching package: 'patchwork'

    ## The following object is masked from 'package:genefilter':
    ## 
    ##     area

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

``` r
(p1 + p2) / p3 +
  plot_annotation(title = "UMAP of bayesdb_out_log (scaled)", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
```

![](MCL_thesis_analysis_files/figure-gfm/scaled%20bayesdb%20umaps-1.png)<!-- -->

``` r
(p4 + p5) +
  plot_annotation(title = "UMAP of bayesdb_out_log (scaled)", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
```

![](MCL_thesis_analysis_files/figure-gfm/scaled%20bayesdb%20umaps-2.png)<!-- -->

``` r
(p1 + p2) / (p4 + p5) +
  plot_annotation(title = "UMAP of bayesdb_out_log (scaled)", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
```

![](MCL_thesis_analysis_files/figure-gfm/scaled%20bayesdb%20umaps-3.png)<!-- -->

## Cola clustering

on bayesdb_out_log

``` r
library(cola)
library(dplyr)

graphics.off()

bayesdb_mat_log <- bayesdb_out_log %>%
  as.matrix() %>%
  t()  

# set seed for reproducibility!
set.seed(123) 

# to prevent the error message "#> Error in dev.off(i2): cannot shut down device 1 (the null device)"
while (!is.null(dev.list()))  dev.off()

#run cola
rl_bayesdb_log <- run_all_consensus_partition_methods(
  bayesdb_mat_log,
  top_value_method = c("SD", "MAD", "ATC"), 
  partition_method = c("hclust", "kmeans", "skmeans"),
  max_k = 6, 
  scale_rows = TRUE,
  cores = 4
)
```

    ## * on a 22x72 matrix.
    ## * calculate top-values.
    ##   - calculate SD score for 22 rows.
    ##   - calculate MAD score for 22 rows.
    ##   - calculate ATC score for 22 rows.
    ## ------------------------------------------------------------
    ## * running partition by SD:skmeans. 1/9
    ## * run SD:skmeans on a 22x72 matrix.
    ## * SD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by SD method
    ##   - skmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).

    ## Loading required package: rngtools

    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * SD:skmeans used 5.064 secs.
    ## ------------------------------------------------------------
    ## * running partition by MAD:skmeans. 2/9
    ## * run MAD:skmeans on a 22x72 matrix.
    ## * MAD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by MAD method
    ##   - skmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * MAD:skmeans used 5.131 secs.
    ## ------------------------------------------------------------
    ## * running partition by ATC:skmeans. 3/9
    ## * run ATC:skmeans on a 22x72 matrix.
    ## * set 4 cores for ATC()
    ## * ATC values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by ATC method
    ##   - skmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * ATC:skmeans used 5.347 secs.
    ## ------------------------------------------------------------
    ## * running partition by SD:kmeans. 4/9
    ## * run SD:kmeans on a 22x72 matrix.
    ## * SD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by SD method
    ##   - kmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * SD:kmeans used 0.5810001 secs.
    ## ------------------------------------------------------------
    ## * running partition by MAD:kmeans. 5/9
    ## * run MAD:kmeans on a 22x72 matrix.
    ## * MAD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by MAD method
    ##   - kmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * MAD:kmeans used 0.552 secs.
    ## ------------------------------------------------------------
    ## * running partition by ATC:kmeans. 6/9
    ## * run ATC:kmeans on a 22x72 matrix.
    ## * ATC values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by ATC method
    ##   - kmeans repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * ATC:kmeans used 0.5569999 secs.
    ## ------------------------------------------------------------
    ## * running partition by SD:hclust. 7/9
    ## * run SD:hclust on a 22x72 matrix.
    ## * SD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by SD method
    ##   - hclust repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * SD:hclust used 0.5350001 secs.
    ## ------------------------------------------------------------
    ## * running partition by MAD:hclust. 8/9
    ## * run MAD:hclust on a 22x72 matrix.
    ## * MAD values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by MAD method
    ##   - hclust repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * MAD:hclust used 0.5250001 secs.
    ## ------------------------------------------------------------
    ## * running partition by ATC:hclust. 9/9
    ## * run ATC:hclust on a 22x72 matrix.
    ## * ATC values have already been calculated. Get from cache.
    ## * rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd
    ## * get top 2 rows by ATC method
    ##   - hclust repeated for 50 times by row-sampling (p = 0.8) from top 2 rows (4 cores).
    ## * wrap results for k = 2
    ## * wrap results for k = 3
    ## * wrap results for k = 4
    ## * wrap results for k = 5
    ## * wrap results for k = 6
    ## * adjust class labels between different k.
    ## * ATC:hclust used 0.5350001 secs.
    ## ------------------------------------------------------------
    ## * adjust class labels according to the consensus classifications from all methods.
    ##   - get reference class labels from all methods, all k.
    ##   - adjust class labels for each single method, each single k.
    ## ------------------------------------------------------------

``` r
best_results_bdb_log <- suggest_best_k(rl_bayesdb_log)
head(best_results_bdb_log, 5)
```

    ##            best_k 1-PAC mean_silhouette concordance    optional_k
    ## SD:hclust       5     1       1.0000000        1.00 **        2,3
    ## SD:kmeans       5     1       0.9915937        0.99 **           
    ## SD:skmeans      6     1       1.0000000        1.00 **        2,3
    ## MAD:hclust      5     1       1.0000000        1.00 **        2,3
    ## MAD:kmeans      5     1       1.0000000        1.00 **          4

``` r
# Generate full report
# eval = FALSE because of knitting problem
cola_report(rl_bayesdb_log, output_dir = "data/processed_data/cola_bayesdb_log", cores = 4)
```
