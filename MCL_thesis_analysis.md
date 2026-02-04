MCL_thesis_analysis
================
Heona
2026-02-04

- [MCL 1-100 Data Analysis](#mcl-1-100-data-analysis)
  - [Import the data](#import-the-data)
  - [Quick raw data check](#quick-raw-data-check)
  - [Sample loading normalization](#sample-loading-normalization)
  - [Principal Component analysis
    (PPCA)](#principal-component-analysis-ppca)
  - [Internal reference scaling](#internal-reference-scaling)
  - [Average replicates and
    log2-transform](#average-replicates-and-log2-transform)

# MCL 1-100 Data Analysis

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
    ## This warning is displayed once every 8 hours.
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

head(mcl_proteome_final, n = 5)
```

    ##                 uniprot_id   gene_id    P7531    P7532    P7533    P7534
    ## 1 P60709;E7EVS6;A0A6Q8PFE4      ACTB 21.25538 21.50400 21.59947 21.59208
    ## 2                   Q16778 HIST2H2BE 20.46697 20.52988 20.45952 20.54794
    ## 3                   P62805  HIST1H4A 20.25392 20.49099 20.36178 20.42972
    ## 4     P08670;B0YJC4;B0YJC5       VIM 20.20285 20.56188 20.41011 20.44556
    ## 5 P68871;A0A2R8Y7R2;F8W6P5       HBB 19.68119 19.55375 18.78971 18.84698
    ##      P7535    P7536    P7537    P7538    P7539   P75310    P7571    P7572
    ## 1 21.29499 21.25602 21.46255 21.26419 21.34793 21.39515 21.41513 21.39643
    ## 2 20.53092 20.70695 20.86797 20.62257 20.44453 20.38031 20.13698 20.41826
    ## 3 20.48301 20.62735 20.51582 20.36780 20.37991 20.36285 20.12024 19.65112
    ## 4 20.17874 20.22294 20.14052 20.03989 20.31954 20.41859 20.13465 20.99085
    ## 5 18.49105 18.66985 18.72256 18.97647 18.76937 21.13533 18.89944 19.32462
    ##      P7573    P7574    P7575    P7576    P7577    P7578    P7579   P75710
    ## 1 21.43436 21.33289 21.47828 21.37624 21.35000 21.47227 20.96708 21.18102
    ## 2 20.30474 20.46835 20.74898 20.73283 20.55075 20.40435 20.49777 20.68563
    ## 3 19.77119 20.18415 20.78352 20.22991 20.40824 20.40575 20.33049 20.09902
    ## 4 20.28441 20.44071 20.37807 20.74960 20.20414 20.48513 20.68957 20.90614
    ## 5 19.17290 20.22663 18.99977 20.39055 19.31895 20.14404 20.90423 19.72342
    ##      P7641    P7642    P7643    P7644    P7645    P7646    P7647    P7648
    ## 1 21.34903 21.59507 21.64406 21.12795 21.56376 21.34022 21.44107 21.32802
    ## 2 20.42230 20.36619 20.54228 20.38481 20.55238 20.78299 20.48386 20.50827
    ## 3 20.40239 20.41916 20.62991 20.38244 20.49320 20.64485 20.51314 20.51041
    ## 4 20.60273 20.66354 20.54610 20.76804 20.35895 20.84403 20.46041 20.43229
    ## 5 19.80044 18.58696 18.64644 19.75790 19.17350 20.59769 18.50594 19.61735
    ##      P7649   P76410    P7721    P7722    P7723    P7724    P7725    P7726
    ## 1 21.52692 21.01937 21.33921 20.92704 21.60750 21.09370 21.09261 21.51267
    ## 2 20.45336 20.03532 20.23541 20.69217 20.47261 20.08648 20.12233 20.55224
    ## 3 20.40663 19.97182 19.66278 20.02035 19.92325 19.46795 19.43286 19.93502
    ## 4 20.43109 20.11469 20.49874 20.73238 20.55397 20.96589 20.08385 20.50914
    ## 5 18.78198 20.15138 18.56865 18.49280 18.32110 19.32107 20.16932 19.32317
    ##      P7727    P7728    P7729   P77210    P7751    P7752    P7753    P7754
    ## 1 21.36028 21.43331 21.21805 21.32893 21.45440 21.42877 21.38309 21.34356
    ## 2 20.69891 20.47119 20.42756 20.44933 20.38400 20.25394 20.17956 20.48505
    ## 3 19.93076 20.44324 19.31492 19.86442 20.39222 20.47930 20.51854 20.75431
    ## 4 20.50586 19.98832 20.05222 20.63264 20.48248 20.41749 20.63893 20.93500
    ## 5 18.57092 18.46474 19.45258 18.84197 18.50902 18.41384 20.43896 19.45596
    ##      P7755    P7756    P7757    P9201    P9202    P9203    P9204    P9205
    ## 1 21.36617 21.36859 21.73759 21.54555 21.21855 21.34537 21.24243 21.21813
    ## 2 20.41600 20.48677 20.53714 21.32035 20.29784 20.67867 20.29872 20.29330
    ## 3 20.52538 20.39171 20.99042 21.45618 20.46075 20.67376 20.47516 20.41941
    ## 4 20.43464 20.38168 20.38959 20.26493 20.50756 20.44750 20.44536 20.39661
    ## 5 19.43949 18.69047 18.22380 18.56686 18.19078 18.25130 19.29689 18.75007
    ##      P9206    P9207    P9208    P9209   P92010    P9281    P9282    P9283
    ## 1 21.61873 21.42747 21.02946 21.52195 21.60911 21.48538 21.44603 21.27240
    ## 2 20.69828 20.58633 20.50871 20.45016 20.55552 20.79433 21.22000 21.16860
    ## 3 20.67824 20.63921 20.67021 20.51547 21.01167 20.37524 20.85611 20.82479
    ## 4 20.76787 20.40386 20.62240 20.32258 20.66715 20.33655 20.32233 20.20510
    ## 5 18.75604 18.86693 19.52339 18.75412 19.08371 18.87944 19.39355 19.06970
    ##      P9284    P9285    P9286    P9287    P9288    P9289   P92810    P9306
    ## 1 20.84597 21.26707 21.03414 21.03364 21.73153 21.04122 21.24085 20.97424
    ## 2 21.12536 20.58400 20.55890 20.93316 20.66639 20.59387 21.58611 20.64435
    ## 3 20.71038 20.22969 20.11753 20.65574 20.33278 20.10565 21.17640 20.62872
    ## 4 20.05512 20.31082 20.60158 20.25692 20.35125 20.48713 19.95409 20.23254
    ## 5 19.36177 19.48795 20.11572 20.83975 18.88002 19.04300 20.07019 20.74820
    ##      P9307    P9308    P9309   P93010    P9351    P9352    P9353    P9354
    ## 1 21.22711 21.49257 21.08876 21.22353 21.31083 20.91156 21.15867 20.92612
    ## 2 20.16151 20.14307 20.56238 21.19937 21.39648 20.51102 20.21887 20.37588
    ## 3 20.47926 20.32799 20.72481 21.13773 21.81802 21.04928 20.65772 20.80578
    ## 4 20.37730 20.33619 20.35912 19.69851 19.85718 20.43477 21.32687 19.99445
    ## 5 21.81922 19.38880 20.23247 19.43306 19.83713 21.58493 20.13333 21.15076
    ##      P9355    P9356    P9357    P9358    P9359   P93510
    ## 1 21.18924 21.13614 21.22950 21.03412 21.14778 21.08201
    ## 2 20.62388 21.17880 21.47880 21.03619 20.72338 20.79518
    ## 3 21.21239 21.52508 21.79645 21.57939 21.14020 20.94360
    ## 4 19.94857 19.99520 19.94448 19.79679 20.16072 19.99660
    ## 5 19.07674 19.92342 20.27656 21.32954 19.87102 20.44066
