MCL_thesis_analysis
================
Heona
2026-02-04

- [MCL 1-100 Data Analysis](#mcl-1-100-data-analysis)
  - [Import the data](#import-the-data)
  - [Quick raw data check](#quick-raw-data-check)

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

### Initial visualiztaion of NAs

``` r
library(tidyverse)

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

``` r
#convert NaN to NA
mcl[mcl == "NaN"] <- NA_integer_

#convert Zero values to NA
mcl[,-c(1:2)][mcl[,-c(1:2)] == "0"] <- NA_integer_

#remove rows with only NA
mcl <- mcl %>%
  filter(rowSums(is.na(mcl[,-c(1:2)])) != ncol(mcl[,-c(1:2)]))

#data check
dim(mcl)
```

    ## [1] 6321  184

``` r
print(paste("Total number of missing values: ", sum(is.na(mcl))))
```

    ## [1] "Total number of missing values:  319128"

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

### Combine replicates & reshape dataframe from wide to long

I combine the replicated to a dataframe with the average of r1 and r2.

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

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
