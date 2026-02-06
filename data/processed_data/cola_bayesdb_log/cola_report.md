cola Report for Consensus Partitioning
==================

**Date**: 2026-02-05 23:03:53 CET, **cola version**: 2.10.0

----------------------------------------------------------------

<style type='text/css'>

body, td, th {
   font-family: Arial,Helvetica,sans-serif;
   background-color: white;
   font-size: 13px;
  max-width: 800px;
  margin: auto;
  margin-left:210px;
  padding: 0px 10px 0px 10px;
  border-left: 1px solid #EEEEEE;
  line-height: 150%;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, 

monospace;
}

h1 {
   font-size:2.2em;
   line-height:150%;
}

h2 {
   font-size:1.8em;
}

h3 {
   font-size:1.4em;
}

h4 {
   font-size:1.0em;
}

h5 {
   font-size:0.9em;
}

h6 {
   font-size:0.8em;
}

a {
  text-decoration: none;
  color: #0366d6;
}

a:hover {
  text-decoration: underline;
}

a:visited {
   color: #0366d6;
}

pre, img {
  max-width: 100%;
}
pre {
  overflow-x: auto;
}
pre code {
   display: block; padding: 0.5em;
}

code {
  font-size: 92%;
  border: 1px solid #ccc;
}

code[class] {
  background-color: #F8F8F8;
}

table, td, th {
  border: 1px solid #ccc;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>




## Summary



First the variable is renamed to `res_list`.


``` r
res_list = rl_bayesdb_log
```



All available functions which can be applied to this `res_list` object:


``` r
res_list
```

```
#> A 'ConsensusPartitionList' object with 9 methods.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows are extracted by 'SD, MAD, ATC' methods.
#>   Subgroups are detected by 'hclust, kmeans, skmeans' method.
#>   Number of partitions are tried for k = 2, 3, 4, 5, 6.
#>   Performed in total 2250 partitions by row resampling.
#> 
#> Following methods can be applied to this 'ConsensusPartitionList' object:
#>  [1] "cola_report"           "collect_classes"       "collect_plots"         "collect_stats"        
#>  [5] "colnames"              "functional_enrichment" "get_anno_col"          "get_anno"             
#>  [9] "get_classes"           "get_matrix"            "get_membership"        "get_stats"            
#> [13] "is_best_k"             "is_stable_k"           "ncol"                  "nrow"                 
#> [17] "rownames"              "show"                  "suggest_best_k"        "test_to_known_factors"
#> [21] "top_rows_heatmap"      "top_rows_overlap"     
#> 
#> You can get result for a single method by, e.g. object["SD", "hclust"] or object["SD:hclust"]
#> or a subset of methods by object[c("SD", "MAD")], c("hclust", "kmeans")]
```

The call of `run_all_consensus_partition_methods()` was:


```
#> run_all_consensus_partition_methods(data = bayesdb_mat_log, top_value_method = c("SD", 
#>     "MAD", "ATC"), partition_method = c("hclust", "kmeans", "skmeans"), max_k = 6, 
#>     cores = 4, scale_rows = TRUE)
```

Dimension of the input matrix:


``` r
mat = get_matrix(res_list)
dim(mat)
```

```
#> [1] 22 72
```

### Density distribution

The density distribution for each sample is visualized as in one column in the
following heatmap. The clustering is based on the distance which is the
Kolmogorov-Smirnov statistic between two distributions.




``` r
library(ComplexHeatmap)
densityHeatmap(mat, ylab = "value", cluster_columns = TRUE, show_column_names = FALSE,
    mc.cores = 4)
```

![](figure_cola/density-heatmap-1.png)<!-- -->





### Suggest the best k



Folowing table shows the best `k` (number of partitions) for each combination
of top-value methods and partitioning methods. Clicking on the method name in
the table goes to the corresponding section for a single combination of methods.

[The cola vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13)
explains the definition of the metrics used for determining the best
number of partitions.


``` r
suggest_best_k(res_list)
```


|            | The best k| 1-PAC| Mean silhouette| Concordance|   |Optional k |
|:-----------|----------:|-----:|---------------:|-----------:|:--|:----------|
|[SD:hclust](#sd-hclust)|          5|     1|           1.000|       1.000|** |2,3        |
|[SD:kmeans](#sd-kmeans)|          5|     1|           0.992|       0.990|** |           |
|[SD:skmeans](#sd-skmeans)|          6|     1|           1.000|       1.000|** |2,3        |
|[MAD:hclust](#mad-hclust)|          5|     1|           1.000|       1.000|** |2,3        |
|[MAD:kmeans](#mad-kmeans)|          5|     1|           1.000|       1.000|** |4          |
|[MAD:skmeans](#mad-skmeans)|          6|     1|           1.000|       1.000|** |2,3        |
|[ATC:hclust](#atc-hclust)|          5|     1|           1.000|       1.000|** |2,3,4      |
|[ATC:kmeans](#atc-kmeans)|          2|     1|           1.000|       1.000|** |           |
|[ATC:skmeans](#atc-skmeans)|          6|     1|           0.962|       0.989|** |2,4        |

\*\*: 1-PAC > 0.95, \*: 1-PAC > 0.9




### CDF of consensus matrices

Cumulative distribution function curves of consensus matrix for all methods.




``` r
collect_plots(res_list, fun = plot_ecdf)
```

![](figure_cola/collect-plots-1.png)<!-- -->



### Consensus heatmap

Consensus heatmaps for all methods. ([What is a consensus heatmap?](https://jokergoo.github.io/cola_vignettes/cola.html#toc_9))





































