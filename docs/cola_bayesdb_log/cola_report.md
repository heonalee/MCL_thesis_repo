cola Report for Consensus Partitioning
==================

**Date**: 2026-02-09 10:10:20 CET, **cola version**: 2.10.0

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

![plot of chunk density-heatmap](figure_cola/density-heatmap-1.png)





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
|[SD:kmeans](#sd-kmeans)|          3|     1|           0.991|       0.987|** |           |
|[SD:skmeans](#sd-skmeans)|          5|     1|           1.000|       1.000|** |2,3        |
|[MAD:hclust](#mad-hclust)|          5|     1|           1.000|       1.000|** |2,3        |
|[MAD:kmeans](#mad-kmeans)|          3|     1|           0.991|       0.987|** |           |
|[MAD:skmeans](#mad-skmeans)|          5|     1|           1.000|       1.000|** |2,3        |
|[ATC:hclust](#atc-hclust)|          5|     1|           1.000|       1.000|** |2,3,4      |
|[ATC:kmeans](#atc-kmeans)|          2|     1|           1.000|       1.000|** |           |
|[ATC:skmeans](#atc-skmeans)|          6|     1|           0.986|       1.000|** |2,4,5      |

\*\*: 1-PAC > 0.95, \*: 1-PAC > 0.9




### CDF of consensus matrices

Cumulative distribution function curves of consensus matrix for all methods.




``` r
collect_plots(res_list, fun = plot_ecdf)
```

![plot of chunk collect-plots](figure_cola/collect-plots-1.png)



### Consensus heatmap

Consensus heatmaps for all methods. ([What is a consensus heatmap?](https://jokergoo.github.io/cola_vignettes/cola.html#toc_9))


<style type='text/css'>



.ui-helper-hidden {
	display: none;
}
.ui-helper-hidden-accessible {
	border: 0;
	clip: rect(0 0 0 0);
	height: 1px;
	margin: -1px;
	overflow: hidden;
	padding: 0;
	position: absolute;
	width: 1px;
}
.ui-helper-reset {
	margin: 0;
	padding: 0;
	border: 0;
	outline: 0;
	line-height: 1.3;
	text-decoration: none;
	font-size: 100%;
	list-style: none;
}
.ui-helper-clearfix:before,
.ui-helper-clearfix:after {
	content: "";
	display: table;
	border-collapse: collapse;
}
.ui-helper-clearfix:after {
	clear: both;
}
.ui-helper-zfix {
	width: 100%;
	height: 100%;
	top: 0;
	left: 0;
	position: absolute;
	opacity: 0;
	filter:Alpha(Opacity=0); 
}

.ui-front {
	z-index: 100;
}



.ui-state-disabled {
	cursor: default !important;
	pointer-events: none;
}



.ui-icon {
	display: inline-block;
	vertical-align: middle;
	margin-top: -.25em;
	position: relative;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
}

.ui-widget-icon-block {
	left: 50%;
	margin-left: -8px;
	display: block;
}




.ui-widget-overlay {
	position: fixed;
	top: 0;
	left: 0;
	width: 100%;
	height: 100%;
}
.ui-accordion .ui-accordion-header {
	display: block;
	cursor: pointer;
	position: relative;
	margin: 2px 0 0 0;
	padding: .5em .5em .5em .7em;
	font-size: 100%;
}
.ui-accordion .ui-accordion-content {
	padding: 1em 2.2em;
	border-top: 0;
	overflow: auto;
}
.ui-autocomplete {
	position: absolute;
	top: 0;
	left: 0;
	cursor: default;
}
.ui-menu {
	list-style: none;
	padding: 0;
	margin: 0;
	display: block;
	outline: 0;
}
.ui-menu .ui-menu {
	position: absolute;
}
.ui-menu .ui-menu-item {
	margin: 0;
	cursor: pointer;
	
	list-style-image: url("data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7");
}
.ui-menu .ui-menu-item-wrapper {
	position: relative;
	padding: 3px 1em 3px .4em;
}
.ui-menu .ui-menu-divider {
	margin: 5px 0;
	height: 0;
	font-size: 0;
	line-height: 0;
	border-width: 1px 0 0 0;
}
.ui-menu .ui-state-focus,
.ui-menu .ui-state-active {
	margin: -1px;
}


.ui-menu-icons {
	position: relative;
}
.ui-menu-icons .ui-menu-item-wrapper {
	padding-left: 2em;
}


.ui-menu .ui-icon {
	position: absolute;
	top: 0;
	bottom: 0;
	left: .2em;
	margin: auto 0;
}


.ui-menu .ui-menu-icon {
	left: auto;
	right: 0;
}
.ui-button {
	padding: .4em 1em;
	display: inline-block;
	position: relative;
	line-height: normal;
	margin-right: .1em;
	cursor: pointer;
	vertical-align: middle;
	text-align: center;
	-webkit-user-select: none;
	-moz-user-select: none;
	-ms-user-select: none;
	user-select: none;

	
	overflow: visible;
}

.ui-button,
.ui-button:link,
.ui-button:visited,
.ui-button:hover,
.ui-button:active {
	text-decoration: none;
}


.ui-button-icon-only {
	width: 2em;
	box-sizing: border-box;
	text-indent: -9999px;
	white-space: nowrap;
}


input.ui-button.ui-button-icon-only {
	text-indent: 0;
}


.ui-button-icon-only .ui-icon {
	position: absolute;
	top: 50%;
	left: 50%;
	margin-top: -8px;
	margin-left: -8px;
}

.ui-button.ui-icon-notext .ui-icon {
	padding: 0;
	width: 2.1em;
	height: 2.1em;
	text-indent: -9999px;
	white-space: nowrap;

}

input.ui-button.ui-icon-notext .ui-icon {
	width: auto;
	height: auto;
	text-indent: 0;
	white-space: normal;
	padding: .4em 1em;
}



input.ui-button::-moz-focus-inner,
button.ui-button::-moz-focus-inner {
	border: 0;
	padding: 0;
}
.ui-controlgroup {
	vertical-align: middle;
	display: inline-block;
}
.ui-controlgroup > .ui-controlgroup-item {
	float: left;
	margin-left: 0;
	margin-right: 0;
}
.ui-controlgroup > .ui-controlgroup-item:focus,
.ui-controlgroup > .ui-controlgroup-item.ui-visual-focus {
	z-index: 9999;
}
.ui-controlgroup-vertical > .ui-controlgroup-item {
	display: block;
	float: none;
	width: 100%;
	margin-top: 0;
	margin-bottom: 0;
	text-align: left;
}
.ui-controlgroup-vertical .ui-controlgroup-item {
	box-sizing: border-box;
}
.ui-controlgroup .ui-controlgroup-label {
	padding: .4em 1em;
}
.ui-controlgroup .ui-controlgroup-label span {
	font-size: 80%;
}
.ui-controlgroup-horizontal .ui-controlgroup-label + .ui-controlgroup-item {
	border-left: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label + .ui-controlgroup-item {
	border-top: none;
}
.ui-controlgroup-horizontal .ui-controlgroup-label.ui-widget-content {
	border-right: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label.ui-widget-content {
	border-bottom: none;
}


.ui-controlgroup-vertical .ui-spinner-input {

	
	width: 75%;
	width: calc( 100% - 2.4em );
}
.ui-controlgroup-vertical .ui-spinner .ui-spinner-up {
	border-top-style: solid;
}

.ui-checkboxradio-label .ui-icon-background {
	box-shadow: inset 1px 1px 1px #ccc;
	border-radius: .12em;
	border: none;
}
.ui-checkboxradio-radio-label .ui-icon-background {
	width: 16px;
	height: 16px;
	border-radius: 1em;
	overflow: visible;
	border: none;
}
.ui-checkboxradio-radio-label.ui-checkboxradio-checked .ui-icon,
.ui-checkboxradio-radio-label.ui-checkboxradio-checked:hover .ui-icon {
	background-image: none;
	width: 8px;
	height: 8px;
	border-width: 4px;
	border-style: solid;
}
.ui-checkboxradio-disabled {
	pointer-events: none;
}
.ui-datepicker {
	width: 17em;
	padding: .2em .2em 0;
	display: none;
}
.ui-datepicker .ui-datepicker-header {
	position: relative;
	padding: .2em 0;
}
.ui-datepicker .ui-datepicker-prev,
.ui-datepicker .ui-datepicker-next {
	position: absolute;
	top: 2px;
	width: 1.8em;
	height: 1.8em;
}
.ui-datepicker .ui-datepicker-prev-hover,
.ui-datepicker .ui-datepicker-next-hover {
	top: 1px;
}
.ui-datepicker .ui-datepicker-prev {
	left: 2px;
}
.ui-datepicker .ui-datepicker-next {
	right: 2px;
}
.ui-datepicker .ui-datepicker-prev-hover {
	left: 1px;
}
.ui-datepicker .ui-datepicker-next-hover {
	right: 1px;
}
.ui-datepicker .ui-datepicker-prev span,
.ui-datepicker .ui-datepicker-next span {
	display: block;
	position: absolute;
	left: 50%;
	margin-left: -8px;
	top: 50%;
	margin-top: -8px;
}
.ui-datepicker .ui-datepicker-title {
	margin: 0 2.3em;
	line-height: 1.8em;
	text-align: center;
}
.ui-datepicker .ui-datepicker-title select {
	font-size: 1em;
	margin: 1px 0;
}
.ui-datepicker select.ui-datepicker-month,
.ui-datepicker select.ui-datepicker-year {
	width: 45%;
}
.ui-datepicker table {
	width: 100%;
	font-size: .9em;
	border-collapse: collapse;
	margin: 0 0 .4em;
}
.ui-datepicker th {
	padding: .7em .3em;
	text-align: center;
	font-weight: bold;
	border: 0;
}
.ui-datepicker td {
	border: 0;
	padding: 1px;
}
.ui-datepicker td span,
.ui-datepicker td a {
	display: block;
	padding: .2em;
	text-align: right;
	text-decoration: none;
}
.ui-datepicker .ui-datepicker-buttonpane {
	background-image: none;
	margin: .7em 0 0 0;
	padding: 0 .2em;
	border-left: 0;
	border-right: 0;
	border-bottom: 0;
}
.ui-datepicker .ui-datepicker-buttonpane button {
	float: right;
	margin: .5em .2em .4em;
	cursor: pointer;
	padding: .2em .6em .3em .6em;
	width: auto;
	overflow: visible;
}
.ui-datepicker .ui-datepicker-buttonpane button.ui-datepicker-current {
	float: left;
}


.ui-datepicker.ui-datepicker-multi {
	width: auto;
}
.ui-datepicker-multi .ui-datepicker-group {
	float: left;
}
.ui-datepicker-multi .ui-datepicker-group table {
	width: 95%;
	margin: 0 auto .4em;
}
.ui-datepicker-multi-2 .ui-datepicker-group {
	width: 50%;
}
.ui-datepicker-multi-3 .ui-datepicker-group {
	width: 33.3%;
}
.ui-datepicker-multi-4 .ui-datepicker-group {
	width: 25%;
}
.ui-datepicker-multi .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-multi .ui-datepicker-group-middle .ui-datepicker-header {
	border-left-width: 0;
}
.ui-datepicker-multi .ui-datepicker-buttonpane {
	clear: left;
}
.ui-datepicker-row-break {
	clear: both;
	width: 100%;
	font-size: 0;
}


.ui-datepicker-rtl {
	direction: rtl;
}
.ui-datepicker-rtl .ui-datepicker-prev {
	right: 2px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next {
	left: 2px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-prev:hover {
	right: 1px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next:hover {
	left: 1px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane {
	clear: right;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button {
	float: left;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button.ui-datepicker-current,
.ui-datepicker-rtl .ui-datepicker-group {
	float: right;
}
.ui-datepicker-rtl .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-rtl .ui-datepicker-group-middle .ui-datepicker-header {
	border-right-width: 0;
	border-left-width: 1px;
}


.ui-datepicker .ui-icon {
	display: block;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
	left: .5em;
	top: .3em;
}
.ui-dialog {
	position: absolute;
	top: 0;
	left: 0;
	padding: .2em;
	outline: 0;
}
.ui-dialog .ui-dialog-titlebar {
	padding: .4em 1em;
	position: relative;
}
.ui-dialog .ui-dialog-title {
	float: left;
	margin: .1em 0;
	white-space: nowrap;
	width: 90%;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-dialog .ui-dialog-titlebar-close {
	position: absolute;
	right: .3em;
	top: 50%;
	width: 20px;
	margin: -10px 0 0 0;
	padding: 1px;
	height: 20px;
}
.ui-dialog .ui-dialog-content {
	position: relative;
	border: 0;
	padding: .5em 1em;
	background: none;
	overflow: auto;
}
.ui-dialog .ui-dialog-buttonpane {
	text-align: left;
	border-width: 1px 0 0 0;
	background-image: none;
	margin-top: .5em;
	padding: .3em 1em .5em .4em;
}
.ui-dialog .ui-dialog-buttonpane .ui-dialog-buttonset {
	float: right;
}
.ui-dialog .ui-dialog-buttonpane button {
	margin: .5em .4em .5em 0;
	cursor: pointer;
}
.ui-dialog .ui-resizable-n {
	height: 2px;
	top: 0;
}
.ui-dialog .ui-resizable-e {
	width: 2px;
	right: 0;
}
.ui-dialog .ui-resizable-s {
	height: 2px;
	bottom: 0;
}
.ui-dialog .ui-resizable-w {
	width: 2px;
	left: 0;
}
.ui-dialog .ui-resizable-se,
.ui-dialog .ui-resizable-sw,
.ui-dialog .ui-resizable-ne,
.ui-dialog .ui-resizable-nw {
	width: 7px;
	height: 7px;
}
.ui-dialog .ui-resizable-se {
	right: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-sw {
	left: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-ne {
	right: 0;
	top: 0;
}
.ui-dialog .ui-resizable-nw {
	left: 0;
	top: 0;
}
.ui-draggable .ui-dialog-titlebar {
	cursor: move;
}
.ui-draggable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable {
	position: relative;
}
.ui-resizable-handle {
	position: absolute;
	font-size: 0.1px;
	display: block;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable-disabled .ui-resizable-handle,
.ui-resizable-autohide .ui-resizable-handle {
	display: none;
}
.ui-resizable-n {
	cursor: n-resize;
	height: 7px;
	width: 100%;
	top: -5px;
	left: 0;
}
.ui-resizable-s {
	cursor: s-resize;
	height: 7px;
	width: 100%;
	bottom: -5px;
	left: 0;
}
.ui-resizable-e {
	cursor: e-resize;
	width: 7px;
	right: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-w {
	cursor: w-resize;
	width: 7px;
	left: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-se {
	cursor: se-resize;
	width: 12px;
	height: 12px;
	right: 1px;
	bottom: 1px;
}
.ui-resizable-sw {
	cursor: sw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	bottom: -5px;
}
.ui-resizable-nw {
	cursor: nw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	top: -5px;
}
.ui-resizable-ne {
	cursor: ne-resize;
	width: 9px;
	height: 9px;
	right: -5px;
	top: -5px;
}
.ui-progressbar {
	height: 2em;
	text-align: left;
	overflow: hidden;
}
.ui-progressbar .ui-progressbar-value {
	margin: -1px;
	height: 100%;
}
.ui-progressbar .ui-progressbar-overlay {
	background: url("data:image/gif;base64,R0lGODlhKAAoAIABAAAAAP///yH/C05FVFNDQVBFMi4wAwEAAAAh+QQJAQABACwAAAAAKAAoAAACkYwNqXrdC52DS06a7MFZI+4FHBCKoDeWKXqymPqGqxvJrXZbMx7Ttc+w9XgU2FB3lOyQRWET2IFGiU9m1frDVpxZZc6bfHwv4c1YXP6k1Vdy292Fb6UkuvFtXpvWSzA+HycXJHUXiGYIiMg2R6W459gnWGfHNdjIqDWVqemH2ekpObkpOlppWUqZiqr6edqqWQAAIfkECQEAAQAsAAAAACgAKAAAApSMgZnGfaqcg1E2uuzDmmHUBR8Qil95hiPKqWn3aqtLsS18y7G1SzNeowWBENtQd+T1JktP05nzPTdJZlR6vUxNWWjV+vUWhWNkWFwxl9VpZRedYcflIOLafaa28XdsH/ynlcc1uPVDZxQIR0K25+cICCmoqCe5mGhZOfeYSUh5yJcJyrkZWWpaR8doJ2o4NYq62lAAACH5BAkBAAEALAAAAAAoACgAAAKVDI4Yy22ZnINRNqosw0Bv7i1gyHUkFj7oSaWlu3ovC8GxNso5fluz3qLVhBVeT/Lz7ZTHyxL5dDalQWPVOsQWtRnuwXaFTj9jVVh8pma9JjZ4zYSj5ZOyma7uuolffh+IR5aW97cHuBUXKGKXlKjn+DiHWMcYJah4N0lYCMlJOXipGRr5qdgoSTrqWSq6WFl2ypoaUAAAIfkECQEAAQAsAAAAACgAKAAAApaEb6HLgd/iO7FNWtcFWe+ufODGjRfoiJ2akShbueb0wtI50zm02pbvwfWEMWBQ1zKGlLIhskiEPm9R6vRXxV4ZzWT2yHOGpWMyorblKlNp8HmHEb/lCXjcW7bmtXP8Xt229OVWR1fod2eWqNfHuMjXCPkIGNileOiImVmCOEmoSfn3yXlJWmoHGhqp6ilYuWYpmTqKUgAAIfkECQEAAQAsAAAAACgAKAAAApiEH6kb58biQ3FNWtMFWW3eNVcojuFGfqnZqSebuS06w5V80/X02pKe8zFwP6EFWOT1lDFk8rGERh1TTNOocQ61Hm4Xm2VexUHpzjymViHrFbiELsefVrn6XKfnt2Q9G/+Xdie499XHd2g4h7ioOGhXGJboGAnXSBnoBwKYyfioubZJ2Hn0RuRZaflZOil56Zp6iioKSXpUAAAh+QQJAQABACwAAAAAKAAoAAACkoQRqRvnxuI7kU1a1UU5bd5tnSeOZXhmn5lWK3qNTWvRdQxP8qvaC+/yaYQzXO7BMvaUEmJRd3TsiMAgswmNYrSgZdYrTX6tSHGZO73ezuAw2uxuQ+BbeZfMxsexY35+/Qe4J1inV0g4x3WHuMhIl2jXOKT2Q+VU5fgoSUI52VfZyfkJGkha6jmY+aaYdirq+lQAACH5BAkBAAEALAAAAAAoACgAAAKWBIKpYe0L3YNKToqswUlvznigd4wiR4KhZrKt9Upqip61i9E3vMvxRdHlbEFiEXfk9YARYxOZZD6VQ2pUunBmtRXo1Lf8hMVVcNl8JafV38aM2/Fu5V16Bn63r6xt97j09+MXSFi4BniGFae3hzbH9+hYBzkpuUh5aZmHuanZOZgIuvbGiNeomCnaxxap2upaCZsq+1kAACH5BAkBAAEALAAAAAAoACgAAAKXjI8By5zf4kOxTVrXNVlv1X0d8IGZGKLnNpYtm8Lr9cqVeuOSvfOW79D9aDHizNhDJidFZhNydEahOaDH6nomtJjp1tutKoNWkvA6JqfRVLHU/QUfau9l2x7G54d1fl995xcIGAdXqMfBNadoYrhH+Mg2KBlpVpbluCiXmMnZ2Sh4GBqJ+ckIOqqJ6LmKSllZmsoq6wpQAAAh+QQJAQABACwAAAAAKAAoAAAClYx/oLvoxuJDkU1a1YUZbJ59nSd2ZXhWqbRa2/gF8Gu2DY3iqs7yrq+xBYEkYvFSM8aSSObE+ZgRl1BHFZNr7pRCavZ5BW2142hY3AN/zWtsmf12p9XxxFl2lpLn1rseztfXZjdIWIf2s5dItwjYKBgo9yg5pHgzJXTEeGlZuenpyPmpGQoKOWkYmSpaSnqKileI2FAAACH5BAkBAAEALAAAAAAoACgAAAKVjB+gu+jG4kORTVrVhRlsnn2dJ3ZleFaptFrb+CXmO9OozeL5VfP99HvAWhpiUdcwkpBH3825AwYdU8xTqlLGhtCosArKMpvfa1mMRae9VvWZfeB2XfPkeLmm18lUcBj+p5dnN8jXZ3YIGEhYuOUn45aoCDkp16hl5IjYJvjWKcnoGQpqyPlpOhr3aElaqrq56Bq7VAAAOw==");
	height: 100%;
	filter: alpha(opacity=25); 
	opacity: 0.25;
}
.ui-progressbar-indeterminate .ui-progressbar-value {
	background-image: none;
}
.ui-selectable {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-selectable-helper {
	position: absolute;
	z-index: 100;
	border: 1px dotted black;
}
.ui-selectmenu-menu {
	padding: 0;
	margin: 0;
	position: absolute;
	top: 0;
	left: 0;
	display: none;
}
.ui-selectmenu-menu .ui-menu {
	overflow: auto;
	overflow-x: hidden;
	padding-bottom: 1px;
}
.ui-selectmenu-menu .ui-menu .ui-selectmenu-optgroup {
	font-size: 1em;
	font-weight: bold;
	line-height: 1.5;
	padding: 2px 0.4em;
	margin: 0.5em 0 0 0;
	height: auto;
	border: 0;
}
.ui-selectmenu-open {
	display: block;
}
.ui-selectmenu-text {
	display: block;
	margin-right: 20px;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-selectmenu-button.ui-button {
	text-align: left;
	white-space: nowrap;
	width: 14em;
}
.ui-selectmenu-icon.ui-icon {
	float: right;
	margin-top: 0;
}
.ui-slider {
	position: relative;
	text-align: left;
}
.ui-slider .ui-slider-handle {
	position: absolute;
	z-index: 2;
	width: 1.2em;
	height: 1.2em;
	cursor: default;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-slider .ui-slider-range {
	position: absolute;
	z-index: 1;
	font-size: .7em;
	display: block;
	border: 0;
	background-position: 0 0;
}


.ui-slider.ui-state-disabled .ui-slider-handle,
.ui-slider.ui-state-disabled .ui-slider-range {
	filter: inherit;
}

.ui-slider-horizontal {
	height: .8em;
}
.ui-slider-horizontal .ui-slider-handle {
	top: -.3em;
	margin-left: -.6em;
}
.ui-slider-horizontal .ui-slider-range {
	top: 0;
	height: 100%;
}
.ui-slider-horizontal .ui-slider-range-min {
	left: 0;
}
.ui-slider-horizontal .ui-slider-range-max {
	right: 0;
}

.ui-slider-vertical {
	width: .8em;
	height: 100px;
}
.ui-slider-vertical .ui-slider-handle {
	left: -.3em;
	margin-left: 0;
	margin-bottom: -.6em;
}
.ui-slider-vertical .ui-slider-range {
	left: 0;
	width: 100%;
}
.ui-slider-vertical .ui-slider-range-min {
	bottom: 0;
}
.ui-slider-vertical .ui-slider-range-max {
	top: 0;
}
.ui-sortable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-spinner {
	position: relative;
	display: inline-block;
	overflow: hidden;
	padding: 0;
	vertical-align: middle;
}
.ui-spinner-input {
	border: none;
	background: none;
	color: inherit;
	padding: .222em 0;
	margin: .2em 0;
	vertical-align: middle;
	margin-left: .4em;
	margin-right: 2em;
}
.ui-spinner-button {
	width: 1.6em;
	height: 50%;
	font-size: .5em;
	padding: 0;
	margin: 0;
	text-align: center;
	position: absolute;
	cursor: default;
	display: block;
	overflow: hidden;
	right: 0;
}

.ui-spinner a.ui-spinner-button {
	border-top-style: none;
	border-bottom-style: none;
	border-right-style: none;
}
.ui-spinner-up {
	top: 0;
}
.ui-spinner-down {
	bottom: 0;
}
.ui-tabs {
	position: relative;
	padding: .2em;
}
.ui-tabs .ui-tabs-nav {
	margin: 0;
	padding: .2em .2em 0;
}
.ui-tabs .ui-tabs-nav li {
	list-style: none;
	float: left;
	position: relative;
	top: 0;
	margin: 1px .2em 0 0;
	border-bottom-width: 0;
	padding: 0;
	white-space: nowrap;
}
.ui-tabs .ui-tabs-nav .ui-tabs-anchor {
	float: left;
	padding: .5em 1em;
	text-decoration: none;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active {
	margin-bottom: -1px;
	padding-bottom: 1px;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-state-disabled .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-tabs-loading .ui-tabs-anchor {
	cursor: text;
}
.ui-tabs-collapsible .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor {
	cursor: pointer;
}
.ui-tabs .ui-tabs-panel {
	display: block;
	border-width: 0;
	padding: 1em 1.4em;
	background: none;
}
.ui-tooltip {
	padding: 8px;
	position: absolute;
	z-index: 9999;
	max-width: 300px;
}
body .ui-tooltip {
	border-width: 2px;
}

.ui-widget {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget .ui-widget {
	font-size: 1em;
}
.ui-widget input,
.ui-widget select,
.ui-widget textarea,
.ui-widget button {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget.ui-widget-content {
	border: 1px solid #c5c5c5;
}
.ui-widget-content {
	border: 1px solid #dddddd;
	background: #ffffff;
	color: #333333;
}
.ui-widget-content a {
	color: #333333;
}
.ui-widget-header {
	border: 1px solid #dddddd;
	background: #e9e9e9;
	color: #333333;
	font-weight: bold;
}
.ui-widget-header a {
	color: #333333;
}


.ui-state-default,
.ui-widget-content .ui-state-default,
.ui-widget-header .ui-state-default,
.ui-button,


html .ui-button.ui-state-disabled:hover,
html .ui-button.ui-state-disabled:active {
	border: 1px solid #c5c5c5;
	background: #f6f6f6;
	font-weight: normal;
	color: #454545;
}
.ui-state-default a,
.ui-state-default a:link,
.ui-state-default a:visited,
a.ui-button,
a:link.ui-button,
a:visited.ui-button,
.ui-button {
	color: #454545;
	text-decoration: none;
}
.ui-state-hover,
.ui-widget-content .ui-state-hover,
.ui-widget-header .ui-state-hover,
.ui-state-focus,
.ui-widget-content .ui-state-focus,
.ui-widget-header .ui-state-focus,
.ui-button:hover,
.ui-button:focus {
	border: 1px solid #cccccc;
	background: #ededed;
	font-weight: normal;
	color: #2b2b2b;
}
.ui-state-hover a,
.ui-state-hover a:hover,
.ui-state-hover a:link,
.ui-state-hover a:visited,
.ui-state-focus a,
.ui-state-focus a:hover,
.ui-state-focus a:link,
.ui-state-focus a:visited,
a.ui-button:hover,
a.ui-button:focus {
	color: #2b2b2b;
	text-decoration: none;
}

.ui-visual-focus {
	box-shadow: 0 0 3px 1px rgb(94, 158, 214);
}
.ui-state-active,
.ui-widget-content .ui-state-active,
.ui-widget-header .ui-state-active,
a.ui-button:active,
.ui-button:active,
.ui-button.ui-state-active:hover {
	border: 1px solid #003eff;
	background: #007fff;
	font-weight: normal;
	color: #ffffff;
}
.ui-icon-background,
.ui-state-active .ui-icon-background {
	border: #003eff;
	background-color: #ffffff;
}
.ui-state-active a,
.ui-state-active a:link,
.ui-state-active a:visited {
	color: #ffffff;
	text-decoration: none;
}


.ui-state-highlight,
.ui-widget-content .ui-state-highlight,
.ui-widget-header .ui-state-highlight {
	border: 1px solid #dad55e;
	background: #fffa90;
	color: #777620;
}
.ui-state-checked {
	border: 1px solid #dad55e;
	background: #fffa90;
}
.ui-state-highlight a,
.ui-widget-content .ui-state-highlight a,
.ui-widget-header .ui-state-highlight a {
	color: #777620;
}
.ui-state-error,
.ui-widget-content .ui-state-error,
.ui-widget-header .ui-state-error {
	border: 1px solid #f1a899;
	background: #fddfdf;
	color: #5f3f3f;
}
.ui-state-error a,
.ui-widget-content .ui-state-error a,
.ui-widget-header .ui-state-error a {
	color: #5f3f3f;
}
.ui-state-error-text,
.ui-widget-content .ui-state-error-text,
.ui-widget-header .ui-state-error-text {
	color: #5f3f3f;
}
.ui-priority-primary,
.ui-widget-content .ui-priority-primary,
.ui-widget-header .ui-priority-primary {
	font-weight: bold;
}
.ui-priority-secondary,
.ui-widget-content .ui-priority-secondary,
.ui-widget-header .ui-priority-secondary {
	opacity: .7;
	filter:Alpha(Opacity=70); 
	font-weight: normal;
}
.ui-state-disabled,
.ui-widget-content .ui-state-disabled,
.ui-widget-header .ui-state-disabled {
	opacity: .35;
	filter:Alpha(Opacity=35); 
	background-image: none;
}
.ui-state-disabled .ui-icon {
	filter:Alpha(Opacity=35); 
}




.ui-icon {
	width: 16px;
	height: 16px;
}
.ui-icon,
.ui-widget-content .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-widget-header .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-state-hover .ui-icon,
.ui-state-focus .ui-icon,
.ui-button:hover .ui-icon,
.ui-button:focus .ui-icon {
	background-image: url("images/ui-icons_555555_256x240.png");
}
.ui-state-active .ui-icon,
.ui-button:active .ui-icon {
	background-image: url("images/ui-icons_ffffff_256x240.png");
}
.ui-state-highlight .ui-icon,
.ui-button .ui-state-highlight.ui-icon {
	background-image: url("images/ui-icons_777620_256x240.png");
}
.ui-state-error .ui-icon,
.ui-state-error-text .ui-icon {
	background-image: url("images/ui-icons_cc0000_256x240.png");
}
.ui-button .ui-icon {
	background-image: url("images/ui-icons_777777_256x240.png");
}


.ui-icon-blank { background-position: 16px 16px; }
.ui-icon-caret-1-n { background-position: 0 0; }
.ui-icon-caret-1-ne { background-position: -16px 0; }
.ui-icon-caret-1-e { background-position: -32px 0; }
.ui-icon-caret-1-se { background-position: -48px 0; }
.ui-icon-caret-1-s { background-position: -65px 0; }
.ui-icon-caret-1-sw { background-position: -80px 0; }
.ui-icon-caret-1-w { background-position: -96px 0; }
.ui-icon-caret-1-nw { background-position: -112px 0; }
.ui-icon-caret-2-n-s { background-position: -128px 0; }
.ui-icon-caret-2-e-w { background-position: -144px 0; }
.ui-icon-triangle-1-n { background-position: 0 -16px; }
.ui-icon-triangle-1-ne { background-position: -16px -16px; }
.ui-icon-triangle-1-e { background-position: -32px -16px; }
.ui-icon-triangle-1-se { background-position: -48px -16px; }
.ui-icon-triangle-1-s { background-position: -65px -16px; }
.ui-icon-triangle-1-sw { background-position: -80px -16px; }
.ui-icon-triangle-1-w { background-position: -96px -16px; }
.ui-icon-triangle-1-nw { background-position: -112px -16px; }
.ui-icon-triangle-2-n-s { background-position: -128px -16px; }
.ui-icon-triangle-2-e-w { background-position: -144px -16px; }
.ui-icon-arrow-1-n { background-position: 0 -32px; }
.ui-icon-arrow-1-ne { background-position: -16px -32px; }
.ui-icon-arrow-1-e { background-position: -32px -32px; }
.ui-icon-arrow-1-se { background-position: -48px -32px; }
.ui-icon-arrow-1-s { background-position: -65px -32px; }
.ui-icon-arrow-1-sw { background-position: -80px -32px; }
.ui-icon-arrow-1-w { background-position: -96px -32px; }
.ui-icon-arrow-1-nw { background-position: -112px -32px; }
.ui-icon-arrow-2-n-s { background-position: -128px -32px; }
.ui-icon-arrow-2-ne-sw { background-position: -144px -32px; }
.ui-icon-arrow-2-e-w { background-position: -160px -32px; }
.ui-icon-arrow-2-se-nw { background-position: -176px -32px; }
.ui-icon-arrowstop-1-n { background-position: -192px -32px; }
.ui-icon-arrowstop-1-e { background-position: -208px -32px; }
.ui-icon-arrowstop-1-s { background-position: -224px -32px; }
.ui-icon-arrowstop-1-w { background-position: -240px -32px; }
.ui-icon-arrowthick-1-n { background-position: 1px -48px; }
.ui-icon-arrowthick-1-ne { background-position: -16px -48px; }
.ui-icon-arrowthick-1-e { background-position: -32px -48px; }
.ui-icon-arrowthick-1-se { background-position: -48px -48px; }
.ui-icon-arrowthick-1-s { background-position: -64px -48px; }
.ui-icon-arrowthick-1-sw { background-position: -80px -48px; }
.ui-icon-arrowthick-1-w { background-position: -96px -48px; }
.ui-icon-arrowthick-1-nw { background-position: -112px -48px; }
.ui-icon-arrowthick-2-n-s { background-position: -128px -48px; }
.ui-icon-arrowthick-2-ne-sw { background-position: -144px -48px; }
.ui-icon-arrowthick-2-e-w { background-position: -160px -48px; }
.ui-icon-arrowthick-2-se-nw { background-position: -176px -48px; }
.ui-icon-arrowthickstop-1-n { background-position: -192px -48px; }
.ui-icon-arrowthickstop-1-e { background-position: -208px -48px; }
.ui-icon-arrowthickstop-1-s { background-position: -224px -48px; }
.ui-icon-arrowthickstop-1-w { background-position: -240px -48px; }
.ui-icon-arrowreturnthick-1-w { background-position: 0 -64px; }
.ui-icon-arrowreturnthick-1-n { background-position: -16px -64px; }
.ui-icon-arrowreturnthick-1-e { background-position: -32px -64px; }
.ui-icon-arrowreturnthick-1-s { background-position: -48px -64px; }
.ui-icon-arrowreturn-1-w { background-position: -64px -64px; }
.ui-icon-arrowreturn-1-n { background-position: -80px -64px; }
.ui-icon-arrowreturn-1-e { background-position: -96px -64px; }
.ui-icon-arrowreturn-1-s { background-position: -112px -64px; }
.ui-icon-arrowrefresh-1-w { background-position: -128px -64px; }
.ui-icon-arrowrefresh-1-n { background-position: -144px -64px; }
.ui-icon-arrowrefresh-1-e { background-position: -160px -64px; }
.ui-icon-arrowrefresh-1-s { background-position: -176px -64px; }
.ui-icon-arrow-4 { background-position: 0 -80px; }
.ui-icon-arrow-4-diag { background-position: -16px -80px; }
.ui-icon-extlink { background-position: -32px -80px; }
.ui-icon-newwin { background-position: -48px -80px; }
.ui-icon-refresh { background-position: -64px -80px; }
.ui-icon-shuffle { background-position: -80px -80px; }
.ui-icon-transfer-e-w { background-position: -96px -80px; }
.ui-icon-transferthick-e-w { background-position: -112px -80px; }
.ui-icon-folder-collapsed { background-position: 0 -96px; }
.ui-icon-folder-open { background-position: -16px -96px; }
.ui-icon-document { background-position: -32px -96px; }
.ui-icon-document-b { background-position: -48px -96px; }
.ui-icon-note { background-position: -64px -96px; }
.ui-icon-mail-closed { background-position: -80px -96px; }
.ui-icon-mail-open { background-position: -96px -96px; }
.ui-icon-suitcase { background-position: -112px -96px; }
.ui-icon-comment { background-position: -128px -96px; }
.ui-icon-person { background-position: -144px -96px; }
.ui-icon-print { background-position: -160px -96px; }
.ui-icon-trash { background-position: -176px -96px; }
.ui-icon-locked { background-position: -192px -96px; }
.ui-icon-unlocked { background-position: -208px -96px; }
.ui-icon-bookmark { background-position: -224px -96px; }
.ui-icon-tag { background-position: -240px -96px; }
.ui-icon-home { background-position: 0 -112px; }
.ui-icon-flag { background-position: -16px -112px; }
.ui-icon-calendar { background-position: -32px -112px; }
.ui-icon-cart { background-position: -48px -112px; }
.ui-icon-pencil { background-position: -64px -112px; }
.ui-icon-clock { background-position: -80px -112px; }
.ui-icon-disk { background-position: -96px -112px; }
.ui-icon-calculator { background-position: -112px -112px; }
.ui-icon-zoomin { background-position: -128px -112px; }
.ui-icon-zoomout { background-position: -144px -112px; }
.ui-icon-search { background-position: -160px -112px; }
.ui-icon-wrench { background-position: -176px -112px; }
.ui-icon-gear { background-position: -192px -112px; }
.ui-icon-heart { background-position: -208px -112px; }
.ui-icon-star { background-position: -224px -112px; }
.ui-icon-link { background-position: -240px -112px; }
.ui-icon-cancel { background-position: 0 -128px; }
.ui-icon-plus { background-position: -16px -128px; }
.ui-icon-plusthick { background-position: -32px -128px; }
.ui-icon-minus { background-position: -48px -128px; }
.ui-icon-minusthick { background-position: -64px -128px; }
.ui-icon-close { background-position: -80px -128px; }
.ui-icon-closethick { background-position: -96px -128px; }
.ui-icon-key { background-position: -112px -128px; }
.ui-icon-lightbulb { background-position: -128px -128px; }
.ui-icon-scissors { background-position: -144px -128px; }
.ui-icon-clipboard { background-position: -160px -128px; }
.ui-icon-copy { background-position: -176px -128px; }
.ui-icon-contact { background-position: -192px -128px; }
.ui-icon-image { background-position: -208px -128px; }
.ui-icon-video { background-position: -224px -128px; }
.ui-icon-script { background-position: -240px -128px; }
.ui-icon-alert { background-position: 0 -144px; }
.ui-icon-info { background-position: -16px -144px; }
.ui-icon-notice { background-position: -32px -144px; }
.ui-icon-help { background-position: -48px -144px; }
.ui-icon-check { background-position: -64px -144px; }
.ui-icon-bullet { background-position: -80px -144px; }
.ui-icon-radio-on { background-position: -96px -144px; }
.ui-icon-radio-off { background-position: -112px -144px; }
.ui-icon-pin-w { background-position: -128px -144px; }
.ui-icon-pin-s { background-position: -144px -144px; }
.ui-icon-play { background-position: 0 -160px; }
.ui-icon-pause { background-position: -16px -160px; }
.ui-icon-seek-next { background-position: -32px -160px; }
.ui-icon-seek-prev { background-position: -48px -160px; }
.ui-icon-seek-end { background-position: -64px -160px; }
.ui-icon-seek-start { background-position: -80px -160px; }

.ui-icon-seek-first { background-position: -80px -160px; }
.ui-icon-stop { background-position: -96px -160px; }
.ui-icon-eject { background-position: -112px -160px; }
.ui-icon-volume-off { background-position: -128px -160px; }
.ui-icon-volume-on { background-position: -144px -160px; }
.ui-icon-power { background-position: 0 -176px; }
.ui-icon-signal-diag { background-position: -16px -176px; }
.ui-icon-signal { background-position: -32px -176px; }
.ui-icon-battery-0 { background-position: -48px -176px; }
.ui-icon-battery-1 { background-position: -64px -176px; }
.ui-icon-battery-2 { background-position: -80px -176px; }
.ui-icon-battery-3 { background-position: -96px -176px; }
.ui-icon-circle-plus { background-position: 0 -192px; }
.ui-icon-circle-minus { background-position: -16px -192px; }
.ui-icon-circle-close { background-position: -32px -192px; }
.ui-icon-circle-triangle-e { background-position: -48px -192px; }
.ui-icon-circle-triangle-s { background-position: -64px -192px; }
.ui-icon-circle-triangle-w { background-position: -80px -192px; }
.ui-icon-circle-triangle-n { background-position: -96px -192px; }
.ui-icon-circle-arrow-e { background-position: -112px -192px; }
.ui-icon-circle-arrow-s { background-position: -128px -192px; }
.ui-icon-circle-arrow-w { background-position: -144px -192px; }
.ui-icon-circle-arrow-n { background-position: -160px -192px; }
.ui-icon-circle-zoomin { background-position: -176px -192px; }
.ui-icon-circle-zoomout { background-position: -192px -192px; }
.ui-icon-circle-check { background-position: -208px -192px; }
.ui-icon-circlesmall-plus { background-position: 0 -208px; }
.ui-icon-circlesmall-minus { background-position: -16px -208px; }
.ui-icon-circlesmall-close { background-position: -32px -208px; }
.ui-icon-squaresmall-plus { background-position: -48px -208px; }
.ui-icon-squaresmall-minus { background-position: -64px -208px; }
.ui-icon-squaresmall-close { background-position: -80px -208px; }
.ui-icon-grip-dotted-vertical { background-position: 0 -224px; }
.ui-icon-grip-dotted-horizontal { background-position: -16px -224px; }
.ui-icon-grip-solid-vertical { background-position: -32px -224px; }
.ui-icon-grip-solid-horizontal { background-position: -48px -224px; }
.ui-icon-gripsmall-diagonal-se { background-position: -64px -224px; }
.ui-icon-grip-diagonal-se { background-position: -80px -224px; }





.ui-corner-all,
.ui-corner-top,
.ui-corner-left,
.ui-corner-tl {
	border-top-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-top,
.ui-corner-right,
.ui-corner-tr {
	border-top-right-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-left,
.ui-corner-bl {
	border-bottom-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-right,
.ui-corner-br {
	border-bottom-right-radius: 3px;
}


.ui-widget-overlay {
	background: #aaaaaa;
	opacity: .3;
	filter: Alpha(Opacity=30); 
}
.ui-widget-shadow {
	-webkit-box-shadow: 0px 0px 5px #666666;
	box-shadow: 0px 0px 5px #666666;
} 
</style>
<script src='js/jquery-1.12.4.js'></script>
<script src='js/jquery-ui.js'></script>

<script>
$( function() {
	$( '#tabs-collect-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-collect-consensus-heatmap'>
<ul>
<li><a href='#tab-collect-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-collect-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-collect-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-collect-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-collect-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-collect-consensus-heatmap-1'>
<pre><code class="language-r">collect_plots(res_list, k = 2, fun = consensus_heatmap, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-consensus-heatmap-2'>
<pre><code class="language-r">collect_plots(res_list, k = 3, fun = consensus_heatmap, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-consensus-heatmap-3'>
<pre><code class="language-r">collect_plots(res_list, k = 4, fun = consensus_heatmap, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-consensus-heatmap-4'>
<pre><code class="language-r">collect_plots(res_list, k = 5, fun = consensus_heatmap, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-consensus-heatmap-5'>
<pre><code class="language-r">collect_plots(res_list, k = 6, fun = consensus_heatmap, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
</div>



### Membership heatmap

Membership heatmaps for all methods. ([What is a membership heatmap?](https://jokergoo.github.io/cola_vignettes/cola.html#toc_12))


<script>
$( function() {
	$( '#tabs-collect-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-collect-membership-heatmap'>
<ul>
<li><a href='#tab-collect-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-collect-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-collect-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-collect-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-collect-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-collect-membership-heatmap-1'>
<pre><code class="language-r">collect_plots(res_list, k = 2, fun = membership_heatmap, cores = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-membership-heatmap-1-1.png" alt="plot of chunk tab-collect-membership-heatmap-1" /></p>

</div>
<div id='tab-collect-membership-heatmap-2'>
<pre><code class="language-r">collect_plots(res_list, k = 3, fun = membership_heatmap, cores = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-membership-heatmap-2-1.png" alt="plot of chunk tab-collect-membership-heatmap-2" /></p>

</div>
<div id='tab-collect-membership-heatmap-3'>
<pre><code class="language-r">collect_plots(res_list, k = 4, fun = membership_heatmap, cores = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-membership-heatmap-3-1.png" alt="plot of chunk tab-collect-membership-heatmap-3" /></p>

</div>
<div id='tab-collect-membership-heatmap-4'>
<pre><code class="language-r">collect_plots(res_list, k = 5, fun = membership_heatmap, cores = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-membership-heatmap-4-1.png" alt="plot of chunk tab-collect-membership-heatmap-4" /></p>

</div>
<div id='tab-collect-membership-heatmap-5'>
<pre><code class="language-r">collect_plots(res_list, k = 6, fun = membership_heatmap, cores = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-membership-heatmap-5-1.png" alt="plot of chunk tab-collect-membership-heatmap-5" /></p>

</div>
</div>



### Signature heatmap

Signature heatmaps for all methods. ([What is a signature heatmap?](https://jokergoo.github.io/cola_vignettes/cola.html#toc_21))


Note in following heatmaps, rows are scaled.



<script>
$( function() {
	$( '#tabs-collect-get-signatures' ).tabs();
} );
</script>
<div id='tabs-collect-get-signatures'>
<ul>
<li><a href='#tab-collect-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-collect-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-collect-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-collect-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-collect-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-collect-get-signatures-1'>
<pre><code class="language-r">collect_plots(res_list, k = 2, fun = get_signatures, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-get-signatures-2'>
<pre><code class="language-r">collect_plots(res_list, k = 3, fun = get_signatures, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-get-signatures-3'>
<pre><code class="language-r">collect_plots(res_list, k = 4, fun = get_signatures, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-get-signatures-4'>
<pre><code class="language-r">collect_plots(res_list, k = 5, fun = get_signatures, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
<div id='tab-collect-get-signatures-5'>
<pre><code class="language-r">collect_plots(res_list, k = 6, fun = get_signatures, cores = 4)
</code></pre>
<pre><code>#&gt; Error in `dev.off()`:
#&gt; ! cannot shut down device 1 (the null device)
</code></pre>

</div>
</div>



### Statistics table

The statistics used for measuring the stability of consensus partitioning.
([How are they
defined?](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13))


<script>
$( function() {
	$( '#tabs-get-stats-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-get-stats-from-consensus-partition-list'>
<ul>
<li><a href='#tab-get-stats-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-get-stats-from-consensus-partition-list-5'>k = 6</a></li>
</ul>
<div id='tab-get-stats-from-consensus-partition-list-1'>
<pre><code class="language-r">get_stats(res_list, k = 2)
</code></pre>
<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; SD:skmeans  2 1.000           1.000       1.000          0.441 0.559   0.559
#&gt; MAD:skmeans 2 1.000           1.000       1.000          0.441 0.559   0.559
#&gt; ATC:skmeans 2 1.000           1.000       1.000          0.504 0.496   0.496
#&gt; SD:kmeans   2 0.421           0.433       0.718          0.383 0.682   0.682
#&gt; MAD:kmeans  2 0.503           0.757       0.799          0.354 0.634   0.634
#&gt; ATC:kmeans  2 1.000           1.000       1.000          0.501 0.499   0.499
#&gt; SD:hclust   2 1.000           1.000       1.000          0.335 0.665   0.665
#&gt; MAD:hclust  2 1.000           1.000       1.000          0.335 0.665   0.665
#&gt; ATC:hclust  2 1.000           1.000       1.000          0.498 0.503   0.503
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-2'>
<pre><code class="language-r">get_stats(res_list, k = 3)
</code></pre>
<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; SD:skmeans  3 1.000           1.000       1.000         0.4285 0.790   0.632
#&gt; MAD:skmeans 3 1.000           1.000       1.000         0.4285 0.790   0.632
#&gt; ATC:skmeans 3 1.000           0.986       1.000         0.0248 0.987   0.975
#&gt; SD:kmeans   3 1.000           0.991       0.987         0.4341 0.666   0.543
#&gt; MAD:kmeans  3 1.000           0.991       0.987         0.5492 0.811   0.701
#&gt; ATC:kmeans  3 0.877           0.865       0.930         0.1924 0.923   0.847
#&gt; SD:hclust   3 1.000           1.000       1.000         0.8816 0.704   0.556
#&gt; MAD:hclust  3 1.000           1.000       1.000         0.8816 0.704   0.556
#&gt; ATC:hclust  3 1.000           1.000       1.000         0.0895 0.955   0.911
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-3'>
<pre><code class="language-r">get_stats(res_list, k = 4)
</code></pre>
<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; SD:skmeans  4 1.000           1.000       1.000         0.0421 0.973   0.928
#&gt; MAD:skmeans 4 1.000           1.000       1.000         0.0421 0.973   0.928
#&gt; ATC:skmeans 4 1.000           0.986       1.000         0.1180 0.939   0.874
#&gt; SD:kmeans   4 0.822           0.788       0.896         0.1698 0.911   0.799
#&gt; MAD:kmeans  4 0.822           0.801       0.899         0.1676 0.911   0.799
#&gt; ATC:kmeans  4 0.738           0.617       0.813         0.1575 0.799   0.564
#&gt; SD:hclust   4 1.000           1.000       1.000         0.0347 0.978   0.941
#&gt; MAD:hclust  4 1.000           1.000       1.000         0.0347 0.978   0.941
#&gt; ATC:hclust  4 1.000           1.000       1.000         0.1643 0.911   0.805
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-4'>
<pre><code class="language-r">get_stats(res_list, k = 5)
</code></pre>
<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; SD:skmeans  5 1.000           1.000       1.000         0.1070 0.930   0.795
#&gt; MAD:skmeans 5 1.000           1.000       1.000         0.1070 0.930   0.795
#&gt; ATC:skmeans 5 1.000           0.985       0.999         0.0789 0.926   0.832
#&gt; SD:kmeans   5 0.796           0.938       0.956         0.0736 0.929   0.807
#&gt; MAD:kmeans  5 0.861           0.928       0.950         0.0801 0.929   0.807
#&gt; ATC:kmeans  5 0.680           0.790       0.676         0.0697 0.904   0.689
#&gt; SD:hclust   5 1.000           1.000       1.000         0.0623 0.959   0.883
#&gt; MAD:hclust  5 1.000           1.000       1.000         0.0623 0.959   0.883
#&gt; ATC:hclust  5 1.000           1.000       1.000         0.1616 0.898   0.723
</code></pre>

</div>
<div id='tab-get-stats-from-consensus-partition-list-5'>
<pre><code class="language-r">get_stats(res_list, k = 6)
</code></pre>
<pre><code>#&gt;             k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#&gt; SD:skmeans  6 0.985           0.974       0.986        0.01699 0.988   0.956
#&gt; MAD:skmeans 6 0.983           0.979       0.989        0.01659 0.988   0.956
#&gt; ATC:skmeans 6 1.000           0.986       1.000        0.04609 0.971   0.923
#&gt; SD:kmeans   6 0.842           0.818       0.897        0.07347 1.000   1.000
#&gt; MAD:kmeans  6 0.843           0.834       0.823        0.08576 0.904   0.700
#&gt; ATC:kmeans  6 0.680           0.718       0.803        0.06427 0.923   0.686
#&gt; SD:hclust   6 1.000           0.986       1.000        0.00677 0.995   0.985
#&gt; MAD:hclust  6 1.000           0.986       1.000        0.00677 0.995   0.985
#&gt; ATC:hclust  6 1.000           1.000       1.000        0.01439 0.989   0.960
</code></pre>

</div>
</div>

Following heatmap plots the partition for each combination of methods and the
lightness correspond to the silhouette scores for samples in each method. On
top the consensus subgroup is inferred from all methods by taking the mean
silhouette scores as weight.


<script>
$( function() {
	$( '#tabs-collect-stats-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-collect-stats-from-consensus-partition-list'>
<ul>
<li><a href='#tab-collect-stats-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-collect-stats-from-consensus-partition-list-5'>k = 6</a></li>
</ul>
<div id='tab-collect-stats-from-consensus-partition-list-1'>
<pre><code class="language-r">collect_stats(res_list, k = 2)
</code></pre>
<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-1-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-1" /></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-2'>
<pre><code class="language-r">collect_stats(res_list, k = 3)
</code></pre>
<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-2-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-2" /></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-3'>
<pre><code class="language-r">collect_stats(res_list, k = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-3-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-3" /></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-4'>
<pre><code class="language-r">collect_stats(res_list, k = 5)
</code></pre>
<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-4-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-4" /></p>

</div>
<div id='tab-collect-stats-from-consensus-partition-list-5'>
<pre><code class="language-r">collect_stats(res_list, k = 6)
</code></pre>
<p><img src="figure_cola/tab-collect-stats-from-consensus-partition-list-5-1.png" alt="plot of chunk tab-collect-stats-from-consensus-partition-list-5" /></p>

</div>
</div>

### Partition from all methods



Collect partitions from all methods:


<script>
$( function() {
	$( '#tabs-collect-classes-from-consensus-partition-list' ).tabs();
} );
</script>
<div id='tabs-collect-classes-from-consensus-partition-list'>
<ul>
<li><a href='#tab-collect-classes-from-consensus-partition-list-1'>k = 2</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-2'>k = 3</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-3'>k = 4</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-4'>k = 5</a></li>
<li><a href='#tab-collect-classes-from-consensus-partition-list-5'>k = 6</a></li>
</ul>
<div id='tab-collect-classes-from-consensus-partition-list-1'>
<pre><code class="language-r">collect_classes(res_list, k = 2)
</code></pre>
<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-1-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-1" /></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-2'>
<pre><code class="language-r">collect_classes(res_list, k = 3)
</code></pre>
<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-2-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-2" /></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-3'>
<pre><code class="language-r">collect_classes(res_list, k = 4)
</code></pre>
<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-3-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-3" /></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-4'>
<pre><code class="language-r">collect_classes(res_list, k = 5)
</code></pre>
<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-4-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-4" /></p>

</div>
<div id='tab-collect-classes-from-consensus-partition-list-5'>
<pre><code class="language-r">collect_classes(res_list, k = 6)
</code></pre>
<p><img src="figure_cola/tab-collect-classes-from-consensus-partition-list-5-1.png" alt="plot of chunk tab-collect-classes-from-consensus-partition-list-5" /></p>

</div>
</div>



### Top rows overlap


Overlap of top rows from different top-row methods:


<script>
$( function() {
	$( '#tabs-top-rows-overlap-by-euler' ).tabs();
} );
</script>
<div id='tabs-top-rows-overlap-by-euler'>
<ul>
<li><a href='#tab-top-rows-overlap-by-euler-1'>top_n = 2</a></li>
</ul>
<div id='tab-top-rows-overlap-by-euler-1'>
<pre><code class="language-r">top_rows_overlap(res_list, top_n = 2, method = &quot;euler&quot;)
</code></pre>
<p><img src="figure_cola/tab-top-rows-overlap-by-euler-1-1.png" alt="plot of chunk tab-top-rows-overlap-by-euler-1" /></p>

</div>
</div>

Also visualize the correspondance of rankings between different top-row methods:


<script>
$( function() {
	$( '#tabs-top-rows-overlap-by-correspondance' ).tabs();
} );
</script>
<div id='tabs-top-rows-overlap-by-correspondance'>
<ul>
<li><a href='#tab-top-rows-overlap-by-correspondance-1'>top_n = 2</a></li>
</ul>
<div id='tab-top-rows-overlap-by-correspondance-1'>
<pre><code class="language-r">top_rows_overlap(res_list, top_n = 2, method = &quot;correspondance&quot;)
</code></pre>
<p><img src="figure_cola/tab-top-rows-overlap-by-correspondance-1-1.png" alt="plot of chunk tab-top-rows-overlap-by-correspondance-1" /></p>

</div>
</div>


Heatmaps of the top rows:


![plot of chunk unnamed-chunk-29](figure_cola/unnamed-chunk-29-1.png)
<script>
$( function() {
	$( '#tabs-top-rows-heatmap' ).tabs();
} );
</script>
<div id='tabs-top-rows-heatmap'>
<ul>
<li><a href='#tab-top-rows-heatmap-1'>top_n = 2</a></li>
</ul>
<div id='tab-top-rows-heatmap-1'>
<pre><code class="language-r">top_rows_heatmap(res_list, top_n = 2)
</code></pre>
<pre><code>#&gt; Error : rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bdf8ba911b93c86026c813bcfdd1888_1_117a8c3605950d.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>
<pre><code>#&gt; Error in `UseMethod()`:
#&gt; ! no applicable method for 'depth' applied to an object of class &quot;NULL&quot;
</code></pre>

</div>
</div>



 
## Results for each method


---------------------------------------------------



### SD:hclust**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["SD", "hclust"]
# you can also extract it by
# res = res_list["SD:hclust"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'SD' method.
#>   Subgroups are detected by 'hclust' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 5.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk SD-hclust-collect-plots](figure_cola/SD-hclust-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk SD-hclust-select-partition-number](figure_cola/SD-hclust-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2     1           1.000           1        0.33517 0.665   0.665
#> 3 3     1           1.000           1        0.88157 0.704   0.556
#> 4 4     1           1.000           1        0.03471 0.978   0.941
#> 5 5     1           1.000           1        0.06229 0.959   0.883
#> 6 6     1           0.986           1        0.00677 0.995   0.985
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 5
#> attr(,"optional")
#> [1] 2 3
```

There is also optional best $k$ = 2 3 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-SD-hclust-get-classes' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-get-classes'>
<ul>
<li><a href='#tab-SD-hclust-get-classes-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-get-classes-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-get-classes-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-get-classes-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-SD-hclust-get-classes-1'>
<p><a id='tab-SD-hclust-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-SD-hclust-get-classes-1-a').parent().next().next().hide();
$('#tab-SD-hclust-get-classes-1-a').click(function(){
  $('#tab-SD-hclust-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-hclust-get-classes-2'>
<p><a id='tab-SD-hclust-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-SD-hclust-get-classes-2-a').parent().next().next().hide();
$('#tab-SD-hclust-get-classes-2-a').click(function(){
  $('#tab-SD-hclust-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-hclust-get-classes-3'>
<p><a id='tab-SD-hclust-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-SD-hclust-get-classes-3-a').parent().next().next().hide();
$('#tab-SD-hclust-get-classes-3-a').click(function(){
  $('#tab-SD-hclust-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-hclust-get-classes-4'>
<p><a id='tab-SD-hclust-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-SD-hclust-get-classes-4-a').parent().next().next().hide();
$('#tab-SD-hclust-get-classes-4-a').click(function(){
  $('#tab-SD-hclust-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-hclust-get-classes-5'>
<p><a id='tab-SD-hclust-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-SD-hclust-get-classes-5-a').parent().next().next().hide();
$('#tab-SD-hclust-get-classes-5-a').click(function(){
  $('#tab-SD-hclust-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-SD-hclust-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-consensus-heatmap'>
<ul>
<li><a href='#tab-SD-hclust-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-hclust-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-consensus-heatmap-1-1.png" alt="plot of chunk tab-SD-hclust-consensus-heatmap-1" /></p>

</div>
<div id='tab-SD-hclust-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-consensus-heatmap-2-1.png" alt="plot of chunk tab-SD-hclust-consensus-heatmap-2" /></p>

</div>
<div id='tab-SD-hclust-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-consensus-heatmap-3-1.png" alt="plot of chunk tab-SD-hclust-consensus-heatmap-3" /></p>

</div>
<div id='tab-SD-hclust-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-consensus-heatmap-4-1.png" alt="plot of chunk tab-SD-hclust-consensus-heatmap-4" /></p>

</div>
<div id='tab-SD-hclust-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-consensus-heatmap-5-1.png" alt="plot of chunk tab-SD-hclust-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-SD-hclust-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-membership-heatmap'>
<ul>
<li><a href='#tab-SD-hclust-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-hclust-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-membership-heatmap-1-1.png" alt="plot of chunk tab-SD-hclust-membership-heatmap-1" /></p>

</div>
<div id='tab-SD-hclust-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-membership-heatmap-2-1.png" alt="plot of chunk tab-SD-hclust-membership-heatmap-2" /></p>

</div>
<div id='tab-SD-hclust-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-membership-heatmap-3-1.png" alt="plot of chunk tab-SD-hclust-membership-heatmap-3" /></p>

</div>
<div id='tab-SD-hclust-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-membership-heatmap-4-1.png" alt="plot of chunk tab-SD-hclust-membership-heatmap-4" /></p>

</div>
<div id='tab-SD-hclust-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-membership-heatmap-5-1.png" alt="plot of chunk tab-SD-hclust-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-SD-hclust-get-signatures' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-get-signatures'>
<ul>
<li><a href='#tab-SD-hclust-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-SD-hclust-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c40433b18.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c1367da90.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c57ac0ab.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c7ce12258.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c92d48dc.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-SD-hclust-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-get-signatures-no-scale'>
<ul>
<li><a href='#tab-SD-hclust-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-SD-hclust-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c6f13e2ea.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4d84def7.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3892d73d.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c46fbe78.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-hclust-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c1d34052.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk SD-hclust-signature_compare](figure_cola/SD-hclust-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-SD-hclust-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-SD-hclust-dimension-reduction'>
<ul>
<li><a href='#tab-SD-hclust-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-SD-hclust-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-SD-hclust-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-SD-hclust-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-SD-hclust-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-SD-hclust-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-dimension-reduction-1-1.png" alt="plot of chunk tab-SD-hclust-dimension-reduction-1" /></p>

</div>
<div id='tab-SD-hclust-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-dimension-reduction-2-1.png" alt="plot of chunk tab-SD-hclust-dimension-reduction-2" /></p>

</div>
<div id='tab-SD-hclust-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-dimension-reduction-3-1.png" alt="plot of chunk tab-SD-hclust-dimension-reduction-3" /></p>

</div>
<div id='tab-SD-hclust-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-dimension-reduction-4-1.png" alt="plot of chunk tab-SD-hclust-dimension-reduction-4" /></p>

</div>
<div id='tab-SD-hclust-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-hclust-dimension-reduction-5-1.png" alt="plot of chunk tab-SD-hclust-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk SD-hclust-collect-classes](figure_cola/SD-hclust-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### SD:kmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["SD", "kmeans"]
# you can also extract it by
# res = res_list["SD:kmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'SD' method.
#>   Subgroups are detected by 'kmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk SD-kmeans-collect-plots](figure_cola/SD-kmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk SD-kmeans-select-partition-number](figure_cola/SD-kmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 0.421           0.433       0.718         0.3825 0.682   0.682
#> 3 3 1.000           0.991       0.987         0.4341 0.666   0.543
#> 4 4 0.822           0.788       0.896         0.1698 0.911   0.799
#> 5 5 0.796           0.938       0.956         0.0736 0.929   0.807
#> 6 6 0.842           0.818       0.897         0.0735 1.000   1.000
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 3
```


Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-SD-kmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-get-classes'>
<ul>
<li><a href='#tab-SD-kmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-SD-kmeans-get-classes-1'>
<p><a id='tab-SD-kmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-SD-kmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-SD-kmeans-get-classes-1-a').click(function(){
  $('#tab-SD-kmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-kmeans-get-classes-2'>
<p><a id='tab-SD-kmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-SD-kmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-SD-kmeans-get-classes-2-a').click(function(){
  $('#tab-SD-kmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-kmeans-get-classes-3'>
<p><a id='tab-SD-kmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-SD-kmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-SD-kmeans-get-classes-3-a').click(function(){
  $('#tab-SD-kmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-kmeans-get-classes-4'>
<p><a id='tab-SD-kmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-SD-kmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-SD-kmeans-get-classes-4-a').click(function(){
  $('#tab-SD-kmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-kmeans-get-classes-5'>
<p><a id='tab-SD-kmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-SD-kmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-SD-kmeans-get-classes-5-a').click(function(){
  $('#tab-SD-kmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-SD-kmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-SD-kmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-kmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-SD-kmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-SD-kmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-SD-kmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-SD-kmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-SD-kmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-SD-kmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-SD-kmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-SD-kmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-SD-kmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-SD-kmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-membership-heatmap'>
<ul>
<li><a href='#tab-SD-kmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-kmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-SD-kmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-SD-kmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-SD-kmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-SD-kmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-SD-kmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-SD-kmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-SD-kmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-SD-kmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-SD-kmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-SD-kmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-get-signatures'>
<ul>
<li><a href='#tab-SD-kmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-SD-kmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-get-signatures-1-1.png" alt="plot of chunk tab-SD-kmeans-get-signatures-1" /></p>

</div>
<div id='tab-SD-kmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c4e51ad4.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c457907fb.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c3ddfde99.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8ca722c42.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-SD-kmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-SD-kmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-SD-kmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-get-signatures-no-scale-1-1.png" alt="plot of chunk tab-SD-kmeans-get-signatures-no-scale-1" /></p>

</div>
<div id='tab-SD-kmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8cea9c1e1.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c7f13dfa1.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c2e34fd36.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-kmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c298ef6d6.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk SD-kmeans-signature_compare](figure_cola/SD-kmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-SD-kmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-SD-kmeans-dimension-reduction'>
<ul>
<li><a href='#tab-SD-kmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-SD-kmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-SD-kmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-SD-kmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-SD-kmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-SD-kmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-SD-kmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-SD-kmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-SD-kmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-SD-kmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-SD-kmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-SD-kmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-SD-kmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-SD-kmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-kmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-SD-kmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk SD-kmeans-collect-classes](figure_cola/SD-kmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### SD:skmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["SD", "skmeans"]
# you can also extract it by
# res = res_list["SD:skmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'SD' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 5.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk SD-skmeans-collect-plots](figure_cola/SD-skmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk SD-skmeans-select-partition-number](figure_cola/SD-skmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000         0.4415 0.559   0.559
#> 3 3 1.000           1.000       1.000         0.4285 0.790   0.632
#> 4 4 1.000           1.000       1.000         0.0421 0.973   0.928
#> 5 5 1.000           1.000       1.000         0.1070 0.930   0.795
#> 6 6 0.985           0.974       0.986         0.0170 0.988   0.956
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 5
#> attr(,"optional")
#> [1] 2 3
```

There is also optional best $k$ = 2 3 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-SD-skmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-get-classes'>
<ul>
<li><a href='#tab-SD-skmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-SD-skmeans-get-classes-1'>
<p><a id='tab-SD-skmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-SD-skmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-SD-skmeans-get-classes-1-a').click(function(){
  $('#tab-SD-skmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-skmeans-get-classes-2'>
<p><a id='tab-SD-skmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-SD-skmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-SD-skmeans-get-classes-2-a').click(function(){
  $('#tab-SD-skmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-skmeans-get-classes-3'>
<p><a id='tab-SD-skmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-SD-skmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-SD-skmeans-get-classes-3-a').click(function(){
  $('#tab-SD-skmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-skmeans-get-classes-4'>
<p><a id='tab-SD-skmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-SD-skmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-SD-skmeans-get-classes-4-a').click(function(){
  $('#tab-SD-skmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-SD-skmeans-get-classes-5'>
<p><a id='tab-SD-skmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-SD-skmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-SD-skmeans-get-classes-5-a').click(function(){
  $('#tab-SD-skmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-SD-skmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-SD-skmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-skmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-SD-skmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-SD-skmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-SD-skmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-SD-skmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-SD-skmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-SD-skmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-SD-skmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-SD-skmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-SD-skmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-SD-skmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-membership-heatmap'>
<ul>
<li><a href='#tab-SD-skmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-SD-skmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-SD-skmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-SD-skmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-SD-skmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-SD-skmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-SD-skmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-SD-skmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-SD-skmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-SD-skmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-SD-skmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-SD-skmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-get-signatures'>
<ul>
<li><a href='#tab-SD-skmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-SD-skmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c410eb872.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c540d55dd.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c357be14c.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c471eb1ba.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c6a20b1e1.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-SD-skmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-SD-skmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-SD-skmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4c5dc26c.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c6f024488.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c701075d3.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3b9fcb16.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-SD-skmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8cad54ef.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk SD-skmeans-signature_compare](figure_cola/SD-skmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-SD-skmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-SD-skmeans-dimension-reduction'>
<ul>
<li><a href='#tab-SD-skmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-SD-skmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-SD-skmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-SD-skmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-SD-skmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-SD-skmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-SD-skmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-SD-skmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-SD-skmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-SD-skmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-SD-skmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-SD-skmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-SD-skmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-SD-skmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-SD-skmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-SD-skmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk SD-skmeans-collect-classes](figure_cola/SD-skmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### MAD:hclust**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["MAD", "hclust"]
# you can also extract it by
# res = res_list["MAD:hclust"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'MAD' method.
#>   Subgroups are detected by 'hclust' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 5.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk MAD-hclust-collect-plots](figure_cola/MAD-hclust-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk MAD-hclust-select-partition-number](figure_cola/MAD-hclust-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2     1           1.000           1        0.33517 0.665   0.665
#> 3 3     1           1.000           1        0.88157 0.704   0.556
#> 4 4     1           1.000           1        0.03471 0.978   0.941
#> 5 5     1           1.000           1        0.06229 0.959   0.883
#> 6 6     1           0.986           1        0.00677 0.995   0.985
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 5
#> attr(,"optional")
#> [1] 2 3
```

There is also optional best $k$ = 2 3 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-MAD-hclust-get-classes' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-get-classes'>
<ul>
<li><a href='#tab-MAD-hclust-get-classes-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-get-classes-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-get-classes-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-get-classes-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-MAD-hclust-get-classes-1'>
<p><a id='tab-MAD-hclust-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-MAD-hclust-get-classes-1-a').parent().next().next().hide();
$('#tab-MAD-hclust-get-classes-1-a').click(function(){
  $('#tab-MAD-hclust-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-hclust-get-classes-2'>
<p><a id='tab-MAD-hclust-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-MAD-hclust-get-classes-2-a').parent().next().next().hide();
$('#tab-MAD-hclust-get-classes-2-a').click(function(){
  $('#tab-MAD-hclust-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-hclust-get-classes-3'>
<p><a id='tab-MAD-hclust-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-MAD-hclust-get-classes-3-a').parent().next().next().hide();
$('#tab-MAD-hclust-get-classes-3-a').click(function(){
  $('#tab-MAD-hclust-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-hclust-get-classes-4'>
<p><a id='tab-MAD-hclust-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-MAD-hclust-get-classes-4-a').parent().next().next().hide();
$('#tab-MAD-hclust-get-classes-4-a').click(function(){
  $('#tab-MAD-hclust-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-hclust-get-classes-5'>
<p><a id='tab-MAD-hclust-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-MAD-hclust-get-classes-5-a').parent().next().next().hide();
$('#tab-MAD-hclust-get-classes-5-a').click(function(){
  $('#tab-MAD-hclust-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-MAD-hclust-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-consensus-heatmap'>
<ul>
<li><a href='#tab-MAD-hclust-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-hclust-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-consensus-heatmap-1-1.png" alt="plot of chunk tab-MAD-hclust-consensus-heatmap-1" /></p>

</div>
<div id='tab-MAD-hclust-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-consensus-heatmap-2-1.png" alt="plot of chunk tab-MAD-hclust-consensus-heatmap-2" /></p>

</div>
<div id='tab-MAD-hclust-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-consensus-heatmap-3-1.png" alt="plot of chunk tab-MAD-hclust-consensus-heatmap-3" /></p>

</div>
<div id='tab-MAD-hclust-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-consensus-heatmap-4-1.png" alt="plot of chunk tab-MAD-hclust-consensus-heatmap-4" /></p>

</div>
<div id='tab-MAD-hclust-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-consensus-heatmap-5-1.png" alt="plot of chunk tab-MAD-hclust-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-MAD-hclust-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-membership-heatmap'>
<ul>
<li><a href='#tab-MAD-hclust-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-hclust-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-membership-heatmap-1-1.png" alt="plot of chunk tab-MAD-hclust-membership-heatmap-1" /></p>

</div>
<div id='tab-MAD-hclust-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-membership-heatmap-2-1.png" alt="plot of chunk tab-MAD-hclust-membership-heatmap-2" /></p>

</div>
<div id='tab-MAD-hclust-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-membership-heatmap-3-1.png" alt="plot of chunk tab-MAD-hclust-membership-heatmap-3" /></p>

</div>
<div id='tab-MAD-hclust-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-membership-heatmap-4-1.png" alt="plot of chunk tab-MAD-hclust-membership-heatmap-4" /></p>

</div>
<div id='tab-MAD-hclust-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-membership-heatmap-5-1.png" alt="plot of chunk tab-MAD-hclust-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-MAD-hclust-get-signatures' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-get-signatures'>
<ul>
<li><a href='#tab-MAD-hclust-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-hclust-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c4f967750.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c1507f0e4.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c6c7483af.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c4b493723.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c4503325c.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-MAD-hclust-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-get-signatures-no-scale'>
<ul>
<li><a href='#tab-MAD-hclust-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-hclust-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8ccffd7a3.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3916d639.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c6a2beeb9.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c5a19f40a.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-hclust-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c47315648.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk MAD-hclust-signature_compare](figure_cola/MAD-hclust-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-MAD-hclust-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-MAD-hclust-dimension-reduction'>
<ul>
<li><a href='#tab-MAD-hclust-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-MAD-hclust-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-MAD-hclust-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-MAD-hclust-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-MAD-hclust-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-hclust-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-dimension-reduction-1-1.png" alt="plot of chunk tab-MAD-hclust-dimension-reduction-1" /></p>

</div>
<div id='tab-MAD-hclust-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-dimension-reduction-2-1.png" alt="plot of chunk tab-MAD-hclust-dimension-reduction-2" /></p>

</div>
<div id='tab-MAD-hclust-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-dimension-reduction-3-1.png" alt="plot of chunk tab-MAD-hclust-dimension-reduction-3" /></p>

</div>
<div id='tab-MAD-hclust-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-dimension-reduction-4-1.png" alt="plot of chunk tab-MAD-hclust-dimension-reduction-4" /></p>

</div>
<div id='tab-MAD-hclust-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-hclust-dimension-reduction-5-1.png" alt="plot of chunk tab-MAD-hclust-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk MAD-hclust-collect-classes](figure_cola/MAD-hclust-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### MAD:kmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["MAD", "kmeans"]
# you can also extract it by
# res = res_list["MAD:kmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'MAD' method.
#>   Subgroups are detected by 'kmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk MAD-kmeans-collect-plots](figure_cola/MAD-kmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk MAD-kmeans-select-partition-number](figure_cola/MAD-kmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 0.503           0.757       0.799         0.3541 0.634   0.634
#> 3 3 1.000           0.991       0.987         0.5492 0.811   0.701
#> 4 4 0.822           0.801       0.899         0.1676 0.911   0.799
#> 5 5 0.861           0.928       0.950         0.0801 0.929   0.807
#> 6 6 0.843           0.834       0.823         0.0858 0.904   0.700
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 3
```


Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-MAD-kmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-get-classes'>
<ul>
<li><a href='#tab-MAD-kmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-MAD-kmeans-get-classes-1'>
<p><a id='tab-MAD-kmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-MAD-kmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-MAD-kmeans-get-classes-1-a').click(function(){
  $('#tab-MAD-kmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-kmeans-get-classes-2'>
<p><a id='tab-MAD-kmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-MAD-kmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-MAD-kmeans-get-classes-2-a').click(function(){
  $('#tab-MAD-kmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-kmeans-get-classes-3'>
<p><a id='tab-MAD-kmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-MAD-kmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-MAD-kmeans-get-classes-3-a').click(function(){
  $('#tab-MAD-kmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-kmeans-get-classes-4'>
<p><a id='tab-MAD-kmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-MAD-kmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-MAD-kmeans-get-classes-4-a').click(function(){
  $('#tab-MAD-kmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-kmeans-get-classes-5'>
<p><a id='tab-MAD-kmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-MAD-kmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-MAD-kmeans-get-classes-5-a').click(function(){
  $('#tab-MAD-kmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-MAD-kmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-MAD-kmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-kmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-MAD-kmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-MAD-kmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-MAD-kmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-MAD-kmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-MAD-kmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-MAD-kmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-MAD-kmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-MAD-kmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-MAD-kmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-MAD-kmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-membership-heatmap'>
<ul>
<li><a href='#tab-MAD-kmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-kmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-MAD-kmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-MAD-kmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-MAD-kmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-MAD-kmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-MAD-kmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-MAD-kmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-MAD-kmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-MAD-kmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-MAD-kmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-MAD-kmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-get-signatures'>
<ul>
<li><a href='#tab-MAD-kmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-kmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c17b79c78.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c6377ff9a.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c32b1dd1f.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c50cc67dd.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c2f14c593.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-MAD-kmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-MAD-kmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-kmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3800bc0f.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3003666f.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c7954da9c.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c7cffdd59.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-kmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c56acb10e.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk MAD-kmeans-signature_compare](figure_cola/MAD-kmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-MAD-kmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-MAD-kmeans-dimension-reduction'>
<ul>
<li><a href='#tab-MAD-kmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-MAD-kmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-MAD-kmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-MAD-kmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-MAD-kmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-kmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-MAD-kmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-MAD-kmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-MAD-kmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-MAD-kmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-MAD-kmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-MAD-kmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-MAD-kmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-MAD-kmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-kmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-MAD-kmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk MAD-kmeans-collect-classes](figure_cola/MAD-kmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### MAD:skmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["MAD", "skmeans"]
# you can also extract it by
# res = res_list["MAD:skmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'MAD' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 5.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk MAD-skmeans-collect-plots](figure_cola/MAD-skmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk MAD-skmeans-select-partition-number](figure_cola/MAD-skmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000         0.4415 0.559   0.559
#> 3 3 1.000           1.000       1.000         0.4285 0.790   0.632
#> 4 4 1.000           1.000       1.000         0.0421 0.973   0.928
#> 5 5 1.000           1.000       1.000         0.1070 0.930   0.795
#> 6 6 0.983           0.979       0.989         0.0166 0.988   0.956
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 5
#> attr(,"optional")
#> [1] 2 3
```

There is also optional best $k$ = 2 3 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-MAD-skmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-get-classes'>
<ul>
<li><a href='#tab-MAD-skmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-MAD-skmeans-get-classes-1'>
<p><a id='tab-MAD-skmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-MAD-skmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-MAD-skmeans-get-classes-1-a').click(function(){
  $('#tab-MAD-skmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-skmeans-get-classes-2'>
<p><a id='tab-MAD-skmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-MAD-skmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-MAD-skmeans-get-classes-2-a').click(function(){
  $('#tab-MAD-skmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-skmeans-get-classes-3'>
<p><a id='tab-MAD-skmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-MAD-skmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-MAD-skmeans-get-classes-3-a').click(function(){
  $('#tab-MAD-skmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-skmeans-get-classes-4'>
<p><a id='tab-MAD-skmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-MAD-skmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-MAD-skmeans-get-classes-4-a').click(function(){
  $('#tab-MAD-skmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-MAD-skmeans-get-classes-5'>
<p><a id='tab-MAD-skmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-MAD-skmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-MAD-skmeans-get-classes-5-a').click(function(){
  $('#tab-MAD-skmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-MAD-skmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-MAD-skmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-skmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-MAD-skmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-MAD-skmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-MAD-skmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-MAD-skmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-MAD-skmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-MAD-skmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-MAD-skmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-MAD-skmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-MAD-skmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-MAD-skmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-membership-heatmap'>
<ul>
<li><a href='#tab-MAD-skmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-skmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-MAD-skmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-MAD-skmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-MAD-skmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-MAD-skmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-MAD-skmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-MAD-skmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-MAD-skmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-MAD-skmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-MAD-skmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-MAD-skmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-get-signatures'>
<ul>
<li><a href='#tab-MAD-skmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-skmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c70fd822a.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c3bbdd23e.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c649ef64b.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c13ddc570.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c9b671e8.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-MAD-skmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-MAD-skmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-skmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4ca76217.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c682621e1.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4d9ef04a.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c54b7fcfe.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-MAD-skmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4785aa7f.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk MAD-skmeans-signature_compare](figure_cola/MAD-skmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-MAD-skmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-MAD-skmeans-dimension-reduction'>
<ul>
<li><a href='#tab-MAD-skmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-MAD-skmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-MAD-skmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-MAD-skmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-MAD-skmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-MAD-skmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-MAD-skmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-MAD-skmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-MAD-skmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-MAD-skmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-MAD-skmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-MAD-skmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-MAD-skmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-MAD-skmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-MAD-skmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-MAD-skmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk MAD-skmeans-collect-classes](figure_cola/MAD-skmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### ATC:hclust**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["ATC", "hclust"]
# you can also extract it by
# res = res_list["ATC:hclust"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'ATC' method.
#>   Subgroups are detected by 'hclust' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 5.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk ATC-hclust-collect-plots](figure_cola/ATC-hclust-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk ATC-hclust-select-partition-number](figure_cola/ATC-hclust-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2     1               1           1         0.4978 0.503   0.503
#> 3 3     1               1           1         0.0895 0.955   0.911
#> 4 4     1               1           1         0.1643 0.911   0.805
#> 5 5     1               1           1         0.1616 0.898   0.723
#> 6 6     1               1           1         0.0144 0.989   0.960
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 5
#> attr(,"optional")
#> [1] 2 3 4
```

There is also optional best $k$ = 2 3 4 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-ATC-hclust-get-classes' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-get-classes'>
<ul>
<li><a href='#tab-ATC-hclust-get-classes-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-get-classes-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-get-classes-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-get-classes-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-ATC-hclust-get-classes-1'>
<p><a id='tab-ATC-hclust-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-ATC-hclust-get-classes-1-a').parent().next().next().hide();
$('#tab-ATC-hclust-get-classes-1-a').click(function(){
  $('#tab-ATC-hclust-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-hclust-get-classes-2'>
<p><a id='tab-ATC-hclust-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-ATC-hclust-get-classes-2-a').parent().next().next().hide();
$('#tab-ATC-hclust-get-classes-2-a').click(function(){
  $('#tab-ATC-hclust-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-hclust-get-classes-3'>
<p><a id='tab-ATC-hclust-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-ATC-hclust-get-classes-3-a').parent().next().next().hide();
$('#tab-ATC-hclust-get-classes-3-a').click(function(){
  $('#tab-ATC-hclust-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-hclust-get-classes-4'>
<p><a id='tab-ATC-hclust-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-ATC-hclust-get-classes-4-a').parent().next().next().hide();
$('#tab-ATC-hclust-get-classes-4-a').click(function(){
  $('#tab-ATC-hclust-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-hclust-get-classes-5'>
<p><a id='tab-ATC-hclust-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-ATC-hclust-get-classes-5-a').parent().next().next().hide();
$('#tab-ATC-hclust-get-classes-5-a').click(function(){
  $('#tab-ATC-hclust-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-ATC-hclust-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-consensus-heatmap'>
<ul>
<li><a href='#tab-ATC-hclust-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-hclust-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-consensus-heatmap-1-1.png" alt="plot of chunk tab-ATC-hclust-consensus-heatmap-1" /></p>

</div>
<div id='tab-ATC-hclust-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-consensus-heatmap-2-1.png" alt="plot of chunk tab-ATC-hclust-consensus-heatmap-2" /></p>

</div>
<div id='tab-ATC-hclust-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-consensus-heatmap-3-1.png" alt="plot of chunk tab-ATC-hclust-consensus-heatmap-3" /></p>

</div>
<div id='tab-ATC-hclust-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-consensus-heatmap-4-1.png" alt="plot of chunk tab-ATC-hclust-consensus-heatmap-4" /></p>

</div>
<div id='tab-ATC-hclust-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-consensus-heatmap-5-1.png" alt="plot of chunk tab-ATC-hclust-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-ATC-hclust-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-membership-heatmap'>
<ul>
<li><a href='#tab-ATC-hclust-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-hclust-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-membership-heatmap-1-1.png" alt="plot of chunk tab-ATC-hclust-membership-heatmap-1" /></p>

</div>
<div id='tab-ATC-hclust-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-membership-heatmap-2-1.png" alt="plot of chunk tab-ATC-hclust-membership-heatmap-2" /></p>

</div>
<div id='tab-ATC-hclust-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-membership-heatmap-3-1.png" alt="plot of chunk tab-ATC-hclust-membership-heatmap-3" /></p>

</div>
<div id='tab-ATC-hclust-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-membership-heatmap-4-1.png" alt="plot of chunk tab-ATC-hclust-membership-heatmap-4" /></p>

</div>
<div id='tab-ATC-hclust-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-membership-heatmap-5-1.png" alt="plot of chunk tab-ATC-hclust-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-ATC-hclust-get-signatures' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-get-signatures'>
<ul>
<li><a href='#tab-ATC-hclust-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-hclust-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c32ddb199.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8cc0c9988.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c772841b6.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c20692905.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c6bc3821d.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-ATC-hclust-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-get-signatures-no-scale'>
<ul>
<li><a href='#tab-ATC-hclust-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-hclust-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c50418630.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c36a74d62.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c76285719.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3edd6ebf.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-hclust-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8ca247c6f.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

```
#> Error in `fit_diagram()`:
#> ! names of elements in `combinations` cannot be duplicated
```

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-ATC-hclust-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-ATC-hclust-dimension-reduction'>
<ul>
<li><a href='#tab-ATC-hclust-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-ATC-hclust-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-ATC-hclust-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-ATC-hclust-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-ATC-hclust-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-hclust-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-dimension-reduction-1-1.png" alt="plot of chunk tab-ATC-hclust-dimension-reduction-1" /></p>

</div>
<div id='tab-ATC-hclust-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-dimension-reduction-2-1.png" alt="plot of chunk tab-ATC-hclust-dimension-reduction-2" /></p>

</div>
<div id='tab-ATC-hclust-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-dimension-reduction-3-1.png" alt="plot of chunk tab-ATC-hclust-dimension-reduction-3" /></p>

</div>
<div id='tab-ATC-hclust-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-dimension-reduction-4-1.png" alt="plot of chunk tab-ATC-hclust-dimension-reduction-4" /></p>

</div>
<div id='tab-ATC-hclust-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-hclust-dimension-reduction-5-1.png" alt="plot of chunk tab-ATC-hclust-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk ATC-hclust-collect-classes](figure_cola/ATC-hclust-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### ATC:kmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["ATC", "kmeans"]
# you can also extract it by
# res = res_list["ATC:kmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'ATC' method.
#>   Subgroups are detected by 'kmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 2.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk ATC-kmeans-collect-plots](figure_cola/ATC-kmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk ATC-kmeans-select-partition-number](figure_cola/ATC-kmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000         0.5013 0.499   0.499
#> 3 3 0.877           0.865       0.930         0.1924 0.923   0.847
#> 4 4 0.738           0.617       0.813         0.1575 0.799   0.564
#> 5 5 0.680           0.790       0.676         0.0697 0.904   0.689
#> 6 6 0.680           0.718       0.803         0.0643 0.923   0.686
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 2
```


Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-ATC-kmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-get-classes'>
<ul>
<li><a href='#tab-ATC-kmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-ATC-kmeans-get-classes-1'>
<p><a id='tab-ATC-kmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-ATC-kmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-ATC-kmeans-get-classes-1-a').click(function(){
  $('#tab-ATC-kmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-kmeans-get-classes-2'>
<p><a id='tab-ATC-kmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-ATC-kmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-ATC-kmeans-get-classes-2-a').click(function(){
  $('#tab-ATC-kmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-kmeans-get-classes-3'>
<p><a id='tab-ATC-kmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-ATC-kmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-ATC-kmeans-get-classes-3-a').click(function(){
  $('#tab-ATC-kmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-kmeans-get-classes-4'>
<p><a id='tab-ATC-kmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-ATC-kmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-ATC-kmeans-get-classes-4-a').click(function(){
  $('#tab-ATC-kmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-kmeans-get-classes-5'>
<p><a id='tab-ATC-kmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-ATC-kmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-ATC-kmeans-get-classes-5-a').click(function(){
  $('#tab-ATC-kmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-ATC-kmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-ATC-kmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-kmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-ATC-kmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-ATC-kmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-ATC-kmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-ATC-kmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-ATC-kmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-ATC-kmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-ATC-kmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-ATC-kmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-ATC-kmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-ATC-kmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-membership-heatmap'>
<ul>
<li><a href='#tab-ATC-kmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-kmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-ATC-kmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-ATC-kmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-ATC-kmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-ATC-kmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-ATC-kmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-ATC-kmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-ATC-kmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-ATC-kmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-ATC-kmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-ATC-kmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-get-signatures'>
<ul>
<li><a href='#tab-ATC-kmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-kmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c222ecb8b.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c249d128b.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8cea173bd.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c7c72c98.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c5053778c.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-ATC-kmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-ATC-kmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-kmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c316f53b7.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c50c858d3.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c71603cfe.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c72cb39bd.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-kmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c55e71213.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk ATC-kmeans-signature_compare](figure_cola/ATC-kmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-ATC-kmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-ATC-kmeans-dimension-reduction'>
<ul>
<li><a href='#tab-ATC-kmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-ATC-kmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-ATC-kmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-ATC-kmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-ATC-kmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-kmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-ATC-kmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-ATC-kmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-ATC-kmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-ATC-kmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-ATC-kmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-ATC-kmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-ATC-kmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-ATC-kmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-kmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-ATC-kmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk ATC-kmeans-collect-classes](figure_cola/ATC-kmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------



### ATC:skmeans**






The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

``` r
res = res_list["ATC", "skmeans"]
# you can also extract it by
# res = res_list["ATC:skmeans"]
```

A summary of `res` and all the functions that can be applied to it:

``` r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
#>   On a matrix with 22 rows and 72 columns.
#>   Top rows (2) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 250 partitions by row resampling.
#>   Best k for subgroups seems to be 6.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

``` r
collect_plots(res)
```

![plot of chunk ATC-skmeans-collect-plots](figure_cola/ATC-skmeans-collect-plots-1.png)

```
#> Error in `UseMethod()`:
#> ! no applicable method for 'depth' applied to an object of class "NULL"
```

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

``` r
select_partition_number(res)
```

![plot of chunk ATC-skmeans-select-partition-number](figure_cola/ATC-skmeans-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

``` r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2     1           1.000       1.000         0.5040 0.496   0.496
#> 3 3     1           0.986       1.000         0.0248 0.987   0.975
#> 4 4     1           0.986       1.000         0.1180 0.939   0.874
#> 5 5     1           0.985       0.999         0.0789 0.926   0.832
#> 6 6     1           0.986       1.000         0.0461 0.971   0.923
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

``` r
suggest_best_k(res)
```

```
#> [1] 6
#> attr(,"optional")
#> [1] 2 4 5
```

There is also optional best $k$ = 2 4 5 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-classes' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-classes'>
<ul>
<li><a href='#tab-ATC-skmeans-get-classes-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-classes-5'>k = 6</a></li>
</ul>

<div id='tab-ATC-skmeans-get-classes-1'>
<p><a id='tab-ATC-skmeans-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-1-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-1-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-2'>
<p><a id='tab-ATC-skmeans-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-2-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-2-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-3'>
<p><a id='tab-ATC-skmeans-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-3-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-3-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-4'>
<p><a id='tab-ATC-skmeans-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-4-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-4-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-ATC-skmeans-get-classes-5'>
<p><a id='tab-ATC-skmeans-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="language-r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<script>
$('#tab-ATC-skmeans-get-classes-5-a').parent().next().next().hide();
$('#tab-ATC-skmeans-get-classes-5-a').click(function(){
  $('#tab-ATC-skmeans-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-ATC-skmeans-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-consensus-heatmap'>
<ul>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-consensus-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-skmeans-consensus-heatmap-1'>
<pre><code class="language-r">consensus_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-1-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-1" /></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-2'>
<pre><code class="language-r">consensus_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-2-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-2" /></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-3'>
<pre><code class="language-r">consensus_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-3-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-3" /></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-4'>
<pre><code class="language-r">consensus_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-4-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-4" /></p>

</div>
<div id='tab-ATC-skmeans-consensus-heatmap-5'>
<pre><code class="language-r">consensus_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-consensus-heatmap-5-1.png" alt="plot of chunk tab-ATC-skmeans-consensus-heatmap-5" /></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-ATC-skmeans-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-membership-heatmap'>
<ul>
<li><a href='#tab-ATC-skmeans-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-membership-heatmap-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-skmeans-membership-heatmap-1'>
<pre><code class="language-r">membership_heatmap(res, k = 2)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-1-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-1" /></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-2'>
<pre><code class="language-r">membership_heatmap(res, k = 3)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-2-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-2" /></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-3'>
<pre><code class="language-r">membership_heatmap(res, k = 4)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-3-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-3" /></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-4'>
<pre><code class="language-r">membership_heatmap(res, k = 5)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-4-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-4" /></p>

</div>
<div id='tab-ATC-skmeans-membership-heatmap-5'>
<pre><code class="language-r">membership_heatmap(res, k = 6)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-membership-heatmap-5-1.png" alt="plot of chunk tab-ATC-skmeans-membership-heatmap-5" /></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-signatures' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-signatures'>
<ul>
<li><a href='#tab-ATC-skmeans-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-skmeans-get-signatures-1'>
<pre><code class="language-r">get_signatures(res, k = 2)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c5d0f9d77.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-2'>
<pre><code class="language-r">get_signatures(res, k = 3)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c566aadf6.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-3'>
<pre><code class="language-r">get_signatures(res, k = 4)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c4974c135.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-4'>
<pre><code class="language-r">get_signatures(res, k = 5)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c6a2b5af2.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-5'>
<pre><code class="language-r">get_signatures(res, k = 6)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_3bce8d10d0e103205e0651295f56585e_1_117a8c33118cb2.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-ATC-skmeans-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-get-signatures-no-scale'>
<ul>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-get-signatures-no-scale-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-skmeans-get-signatures-no-scale-1'>
<pre><code class="language-r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c4bc3d911.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-2'>
<pre><code class="language-r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c13410686.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-3'>
<pre><code class="language-r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c7ad6f571.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-4'>
<pre><code class="language-r">get_signatures(res, k = 5, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3db67612.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
<div id='tab-ATC-skmeans-get-signatures-no-scale-5'>
<pre><code class="language-r">get_signatures(res, k = 6, scale_rows = FALSE)
</code></pre>
<pre><code>#&gt; Error:
#&gt; ! rsession-arm64: UnableToOpenBlob `/var/folders/6d/7_7gcqq502dbbjl82c6tz9r40000gn/T//RtmpRGlT7N/.heatmap_body_18d54481c36ea4f833abe59ded88edd4_1_117a8c3c7d9402.png': No such file or directory @ error/blob.c/OpenBlob/2967
</code></pre>

</div>
</div>



Compare the overlap of signatures from different k:

``` r
compare_signatures(res)
```

![plot of chunk ATC-skmeans-signature_compare](figure_cola/ATC-skmeans-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

``` r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

``` r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

``` r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-ATC-skmeans-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-ATC-skmeans-dimension-reduction'>
<ul>
<li><a href='#tab-ATC-skmeans-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-ATC-skmeans-dimension-reduction-5'>k = 6</a></li>
</ul>
<div id='tab-ATC-skmeans-dimension-reduction-1'>
<pre><code class="language-r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-dimension-reduction-1-1.png" alt="plot of chunk tab-ATC-skmeans-dimension-reduction-1" /></p>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-2'>
<pre><code class="language-r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-dimension-reduction-2-1.png" alt="plot of chunk tab-ATC-skmeans-dimension-reduction-2" /></p>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-3'>
<pre><code class="language-r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-dimension-reduction-3-1.png" alt="plot of chunk tab-ATC-skmeans-dimension-reduction-3" /></p>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-4'>
<pre><code class="language-r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-dimension-reduction-4-1.png" alt="plot of chunk tab-ATC-skmeans-dimension-reduction-4" /></p>

</div>
<div id='tab-ATC-skmeans-dimension-reduction-5'>
<pre><code class="language-r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>
<p><img src="figure_cola/tab-ATC-skmeans-dimension-reduction-5-1.png" alt="plot of chunk tab-ATC-skmeans-dimension-reduction-5" /></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

``` r
collect_classes(res)
```

![plot of chunk ATC-skmeans-collect-classes](figure_cola/ATC-skmeans-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

## Session info


``` r
sessionInfo()
```

```
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sonoma 14.6
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> Random number generation:
#>  RNG:     L'Ecuyer-CMRG 
#>  Normal:  Inversion 
#>  Sample:  Rejection 
#>  
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] fastDummies_1.7.5     GSVA_1.52.3           matrixStats_1.5.0     fgsea_1.30.0         
#>  [5] plotly_4.12.0         ggrepel_0.9.6         gridExtra_2.3         markdown_2.0         
#>  [9] knitr_1.51            doRNG_1.8.6.2         rngtools_1.5.2        patchwork_1.3.2      
#> [13] impute_1.78.0         devtools_2.4.6        usethis_3.2.1         Rcpp_1.1.1           
#> [17] cola_2.10.0           umap_0.2.10.0         BayesDeBulk_1.0       DreamAI_0.1.0        
#> [21] itertools_0.1-3       iterators_1.0.14      foreach_1.5.2         glmnet_4.1-10        
#> [25] Matrix_1.7-4          missForest_1.6.1      randomForest_4.7-1.2  survival_3.8-6       
#> [29] cluster_2.1.8.1       visdat_0.6.0          sva_3.52.0            BiocParallel_1.38.0  
#> [33] genefilter_1.86.0     mgcv_1.9-4            nlme_3.1-168          limma_3.60.6         
#> [37] readxl_1.4.5          HarmonizR_1.2.0       lubridate_1.9.4       stringr_1.6.0        
#> [41] purrr_1.2.1           readr_2.1.6           tidyverse_2.0.0       scrime_1.3.7         
#> [45] psych_2.6.1           tibble_3.3.1          circlize_0.4.17       ComplexHeatmap_2.20.0
#> [49] randomcoloR_1.1.0.1   cowplot_1.2.0         ggplot2_4.0.2         forcats_1.0.1        
#> [53] tidyr_1.3.2           dplyr_1.2.0           janitor_2.2.1         pcaMethods_1.96.0    
#> [57] Biobase_2.64.0        BiocGenerics_0.50.0  
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_1.6.6                    httr_1.4.7                  RColorBrewer_1.1-3         
#>   [4] doParallel_1.0.17           tools_4.4.1                 R6_2.6.1                   
#>   [7] HDF5Array_1.32.1            lazyeval_0.2.2              rhdf5filters_1.16.0        
#>  [10] GetoptLong_1.1.0            litedown_0.9                withr_3.0.2                
#>  [13] quantreg_6.1                cli_3.6.5                   microbenchmark_1.5.0       
#>  [16] labeling_0.4.3              slam_0.1-55                 S7_0.2.1                   
#>  [19] askpass_1.2.1               commonmark_2.0.0            MCMCpack_1.7-1             
#>  [22] sessioninfo_1.2.3           rstudioapi_0.18.0           RSQLite_2.4.5              
#>  [25] generics_0.1.4              shape_1.4.6.1               crosstalk_1.2.2            
#>  [28] S4Vectors_0.42.1            abind_1.4-8                 lifecycle_1.0.5            
#>  [31] yaml_2.3.12                 edgeR_4.2.2                 snakecase_0.11.1           
#>  [34] SummarizedExperiment_1.34.0 rhdf5_2.48.0                SparseArray_1.4.8          
#>  [37] Rtsne_0.17                  blob_1.3.0                  crayon_1.5.3               
#>  [40] lattice_0.22-7              beachmat_2.20.0             annotate_1.82.0            
#>  [43] KEGGREST_1.44.1             magick_2.9.0                pillar_1.11.1              
#>  [46] GenomicRanges_1.56.2        rjson_0.2.23                codetools_0.2-20           
#>  [49] fastmatch_1.1-8             glue_1.8.0                  V8_8.0.1                   
#>  [52] data.table_1.18.2.1         remotes_2.5.0               vctrs_0.7.1                
#>  [55] png_0.1-8                   Rdpack_2.6.5                cellranger_1.1.0           
#>  [58] gtable_0.3.6                cachem_1.1.0                xfun_0.56                  
#>  [61] rbibutils_2.4.1             S4Arrays_1.4.1              skmeans_0.2-18             
#>  [64] coda_0.19-4.1               SingleCellExperiment_1.26.0 statmod_1.5.1              
#>  [67] ellipsis_0.3.2              brew_1.0-10                 mcmc_0.9-8                 
#>  [70] bit64_4.6.0-1               GenomeInfoDb_1.40.1         irlba_2.3.7                
#>  [73] otel_0.2.0                  colorspace_2.1-2            DBI_1.2.3                  
#>  [76] mnormt_2.1.2                tidyselect_1.2.1            bit_4.6.0                  
#>  [79] compiler_4.4.1              curl_7.0.0                  graph_1.82.0               
#>  [82] SparseM_1.84-2              xml2_1.5.2                  DelayedArray_0.30.1        
#>  [85] scales_1.4.0                SpatialExperiment_1.14.0    digest_0.6.39              
#>  [88] rmarkdown_2.30              XVector_0.44.0              htmltools_0.5.9            
#>  [91] pkgconfig_2.0.3             sparseMatrixStats_1.16.0    MatrixGenerics_1.16.0      
#>  [94] fastmap_1.2.0               rlang_1.1.7                 GlobalOptions_0.1.3        
#>  [97] htmlwidgets_1.6.4           UCSC.utils_1.0.0            farver_2.1.2               
#> [100] jsonlite_2.0.0              mclust_6.1.2                BiocSingular_1.20.0        
#> [103] magrittr_2.0.4              GenomeInfoDbData_1.2.12     Rhdf5lib_1.26.0            
#> [106] reticulate_1.44.1           stringi_1.8.7               zlibbioc_1.50.0            
#> [109] MASS_7.3-65                 plyr_1.8.9                  pkgbuild_1.4.8             
#> [112] parallel_4.4.1              Biostrings_2.72.1           splines_4.4.1              
#> [115] hms_1.1.4                   locfit_1.5-9.12             polylabelr_1.0.0           
#> [118] ranger_0.18.0               stats4_4.4.1                ScaledMatrix_1.12.0        
#> [121] pkgload_1.5.0               XML_3.99-0.20               evaluate_1.0.5             
#> [124] BiocManager_1.30.27         tzdb_0.5.0                  MatrixModels_0.5-4         
#> [127] openssl_2.3.4               polyclip_1.10-7             clue_0.3-66                
#> [130] rsvd_1.0.5                  eulerr_7.0.4                xtable_1.8-4               
#> [133] RSpectra_0.16-2             viridisLite_0.4.2           truncnorm_1.0-9            
#> [136] memoise_2.0.1               AnnotationDbi_1.66.0        IRanges_2.38.1             
#> [139] timechange_0.4.0            GSEABase_1.66.0
```




<script type="text/javascript">
$(function() {
    $("#TOC > ul > li > a").remove();
}); 
</script>
