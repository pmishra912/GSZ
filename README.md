# GSZ
Robust gene set analysis methods using Gene Set Z-score



## This is an example analysis of airway RNAseq gene expression data using mGSZm

Load RNAseq gene expression data


``` r
# source("https://bioconductor.org/biocLite.R")
# biocLite("airway")
suppressWarnings(suppressMessages(library(airway)))
data(airway)
```

Count data


``` r
countdata <- assay(airway)
```

Metadata

``` r
coldata <- colData(airway)
```

Sample labels

``` r
sample.labels <- coldata$dex
```

Filter out non-expressed genes (criteria - must have more than 5 reads in at least two samples for each gene)


``` r
filter <- apply(countdata, 1, function(x) length(x[x>5])>=2)
countdata <- countdata[filter,]
```

Install (attached source code) and load mGSZ

``` r
library(mGSZm)
```

Load example gene set data (attached) and take a subset of it to make runs faster


``` r
load("airway.gene.sets.RData")
gene.sets <- gene.sets.airway.go[1:50]
```

Create mGSZm object

``` r
# exprdat:  raw count data
# genesets: gene set data (list or binary matrix)
# sampleLabels: sample labels (example, c("treatment","treatment","treatment","control","control","control"))
gsz_obj <- GSZ(exprdat=countdata,genesets=gene.sets,sampleLabels=sample.labels)
```

Run mGSZm (we have only two treatments/groups)

``` r
# rnaseq: logical. TRUE for RNAseq data
# perm.number: number of permutations for P-value estimation
res1 <- mGSZm(gsz_obj,rnaseq=TRUE,perm.number=10)
```

```
## [1] "Gene set analysis of two-group gene expression data"
## [1] "RNAseq data inputed. This software expects raw counts!"
```


Table with results

``` r
# n : number of most significant genes sets
# order.by: criteria for ordering the table (P-value, FDR)
top <- topResults(res1,n=10,order.by = "FDR")
top
```

```
##     Gene sets No. of member genes upper tail lower tail GSZ scores direction      P-value         FDR
## 11 GO:0033194                   6   2.618348  2.4107132   2.618348        up 3.137923e-05 0.001286548
## 38 GO:0055093                  15   2.984119  1.1922521   2.984119        up 7.303299e-04 0.014971763
## 40 GO:0071385                   7   1.905875  1.0603899   1.905875        up 1.976136e-03 0.027007189
## 1  GO:0031966                  71   2.647247  1.5935251   2.647247        up 8.253267e-03 0.084595987
## 7  GO:0042493                 234   2.287894  2.6675029   2.667503      down 2.514356e-02 0.206177174
## 14 GO:0019901                 331   2.110277  2.6889365   2.688936      down 6.237142e-02 0.370594036
## 36 GO:0015078                  16   2.145527  1.4790183   2.145527        up 6.327215e-02 0.370594036
## 39 GO:0008535                   9   1.430228  0.7273896   1.430228        up 1.462364e-01 0.749461745
## 21 GO:0005751                   6   1.687206  0.8951459   1.687206        up 1.827455e-01 0.753139552
## 41 GO:0009725                  45   1.701114  2.0553979   2.055398      down 1.836926e-01 0.753139552
```

In order to make an example run with multi-group (> 2 groups) data, lets create multi-group (>2 groups) gene expression data by modifying the sample labels as follows,


``` r
sample.labels.multi <- c("treatment1","treatment1","treatment2","treatment2","treatment2","control","control","control")
```

Update the sample.labels of gsz_obj


``` r
sampleLabels(gsz_obj) <- sample.labels.multi
```

Run mGSZm for following comparisons: group1-group2 and group1-group2


``` r
comparisons <- c("treatment1-control","treatment2-control")
res1 <- mGSZm(gsz_obj,rnaseq=TRUE,comparison=comparisons,perm.number=10)
```

```
## [1] "Gene set analysis with > 2 groups. Advanced permutation technique will be used if the number of replicates is < 6."
## [1] "RNAseq data inputed. This software expects raw counts!"
```

Table with top 10 gene sets


``` r
# direction: up regulated, down regulated or mixed. 
# In Control vs Treatment setting, up represents that a geneset is upregulated in Treatment. The program orders the groups and the first group is used as reference.

# comparison: comparison of interest
top_results <- topResults(res1,comparison=comparisons[1], order.by="FDR",direction="up",n=10)
top_results
```

```
##     Gene sets No. of member genes upper tail lower tail GSZ scores direction    P-value       FDR
## 15 GO:0035255                  18   1.481503  1.1859946   1.481503        up 0.12889862 0.4962522
## 20 GO:0009060                  28   2.157680  2.1359924   2.157680        up 0.11629276 0.4962522
## 22 GO:0007568                 133   2.185555  2.1367136   2.185555        up 0.09569033 0.4962522
## 24 GO:0021549                  30   1.964561  1.7841433   1.964561        up 0.06229226 0.4962522
## 17 GO:0020037                  69   2.271512  1.3018626   2.271512        up 0.19238838 0.5634231
## 35 GO:0005753                  21   3.168748  2.9606590   3.168748        up 0.21415379 0.5853537
## 30 GO:0009409                  27   1.859512  0.9640503   1.859512        up 0.28665777 0.5882292
## 33 GO:0006754                  27   2.580623  2.1095445   2.580623        up 0.28694106 0.5882292
## 16 GO:0005506                  78   2.033882  1.7918842   2.033882        up 0.38553639 0.6124060
## 25 GO:0051602                  18   1.417770  1.3612023   1.417770        up 0.38835502 0.6124060
```


``` r
# generate result table with all comparisons
top_results1 <- topResults(res1,comparison=comparisons, order.by="FDR",direction="up",n=10)
top_results1
```

```
##              Contrast  Gene sets No. of member genes upper tail lower tail GSZ scores direction
## 1  treatment1-control GO:0035255                  18   1.481503  1.1859946   1.481503        up
## 2  treatment1-control GO:0009060                  28   2.157680  2.1359924   2.157680        up
## 3  treatment1-control GO:0007568                 133   2.185555  2.1367136   2.185555        up
## 4  treatment1-control GO:0021549                  30   1.964561  1.7841433   1.964561        up
## 5  treatment1-control GO:0020037                  69   2.271512  1.3018626   2.271512        up
## 6  treatment1-control GO:0005753                  21   3.168748  2.9606590   3.168748        up
## 7  treatment1-control GO:0009409                  27   1.859512  0.9640503   1.859512        up
## 8  treatment1-control GO:0006754                  27   2.580623  2.1095445   2.580623        up
## 9  treatment1-control GO:0005506                  78   2.033882  1.7918842   2.033882        up
## 10 treatment1-control GO:0051602                  18   1.417770  1.3612023   1.417770        up
##       P-value       FDR           Contrast  Gene sets No. of member genes upper tail lower tail
## 1  0.12889862 0.4962522 treatment2-control GO:0020037                  69   2.762692  2.1281358
## 2  0.11629276 0.4962522 treatment2-control GO:0016887                 129   2.464663  2.1588607
## 3  0.09569033 0.4962522 treatment2-control GO:0055093                  15   2.098075  2.0814667
## 4  0.06229226 0.4962522 treatment2-control GO:0005506                  78   2.551349  1.6590735
## 5  0.19238838 0.5634231 treatment2-control GO:0021549                  30   1.675350  1.1937824
## 6  0.21415379 0.5853537 treatment2-control GO:0031966                  71   1.853553  1.4361492
## 7  0.28665777 0.5882292 treatment2-control GO:0043025                 246   1.651968  1.4823446
## 8  0.28694106 0.5882292 treatment2-control GO:0072593                  33   1.441924  1.0076704
## 9  0.38553639 0.6124060 treatment2-control GO:0035255                  18   1.073631  0.7170352
## 10 0.38835502 0.6124060 treatment2-control GO:0015078                  16   1.933010  1.6485572
##    GSZ scores direction    P-value       FDR
## 1    2.762692        up 0.02878150 0.5378274
## 2    2.464663        up 0.06995491 0.5378274
## 3    2.098075        up 0.06531211 0.5378274
## 4    2.551349        up 0.25053153 0.7102841
## 5    1.675350        up 0.35080690 0.7990602
## 6    1.853553        up 0.48898845 0.8616139
## 7    1.651968        up 0.55630267 0.8616139
## 8    1.441924        up 0.42911666 0.8616139
## 9    1.073631        up 0.58841926 0.8616139
## 10   1.933010        up 0.55820749 0.8616139
```

Generate result table with all the analyzed gene sets without any ordering


``` r
all_results <- topResults(res1,comparison=comparisons, order.by="none",n=Inf)
head(all_results)
```

```
##             Contrast  Gene sets No. of member genes upper tail lower tail GSZ scores direction    P-value
## 1 treatment1-control GO:0031966                  71   1.017645   1.130806   1.130806      down 0.85394013
## 2 treatment1-control GO:0005743                 401   6.148015   6.631024   6.631024      down 0.17211280
## 3 treatment1-control GO:0006120                  48   3.242584   3.970406   3.970406      down 0.18464652
## 4 treatment1-control GO:0008137                  44   2.943147   3.624178   3.624178      down 0.27089501
## 5 treatment1-control GO:0005747                  48   3.776392   4.524348   4.524348      down 0.07499684
## 6 treatment1-control GO:0032981                  61   3.794903   4.281379   4.281379      down 0.13314084
##         FDR           Contrast  Gene sets No. of member genes upper tail lower tail GSZ scores direction
## 1 0.8752886 treatment2-control GO:0031966                  71   1.853553   1.436149   1.853553        up
## 2 0.5634231 treatment2-control GO:0005743                 401   3.947784   4.092942   4.092942      down
## 3 0.5634231 treatment2-control GO:0006120                  48   4.166147   4.381710   4.381710      down
## 4 0.5882292 treatment2-control GO:0008137                  44   4.001425   4.120358   4.120358      down
## 5 0.4962522 treatment2-control GO:0005747                  48   3.727895   3.959048   3.959048      down
## 6 0.4962522 treatment2-control GO:0032981                  61   4.108161   4.405612   4.405612      down
##      P-value       FDR
## 1 0.48898845 0.8616139
## 2 0.55441683 0.8616139
## 3 0.09718069 0.5378274
## 4 0.14421166 0.5912678
## 5 0.11805967 0.5378274
## 6 0.04827556 0.5378274
```

