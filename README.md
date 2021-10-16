**Direct RNA nanopore sequencing of Pseudomonas aeruginosa clone C transcriptomes**<br/>
Marie-Madlen Pust<sup>1,2</sup>, Colin Davenport<sup>3</sup>, Lutz Wiehlmann<sup>3</sup>, Burkhard Tümmler<sup>1,2*</sup> 
<br/><br/>
<sup>1</sup>Department of Paediatric Pneumology, Allergology, and Neonatology, Hannover Medical School, Germany <br/>
<sup>2</sup>Biomedical Research in Endstage and Obstructive Lung Disease Hannover (BREATH), German Center for Lung Research, Hannover Medical School, Germany <br/>
<sup>3</sup>Research Core Unit Genomics, Hannover Medical School, Germany <br/>
<br/>

**Running title** <br/>
Direct RNA nanopore sequencing of SG17M and NN2
<br/>
<br/>
R Session Info


```{r}
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggcorrplot_0.1.3            ggdendro_0.1.22             forcats_0.5.1               purrr_0.3.4                
 [5] tibble_3.1.5                tidyverse_1.3.1             lemon_0.4.5                 RColorBrewer_1.1-2         
 [9] factoextra_1.0.7            scales_1.1.1                sjmisc_2.8.7                ggplotify_0.1.0            
[13] ggrepel_0.9.1               ashr_2.2-47                 DESeq2_1.32.0               plyr_1.8.6                 
[17] seqinr_4.2-8                TmCalculator_1.0.1          pheatmap_1.0.12             corrplot_0.90              
[21] Hmisc_4.5-0                 Formula_1.2-4               survival_3.2-13             rcompanion_2.4.1           
[25] vegan_2.5-7                 lattice_0.20-44             permute_0.9-5               GenomicAlignments_1.28.0   
[29] SummarizedExperiment_1.22.0 Biobase_2.52.0              MatrixGenerics_1.4.3        matrixStats_0.61.0         
[33] Rsamtools_2.8.0             Biostrings_2.60.2           XVector_0.32.0              GenomicRanges_1.44.0       
[37] GenomeInfoDb_1.28.4         IRanges_2.26.0              S4Vectors_0.30.1            BiocGenerics_0.38.0        
[41] dplyr_1.0.7                 circlize_0.4.13             ggpubr_0.4.0                tidyr_1.1.4                
[45] stringr_1.4.0               edgeR_3.34.1                limma_3.48.3                Rsubread_2.6.4             
[49] readxl_1.3.1                readr_2.0.2                 ggplot2_3.3.5               rlist_0.4.6.2              

loaded via a namespace (and not attached):
  [1] utf8_1.2.2             tidyselect_1.1.1       RSQLite_2.2.8          AnnotationDbi_1.54.1   htmlwidgets_1.5.4     
  [6] grid_4.1.1             BiocParallel_1.26.2    munsell_0.5.0          codetools_0.2-18       withr_2.4.2           
 [11] colorspace_2.0-2       knitr_1.36             rstudioapi_0.13        DescTools_0.99.43      ggsignif_0.6.3        
 [16] labeling_0.4.2         GenomeInfoDbData_1.2.6 mixsqp_0.3-43          farver_2.1.0           bit64_4.0.5           
 [21] vctrs_0.3.8            generics_0.1.0         TH.data_1.1-0          xfun_0.26              R6_2.5.1              
 [26] invgamma_1.1           locfit_1.5-9.4         bitops_1.0-7           cachem_1.0.6           gridGraphics_0.5-1    
 [31] DelayedArray_0.18.0    assertthat_0.2.1       vroom_1.5.5            multcomp_1.4-17        nnet_7.3-16           
 [36] rootSolve_1.8.2.3      gtable_0.3.0           multcompView_0.1-8     lmom_2.8               sandwich_3.0-1        
 [41] rlang_0.4.11           genefilter_1.74.0      GlobalOptions_0.1.2    splines_4.1.1          rstatix_0.7.0         
 [46] broom_0.7.9            checkmate_2.0.0        modelr_0.1.8           yaml_2.2.1             abind_1.4-5           
 [51] backports_1.2.1        tools_4.1.1            ellipsis_0.3.2         proxy_0.4-26           Rcpp_1.0.7            
 [56] base64enc_0.1-3        zlibbioc_1.38.0        RCurl_1.98-1.5         rpart_4.1-15           cowplot_1.1.1         
 [61] zoo_1.8-9              haven_2.4.3            cluster_2.1.2          fs_1.5.0               magrittr_2.0.1        
 [66] data.table_1.14.2      openxlsx_4.2.4         reprex_2.0.1           lmtest_0.9-38          truncnorm_1.0-8       
 [71] mvtnorm_1.1-2          SQUAREM_2021.1         hms_1.1.1              evaluate_0.14          xtable_1.8-4          
 [76] XML_3.99-0.8           rio_0.5.27             jpeg_0.1-9             gridExtra_2.3          shape_1.4.6           
 [81] compiler_4.1.1         crayon_1.4.1           htmltools_0.5.2        mgcv_1.8-36            tzdb_0.1.2            
 [86] geneplotter_1.70.0     libcoin_1.0-9          expm_0.999-6           Exact_3.0              lubridate_1.8.0       
 [91] DBI_1.1.1              sjlabelled_1.1.8       dbplyr_2.1.1           MASS_7.3-54            boot_1.3-28           
 [96] Matrix_1.3-4           ade4_1.7-18            car_3.0-11             cli_3.0.1              insight_0.14.4        
[101] pkgconfig_2.0.3        coin_1.4-1             foreign_0.8-81         xml2_1.3.2             annotate_1.70.0       
[106] rvest_1.0.1            yulab.utils_0.0.2      digest_0.6.28          rmarkdown_2.11         cellranger_1.1.0      
[111] htmlTable_2.2.1        nortest_1.0-4          gld_2.6.2              curl_4.3.2             modeltools_0.2-23     
[116] jsonlite_1.7.2         lifecycle_1.0.1        nlme_3.1-152           carData_3.0-4          fansi_0.5.0           
[121] pillar_1.6.3           KEGGREST_1.32.0        fastmap_1.1.0          httr_1.4.2             glue_1.4.2            
[126] zip_2.2.0              png_0.1-7              bit_4.0.4              class_7.3-19           stringi_1.7.4         
[131] blob_1.2.2             latticeExtra_0.6-29    memoise_2.0.0          irlba_2.3.3            e1071_1.7-9 
```
