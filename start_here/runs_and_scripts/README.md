
# Read before running scripts

If you are interested in following what we did to run Studies 1 and 2, and Auxiliary Study 2e, or if you are interested in applying our R scripts to your own data, this is the right place to start. A few things to know before starting:

* Make sure to read the top-level README file before continuing to read this file.
* Start with the script for Study 1, as it has the most comments in it which document what the code is doing. Similar structures and code is used in the other R scripts.
* In the code 'SETAL' (Stack Harrington et al., 2018) refers to what we call 'HS18' in the main text and 'FETAL' (Fine et al., 2013) refers to 'F13' in the main text. The third data set is referred to as NSC (Natural Story Corpus, Futrell et al., 2018) in both the R code and the main corpus.

If you discover a bug in our code, or think that additional comments would be helpful, please feel free to submit a pull request to the github repo at https://github.com/hlplab/BurchillJaegerRTAnalysis/. If you have questionsabout the code, please reach out to the code developer, Zach Burchill (zach.burchill@gmail.com). For conceptual questions about the purpose of the code or other questions about the paper that this code was used for, please Florian Jaeger (fjaeger@ur.rochester.edu).

# Session info

We have tested the R scripts in this folder on two environments. On the first environment, Studies 1 and 2 each took about 1 week of compute time to complete. 

1) 2013 MacPro with 24 cores and 64GB RAM. 

```
R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Big Sur 11.6.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Toronto
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base    

other attached packages:
 [1] lme4_1.1-33         Matrix_1.5-4.1      broom.mixed_0.2.9.4 beepr_1.3           cs_0.1.7          
 [6] furrr_0.3.1         future_1.32.0       lubridate_1.9.2     forcats_1.0.0       stringr_1.5.0      
[11] dplyr_1.1.2         purrr_1.0.1         readr_2.1.4         tidyr_1.3.0         tibble_3.2.1      
[16] ggplot2_3.4.2       tidyverse_2.0.0    

loaded via a namespace (and not attached):
 [1] utf8_1.2.3        generics_0.1.3    stringi_1.7.12    lattice_0.21-8    listenv_0.9.0     hms_1.1.3        
 [7] digest_0.6.31     magrittr_2.0.3    grid_4.3.0        timechange_0.2.0  zplyr_0.2.3.1     backports_1.4.1  
[13] audio_0.1-10      fansi_1.0.4       scales_1.2.1      codetools_0.2-19  cli_3.6.1         rlang_1.1.1      
[19] parallelly_1.36.0 splines_4.3.0     munsell_0.5.0     withr_2.5.0       tools_4.3.0       parallel_4.3.0  
[25] tzdb_0.4.0        nloptr_2.0.3      minqa_1.2.5       colorspace_2.1-0  boot_1.3-28.1     globals_0.16.2  
[31] broom_1.0.5       vctrs_0.6.2       R6_2.5.1          lifecycle_1.0.3   MASS_7.3-60       pkgconfig_2.0.3  
[37] pillar_1.9.0      gtable_0.3.3      Rcpp_1.0.10       glue_1.6.2        tidyselect_1.2.0  rstudioapi_0.14  
[43] nlme_3.1-162      compiler_4.3.0 
```

3) 2020 Macbook Pro

```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lme4_1.1-25       Matrix_1.2-18     broom.mixed_0.2.6 beepr_1.3        
 [5] cs_0.1.7          furrr_0.1.0       future_1.19.1     forcats_0.5.0    
 [9] stringr_1.4.0     dplyr_1.0.5       purrr_0.3.4       readr_1.3.1      
[13] tidyr_1.1.2       tibble_3.0.3.9000 ggplot2_3.4.3     tidyverse_1.3.0  

loaded via a namespace (and not attached):
 [1] httr_1.4.2        jsonlite_1.7.1    splines_3.6.1     modelr_0.1.8     
 [5] assertthat_0.2.1  statmod_1.4.34    blob_1.2.1        cellranger_1.1.0 
 [9] yaml_2.2.1        remotes_2.2.0     globals_0.13.0    pillar_1.4.6     
[13] backports_1.1.10  lattice_0.20-41   glue_1.4.2        digest_0.6.25    
[17] rvest_0.3.6       minqa_1.2.4       colorspace_1.4-1  plyr_1.8.6       
[21] pkgconfig_2.0.3   broom_0.7.2       listenv_0.8.0     haven_2.3.1      
[25] scales_1.2.1      processx_3.4.4    generics_0.0.2    ellipsis_0.3.1   
[29] withr_2.5.0       TMB_1.7.18        cli_3.6.1         magrittr_1.5     
[33] crayon_1.3.4      readxl_1.3.1      ps_1.3.4          fs_1.5.0         
[37] nlme_3.1-149      MASS_7.3-53       xml2_1.3.2        pkgbuild_1.1.0   
[41] tools_3.6.1       prettyunits_1.1.1 hms_0.5.3         lifecycle_1.0.3  
[45] munsell_0.5.0     reprex_0.3.0      callr_3.4.4       compiler_3.6.1   
[49] rlang_1.1.1       nloptr_1.2.2.2    rstudioapi_0.11   zplyr_0.2.3.1    
[53] boot_1.3-25       gtable_0.3.0      codetools_0.2-16  DBI_1.1.0        
[57] curl_4.3          reshape2_1.4.4    R6_2.4.1          lubridate_1.7.9  
[61] rprojroot_1.3-2   stringi_1.5.3     parallel_3.6.1    Rcpp_1.0.5       
[65] vctrs_0.6.3       dbplyr_1.4.4      audio_0.1-7       tidyselect_1.1.0 
[69] coda_0.19-4   

```
