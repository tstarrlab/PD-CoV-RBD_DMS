Collapse barcodes to final per-RBD/mutant phenotype scores
================
Tyler Starr
01/02/2023

- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#calculate-per-variant-mean-scores-within-replicates"
  id="toc-calculate-per-variant-mean-scores-within-replicates">Calculate
  per-variant mean scores within replicates</a>
- <a href="#calculate-per-mutant-score-across-libraries"
  id="toc-calculate-per-mutant-score-across-libraries">Calculate
  per-mutant score across libraries</a>
- <a
  href="#correlations-among-backgrounds-and-to-prior-wuhan-hu-1-dms-data"
  id="toc-correlations-among-backgrounds-and-to-prior-wuhan-hu-1-dms-data">Correlations
  among backgrounds and to prior Wuhan-Hu-1 DMS data</a>
- <a href="#heatmaps" id="toc-heatmaps">Heatmaps!</a>

This notebook reads in the per-barcode Kd and MFI measurements from the
`compute_xAPN` scripts. It synthesizes these two sets of results and
calculates the final ‘mean’ phenotypes for each variant, and generates
some coverage and QC analyses.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$final_variant_scores_dir)){
  dir.create(file.path(config$final_variant_scores_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [5] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0       tibble_3.0.2     
    ##  [9] ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0       
    ## [13] knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.5     compiler_3.6.2  
    ##  [5] dbplyr_1.4.2     tools_3.6.2      digest_0.6.23    lubridate_1.7.4 
    ##  [9] jsonlite_1.6     evaluate_0.14    lifecycle_0.2.0  gtable_0.3.0    
    ## [13] pkgconfig_2.0.3  rlang_0.4.7      reprex_0.3.0     cli_2.0.0       
    ## [17] rstudioapi_0.10  DBI_1.1.0        haven_2.2.0      xfun_0.11       
    ## [21] withr_2.1.2      xml2_1.3.3       httr_1.4.1       fs_1.3.1        
    ## [25] hms_0.5.2        generics_0.0.2   vctrs_0.3.1      grid_3.6.2      
    ## [29] tidyselect_1.1.0 glue_1.3.1       R6_2.4.1         fansi_0.4.0     
    ## [33] readxl_1.3.1     rmarkdown_2.0    modelr_0.1.5     magrittr_1.5    
    ## [37] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   htmltools_0.4.0 
    ## [41] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
    ## [45] munsell_0.5.0    broom_0.7.0      crayon_1.3.4

## Setup

Read in tables of per-barcode expression and binding Kd measurements and
combine.

``` r
dt_gAPN <- data.table(read.csv(config$gAPN_Kds_file),stringsAsFactors=F)
dt_gAPN[library=="lib55A",library:="lib55"]
dt_gAPN[library=="lib56A",library:="lib56"]
setnames(dt_gAPN,"TiteSeq_avgcount","gAPN_count")
dt_hAPN <- data.table(read.csv(config$hAPN_meanF_file),stringsAsFactors=F)
dt_pAPN <- data.table(read.csv(config$pAPN_meanF_file),stringsAsFactors=F)
```

## Calculate per-variant mean scores within replicates

Calculate the mean binding and expression score collapsed by genotype.
Also output the number of barcodes across which a variant score was
determined in each library.

``` r
dt_gAPN[is.na(log10Ka_gAPN),gAPN_count:=NA]
dt_hAPN[is.na(hAPN_meanF),hAPN_count:=NA]
dt_pAPN[is.na(pAPN_meanF),pAPN_count:=NA]

dt_gAPN[,mean_bind_gAPN:=mean(log10Ka_gAPN,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_gAPN[,sd_bind_gAPN:=sd(log10Ka_gAPN,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_gAPN[,n_bc_bind_gAPN:=sum(!is.na(log10Ka_gAPN)),by=c("library","target","variant_class","aa_substitutions")]
dt_gAPN[,avg_count_bind_gAPN:=mean(gAPN_count,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_gAPN <- unique(dt_gAPN[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_bind_gAPN,sd_bind_gAPN,n_bc_bind_gAPN,avg_count_bind_gAPN)])

dt_hAPN[,mean_bind_hAPN:=mean(hAPN_meanF,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_hAPN[,sd_bind_hAPN:=sd(hAPN_meanF,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_hAPN[,n_bc_bind_hAPN:=sum(!is.na(hAPN_meanF)),by=c("library","target","variant_class","aa_substitutions")]
dt_hAPN[,avg_count_bind_hAPN:=mean(hAPN_count,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_hAPN <- unique(dt_hAPN[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_bind_hAPN,sd_bind_hAPN,n_bc_bind_hAPN,avg_count_bind_hAPN)])

dt_pAPN[,mean_bind_pAPN:=mean(pAPN_meanF,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_pAPN[,sd_bind_pAPN:=sd(pAPN_meanF,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_pAPN[,n_bc_bind_pAPN:=sum(!is.na(pAPN_meanF)),by=c("library","target","variant_class","aa_substitutions")]
dt_pAPN[,avg_count_bind_pAPN:=mean(pAPN_count,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_pAPN <- unique(dt_pAPN[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_bind_pAPN,sd_bind_pAPN,n_bc_bind_pAPN,avg_count_bind_pAPN)])
```

Some QC plots. First, look at distribution of number barcodes for
binding measurements for single mutant detemrinations. These are
‘left-justified’ histograms, so the leftmost bar represents the number
of genotypes for which no barcodes were collapsed to final measurement
in a pool.

``` r
par(mfrow=c(3,2))
hist(dt_gAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],main="lib55, gAPN bind",right=F,breaks=max(dt_gAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],na.rm=T),xlab="")
hist(dt_gAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],main="lib56, gAPN bind",right=F,breaks=max(dt_gAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],na.rm=T),xlab="")

hist(dt_hAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],main="lib55, hAPN bind",right=F,breaks=max(dt_hAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],na.rm=T),xlab="")
hist(dt_hAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],main="lib56, hAPN bind",right=F,breaks=max(dt_hAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],na.rm=T),xlab="")

hist(dt_pAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],main="lib55, pAPN bind",right=F,breaks=max(dt_pAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],na.rm=T),xlab="")
hist(dt_pAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],main="lib56, pAPN bind",right=F,breaks=max(dt_pAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],na.rm=T),xlab="")
```

<img src="collapse_scores_files/figure-gfm/hist_n_bc_per_mutant-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_sep-libs.pdf",sep=""),useDingbats=F))
```

What about how SEM tracks with number of barcodes collapsed? This could
help for choosing a minimum number of barcodes to use.

``` r
par(mfrow=c(3,2))
plot(dt_gAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],
     dt_gAPN[library=="lib55" & variant_class=="1 nonsynonymous",sd_bind_gAPN/sqrt(n_bc_bind_gAPN)],
     pch=16,col="#00000005",main="lib55, gAPN bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_gAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_gAPN],
     dt_gAPN[library=="lib56" & variant_class=="1 nonsynonymous",sd_bind_gAPN/sqrt(n_bc_bind_gAPN)],
     pch=16,col="#00000005",main="lib56, gAPN bind",ylab="SEM",xlab="number barcodes collapsed")

plot(dt_hAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],
     dt_hAPN[library=="lib55" & variant_class=="1 nonsynonymous",sd_bind_hAPN/sqrt(n_bc_bind_hAPN)],
     pch=16,col="#00000005",main="lib55, hAPN bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_hAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_hAPN],
     dt_hAPN[library=="lib56" & variant_class=="1 nonsynonymous",sd_bind_hAPN/sqrt(n_bc_bind_hAPN)],
     pch=16,col="#00000005",main="lib56, hAPN bind",ylab="SEM",xlab="number barcodes collapsed")

plot(dt_pAPN[library=="lib55" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],
     dt_pAPN[library=="lib55" & variant_class=="1 nonsynonymous",sd_bind_pAPN/sqrt(n_bc_bind_pAPN)],
     pch=16,col="#00000005",main="lib55, pAPN bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_pAPN[library=="lib56" & variant_class=="1 nonsynonymous",n_bc_bind_pAPN],
     dt_pAPN[library=="lib56" & variant_class=="1 nonsynonymous",sd_bind_pAPN/sqrt(n_bc_bind_pAPN)],
     pch=16,col="#00000005",main="lib56, pAPN bind",ylab="SEM",xlab="number barcodes collapsed")
```

<img src="collapse_scores_files/figure-gfm/sem_v_n-bc-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/sem_v_n-bc.pdf",sep=""),useDingbats=F))
```

Format into a ‘mutation lookup table’, where we focus just on the single
mutants (and wildtype), breakup the string of mutations, and fill in the
table to also include any missing mutants.

``` r
dt <- merge(merge(dt_gAPN, dt_hAPN), dt_pAPN)

dt_mutant <- dt[variant_class %in% "1 nonsynonymous",]

#split mutation string
#define function to apply
split_mut <- function(x){
  split <- strsplit(x,split="")[[1]]
  return(list(split[1],as.numeric(paste(split[2:(length(split)-1)],collapse="")),split[length(split)]))
}
dt_mutant[,c("wildtype","position","mutant"):=split_mut(as.character(aa_substitutions)),by=aa_substitutions]

dt_mutant <- dt_mutant[,.(library,target,wildtype,position,mutant,
                          mean_bind_gAPN,sd_bind_gAPN,n_bc_bind_gAPN,avg_count_bind_gAPN,
                          mean_bind_hAPN,sd_bind_hAPN,n_bc_bind_hAPN,avg_count_bind_hAPN,
                          mean_bind_pAPN,sd_bind_pAPN,n_bc_bind_pAPN,avg_count_bind_pAPN)]

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#fill out missing values in table with a hideous loop, so the table is complete for all mutaitons (including those that are missing). If you are somebody who is reading this code, I apologize.
for(lib in c("lib55","lib56")){
  for(bg in as.character(unique(dt_mutant$target))){
    for(pos in 1:max(dt_mutant$position)){
      for(aa in aas){
        if(!(aa %in% as.character(dt_mutant[library==lib & target==bg & position==pos,mutant]))){
          dt_mutant <- rbind(dt_mutant,list(lib, bg, dt_mutant[target==bg & position==pos,wildtype][1],pos,aa),fill=T) #note this will leave NA for wildtype if a position is completely missing in both libraries
        }
      }
    }
  }
}
setkey(dt_mutant,library,target,position,mutant)

#fill in wildtype values -- should vectorize in data table but being so stupid so just going to write for loop
for(bg in c("PDCoV")){
  for(lib in c("lib55","lib56")){
    dt_mutant[library==lib & target==bg & wildtype==mutant, c("mean_bind_gAPN","sd_bind_gAPN","n_bc_bind_gAPN","avg_count_bind_gAPN",
                                                              "mean_bind_hAPN","sd_bind_hAPN","n_bc_bind_hAPN","avg_count_bind_hAPN",
                                                              "mean_bind_pAPN","sd_bind_pAPN","n_bc_bind_pAPN","avg_count_bind_pAPN") := 
                dt[library==lib & target==bg & variant_class=="wildtype",
                        .(mean_bind_gAPN,sd_bind_gAPN,n_bc_bind_gAPN,avg_count_bind_gAPN,
                          mean_bind_hAPN,sd_bind_hAPN,n_bc_bind_hAPN,avg_count_bind_hAPN,
                          mean_bind_pAPN,sd_bind_pAPN,n_bc_bind_pAPN,avg_count_bind_pAPN)]]
  }
}

#add delta bind measures
for(bg in c("PDCoV")){
  for(lib in c("lib55","lib56")){
    ref_bind_gAPN <- dt[library==lib & target==bg & variant_class=="wildtype",mean_bind_gAPN]
    dt_mutant[library==lib & target==bg,delta_bind_gAPN := mean_bind_gAPN - ref_bind_gAPN]
    ref_bind_hAPN <- dt[library==lib & target==bg & variant_class=="wildtype",mean_bind_hAPN]
    dt_mutant[library==lib & target==bg,delta_bind_hAPN := mean_bind_hAPN - ref_bind_hAPN]
    ref_bind_pAPN <- dt[library==lib & target==bg & variant_class=="wildtype",mean_bind_pAPN]
    dt_mutant[library==lib & target==bg,delta_bind_pAPN := mean_bind_pAPN - ref_bind_pAPN]
  }
}
```

We have duplicates for each measurement. Let’s look at correlations!

``` r
par(mfrow=c(1,3))
x <- dt_mutant[library=="lib55" & wildtype!=mutant,mean_bind_gAPN]; y <- dt_mutant[library=="lib56" & wildtype!=mutant,mean_bind_gAPN]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="gAPN binding");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_mutant[library=="lib55" & wildtype!=mutant,mean_bind_hAPN]; y <- dt_mutant[library=="lib56" & wildtype!=mutant,mean_bind_hAPN]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="hAPN binding");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_mutant[library=="lib55" & wildtype!=mutant,mean_bind_pAPN]; y <- dt_mutant[library=="lib56" & wildtype!=mutant,mean_bind_pAPN]; plot(x,y,pch=16,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="pAPN binding");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/replicate_correlations.pdf",sep=""),useDingbats=F))
```

## Calculate per-mutant score across libraries

Collapse down to mean from both replicates, and total n barcodes between
the two/three replicates. Also record the number of the replicates the
variant was quantified within. Note, we are currently keeping a value
even if it’s determined from a single bc fit in a single pool. Later on,
we may want to require some combination of minimum number of bcs within
or between libraries for retention.

``` r
dt_final <- copy(dt_mutant)

dt_final[ ,bind_gAPN:=mean(mean_bind_gAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,delta_bind_gAPN:=mean(delta_bind_gAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_bc_gAPN:=sum(n_bc_bind_gAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_libs_gAPN:=sum(!is.na(mean_bind_gAPN)),by=c("target","position","mutant")]

dt_final[ ,bind_hAPN:=mean(mean_bind_hAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,delta_bind_hAPN:=mean(delta_bind_hAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_bc_hAPN:=sum(n_bc_bind_hAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_libs_hAPN:=sum(!is.na(mean_bind_hAPN)),by=c("target","position","mutant")]

dt_final[ ,bind_pAPN:=mean(mean_bind_pAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,delta_bind_pAPN:=mean(delta_bind_pAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_bc_pAPN:=sum(n_bc_bind_pAPN,na.rm=T),by=c("target","position","mutant")]
dt_final[ ,n_libs_pAPN:=sum(!is.na(mean_bind_pAPN)),by=c("target","position","mutant")]

#switch to spike indexing of postitions
dt_final$position <- dt_final$position + config$site_number_offset

#add single mutation string
dt_final[,mutation:=paste(wildtype,position,mutant,sep=""),by=c("wildtype","position","mutant")]

dt_final <- unique(dt_final[,.(target,wildtype,position,mutant,mutation,
                               bind_gAPN,delta_bind_gAPN,n_bc_gAPN,n_libs_gAPN,
                               bind_hAPN,delta_bind_hAPN,n_bc_hAPN,n_libs_hAPN,
                               bind_pAPN,delta_bind_pAPN,n_bc_pAPN,n_libs_pAPN)])

setkey(dt_final,target,position,mutant)

#add the rep1 and rep2 sub-averages
dt_final[,gAPN_rep1 := dt_mutant[library=="lib55", mean_bind_gAPN]]
dt_final[,gAPN_rep2 := dt_mutant[library=="lib56", mean_bind_gAPN]]

dt_final[,hAPN_rep1 := dt_mutant[library=="lib55", mean_bind_hAPN]]
dt_final[,hAPN_rep2 := dt_mutant[library=="lib56", mean_bind_hAPN]]

dt_final[,pAPN_rep1 := dt_mutant[library=="lib55", mean_bind_pAPN]]
dt_final[,pAPN_rep2 := dt_mutant[library=="lib56", mean_bind_pAPN]]
```

Censor any measurements that are from \<3 bc or only sampled in a single
replicate? Don’t do this for now.

``` r
# min_bc <- 2
# min_lib <- 2
# 
# dt_final[n_bc_bind < min_bc & n_libs_bind < min_lib, c("bind","delta_bind","n_bc_bind","n_libs_bind") := list(NA,NA,NA,NA)]
# dt_final[n_bc_expr < min_bc & n_libs_expr < min_lib, c("expr","delta_expr","n_bc_expr","n_libs_expr") := list(NA,NA,NA,NA)]
```

Coverage stats on n_barcodes for different measurements in the final
pooled measurements.

``` r
par(mfrow=c(1,3))
hist(dt_final[wildtype!=mutant, n_bc_gAPN],col="gray50",main=paste("mutant gAPN bind score,\nmedian ",median(dt_final[wildtype!=mutant, n_bc_gAPN],na.rm=T),sep=""),right=F,breaks=max(c(dt_final[wildtype!=mutant, n_bc_gAPN],dt_final[wildtype!=mutant, n_bc_hAPN],dt_final[wildtype!=mutant, n_bc_pAPN]))/2,xlab="number barcodes", xlim=c(0,100))

hist(dt_final[wildtype!=mutant, n_bc_hAPN],col="gray50",main=paste("mutant hAPN bind score,\nmedian ",median(dt_final[wildtype!=mutant, n_bc_hAPN],na.rm=T),sep=""),right=F,breaks=max(c(dt_final[wildtype!=mutant, n_bc_gAPN],dt_final[wildtype!=mutant, n_bc_hAPN],dt_final[wildtype!=mutant, n_bc_pAPN]))/2,xlab="number barcodes", xlim=c(0,100))

hist(dt_final[wildtype!=mutant, n_bc_pAPN],col="gray50",main=paste("mutant pAPN bind score,\nmedian ",median(dt_final[wildtype!=mutant, n_bc_pAPN],na.rm=T),sep=""),right=F,breaks=max(c(dt_final[wildtype!=mutant, n_bc_gAPN],dt_final[wildtype!=mutant, n_bc_hAPN],dt_final[wildtype!=mutant, n_bc_pAPN]))/2,xlab="number barcodes", xlim=c(0,100))
```

<img src="collapse_scores_files/figure-gfm/n_barcode_plots-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_pooled-libs.pdf",sep="")))
```

## Correlations among backgrounds and to prior Wuhan-Hu-1 DMS data

Look at correlations in mutation effects between each background, for
bind phenotype

``` r
par(mfrow=c(2,2))

x <- dt_final[,bind_gAPN]; y <- dt_final[,bind_hAPN]; plot(x,y,pch=16,col="#00000020",xlab="galline APN (log10Ka)",ylab="human APN (MFI)",main="", xlim=c(4,10.5),ylim=c(4,10.5));model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[,bind_gAPN]; y <- dt_final[,bind_pAPN]; plot(x,y,pch=16,col="#00000020",xlab="galline APN (log10Ka)",ylab="porcine APN (MFI)",main="", xlim=c(4,10.5),ylim=c(4,10.5));model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[,bind_hAPN]; y <- dt_final[,bind_pAPN]; plot(x,y,pch=16,col="#00000020",xlab="human APN (MFI)",ylab="porcine APN (MFI)",main="", xlim=c(4,10.5),ylim=c(4,10.5));model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_by_bg_bind-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/background_correlations_bind.pdf",sep=""),useDingbats=F))
```

## Heatmaps!

Order factor variables for plotting

``` r
#order mutant as a factor for grouping by rough biochemical grouping
dt_final$mutant <- factor(dt_final$mutant, levels=c("C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
#add character vector indicating wildtype to use as plotting symbols for wt
dt_final[,wildtype_indicator := ""]
dt_final[as.character(mutant)==as.character(wildtype),wildtype_indicator := "x"]

#make temp long-form data frame
temp <- data.table::melt(dt_final[, .(target,position,mutant,
                                      bind_gAPN,delta_bind_gAPN,
                                      bind_hAPN,delta_bind_hAPN,
                                      bind_pAPN,delta_bind_pAPN,
                                      wildtype_indicator)],
                         id.vars=c("target","position","mutant","wildtype_indicator"),
                         measure.vars=c("bind_gAPN","delta_bind_gAPN",
                                        "bind_hAPN","delta_bind_hAPN",
                                        "bind_pAPN","delta_bind_pAPN"),
                         variable.name="measurement",value.name="value")

#for method to duplicate aa labels on right side of plot https://github.com/tidyverse/ggplot2/issues/3171
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}
```

Make heatmaps showing raw affinity and delta-affinity of muts relative
to wildtype for each of the three APNs. First, log10Ka for galline APN

``` r
p1 <- ggplot(temp[measurement=="bind_gAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,10.5),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_log10Ka_gAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_gAPN_log10Ka.pdf",sep="")))
```

Second, illustrating delta_log10Ka for gAPN

``` r
p1 <- ggplot(temp[measurement=="delta_bind_gAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-3,2.05),values=c(0/5.05,1/5.05,2/5.05,3/5.05,4/5.05,5.05/5.05),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-log10Ka_gAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_gAPN_delta-log10Ka.pdf",sep="")))
```

Next, MFI metrics for hAPN

``` r
p1 <- ggplot(temp[measurement=="bind_hAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(4,9.5),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_MFI_hAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_hAPN_MFI.pdf",sep="")))
```

``` r
p1 <- ggplot(temp[measurement=="delta_bind_hAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-2,3),values=c(0/5,1/5,2/5,3.5/5,5/5),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-MFI-hAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_hAPN_delta-MFI.pdf",sep="")))
```

Last, MFI metrics for pAPN

``` r
p1 <- ggplot(temp[measurement=="bind_pAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(4,9.5),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_MFI_pAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_pAPN_MFI.pdf",sep="")))
```

``` r
p1 <- ggplot(temp[measurement=="delta_bind_pAPN",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-1,4.5),values=c(0/5.5,0.5/5.5,1/5.5,3.25/5.5,5.5/5.5),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(303,seq(305,416,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-MFI-pAPN-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_pAPN_delta-MFI.pdf",sep="")))
```

Save output files.

``` r
dt_final[,.(target,wildtype,position,mutant,mutation,
            bind_gAPN,delta_bind_gAPN,n_bc_gAPN,n_libs_gAPN,gAPN_rep1,gAPN_rep2,
            bind_hAPN,delta_bind_hAPN,n_bc_hAPN,n_libs_hAPN,hAPN_rep1,hAPN_rep2,
            bind_pAPN,delta_bind_pAPN,n_bc_pAPN,n_libs_pAPN,pAPN_rep1,pAPN_rep2)] %>%
  mutate_if(is.numeric, round, digits=5) %>%
  write.csv(file=config$final_variant_scores_mut_file, row.names=F,quote=F)
```
