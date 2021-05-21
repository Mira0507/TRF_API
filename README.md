

## Utilization of Tandem Repeats Finder (TRF) - Mira Sohn 




### Aims

#### - This workflow is aimed at finding and counting repetitive elements based on ensembl gene id using the Tandem Repeats Finder (Unix/Linux ver.)

#### - References: [TRF Web](https://tandem.bu.edu/trf/trf.html), [TRF linux](https://tandem.bu.edu/trf/trf.unix.help.html), [TRF paper](https://pubmed.ncbi.nlm.nih.gov/9862982/), [Ensembl REST API](http://useast.ensembl.org/info/docs/api/index.html)


### Conda environment

#### - This analysis was performed under [conda](https://conda.io/projects/conda/en/latest/index.html) environment (see the [config](https://github.com/Mira0507/TRF_API/blob/master/conda.yaml))






### Loading packages


```r
library(data.table)
library(tidyverse)
library(httr)
library(jsonlite)
library(xml2)
```



### Importing gene ids 

#### - This demonstration uses ensembl human gene id 



```r
# Importing my genes to a vector 
geneid <- scan("geneid.txt", character(), quote="")
# Exploring the imported object 
geneid[1:10]
```

```
##  [1] "ENSG00000173267" "ENSG00000182985" "ENSG00000151136" "ENSG00000238045"
##  [5] "ENSG00000262814" "ENSG00000227619" "ENSG00000164062" "ENSG00000235142"
##  [9] "ENSG00000116194" "ENSG00000077380"
```

```r
class(geneid)
```

```
## [1] "character"
```

```r
length(geneid)
```

```
## [1] 12
```

### Retrieving cDNA sequences using API



#### - Tandem Repeats Finder (TRF): [Web](https://tandem.bu.edu/trf/trf.html), [Linux Version](https://tandem.bu.edu/trf/trf.unix.help.html), [Reference paper](https://pubmed.ncbi.nlm.nih.gov/9862982)


#### - Genomic sequences were retrieved by using the Ensembl  REST API [GET sequence/id/:id](https://rest.ensembl.org/documentation/info/sequence_id), [GET archive/id/:id](https://rest.ensembl.org/documentation/info/archive_id_get)



```r
# Assign the server for API 
server <- "https://rest.ensembl.org"



for (i in 1:length(geneid)) {

    # Test whether the gene is currently in the database
    ext.id <- paste0("/archive/id/", geneid[i], "?")

    r.id <- GET(paste(server, ext.id, sep = ""), content_type("application/json"))

    stop_for_status(r.id)

    # If the gene is currently in the database
    if (content(r.id)$is_current == "1") {

        # Retrieve data using ensembl gene ids 
        ext.seq <- paste0("/sequence/id/",
                      geneid[i],
                      "?")

        r.seq <- GET(paste(server, ext.seq, sep = ""), content_type("text/plain"))
         
        stop_for_status(r.seq)

        line1 <- paste0(">", geneid[i]) 
        line2 <- content(r.seq)

        # Save the gene id and sequences as fasta format 
        write(line1, file="trf_input.fa", append=T)
        write(line2, file="trf_input.fa", append=T)
    }
}

# Running tandem repeats finder (trf)
# Command: 
# trf <File> <Match> <Mismatch> <Delta> <PM> <PI> <Minscore> <MaxPeriod>
# -h: avoids to create html output
# -d: creates a .dat file 
system("trf trf_input.fa 2 7 7 80 10 50 500 -d -h")  # Returns a .dat file 
# Done.> pops up when it's done  
# "trf_input.fa.2.7.7.80.10.50.500.dat" (=output) is created in my working directory 
```



### Exploring analyzed repetitive elements 



```r
# Converting the TRF results from the .dat file to a data frame
# 
# 
# Isolating every single line of info to the trflines object
trflines <- readLines("trf_input.fa.2.7.7.80.10.50.500.dat")  

# Exploring the trflines object
class(trflines)
```

```
## [1] "character"
```

```r
length(trflines)
```

```
## [1] 374
```

```r
trflines[1:20]
```

```
##  [1] "Tandem Repeats Finder Program written by:"
##  [2] ""                                         
##  [3] "Gary Benson"                              
##  [4] "Program in Bioinformatics"                
##  [5] "Boston University"                        
##  [6] "Version 4.09"                             
##  [7] ""                                         
##  [8] ""                                         
##  [9] "Sequence: ENSG00000173267"                
## [10] ""                                         
## [11] ""                                         
## [12] ""                                         
## [13] "Parameters: 2 7 7 80 10 50 500"           
## [14] ""                                         
## [15] ""                                         
## [16] ""                                         
## [17] ""                                         
## [18] "Sequence: ENSG00000182985"                
## [19] ""                                         
## [20] ""
```

```r
# Creating an empty data frame
trfTable <- data.frame()

# Storing every single line to the data frame except empty ("") and 
# unnecessary (1 to 6) lines 
# (The line 1 to 6 of trflines have title & description about this run)
for (i in 7:length(trflines)) {

    if (trflines[i] != "" & !grepl("Parameters", trflines[i], fixed=T)) {

        trflines[i] <- str_replace(trflines[i], "Sequence: ", "")

        df <- data.frame(Line=trflines[i])

        trfTable <- rbind(trfTable, df)
    }
}

# Exploring the output data frame
head(trfTable)
```

```
##                                                                                                                                                   Line
## 1                                                                                                                                      ENSG00000173267
## 2                                                                                                                                      ENSG00000182985
## 3                                                                            2603 2630 1 28.0 1 100 0 56 100 0 0 0 0.00 A AAAAAAAAAAAAAAAAAAAAAAAAAAAA
## 4                                                       4990 5028 10 3.9 10 86 0 51 51 15 0 33 1.44 ATACATATAC ATACAAATACATACATATACATACTTATATATACATATA
## 5                                                                       11996 12025 2 15.0 2 100 0 60 0 0 50 50 1.00 GT GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
## 6 18689 18753 34 1.9 34 90 0 103 38 10 24 26 1.88 ATGAATCCTACAATGATAAGAAATTGAGAGGTTC ATGAATCCTGCAATGATAAGAAATTGAGAGGTTCATGAGTCCTACAATGATTAGAAATTGAGAGG
```

```r
dim(trfTable)
```

```
## [1] 280   1
```

```r
class(trfTable)
```

```
## [1] "data.frame"
```

```r
# Clean the TRF data frame 
# (removing rows containing parameters and adding line codes)
trfTable <- trfTable %>%
    mutate(code=ifelse(grepl("ENSG", Line, fixed=T), 
                       "ID", 
                       "CDNA")) 

# Exploring the output data frame
head(trfTable)
```

```
##                                                                                                                                                   Line
## 1                                                                                                                                      ENSG00000173267
## 2                                                                                                                                      ENSG00000182985
## 3                                                                            2603 2630 1 28.0 1 100 0 56 100 0 0 0 0.00 A AAAAAAAAAAAAAAAAAAAAAAAAAAAA
## 4                                                       4990 5028 10 3.9 10 86 0 51 51 15 0 33 1.44 ATACATATAC ATACAAATACATACATATACATACTTATATATACATATA
## 5                                                                       11996 12025 2 15.0 2 100 0 60 0 0 50 50 1.00 GT GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
## 6 18689 18753 34 1.9 34 90 0 103 38 10 24 26 1.88 ATGAATCCTACAATGATAAGAAATTGAGAGGTTC ATGAATCCTGCAATGATAAGAAATTGAGAGGTTCATGAGTCCTACAATGATTAGAAATTGAGAGG
##   code
## 1   ID
## 2   ID
## 3 CDNA
## 4 CDNA
## 5 CDNA
## 6 CDNA
```

```r
dim(trfTable)
```

```
## [1] 280   2
```

```r
# Extracting the number of repetitive elements per GENEID
#
# Creating an empty data frame 
trfnumTable <- data.frame()
code.vector <- trfTable$code

GENE <- c()
COUNT <- c()

index.id <- 0
index.cdna <- 0

for (i in 1:nrow(trfTable)) {

    if (code.vector[i] == "ID") {

        index.id <- index.id + 1
        index.cdna <- 0

    } else {

        index.id <- 0
        index.cdna <- index.cdna + 1
    }

    GENE[i] <- index.id # Indicates GENEID coded "ID"
    COUNT[i] <- index.cdna  # Indicates repetitive elements coded "CDNA" 
}

# Exploring the GENE and Count vectors
GENE
```

```
##   [1] 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
##  [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
##  [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [149] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [186] 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0
## [223] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1
## [260] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
```

```r
COUNT
```

```
##   [1]  0  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
##  [26] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
##  [51] 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73
##  [76] 74 75 76 77 78 79 80 81 82 83 84 85 86 87  0  1  2  3  4  5  6  7  8  9 10
## [101] 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
## [126] 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
## [151] 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85
## [176] 86 87 88 89 90 91 92 93 94 95 96 97  0  1  2  3  4  5  6  0  1  2  3  0  1
## [201]  2  3  4  5  6  7  8  9 10 11 12  0  1  0  1  2  3  4  5  6  7  8  9 10 11
## [226] 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
## [251] 37  0  1  2  3  4  5  6  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
## [276] 17 18 19 20  0
```

```r
# GENE[i] = 0: repetitive element line 
# GENE[i] > 0: gene id line
#   1  2  3  0  1  2  3  4  0  0  1
#         |  |           |  |  |
#     geneid |       geneid |  |
#           repeat1     repeat1|
#                            repeat2
#
#
#
# COUNT[i] = 0: gene id line
# COUNT[i] > 0: repetitive element line 
#   0  0  0  1  0  0  0  0  1  2  0 
#         |  |           |  |  |  
#     geneid |       geneid |  |
#          repeat1      repeat1|
#                           repeat2
# Assign vectors storing line indices of 
# gene id (repetitive.geneid), the first repetitive elements (repetitive.first), 
# and the last repetitive elements (repetitive.last)
repetitive.geneid <- which(COUNT == 1) - 1  # right before the first repetitive elements
repetitive.first <- which(COUNT == 1)  # the first repetitive elements 
repetitive.last <- which(GENE == 1) - 1  # the last repetitive elements per gene id

# Exploring the indices
repetitive.geneid
```

```
## [1]   2  90 188 195 199 212 214 252 259
```

```r
repetitive.first
```

```
## [1]   3  91 189 196 200 213 215 253 260
```

```r
repetitive.last
```

```
##  [1]   0  89 187 194 198 211 213 251 258 279
```

```r
# In case that the first and/or the last gene had one repetitive element: 
if (repetitive.last[1] < repetitive.first[1]) {

    repetitive.last <- repetitive.last[-1]

    if (repetitive.last[length(repetitive.last)] < 
        repetitive.first[length(repetitive.first)]) {

        repetitive.last <- c(repetitive.last, 
                             repetitive.first[length(repetitive.first)])
    }
    
}

# Explore the indices (all three vectors have to have the same length!)
repetitive.geneid
```

```
## [1]   2  90 188 195 199 212 214 252 259
```

```r
repetitive.first
```

```
## [1]   3  91 189 196 200 213 215 253 260
```

```r
repetitive.last
```

```
## [1]  89 187 194 198 211 213 251 258 279
```

```r
# Summarize the gene id (GENEID), indices of the first (Repeats_First) and 
# the last (Repeats_Last) repetitive elements in each gene id
trfnumTable <- data.frame(GENEID=trfTable$Line[repetitive.geneid],  
                          Repeats_First=repetitive.first, 
                          Repeats_Last=repetitive.last) %>% mutate(Repeats_Count=Repeats_Last - Repeats_First + 1,
                                   GENEID=str_replace(GENEID, "Sequence: ", "")) %>% dplyr::select(GENEID, Repeats_Count)

# Exploring the output data frame
head(trfnumTable)
```

```
##            GENEID Repeats_Count
## 1 ENSG00000182985            87
## 2 ENSG00000151136            97
## 3 ENSG00000238045             6
## 4 ENSG00000262814             3
## 5 ENSG00000227619            12
## 6 ENSG00000164062             1
```

```r
dim(trfnumTable)
```

```
## [1] 9 2
```

```r
# Explore the output data frame
glimpse(trfnumTable)
```

```
## Rows: 9
## Columns: 2
## $ GENEID        <chr> "ENSG00000182985", "ENSG00000151136", "ENSG00000238045",…
## $ Repeats_Count <dbl> 87, 97, 6, 3, 12, 1, 37, 6, 20
```

```r
# Export as a csv file
write.csv(trfnumTable, "TRF_demo.csv")
```




### Session Info


```r
sessionInfo()
```

```
## R version 4.0.5 (2021-03-31)
## Platform: x86_64-conda-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.2 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /home/mira/miniconda3/envs/trf_api/lib/libopenblasp-r0.3.15.so
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
##  [1] xml2_1.3.2        jsonlite_1.7.2    httr_1.4.2        forcats_0.5.1    
##  [5] stringr_1.4.0     dplyr_1.0.6       purrr_0.3.4       readr_1.4.0      
##  [9] tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.3     tidyverse_1.3.1  
## [13] data.table_1.14.0
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.1  xfun_0.23         haven_2.4.1       colorspace_2.0-1 
##  [5] vctrs_0.3.8       generics_0.1.0    htmltools_0.5.1.1 yaml_2.2.1       
##  [9] utf8_1.2.1        rlang_0.4.11      pillar_1.6.1      glue_1.4.2       
## [13] withr_2.4.2       DBI_1.1.1         dbplyr_2.1.1      modelr_0.1.8     
## [17] readxl_1.3.1      lifecycle_1.0.0   munsell_0.5.0     gtable_0.3.0     
## [21] cellranger_1.1.0  rvest_1.0.0       evaluate_0.14     knitr_1.33       
## [25] ps_1.6.0          curl_4.3.1        fansi_0.4.2       broom_0.7.6      
## [29] Rcpp_1.0.6        scales_1.1.1      backports_1.2.1   fs_1.5.0         
## [33] hms_1.1.0         digest_0.6.27     stringi_1.6.2     grid_4.0.5       
## [37] cli_2.5.0         tools_4.0.5       magrittr_2.0.1    crayon_1.4.1     
## [41] pkgconfig_2.0.3   ellipsis_0.3.2    reprex_2.0.0      lubridate_1.7.10 
## [45] assertthat_0.2.1  rmarkdown_2.8     rstudioapi_0.13   R6_2.5.0         
## [49] compiler_4.0.5
```
