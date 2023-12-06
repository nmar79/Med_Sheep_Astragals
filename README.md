# Med_Sheep_Astragals
R script and data input files, supplement to Harding et al.

"gmm_pipeline()" is a wrapper function that executes the geometric morphometric analysis presented in the manuscript "The contribution of Mediterranean connectivity to morphological variability in Iron Age sheep of the Eastern Mediterranean" by Sierra Harding, Angelos Hadjikoumis, Shyama Vermeersch & Nimrod Marom. The code was written by NM (nimrod.arch@gmail.com).

# Instructions

The function requires the packages 'geomorph', 'Morpho', 'ape', 'ggplot' and 'ggsci'. 
run the gmm_pipeline() function, and then

gmm_pipeline("sheep_specinfo_20230824.csv", "sheep_sliders_202212180739.csv", "sheep_lmrks_20230824.TPS", c(1:25), library_check = FALSE)

(make sure that the csv and TPS files are in the working directory...)

The function will output results to the console, to the "results" folder, and to the "plots" tab in R-studio. 

# Session Info

The integrity of the code and the output was checked for the last time on 6-12-2023. The session info was:

R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 23.10

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.11.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Jerusalem
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Morpho_2.11    ape_5.7-1      ggsci_3.0.0    ggplot2_3.4.4 
[5] geomorph_4.0.6 Matrix_1.6-1   rgl_1.2.1      RRPP_1.4.0    
