
#############################################################################
## Script to plot the results from the differential accessibility analysis ##
#############################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

source("/Users/ricard/human_embryo_multiomics/acc/differential/feature_level/analysis/utils.R")

#####################
## Define settings ##
#####################

## I/O ##
io <- list()
io$input.dir <- "/Users/ricard/data/human_embryo_multiomics/acc/differential/feature_level"
io$outdir <- "/Users/ricard/data/human_embryo_multiomics/acc/differential/feature_level/pdf"
dir.create(io$outdir, showWarnings = F)

## Options ##
opts <- list()

opts$comparisons <- c(
  "E6ICMEPI_vs_E6TE"
)

# Select genomic contexts
opts$annos <- c(
  "2cell_H3K27me3",
  "4cell_3PN_H3K4me3",
  "4cell_H3K27me3",
  "8cell_3PN_H3K27me3",
  "8cell_H3K27ac",
  "CGI",
  "H1_H3K27ac",
  "H1_H3K27me3",
  "H1_H3K4me1",
  "H1_H3K4me3",
  "H1_distal_H3K27ac",
  "ICM_H3K27ac",
  "ICM_H3K4me3",
  "TE_H3K27me3",
  "atac_peaks_2cell_2pn",
  "atac_peaks_2cell_3pn",
  "atac_peaks_4cell_3pn",
  "atac_peaks_8cell_2pn",
  "atac_peaks_8cell_3pn",
  "atac_peaks_icm_2pn",
  "distal_4cell_3PN_H3K4me3",
  "distal_8cell_H3K27ac",
  "distal_ICM_H3K27ac",
  "distal_ICM_H3K4me3",
  "genebody",
  "prom_2000_2000",
  "prom_2000_2000_cgi",
  "prom_2000_2000_noncgi",
  "prom_200_200",
  "prom_200_200_cgi",
  "prom_200_200_noncgi",
  "proximal_4cell_3PN_H3K4me3",
  "proximal_8cell_H3K27ac",
  "proximal_ICM_H3K27ac",
  "proximal_ICM_H3K4me3"
)

# Minimum differential levels (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j) {
    file <- sprintf("%s/%s_%s.txt.gz",io$input.dir,i,j)
    if (file.exists(file)) {
      fread(file) %>% .[,anno:=as.factor(j)]
    } else {
      cat(sprintf("%s does not exist",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[complete.cases(.)]



# Select statistically significant hits
diff.results.[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)]

###################
## Volcano plots ##
###################

for (i in unique(diff.results$comparison)) {
  for (j in unique(diff.results$anno)) {
  
    p <- gg_volcano_plot(diff.results[comparison==i & anno==j], top_genes = 0)
    
    # pdf(sprintf("%s/volcano_%s_%s.pdf",io$outdir,i,j), width=7, height=5)
    png(sprintf("%s/volcano_%s_%s.png",io$outdir,i,j), width=7, height=5, units="in", res=400)
    print(p)
    dev.off()
  }
}
