
#############################################################################
## Script to plot the results from the differential methylation analysis ##
#############################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

source("/Users/ricard/human_embryo_multiomics/met/differential/feature_level/analysis/utils.R")

#####################
## Define settings ##
#####################

## I/O ##
io <- list()
io$input.dir <- "/Users/ricard/data/human_embryo_multiomics/met/differential/feature_level"
io$outdir <- "/Users/ricard/data/human_embryo_multiomics/met/differential/feature_level/pdf"
dir.create(io$outdir, showWarnings = F)

## Options ##
opts <- list()

opts$comparisons <- c(
  "E56ICM_vs_E56TE"
)

# Select genomic contexts
opts$annos <- c(
  "genebody",
  "prom_200_200_cgi",
  "prom_200_200_noncgi",
  "prom_200_200",
  "prom_2000_2000",
  "prom_2000_2000_cgi",
  "prom_2000_2000_noncgi",
  "H1_H3K27ac",
  "H1_H3K27me3",
  "H1_H3K4me1",
  "H1_H3K4me3"
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
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$input.dir,i,j)) %>% .[,anno:=as.factor(j)]
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
