
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

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$input.dir,i,j)) %>% .[,anno:=as.factor(j)]
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[complete.cases(.)]

##############
## Barplots ##
##############

tmp <- diff.results %>%
  .[,.(number_positive_hits=sum(sig==T & diff>0), number_negative_hits=sum(sig==T & diff<0)), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(tmp$value), max(tmp$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim)
  
  # pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=6, height=4)
  png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
