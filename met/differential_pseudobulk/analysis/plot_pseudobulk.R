
########################################################################################
## Script to plot the results from the pseudobulked differential methylation analysis ##
########################################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

source("/Users/ricard/Guo_2017/differential/differential_pseudobulk/utils.R")

#####################
## Define settings ##
#####################

## I/O ##
io <- list()
io$input.dir <- "/Users/ricard/data/Guo_2017/met/differential/pseudobulk"
io$outdir <- "/Users/ricard/data/Guo_2017/met/differential/pseudobulk/pdf"

## Options ##
opts <- list()
opts$comparisons <- c(
  "ICM_vs_TE"
)

# Select genomic contexts
opts$annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "prom_2000_2000",
  "window2000_step1000",
  "CGI",
  "genebody"
)

# opts$colors <- c(
#   "H3K27ac_distal_E7.5_Ect_intersect12_500" = "steelblue",
#   "H3K27ac_distal_E7.5_End_intersect12_500" = "#43CD80",
#   "H3K27ac_distal_E7.5_Mes_intersect12_500" = "violetred",
#   "prom_2000_2000" = "grey60"
# )


###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$input.dir,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[complete.cases(.)]

# diff.results[,c("rate_A","rate_B","diff"):=list(rate_A*100,rate_B*100,diff*100)]

##############
## Barplots ##
##############

opts$min.diff <- 15

tmp <- diff.results %>%
  .[,.(number_positive_hits=mean(diff>opts$min.diff), 
       number_negative_hits=mean(diff<(-opts$min.diff))), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(tmp$value), max(tmp$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim)
  
  # p <- p + scale_fill_manual(values=opts$colors)
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=6, height=4)
  print(p)
  # dev.off()
}
