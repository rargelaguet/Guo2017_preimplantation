
######################################################################################
## Script to quantify the fraction of differential hits that are distal vs proximal ##
######################################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

#####################
## Define settings ##
#####################

## I/O ##
io <- list()
io$tss <- "/Users/ricard/data/hg38_regulation/promoters/TSS.bed"
io$input.dir <- "/Users/ricard/data/human_embryo_multiomics/acc/differential/feature_level"
io$outdir <- "/Users/ricard/data/human_embryo_multiomics/acc/differential/feature_level/link_to_genes"
dir.create(io$outdir, showWarnings = F)

## Options ##
opts <- list()
opts$comparisons <- c(
  "E6ICM_vs_E6TE"
)

# Select genomic contexts
opts$annos <- c(
  "window500_step250"
)

# Window for the overlap
opts$window.overlap <- 15000

# opts$colors <- c(
#   "H3K27ac_distal_E7.5_Ect_intersect12_500" = "steelblue",
#   "H3K27ac_distal_E7.5_End_intersect12_500" = "#43CD80",
#   "H3K27ac_distal_E7.5_Mes_intersect12_500" = "violetred",
#   "prom_2000_2000" = "grey60"
# )


# Minimum differential accessibility (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$input.dir,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[complete.cases(.)]

diff.results %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)]

# Load TSS
tss <- fread(io$tss, header=F,
  select=c("V1"="character", "V2"="integer", "V3"="integer","V4"="character")
) %>% setnames(c("chr","start","end","gene")) %>%
  .[,chr:=as.factor(gsub("chr","",chr))] %>%
  setkey(chr,start,end)

################
## Parse data ##
################

# Extend TSS
tss.extended <- tss %>% copy %>%
  .[,c("start","end"):=list(start-opts$window.overlap,end+opts$window.overlap)] %>% 
  .[start<0,start:=0] %>%
  setkey(chr,start,end)

#############
## Overlap ##
#############

# TO-DO: ONLY OVERLAP IF ITS NOT A GENE-CENTRIC ANNOTATION

dt.ov <- diff.results[,c("chr","start","end","id","comparison","sig")] %>%
  setkey(chr,start,end) %>%
  foverlaps(tss.extended, mult="all", nomatch = 0) %>%
  setnames(c("start","end"), c("promoter_start","promoter_end")) %>%
  setnames(c("i.start","i.end"), c("feature_start","feature_end")) %>%
  # FIX TO MINIMUM DISTANCE
  .[,dist:=round(abs(mean(abs(promoter_start-promoter_end))-mean(abs(feature_start-feature_end))))] %>%
  .[,c("chr","gene","promoter_start","promoter_end","id","feature_start","feature_end","dist","sig","comparison")]

##########
## Save ##
##########

for (i in opts$comparisons) {
  tmp <- dt.ov[comparison==i] %>% .[,comparison:=NULL]
  fwrite(tmp, sprintf("%s/%s.txt.gz",io$outdir,i), na = "NA", sep="\t", quote=F)
}

