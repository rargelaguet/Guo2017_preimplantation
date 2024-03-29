
#############################################################################
## Script to plot the results from the differential methylation analysis ##
#############################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

#####################
## Define settings ##
#####################

## I/O ##
io <- list()
io$input.dir <- "/Users/ricard/data/gastrulation/met/differential/feature_level/test"
io$outdir <- "/Users/ricard/data/gastrulation/met/differential/feature_level/pdf"

## Options ##
opts <- list()
opts$comparisons <- c(
  # "E6.5Primitive_Streak_vs_E7.5Nascent_mesoderm",
  # "E6.5Primitive_Streak_vs_E7.5Primitive_Streak",
  # "E6.5E7.5Nascent_mesoderm_vs_E7.5Mature_mesoderm",
  # "E6.5E7.5Primitive_Streak_vs_E6.5E7.5Nascent_mesoderm",
  # "E6.5E7.5Primitive_Streak_vs_E7.5Mature_mesoderm"
  
  "E6.5Epiblast_vs_E6.5Primitive_Streak",
  "E7.5Ectoderm_vs_E7.5MesodermEndoderm",
  "E7.5EpiblastEctoderm_vs_E7.5MesodermEndoderm",
  "E7.5Epiblast_vs_E7.5MesodermEndoderm"
)

# Select genomic contexts
opts$annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "prom_2000_2000"
)

opts$colors <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12_500" = "steelblue",
  "H3K27ac_distal_E7.5_End_intersect12_500" = "#43CD80",
  "H3K27ac_distal_E7.5_Mes_intersect12_500" = "violetred",
  "prom_2000_2000" = "grey60"
)


###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(cmd=sprintf("zcat < %s/%s_%s.txt.gz",io$input.dir,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[complete.cases(.)]

# diff.results[,c("prop1","prop2","diff"):=list(rateA*100,rateB*100,diff*100)]

##############
## Barplots ##
##############

gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value, group=anno)) +
    geom_bar(aes(fill=anno), color="black", stat="identity", position="dodge") +
    scale_fill_manual(values=opts$colors) +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title="", x="", y="Number of hits") +
    theme(
      plot.title = element_text(size=11, face='bold', hjust=0.5),
      axis.text = element_text(size=rel(1.0), color='black'),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position="none",
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  return(p)
}


tmp <- diff.results %>%
  .[,.(number_positive_hits=sum(sig==T & diff>0), number_negative_hits=sum(sig==T & diff<0)), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(tmp$value), max(tmp$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim)
  
  pdf(sprintf("%s/%s.pdf",io$outdir,i), width=6, height=4)
  print(p)
  dev.off()
}
