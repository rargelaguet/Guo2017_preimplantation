suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
  source("/Users/ricard/Guo_2017/metacc/plot_examples/load_data.R")
} else {
  stop("Computer not recognised")
}

## I/O ##
io$outdir <- paste0(io$basedir,"/metacc/plot_examples")
dir.create(io$outdir, showWarnings = F)


## Options ##

# Define genomic feature
id <- c("ENSMUSG00000025907","ENSMUSG00000049775")
anno <- "prom_2000_2000"


# Define stages
opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

# Define colors for the omics
opts$color <- c(
  "Chromatin accessibility" = "#00BFC4",
  "DNA methylation" = "#F37A71"
)

# Define minimum coverage
opts$min.cpg <- 1
opts$min.gpc <- 3

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] 
opts$met_cells <- sample_metadata %>% .[,id_met]
opts$acc_cells <- sample_metadata %>% .[,id_acc]

###############
## Load data ##
###############

# Load the three omics
dt <- load_data(id, anno, opts$met_cells, opts$acc_cells, opts$min.cpg, opts$min.gpc)

# Merge data with sample metadata
dt$acc <- merge(dt$acc, sample_metadata[,c("sample","id_acc","stage")], by="id_acc")
dt$met <- merge(dt$met, sample_metadata[,c("sample","id_met","stage")], by="id_met")

# bind in a single data table
dt_all <- do.call("rbind",list(
  dt$met[,c("sample","id","stage","rate")] %>% .[,assay:="DNA methylation"],
  dt$acc[,c("sample","id","stage","rate")] %>% .[,assay:="Chromatin accessibility"]
))

dt_all[,stage:=gsub("_"," ",stage)]

dt_all[,assay:=factor(assay,levels=c("DNA methylation","Chromatin accessibility"))]
  
#######################
## Generate Boxplots ##
#######################

for (i in unique(dt_all$id)) {
  
  p <- ggplot(dt_all[id==i], aes(x=stage, y=rate)) +
    facet_wrap(~assay, ncol=1, scales="free_y") +
    geom_jitter(aes(color=assay), size=0.5) +
    geom_violin(aes(fill=assay), alpha=0.5, size=0.25) +
    geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
    scale_fill_manual(values=opts$color) +
    scale_color_manual(values=opts$color) +
    coord_cartesian(ylim=c(0,100)) +
    labs(x="", y="Levels (%)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(color="black", hjust = 0.5),
      axis.title.y = element_text(colour="black", size=rel(1.1), vjust=1.5),
      axis.text.x = element_text(size=rel(1.2), color="black"),
      axis.text.y = element_text(colour="black",size=rel(1.3)),
      axis.line = element_line(colour="black", size=rel(0.7)),
      axis.ticks.y = element_line(colour="black"),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      legend.text = element_text(size=15),
      legend.title = element_blank()
    )

  # pdf(sprintf("%s/boxplot_%s_%s.pdf",io$outdir,id,anno), useDingbats=F, width=6, height=5)
  print(p)
  # dev.off()
}
