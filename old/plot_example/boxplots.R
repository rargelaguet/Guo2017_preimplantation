library(ggplot2)
library(RColorBrewer)
source("/Users/ricard/gastrulation/metaccrna/dynamics/load_data.R")

scale <- function(value, min.data, max.data, min.scaled, max.scaled) {
  stopifnot(is.numeric(value))
  stopifnot(value<=max.data & value>=min.data)
  return ((max.scaled - min.scaled) * (value - min.data) / (max.data - min.data)) + min.scaled
}

################
## Define I/O ##
################

io <- list()
io$outdir <- "/Users/ricard/gastrulation/metaccrna/dynamics/out"
io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata_scNMT.txt"

####################
## Define options ##
####################

opts <- list()

# Define gene identity
rna.id <- "ENSMUSG00000023043"
met.id <- "H3K27ac_distal_E7.5_End_intersect12_500_766"
met.anno <- "H3K27ac_distal_E7.5_End_intersect12_500"
acc.id <- "H3K27ac_distal_E7.5_End_intersect12_766"
acc.anno <- "H3K27ac_distal_E7.5_End_intersect12"

# Define stages to plot
opts$stage_lineages <- c(
  "E4.5_EPI",
  "E5.5_EPI",
  "E6.5_EPI","E6.5_PS",
  "E7.5_Ectoderm","E7.5_Mesoderm","E7.5_Endoderm"
)

# Define colors for the omics
# opts$color <- c("rna"="forestgreen", "acc"="royalblue4", "met"="red4")
opts$color <- c("RNA expression"="#3CB54E", "Chromatin accessibility"="#00BFC4", "DNA methylation"="#F37A71")

# Define minimum coverage
opts$min.cpg <- 2
opts$min.gpc <- 3

# Define cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[!is.na(id_met) & !is.na(id_acc)]
opts$met_cells <- tmp %>% .[pass_metQC==T & outlier==F & stage_lineage%in%opts$stage_lineage,id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T & outlier==F & stage_lineage%in%opts$stage_lineage,id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T & outlier==F & stage_lineage%in%opts$stage_lineage,id_acc]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample.metadata,stringsAsFactors=F) %>%
  .[,c("sample","id_met","id_rna","id_acc","stage","lineage")] %>%
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ] %>%
  droplevels()

# Load the three omics
dt <- load_data(rna.id, met.id, met.anno, acc.id, acc.anno, opts$min.cpg, opts$min.gpc)

# Merge data with sample metadata
dt$acc <- merge(dt$acc, sample_metadata[id_acc%in%opts$acc_cells,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
dt$met <- merge(dt$met, sample_metadata[id_met%in%opts$met_cells,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()
dt$rna <- merge(dt$rna, sample_metadata[id_rna%in%opts$rna_cells,c("sample","id_rna","stage","stage_lineage")], by="id_rna") %>% droplevels()

# bind in a single data table
dt_all <- do.call("rbind",list(
  dt$rna[,c("sample","gene","stage","stage_lineage","value")] %>% .[,assay:="RNA expression"] %>% setnames("gene","id"),
  dt$met[,c("sample","id","stage","stage_lineage","value")] %>% .[,assay:="DNA methylation"],
  dt$acc[,c("sample","id","stage","stage_lineage","value")] %>% .[,assay:="Chromatin accessibility"]
))

dt_all[,stage_lineage:=gsub("_"," ",stage_lineage)]

dt_all[,assay:=factor(assay,levels=c("Chromatin accessibility","RNA expression","DNA methylation"))]
  
#######################
## Generate Boxplots ##
#######################

p <- ggplot(dt_all, aes(x=stage_lineage, y=value)) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  geom_jitter(aes(color=assay), size=1.0) +
  geom_violin(aes(fill=assay), alpha=0.5) +
  scale_fill_manual(values=opts$color) +
  scale_color_manual(values=opts$color) +
  geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15) +
  labs(x="", y="", title="") +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(colour="black", size=rel(1.2), vjust=1.5),
    axis.title.x = element_text(colour="black", size=rel(1.2), vjust=1.5),
    axis.text.x = element_text(size=rel(1.4), angle=90, hjust=1, vjust=0.5, color="black"),
    axis.text.y = element_text(colour="black",size=rel(1.3)),
    axis.line = element_line(colour="black"),
    axis.ticks.y = element_line(colour="black"),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size=rel(1.3), color="black"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    legend.text=element_text(size=15),
    legend.title=element_blank(),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
print(p)

outfile <- sprintf("%s/boxplot_rna%s_met%s_acc%s.pdf",io$outdir,rna.id,met.id,acc.id)
pdf(outfile, useDingbats = F, width=6, height=6)
print(p)
dev.off()