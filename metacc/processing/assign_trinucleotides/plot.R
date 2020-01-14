
#############################################################################
## Script to plot the results from the differential methylation analysis ##
#############################################################################

library(data.table)
library(purrr)
library(ggplot2)

#####################
## Define settings ##
#####################

source("/Users/ricard/Guo_2017/settings.R")

io$input.file <- "/Users/ricard/data/Guo_2017/met/trinucleotides/trinucleotide_stats.txt.gz"
io$outdir <- "/Users/ricard/data/Guo_2017/met/trinucleotides"

# opts$stages <- c(
#   "Zygote",
#   "2-cell",
#   "4-cell",
#   "8-cell",
#   "Morula",
#   "ICM",
#   "TE"
# )

###############
## Load data ##
###############

data <- fread(io$input.file) %>%
  setnames("sample","id_met") %>%
  merge(sample_metadata, by="id_met")

##########
## Plot ##
##########

to.plot <- data[,.(N=sum(N)),by=c("id_met","stage")]
# p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim)

# pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=6, height=4)
png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
print(p)
dev.off()
