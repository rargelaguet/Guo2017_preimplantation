
#############################################################################
## Script to plot the results from the differential methylation analysis ##
#############################################################################

library(data.table)
library(purrr)
library(ggplot2)

#####################
## Define settings ##
#####################

source("/Users/ricard/Guo2017_preimplantation/settings.R")

io$input.file <- "/Users/ricard/data/Guo2017_preimplantation/met/trinucleotides/trinucleotide_stats.txt.gz"
io$outdir <- "/Users/ricard/data/Guo2017_preimplantation/met/trinucleotides"

###############
## Load data ##
###############

data <- fread(io$input.file) %>%
  setnames("sample","id_met") %>%
  merge(sample_metadata, by="id_met")

##########
## Plot ##
##########

to.plot <- data %>%
  .[,.(N=sum(N)),by=c("id_met","stage","trinucleotide")] %>%
  .[,N_total:=sum(N),by="id_met"] %>%
  .[,value:=N/N_total]

ggscatter(to.plot, x="trinucleotide", y="value", size=0.5, position = position_jitter(0.1)) +
  labs(x="Trinucleotide", y="Fraction of CpG sites")

# pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=6, height=4)
png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
print(p)
dev.off()
