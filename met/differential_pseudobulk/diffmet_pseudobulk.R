#####################################################################
## Script to compute differential methylation in pseudobulked data ##
#####################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(argparse)

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-a',  '--anno',           type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('-s1', '--stage1', type="character",  nargs='+',  help='stage 1 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-s2', '--stage2', type="character",  nargs='+',  help='stage 2 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-o',  '--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## Define I/O ##
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
  source("/Users/ricard/Guo_2017/met/differential_pseudobulk/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/Guo_2017/settings.R")
  source("/homes/ricard/Guo_2017/met/differential_pseudobulk/utils.R")
} else {
  stop("Computer not recognised")
}
io$outfile <- args$outfile

io$met_data_parsed <- paste0(io$met_data_parsed,"/pseudobulk")
io$acc_data_parsed <- paste0(io$acc_data_parsed,"/pseudobulk")
io$metadata <- paste0(io$basedir,"/sample_metadata_pseudobulk.txt") 

## Define options ##
opts <- list()

# Define genomic contexts
opts$annos <- args$anno

# Define stage and lineage
opts$groupA <- args$stage1
opts$groupB <- args$stage2

# Filter by coverage
opts$min.observations <- 5  # Minimum number of observations (reads) per feature

## START TESTING ##
# opts$anno <- "H3K27ac_distal_E7.5_End_intersect12"
# opts$groupA <- "Zygote"; opts$groupB <- c("ICM","TE")
## END TESTING ##

# Load sample metadata
sample_metadata <- fread(io$metadata) %>% 
  .[stage%in%c(opts$groupA,opts$groupB)] %>%
  .[,c("id_met","stage")]


###############
## Load data ##
###############

# Load methylation data
# data <- fread(sprintf("%s/%s.tsv.gz",io$data.dir,opts$anno)) %>% 
#   setnames(c("id_met","id","anno","Nmet","Ntotal","rate"))

data <- fread(
  file = sprintf("%s/%s.tsv.gz",io$met_data_parsed,opts$anno)
) %>% .[V1%in%sample_metadata$id_met] %>% setnames(c("id_met","id","anno","Nmet","Ntotal","rate"))


############################
## Parse methylation data ##
############################

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata[,c("id_met","stage")], by="id_met")

# Define the two exclusive groups
data[,group:=as.factor( c("A","B")[as.numeric(stage%in%opts$groupB)+1] )]
sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage%in%opts$groupB)+1] )]

#############################
## Filter methylation data ##
#############################

# Filter features by coverage
data <- data[Ntotal>=opts$min.observations]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by=c("id","anno")] %>% .[Ngroup==2] %>% .[,Ngroup:=NULL]

#######################################
## Differential methylation analysis ##
#######################################

diff <- data %>%
  dcast(id+anno~group, value.var=c("rate","Ntotal"), fun.aggregate=mean) %>%
  .[,diff:=rate_B-rate_A] %>%
  .[,c("rate_A","rate_B","diff"):=list(round(rate_A,2),round(rate_B,2),round(diff,2))]

#############################
## Add genomic coordinates ##
#############################

feature_metadata <- fread(
  sprintf("%s/%s.bed.gz",io$features.dir,opts$anno), header=FALSE, 
  select=c("V1"="factor", "V2"="integer", "V3"="integer","V5"="character")
) %>% setnames(c("chr","start","end","id"))

# stopifnot(all(unique(diff$id) %in% unique(feature_metadata$id)))
diff <- merge(feature_metadata,diff, by="id")

# Sort results by p-value
diff %>% setorderv("diff")

##################
## Save results ##
##################

fwrite(diff, file=io$outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="NA")
