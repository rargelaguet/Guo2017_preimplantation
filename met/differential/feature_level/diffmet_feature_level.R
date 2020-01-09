#####################################################################
## Script to compute differential methylation at the feature level ##
#####################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-a',     '--anno',       type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('-g1',    '--group1',     type="character",  nargs='+',  help='group 1 (stage)')
p$add_argument('-g2',    '--group2',     type="character",  nargs='+',  help='group 2 (stage)')
p$add_argument('-cells', '--min.cells',  type="integer",                help='Minimum number of cells per group')
p$add_argument('-o',     '--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## Define I/O ##
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
  source("/Users/ricard/Guo_2017/met/differential/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/Guo_2017/settings.R")
  source("/homes/ricard/Guo_2017/met/differential/utils.R")
} else {
  stop("Computer not recognised")
}
io$outfile <- args$outfile

## Define options ##

# Define genomic contexts
stopifnot(length(args$anno)==1)
opts$anno <- args$anno

# Define groups
opts$groupA <- args$group1
opts$groupB <- args$group2

# Subset top most variable sites
opts$number_features <- NA

# Filter by coverage
opts$min.CpGs <- 1                # Minimum number of CpG per feature in each cell
opts$min.cells <- args$min.cells  # Minimum number of cells per feature in each group

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"

# Minimum differential methylation (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

## START TESTING ##
# opts$anno <- "prom_2000_2000"
# opts$groupA <- c("ICM")
# opts$groupB <- c("Morula")
# opts$min.cells <- 5
## END TESTING ##

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[,c("id_met","stage")] %>%
  .[stage%in%c(opts$groupA,opts$groupB)]

table(sample_metadata$stage)

###############
## Load data ##
###############

print("Loading data...")

# Load methylation data
data <- fread(
  file = sprintf("%s/%s.tsv.gz",io$met_data_parsed,opts$anno), 
  showProgress=F, header=F,
  select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer")
) %>% .[V1%in%sample_metadata$id_met] %>% droplevels %>% setnames(c("id_met","id","Nmet","N","rate"))

# Load methylation data (option 2, cell by cell)
# data <- sample_metadata$id_met %>% map(function(i) {
#   fread(
#     file = sprintf("%s/tmp/%s_%s.gz",io$met_data_parsed,i,opts$anno), 
#     showProgress = FALSE, header = FALSE,
#     select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer"))
# }) %>% rbindlist %>% setnames(c("id_met","id","Nmet","N","rate"))

print(object.size(data), units='auto')

############################
## Parse methylation data ##
############################

print("Parsing data...")

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata[,c("id_met","stage")], by="id_met")

# Define the two exclusive groups
data[,group:=as.factor( c("A","B")[as.numeric(stage%in%opts$groupB)+1] )]
sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage%in%opts$groupB)+1] )]
data[,stage:=NULL]

# Convert beta value to M value
if (opts$statistical.test == "t.test")
  data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

print(object.size(data), units='auto')

#############################
## Filter methylation data ##
#############################

print("Filtering data...")

# Filter features by coverage
data <- data[N>=opts$min.CpGs]

# Filter features by minimum number of cells per group
data <- data[,Ncells:=min(.N), by=c("id","group")] %>% .[Ncells>=opts$min.cells] %>% .[,Ncells:=NULL]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by=c("id")] %>% .[Ngroup==2] %>% .[,Ngroup:=NULL]

# Filter by variance
if (!is.na(opts$number_features)) {
  data[,var := var(rate), by="id"] %>% setorder(-var) %>% head(n=opts$number_features) %>% .[,var:=NULL]
}

#########################################
## Differential methylation analysis ##
#########################################

print("Doing differential testing...")

# Binomial assumption: test of equal proportions using Fisher exact test
if (opts$statistical.test == "binomial") {
  diff <- data %>% dcast(id~group, value.var=c("Nmet","N"), fun.aggregate=sum) %>%
    setnames(c("Nmet_A","Nmet_B"),c("A_met","B_met")) %>%
    .[,c("A_unmet","B_unmet"):=list(N_A-A_met,N_B-B_met)] %>% 
    .[,c("N_A","N_B"):=NULL] %>%
    .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id")] %>%
    .[,c("rateA","rateB"):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))]
  
# T-test under normality assumption
} else if (opts$statistical.test == "t.test") {
  warning("This data.table code is very slow, requires optimisation....")
  diff <- data[, .(
    N_A = .SD[group=="A",.N], N_B = .SD[group=="B",.N],
    rateA = mean(.SD[group=="A",rate]), rateB = mean(.SD[group=="B",rate]),
    p.value = t.test(x=.SD[group=="B",m], y=.SD[group=="A",m], var.equal=FALSE)[["p.value"]]), by = c("id")]
}

# Multiple testing correction and define significant hits
diff %>%
  .[,diff:=rateB-rateA] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr"))] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,(c("rateA","rateB","diff","log_padj_fdr")) := lapply(.SD,round,2), .SDcols = (c("rateA","rateB","diff","log_padj_fdr"))]

rm(data); gc(reset=TRUE)

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
diff %>% setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff, file=io$outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="NA")
