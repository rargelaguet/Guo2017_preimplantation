########################################################################
## Script to assign strand information to CpG sites from bismark file ##
########################################################################

# Input: 
# single-cell methylation files output from Bismark. In either one of the two following formats:
# input_format=1:
  # chr     pos     rate
  # 1       3019021 0
  # 1       3027398 100
  # 1       3052955 100

# input_format=2:
  # chr     pos     met_reads non nomet_reads
  # 1       3019021 0 1
  # 1       3027398 1 1
  # 1       3052955 1 0

# input_format=3:
# chr     pos     met_reads non nomet_reads rate
# 1       3019021 0 1 0
# 1       3027398 1 2 0 
# 1       3052955 1 0 1

# Output:
# output_format=1: we include the strand information in a new column but we don't modify the coordinates
# output_format=2: we modify the coordinates to map negative CpGs to the positive strand 

# Load libraries
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
# suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores', default=1)
  
opts <- p$parse_args(commandArgs(TRUE))
# opts$context <- toupper(opts$context); stopifnot(opts$context %in% c("CG","GC"))
opts$context <- "CG"

# Define options
opts$input_format <- 1
opts$output_format <- 1

# Define I/0
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  if (opts$context=="CG") {
    io$indir <- "/Users/ricard/data/Guo_2017/met/cpg_level"
    io$outdir <- "/Users/ricard/data/Guo_2017/met/strand"
  } else {
    stop()
  }
} else {
  if (opts$context=="CG") {
    # io$indir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/met/liftover"
    # io$outdir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/met/liftover/strand"
  } else {
    # io$indir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/acc/liftover"
    # io$outdir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/acc/liftover/strand"
  }
}
dir.create(io$outdir)

# List of valid chromosomes
# opts$chr <- c(1:19,"X","Y","M")
opts$chr <- c(1,"X","M")

# Load samples
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))

# Create stats data.frame
stats <- data.table(expand.grid(samples,opts$chr)) %>% setnames(c("sample","chr")) %>%
  .[,c("rate_forward","rate_reverse","coverage_forward","coverage_reverse"):=as.numeric(NA)]

for (i in 1:length(samples)) {
  cat(sprintf("Processing %s...\n",samples[i]))
  
  # Load data
  data <- fread(sprintf("%s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F) 
  
  # Input format 1 (chr,pos,rate)
  if (opts$input_format == 1) {
    colnames(data) <- c("chr","pos","rate")
    # Input format 2 (chr,pos,met_reads,nonmet_reads)
  } else if (opts$input_format == 2) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
    # Input format 2 (chr,pos,met_reads,nonmet_reads,rate)
  }  else if (opts$input_format == 3) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads","rate")
  }
  
  data <- data %>% .[chr%in%opts$chr]
  
  # Add one to the coordinates do to 0-based index
  data[,pos:=pos+1]
  
  # Get sequence
  seq <- getSeq(Mmusculus, names=paste0("chr",data$chr), start=data$pos-1, end=data$pos+1)
  data[,c("base_up","base","base_down") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]
  

  # Do sanity checks
  # if (opts$context=="CG") {
  #   data[base=="C",table(base_up)] # for CG methylation, this should contain only A and Ts
  #   data[base=="C",table(base_down)] # for CG methylation, this should contain only Gs
  #   data[base=="G",table(base_up)] # for methylation, this should only be Cs
  #   data[base=="G",table(base_down)] # for methylation, this should only be As and Ts
  # } else {
  #   stop()
  # }
  
  # Remove sites which are not CG in the reference
  data <- data[base%in%c("C","G")]
  stopifnot(unique(data$base) %in% c("G","C"))
  
  # Add strand information and remove base information
  if (opts$context=="CG") {
    data[,strand:=ifelse(base=="C","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
  } else {
    data[,strand:=ifelse(base=="G","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
  }
  
  # "Positivise" the dinucleotides and remove strand information
  # if (opts$output_format == 2) {
  #   data[,pos:=ifelse(strand=="+",pos,pos-1)] %>% .[,strand:=NULL]
  #   # Sometimes we collect the same CpG site in both strands, in that case we can merge the information
  #   if (opts$input_format == 1) {
  #     data <- data[,.(rate=as.integer(round(mean(rate)))), keyby=c("chr","pos")]
  #   } else if (opts$input_format == 2) {
  #     data <- data[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), keyby=c("chr","pos")]
  #   } else if (opts$input_format == 3) {
  #     data <- data[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), keyby=c("chr","pos")] %>%
  #       .[,rate:=round(met_reads/(met_reads+nonmet_reads))]
  #   }
  # }
  for (j in unique(stats$chr)) {
    stats[chr==j & sample==samples[i],rate_forward:=data[chr==j & strand=="+",round(mean(rate),4)]]
    stats[chr==j & sample==samples[i],coverage_forward:=data[chr==j & strand=="+",.N]]
    stats[chr==j & sample==samples[i],rate_reverse:=data[chr==j & strand=="-",round(mean(rate),4)]]
    stats[chr==j & sample==samples[i],coverage_reverse:=data[chr==j & strand=="-",.N]]
  }
  
  # Save results
  # fwrite(data, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=FALSE)
}

fwrite(stats, sprintf("%s/strand_stats.txt",io$outdir), quote=F, col.names=F, row.names=F, na="NA")
