suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(doParallel))

#####################
## Define settings ##
#####################

## I/O ##

io <- list()
io$sample.metadata <- "/Users/ricard/data/Guo_2017/sample_metadata.txt"
# io$indir <- "/Users/ricard/data/Guo_2017/met/cpg_level"
# io$outdir <- "/Users/ricard/data/Guo_2017/met/cpg_level/seqmonk"
io$indir <- "/Users/ricard/data/Guo_2017/acc/gpc_level"
io$outdir <- "/Users/ricard/data/Guo_2017/acc/gpc_level/seqmonk"; dir.create(io$outdir)

## options ##
opts <- list()

# Number of cores
opts$cores <- 2

# Stages to use
opts$stage <- c(
  # "ESC",
  # "TE",
  # "ICM",
  # "16-cell"
  # "8-cell",
  "4-cell",
  "2-cell",
  "Zygote"
)

###############
## Load data ##
###############

sample_metadata <- fread(io$sample.metadata) %>%
    .[stage%in%opts$stage]

# samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))
# samples <- sample_metadata$id_met
samples <- sample_metadata$id_acc

################
## Pseudobulk ##
################

registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples)) %dopar% {
# for (i in 1:length(samples)) {
  outfile <- sprintf("%s/%s.cov.gz",io$outdir,samples[i])
  if (file.exists(outfile)) {
    cat(sprintf("File %s already exists, skipping...\n",outfile))
  } else {
    
    # Load data
    cat(sprintf("Processing %s...\n",samples[i]))
    data <- fread(sprintf("%s/%s.tsv.gz",io$indir,samples[i])) %>% setnames(c("chr","pos","rate"))
    
    # Add columns for met_reads and nonmet_reads
    data %>% 
    .[rate==0.50,c("met_reads","nonmet_reads"):=.(1L,1L)] %>%
    .[rate>0.50,c("met_reads","nonmet_reads"):=.(1L,0L)] %>%
    .[rate<0.50,c("met_reads","nonmet_reads"):=.(0L,1L)]

    # Convert rate from 0-1 to 0-100
    data[,rate:=rate*100]
    
    # Add columns 'start' and 'end'
    data[,c("start","end"):=pos] %>% .[,pos:=NULL]
    
    # Reorder columns
    data <- data %>% setcolorder(c("chr","start","end","rate","met_reads","nonmet_reads"))
    
    # Save results
    fwrite(data, outfile, sep="\t", col.names=FALSE)
  }
})
