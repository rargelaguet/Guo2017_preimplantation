suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(doParallel))

#####################
## Define settings ##
#####################

## I/O ##

io <- list()
io$indir <- "/Users/ricard/data/Guo2017_preimplantation/acc/gpc_level/pseudobulk" 
io$outdir <- "/Users/ricard/data/Guo2017_preimplantation/acc/gpc_level/pseudobulk/seqmonk"
dir.create(io$outdir, showWarnings=F)

## options ##
opts <- list()
opts$cores <- 1

################
## Define I/0 ##
################

# Load samples
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))

registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples)) %dopar% {
# for (i in 1:length(samples)) {
  outfile <- sprintf("%s/%s.cov.gz",io$outdir,samples[i])
  if (file.exists(outfile)) {
    cat(sprintf("File %s already exists, skipping...\n",outfile))
  } else {
    
    # Load data
    cat(sprintf("Processing %s...\n",samples[i]))
    data <- fread(sprintf("%s/%s.tsv.gz",io$indir,samples[i])) %>%
      
      # Add columns 'start' and 'end'
      .[,c("start","end"):=pos] %>% .[,pos:=NULL] %>%
      
      # Reorder columns
      setcolorder(c("chr","start","end","rate","met_reads","nonmet_reads"))
      
    # Save results
    fwrite(data, outfile, sep="\t", col.names=FALSE)
  }
})
