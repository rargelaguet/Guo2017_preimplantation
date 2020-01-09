library(data.table)
library(purrr)
library(furrr)

######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_reads.x), met_reads.x := 0] %>%
    .[is.na(met_reads.y), met_reads.y := 0] %>%
    .[is.na(nonmet_reads.x), nonmet_reads.x := 0] %>%
    .[is.na(nonmet_reads.y), nonmet_reads.y := 0] %>%
    .[,.(chr=chr, pos=pos, met_reads=met_reads.x+met_reads.y, nonmet_reads=nonmet_reads.x+nonmet_reads.y)]
}


fread_and_merge <- function(dt, file){
  setkey(dt, chr, pos)
  fread(file, header = F, select=1:3, colClasses=list(factor=1L)) %>% 
    setnames(c("chr","pos","rate")) %>%
    .[rate==0.50,c("met_reads","nonmet_reads"):=.(1L,1L)] %>%
    .[rate>0.50,c("met_reads","nonmet_reads"):=.(1L,0L)] %>%
    .[rate<0.50,c("met_reads","nonmet_reads"):=.(0L,1L)] %>%
    .[,rate:=NULL] %>%
    merge_and_sum(dt)
}



################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/acc/gpc_level/pseudobulk")
dir.create(io$outdir, showWarnings = F)

####################
## Define options ##
####################

# Define groups to pseudobulk (stages)
opts$groups <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

sample_metadata <- sample_metadata %>%
  .[,group:=stage] %>%
  .[group%in%opts$groups]
opts$cells <- sample_metadata$id_met

table(sample_metadata$group)

# Parallel processing
opts$parallel <- TRUE    # do parallel processing?
opts$ncores <- 2         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????

##############################
## Load data and pseudobulk ##
##############################

# Parallel processing options
if (opts$parallel){
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}

for (i in opts$groups) {
  
  # Define input files 
  cells <- sample_metadata[group%in%i,id_met]
  # cells <- head(cells,n=20)
  
  files <- paste0(io$met_data_raw, "/", cells, ".tsv.gz")
  
  # split into chunks for parallel processing
  if (opts$parallel) {
    chunks <- ceiling(seq_along(files)/opts$chunk_size)
    file_list <- split(files, chunks)
  } else {
    file_list <- list(files)
  }
  
  # pseudobulk
  init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_reads=as.integer(NA), nonmet_reads=as.integer(NA))
  data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
  # data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
    purrr::reduce(merge_and_sum) %>%
    .[,rate:=round(100*met_reads/(met_reads+nonmet_reads))] %>%
    .[complete.cases(.)] %>%
    .[,.(chr,pos,met_reads,nonmet_reads,rate)]
  
  # Save
  outfile = sprintf("%s/%s.tsv.gz",io$outdir,i)
  fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
}
