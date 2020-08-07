#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/human_embryo_multiomics/settings.R")
} else {
  stop("Computer not recognised")
}

io$outdir <- paste0(io$basedir,"/met/stats")

# Define cells
# opts$cells <- metadata %>% .[!is.na(id_met),id_met]# %>% head(n=5)
opts$cells <- list.files(io$met_data_raw, pattern = "*.tsv.gz") %>% 
  stringr::str_replace_all(".tsv.gz","")

####################################################
## Load methylation data and calculate statistics ##
####################################################

stats <- data.table(id_met=opts$cells) %>% 
  .[,c("nreads","nCG","rate"):=as.numeric(NA)]

for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",io$met_data_raw,i))) {
    print(i)

    data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), sep="\t", verbose=F, showProgress=F)

    # Compute genome-wide statistics
    stats[id_met==i, c("nreads","nCG","rate"):=list(sum(data$met_reads+data$nonmet_reads), nrow(data), round(100*mean(data$rate),2))]

  } else {
    print(sprintf("Sample %s not found for methylation",i))
  }
}

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/sample_stats.txt")
fwrite(stats, outfile, sep="\t", na = "NA")
