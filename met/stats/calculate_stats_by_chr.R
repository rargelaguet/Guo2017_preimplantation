#####################
## Define settings ##
#####################

source("/homes/ricard/Guo2017_preimplantation/settings.R")

# Define I/O
io$outdir <- paste0(io$basedir,"/met/results/stats")

# Define chromosome
# opts$chr <- c(1:19,"X","Y")

# Update sample metadata
# sample_metadata <- sample_metadata %>% .[pass_metQC==TRUE]
sample_metadata <- sample_metadata %>% .[!is.na(id_met)]

# Test mode
# sample_metadata <- head(sample_metadata,n=3)

#####################
## Calculate stats ##
#####################

stats <- data.table(expand.grid(as.character(sample_metadata$id_met),opts$chr)) %>% 
  setnames(c("id_met","chr")) %>%
  .[,c("coverage","mean"):=as.numeric(NA)]

for (i in sample_metadata$id_met) {
  if (file.exists(sprintf("%s/%s.tsv.gz",io$met_data_raw,i))) {
    print(i)

    # Load sample methylation data
    data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i)) %>%
      setnames(c("chr","pos","rate"))
      # .[,chr:=paste0("chr",chr)]

    # Compute methylation statistics per chromosome
    for (j in opts$chr) {
      data_j <- data[chr==j]
      # stats[id_met==i & chr==j, c("nreads","coverage","mean"):=list(sum(data_j$met_reads+data_j$nonmet_reads), nrow(data_j),round(mean(data_j$rate)*100,2))]  
      stats[id_met==i & chr==j, c("coverage","mean"):=list(nrow(data_j),round(mean(data_j$rate)*100,2))]
    }

  } else {
    print(sprintf("Sample %s not found for methylation",i))
  }
}

fwrite(stats, paste0(io$outdir,"/stats_per_chromosome.txt.gz"), sep="\t")
