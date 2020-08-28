##  Script to overlap bismark output files with genomic features ##
###################################################################

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",  help='cg/CG or gc/GC')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

# args <- list()
# args$context <- "CG"

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
args$context <- toupper(args$context)
stopifnot(args$context %in% c("CG","GC"))

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
} else if (grepl("yoda",Sys.info()['nodename'])) {
  source("/homes/ricard/Guo2017_preimplantation/settings.R")
} else {
  stop("Computer not recognised")
}

if (args$context == "GC") {
  io$in.folder <- io$acc_data_raw
  io$out.folder <- io$acc_data_parsed
} else if (args$context=="CG") {
  io$in.folder <- io$met_data_raw
  io$out.folder <- io$met_data_parsed
}

## Options ##

# Annotations to analyse
# opts$annos <- "all"
# opts$annos <- c("window500_step250")
opts$annos <- c(
  "CGI"
  # "genebody",
  # "prom_2000_2000",
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi",
  # "prom_200_200",
  # "prom_200_200_cgi",
  # "prom_200_200_noncgi",
  # "LINE",
  # "LTR",
  # "ESC_CTCF",
  # "ESC_DHS",
  # "ESC_p300",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Mes_intersect12"
  # "window2000_step250",
  # "window2000_step1000"
)


if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$features.dir, pattern = "\\.bed.gz$"),"\\.bed.gz"),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",args$context))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
if (args$context=="CG") {
  samples_keep <- fread(io$metadata, header=T) %>% .[!is.na(id_met),id_met]
} else if (args$context=="GC") {
  samples_keep <- fread(io$metadata, header=T) %>% .[!is.na(id_acc),id_acc]
} else{
  stop()
}
stopifnot(all(!duplicated(samples_keep)))


############################
## Preprocess annotations ##
############################

anno_list <- list()
for (i in 1:length(opts$annos)) {
  anno_list[[i]] <- fread(sprintf("%s/%s.bed.gz",io$features.dir,opts$annos[i]), header=F, select=c(1,2,3,4,5)) %>% 
    setnames(c("chr","start","end","strand","id")) %>%
    setkey(chr,start,end)
}
names(anno_list) <- opts$annos


#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",sample,opts$annos) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    filename <- sprintf("%s/%s.tsv.gz",io$in.folder,sample)
    print(filename)
    if (file.exists(filename)) {
      dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) 
      if (nrow(dat_sample)>1) {
        dat_sample <- dat_sample %>%
          setnames(c("chr","pos","rate")) %>%
          .[,rate:=rate*100] %>%
          .[,c("start","end") := list(pos,pos)] %>%
          .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
          setkey(chr,start,end)

        # stopifnot(all(dat_sample$rate %in% c(0,100)))
        
        # Overlap data with annotations
        for (anno in opts$annos) {
          # fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
          fname.out <- sprintf("%s/tmp/%s_%s.gz",io$out.folder,sample,anno)
          if (file.exists(fname.out)) {
            cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
          } else {
            cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
            
            # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
            ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
              .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
              .[,c("sample","anno") := list(sample,anno)] %>%
              .[,.(rate=round(mean(rate)), Nmet=sum(rate>50), N=.N), keyby=.(sample,id,anno)] %>%
              .[,c("sample","id","anno","Nmet","N","rate")]
            

            # Store and save results
            fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
          }
        }
        rm(dat_sample)
      }
    }
  }
}



# Concatenate everything and save it
for (i in opts$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    # files <- list.files(paste(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample' | gzip > %s", io$out.folder, i, outfile))
    cat("\n")
  }
}
