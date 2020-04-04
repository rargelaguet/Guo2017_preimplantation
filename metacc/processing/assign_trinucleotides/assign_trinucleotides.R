#####################################
## Script to assign trinucleotides ##
#####################################

suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('-c', '--context', type="character",help='cg/CG or gc/GC')
  
args <- p$parse_args(commandArgs(TRUE))

args$context <- toupper(args$context)
stopifnot(args$context %in% c("CG","GC"))

################
## Define I/0 ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
} else if (grepl("yoda",Sys.info()['nodename'])) {
  source("/homes/ricard/Guo2017_preimplantation/settings.R")
} else {
  stop("Computer not recognised")
}

if (args$context=="CG") {
  io$indir <- paste0(io$basedir,"/met/cpg_level")
  io$outdir <- paste0(io$basedir,"/met/trinucleotides")
} else if (args$context=="GC") {
  io$indir <- paste0(io$basedir,"/acc/gpc_level")
  io$outdir <- paste0(io$basedir,"/acc/trinucleotides")
}
dir.create(io$outdir, showWarnings = F)

####################
## Define options ##
####################

# List of chromosomes
opts$chr <- c(1:19,"X","Y","M")

if (args$context=="CG") {
  opts$trinucleotides <- c("ACG","TCG","CCG","GCG")
} else if (args$context=="GC") {
  opts$trinucleotides <- c("AGC","TGC","CGC","GGC")
}

########################################
## Load data and fetch trinucleotides ##
########################################

# Define samples
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))

# Create data.frame to store the output
stats <- data.table(expand.grid(samples,opts$chr)) %>% setnames(c("sample","chr"))
if (args$context=="CG") {
  stats %>% .[,c("ACG","TCG","CCG","GCG"):=as.numeric(NA)]
} else if (args$context=="GC") {
  stats %>% .[,c("AGC","TGC","CGC","GGC"):=as.numeric(NA)]
}
stats <- stats %>% melt(id.vars=c("sample","chr"), variable.name="trinucleotide", value.name="N")


# for (i in 1:5) {
for (i in 1:length(samples)) {
  cat(sprintf("\nProcessing %s...\n",samples[i]))
  
  # Load data
  data <- fread(sprintf("%s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F) 
  
  # Assign column names
  if (ncol(data)==3) {
    colnames(data) <- c("chr","pos","rate")
  } else if (ncol(data)==5) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
  } else if (ncol(data)==6) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads","rate")
  } else {
    stop("Format not recognised")
  }
  
  # Filter chromosomes
  data <- data %>% .[chr%in%opts$chr]
  
  # Add one to the coordinates do to 0-based index
  TO-DO: CHECK IF THIS IS NECESSARY.
  data[,pos:=pos+1]
    
  # Get sequence for the central base
  seq <- getSeq(Mmusculus, names=paste0("chr",data$chr), start=data$pos, end=data$pos)
  data[,base:=as.character(seq)]
  
  # Remove sites which are not C/G in the reference
  data <- data[base%in%c("C","G")]
  
  # Add strand information
  if (args$context=="CG") {
    data[,strand:=ifelse(base=="C","+","-")]
    data[,base:="C"]
  } else if (args$context=="GC") {
    data[,strand:=ifelse(base=="G","+","-")]
    data[,base:="G"]
  }
  
  # Get sequence of sorrounding nucleotides
  seq <- getSeq(Mmusculus, names=paste0("chr",data$chr), start=data$pos-1, end=data$pos+1, strand=data$strand)
  data[,c("base_up","base_down") := list(substr(as.character(seq),1,1), substr(as.character(seq),3,3))]
  
  # Assign trinucleotide
  data[,trinucleotide:=paste(base_up,base,base_down,sep="")]
  print(table(data$trinucleotide))
  
  # Store statistics
  for (j in unique(stats$chr)) {
    for (k in opts$trinucleotides) {
      stats[chr==j & sample==samples[i] & trinucleotide==k, N:=data[chr==j & trinucleotide==k,.N]]
    }
  }
}

fwrite(stats, sprintf("%s/trinucleotide_stats.txt.gz",io$outdir), quote=F, na="NA", sep="\t")
