suppressPackageStartupMessages(library(SingleCellExperiment))

# Define functions
mean_sd <- function(x) data.frame(y=mean(x), ymin=mean(x)-sd(x)/2, ymax=mean(x)+sd(x)/2) 

######################
## Define settings  ##
######################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
  io$tss <- "/Users/ricard/data/hg38_regulation/promoters/TSS_mRNA.bed"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/Guo2017_preimplantation/settings.R")
  io$tss <- "/hps/nobackup2/research/stegle/users/ricard/hg38_regulation/promoters/TSS_mRNA.bed"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/metacc/coupling/local_coupling")

# Define window options
opts$up <- 3000
opts$down <- 3000
opts$window <- 150
opts$slide <- 50

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  # .[pass_metQC==TRUE & pass_accQC==TRUE]
  .[!is.na(id_met) & !is.na(id_acc)] %>% head(n=3)

# [TESTING MODE] Subset cells to reduce memory burden  ##
# sample_metadata <- sample_metadata[,head(.SD,n=3),by="lineage"]

#########################
## Load gene metadata  ##
#########################

tss <- fread(io$tss, colClasses=c("factor","integer","integer","factor","factor","factor")) %>%
  setnames(c("chr","start","end","id","score","strand")) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,score:=NULL] %>%
  .[,tss:=start] %>%
  .[,c("start","end"):=.(start-opts$up,end+opts$down)] %>%
  setkey(chr,start,end)

# Split between CGI or non-CGI promoters
# non_cgi_promoters <- fread("/Users/ricard/data/gastrulation/features/filt/prom_2000_0_noncgi.bed",stringsAsFactors=TRUE)[,V5]
# cgi_promoters <- fread("/Users/ricard/data/gastrulation/features/filt/prom_2000_0_cgi.bed",stringsAsFactors=TRUE)[,V5]
# tss <- tss[,type:=ifelse(id%in%non_cgi_promoters,'non_cgi','cgi')]
# tss <- tss[type=="non_cgi"]

##########
## Run  ##
##########

# Define genomic window around TSS
tmp <- seq(from=0-opts$up, to=0+opts$down-opts$window, by=opts$slide)
foo <- data.table(window_center=tmp+(opts$window/2), rel_start=tmp, rel_end=tmp+opts$window) %>%
  .[,(c("window_center","rel_start","rel_end")):=lapply(.SD, as.integer),.SDcols=(c("window_center","rel_start","rel_end"))] %>%
  setkey(rel_start,rel_end)

dt <- list()
for (i in sample_metadata$sample) {
  print(i)
  
  # Define chromatin accessibility data.table
  id.acc <- sample_metadata[sample==i,id_acc]
  acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,id.acc), sep=" ", showProgress = F, header = T) %>%
    setnames(c("chr","pos","rate")) %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>%  
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(opts$slide*round(dist/opts$slide))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%
    
    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,.(acc_rate=mean(rate)), by=c("id","window_center")] %>%
    .[,sample:=i]
  
  # Define DNA methylation data.table
  id.met <- sample_metadata[sample==i,id_met]
  met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,id.met), sep=" ", showProgress = F, header = T) %>%
    setnames(c("chr","pos","rate")) %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>%  
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(opts$slide*round(dist/opts$slide))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%

    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,.(met_rate=mean(rate)), by=c("id","window_center")] %>%
    .[,sample:=i]
    
    
    # Calculate Met vs Acc correlation coefficient per window
    dt[[i]] <- merge(met_dt, acc_dt, by = c("window_center","id","sample")) %>%
      .[,.(r = cor(x=met_rate, y=acc_rate, method="pearson")), by = c("sample","window_center")]
}  
  
dt <- rbindlist(dt)

##########
## Save ##
##########

fwrite(dt, paste0(io$outdir,"/data.txt.gz"))
