
#######################################
## Load modules and define functions ##
#######################################


library(data.table)
library(purrr)
library(scater)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

source("/Users/ricard/Guo2017_preimplantation/mofa/load_settings.R")

###########################
## Load Methylation data ##
###########################

met_dt <- lapply(opts$met.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met.dir,n)) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist
colnames(met_dt) <- c("id_met","id","anno","Nmet","Ntotal","rate")

#############################
## Load Accessibility data ##
#############################

acc_dt <- lapply(opts$acc.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc.dir,n)) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist
colnames(acc_dt) <- c("id_acc","id","anno","Nmet","Ntotal","rate")

##############################
## Merge data with metadata ##
##############################

met_dt <- merge(met_dt, sample_metadata, by="id_met")
acc_dt <- merge(acc_dt, sample_metadata, by="id_acc")

################
## Parse data ##
################

# Parse accessibility data
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


# Parse methylation data
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[Ntotal>=opts$met_min.CpGs]

# Filter features by  minimum number of cells
met_dt <- met_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by  minimum number of cells (per stage)
# for (i in unique(met_dt$stage)) {
#   met_dt[stage==i,Ntotal:=sample_metadata[id_met%in%opts$met_cells & stage==i,.N]]
# }
# keep_cov_sites <- met_dt %>% split(.$stage) %>% 
#   map(~ .[, ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$met_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# met_dt <- met_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n = opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

# Filter features by variance (per stage)
# met_dt %>% .[,id_anno:=paste(id,anno,sep="_")]
# keep_hv_genes <- met_dt %>% split(.$stage) %>%
#   map(~ split(.,.$anno) %>%
#     map(~ .[,.(var=var(rate)), by=c("id_anno")] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id_anno) %>% unlist
#   )
# met_dt <- met_dt %>% split(.$stage) %>% map2(.,names(.), function(x,y) x[id_anno%in%keep_hv_genes[[y]]]) %>% rbindlist

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[Ntotal>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt <- acc_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by  minimum number of cells (per stage)
# for (i in unique(acc_dt$stage)) {
#   acc_dt[stage==i,Ntotal:=sample_metadata[id_acc%in%opts$acc_cells & stage==i,.N]]
# }
# keep_cov_sites <- acc_dt %>% split(.$stage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$acc_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# acc_dt <- acc_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n = opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

# Filter features by variance (per stage)
# acc_dt %>% .[,id_anno:=paste(id,anno,sep="_")]
# keep_hv_genes <- acc_dt %>% split(.$stage) %>%
#   map(~ split(.,.$anno) %>%
#     map(~ .[,.(var=var(rate)), by=c("id_anno")] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id_anno) %>% unlist
#   )
# acc_dt <- acc_dt %>% split(.$stage) %>% map2(.,names(.), function(x,y) x[id_anno%in%keep_hv_genes[[y]]]) %>% rbindlist
# acc_dt <- acc_dt %>% droplevels()


########################
## Join the data sets ##
########################

# Don't group
data1 <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","view")) %>% 
  .[,c("feature","view","group"):=list(paste0("met_",feature), paste0("met_",view), "group")]
data2 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","view")) %>% 
  .[,c("feature","view","group"):=list(paste0("acc_",feature), paste0("acc_",view), "group")]

# Group by stage
# data1 <- met_dt %>% .[,c("sample","stage","id","m","anno")] %>%  
#   setnames(c("sample","group","feature","value","view")) %>% .[,c("feature","view"):=list(paste0("met_",feature), paste0("met_",view))]
# data2 <- acc_dt %>% .[,c("sample","stage","id","m","anno")] %>%  
#   setnames(c("sample","group","feature","value","view")) %>% .[,c("feature","view"):=list(paste0("acc_",feature), paste0("acc_",view))]

# Concatenate
data <- rbind(data1,data2)

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/data.txt")
fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
system(sprintf("pigz -f %s",outfile))
