###########################
## Load Methylation data ##
###########################

met_dt1 <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir1,n), showProgress=F, stringsAsFactors=T, quote="") %>%
    .[V1%in%opts$met_cells1]
}) %>% rbindlist
colnames(met_dt1) <- c("id_met","id","anno","Nmet","N","rate")

met_dt2 <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir2,n), showProgress=F, stringsAsFactors=T, quote="") %>%
    .[V1%in%opts$met_cells2]
}) %>% rbindlist
colnames(met_dt2) <- c("id_met","id","anno","Nmet","N","rate")

#############################
## Load Accessibility data ##
#############################

acc_dt1 <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir1,n), showProgress=F, stringsAsFactors=T, quote="") %>%
    .[V1%in%opts$acc_cells1]
}) %>% rbindlist
colnames(acc_dt1) <- c("id_acc","id","anno","Nmet","N","rate")

acc_dt2 <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir2,n), showProgress=F, stringsAsFactors=T, quote="") %>%
    .[V1%in%opts$acc_cells2]
}) %>% rbindlist
colnames(acc_dt2) <- c("id_acc","id","anno","Nmet","N","rate")


###################
## Load RNA data ##
###################

sce <- readRDS(io$rna.file) %>% .[,opts$rna_cells1]

# Convert to data.table
rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))
rna_dt[,c("id_rna","gene","ens_id"):=list(as.factor(id_rna),as.factor(gene),as.factor(ens_id))]

################
## Merge data ##
################

met_dt <- rbind(met_dt1, met_dt2)
acc_dt <- rbind(acc_dt1, acc_dt2)

met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()
acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage","stage_lineage")], by="id_rna") %>% droplevels()

################
## Parse data ##
################

## Parse gene and feature metadata ##

feature_metadata_filt.met <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% met_dt[anno==y,id]] ) %>%
  rbindlist

feature_metadata_filt.acc <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% acc_dt[anno==y,id]] ) %>%
  rbindlist

gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","gene")] %>% 
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)

## Parse accessibility data ##
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


#####################################################################
## Associate the non-genic genomic contexts with overlapping genes ##
#####################################################################
  
# Methylation
if (opts$overlapGenes) {
  met_list <- list()
  for (i in unique(met_dt$anno)){
    
    # Subset corresponding anno
    met_tmp <- met_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(met_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.met[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      met_list[[i]] <- merge(met_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(met_tmp$id)) == TRUE)) {
      met_list[[i]] <- merge(met_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  met_dt <- rbindlist(met_list)
  rm(met_list, met_tmp,feature_metadata_tmp,ov)
} else {
  met_dt[,gene:="NA"]
}

# Accessibility
if (opts$overlapGenes) {
  acc_list <- list()
  for (i in unique(acc_dt$anno)){
    
    # Subset corresponding anno
    acc_tmp <- acc_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.acc[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      acc_list[[i]] <- merge(acc_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == TRUE)) {
      acc_list[[i]] <- merge(acc_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  acc_dt <- rbindlist(acc_list)
  rm(acc_list, acc_tmp,feature_metadata_tmp,ov)
} else {
  acc_dt[,gene:="NA"]
}

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Remove Y chromosome because they might have inherent differences due to gender
# met_dt <- merge(met_dt,feature_stats[,c("id","chr")], by="id", allow.cartesian=T) %>%  .[!chr %in% c("chrY")] %>% .[,chr:=NULL]
# met_dt <- met_dt[!id%in%feature_stats[chr=="Y",id]]

# Filter features by  minimum number of cells (by stage_lineage)
for (i in unique(met_dt$stage)) {
  met_dt[stage==i,Ntotal:=sample_metadata[id_met%in%opts$met_cells & stage==i,.N]]
}
keep_cov_sites <- met_dt %>% split(.$stage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$met_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
met_dt <- met_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
  .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
# keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>1] %>% setorder(-var) %>% head(n = opts$met_nfeatures) %>% .$id)
# met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

# Filter features by variance (per stage separately, otherwise you enrich for differences between stages)
met_dt %>% .[,id_anno:=paste(id,anno,sep="_")]
keep_hv_genes <- met_dt %>% split(.$stage) %>%
  map(~ split(.,.$anno) %>%
    map(~ .[,.(var=var(rate)), by=c("id_anno")] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id_anno) %>% unlist
  )
met_dt <- met_dt %>% split(.$stage) %>% map2(.,names(.), function(x,y) x[id_anno%in%keep_hv_genes[[y]]]) %>% rbindlist

met_dt <- met_dt %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells (by stage)
for (i in unique(acc_dt$stage)) {
  acc_dt[stage==i,Ntotal:=sample_metadata[id_acc%in%opts$acc_cells & stage==i,.N]]
}
keep_cov_sites <- acc_dt %>% split(.$stage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$acc_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
acc_dt <- acc_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
  .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
# keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>1] %>% setorder(-var) %>% head(n = opts$acc_nfeatures) %>% .$id)
# acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

# Filter features by variance (per stage separately, otherwise you enrich for differences between stages)
acc_dt %>% .[,id_anno:=paste(id,anno,sep="_")]
keep_hv_genes <- acc_dt %>% split(.$stage) %>%
  map(~ split(.,.$anno) %>%
    map(~ .[,.(var=var(rate)), by=c("id_anno")] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id_anno) %>% unlist
  )
acc_dt <- acc_dt %>% split(.$stage) %>% map2(.,names(.), function(x,y) x[id_anno%in%keep_hv_genes[[y]]]) %>% rbindlist

acc_dt <- acc_dt %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Remove lowly expressed genes
rna_dt <- rna_dt[,mean:=mean(expr),by="ens_id"] %>% .[mean>=0.5] %>% .[,mean:=NULL]

# Remove genes with constant expression levels
rna_dt <- rna_dt[,var:=var(expr),by="ens_id"] %>% .[var>0.01] %>% .[,var:=NULL]

# Filter genes with low cellular detection rate and sites with low coverage across samples
rna_dt <- rna_dt[,cdr:=sum(expr>0)/length(opts$rna_cells), by="ens_id"] %>% .[cdr>=opts$rna_min.cdr] %>% .[,cdr:=NULL]

# Extract top N highly variable genes
# keep_hv_genes <- rna_dt[,.(var=var(expr)), by="ens_id"] %>% setorder(-var)  %>% head(n = opts$rna_ngenes) %>% .$ens_id
# rna_dt <- rna_dt[ens_id%in%as.character(keep_hv_genes)]

# Extract top N highly variable genes (per stage)
keep_hv_genes <- rna_dt %>% split(.$stage) %>% 
  map(~ .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var)  %>% head(n=opts$rna_ngenes) %>% .$ens_id) %>%
  unlist %>% unname %>% as.character %>% unique
rna_dt <- rna_dt[ens_id%in%keep_hv_genes]

rna_dt <- rna_dt %>% droplevels()

# Join the data set and save it for biofam

# Prepare data for biofam
data1 <- rna_dt %>% .[,c("sample","stage","gene","expr")] %>%  
  setnames(c("sample","sample_group","feature","value")) %>% .[,c("feature_group"):="RNA"]
data2 <- met_dt %>% .[,c("sample","stage","id","m","anno")] %>%  
  setnames(c("sample","sample_group","feature","value","feature_group")) %>% .[,c("feature","feature_group"):=list(paste0("met_",feature), paste0("met_",feature_group))]
data3 <- acc_dt %>% .[,c("sample","stage","id","m","anno")] %>%  
  setnames(c("sample","sample_group","feature","value","feature_group")) %>% .[,c("feature","feature_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group))]

data <- rbind(data1,data2,data3)

fwrite(data, file=paste0(io$outdir,"/data_test.txt"), col.names=T, quote=F, sep="\t")
# fwrite(data1, file=paste0(io$outdir,"/rna.txt"), col.names=T, quote=F, sep="\t")
# fwrite(data2[sample_group%in%c("Morula","E4.5")], file=paste0(io$outdir,"/met.txt"), col.names=T, quote=F, sep="\t")





# STOP HERE, NO NEED TO CREATE THE MATRICES

# ## Create matrix from the data.table ##
# met_cells <- as.character(unique(met_dt$sample))
# rna_cells <- as.character(unique(rna_dt$sample))
# acc_cells <- as.character(unique(acc_dt$sample))
# 
# rna_matrix <- rna_dt[,c("gene","expr","sample")] %>%
#   .[,c("sample","gene"):=list(as.character(sample),as.character(gene))] %>%
#   .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
#   dcast(sample~gene, value.var="expr", drop=F) %>% matrix.please() %>% t
# 
# met_matrix_list <- list()
# for (n in unique(met_dt$anno)) {
#   met_matrix_list[[paste("met",n,sep="_")]] <- met_dt[anno==n,c("id","gene","m","sample")] %>%
#     .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
#     # .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
#     .[,sample:=factor(sample,levels=Reduce(union,list(acc_cells,met_cells)))] %>%
#     .[,id_gene:=paste(id,gene,sep="_")] %>%
#     dcast(sample~id_gene, value.var="m", drop=F) %>% matrix.please() %>% t
#   
#   
#   cat(sprintf("%s methylation matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
#               nrow(met_matrix_list[[paste("met",n,sep="_")]]), ncol(met_matrix_list[[paste("met",n,sep="_")]]),
#               100*mean(is.na(met_matrix_list[[paste("met",n,sep="_")]]))))
# }
# 
# cat("\n")
# 
# acc_matrix_list <- list()
# for (n in unique(acc_dt$anno)) {
#   acc_matrix_list[[paste("acc",n,sep="_")]] <- acc_dt[anno==n,c("id","gene","m","sample")] %>%
#     .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
#     # .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
#     .[,sample:=factor(sample,levels=Reduce(union,list(acc_cells,met_cells)))] %>%
#     .[,id_gene:=paste(id,gene,sep="_")] %>%
#     dcast(sample~id_gene, value.var="m", drop=F) %>% matrix.please() %>% t
#   
#   
#   cat(sprintf("%s accessibility matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
#               nrow(acc_matrix_list[[paste("acc",n,sep="_")]]), ncol(acc_matrix_list[[paste("acc",n,sep="_")]]),
#               100*mean(is.na(acc_matrix_list[[paste("acc",n,sep="_")]]))))
# }
# 
# # join everything filling with missing values
# all_matrix_list <- c(rna=list(rna_matrix),met_matrix_list,acc_matrix_list)
# # all_matrix_list <- c(met_matrix_list,acc_matrix_list)
