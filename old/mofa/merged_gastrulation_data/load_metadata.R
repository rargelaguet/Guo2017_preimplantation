
###################
## Load metadata ##
###################

# Load sample metadata
sample_metadata1 <- fread(io$sample.metadata1,stringsAsFactors=T) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage,sep="_"))] %>%
  .[id_met%in%opts$met_cells1 | id_acc %in% opts$acc_cells1 | id_rna %in% opts$rna_cells1 ] %>%
  droplevels()
sample_metadata2 <- fread(io$sample.metadata2,stringsAsFactors=T) %>%
  .[id_met%in%opts$met_cells2 | id_acc %in% opts$acc_cells2 ] %>%
  .[,id_rna:=NA] %>%
  droplevels()

# Load annotation metadata
feature_metadata <- lapply(unique(c(opts$met.annos,opts$acc.annos)), function(i) 
  fread(sprintf("%s/%s.bed",io$annos_dir,i), stringsAsFactors=T)[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

# Load gene metadata 
gene_metadata <- fread(io$gene_metadata,stringsAsFactors=T) %>% 
  setnames(c("ens_id","symbol"),c("id","gene")) %>% 
  .[,chr:=as.factor(stringr::str_replace_all(chr,"chr",""))]


####################
## Parse metadata ##
####################

# Merge all pre-implantation stages
sample_metadata2[stage%in%opts$stage_lineage2,stage:="Preimplantation"]

# Merge ICM and Morula
# sample_metadata2[stage%in%c("ICM","Morula"),stage:="Morula_ICM"]

# 
# sample_metadata2[stage%in%c("ICM","Morula"),stage:="Morula_ICM"]

# Merge

# Concatenate the two data sets
# sample_metadata <- rbind(
#   sample_metadata1[,c("id_rna","id_acc","id_met","sample","stage_lineage")],
#   sample_metadata2[,c("id_rna","id_acc","id_met","sample","stage")] %>% setnames(c("stage"),c("stage_lineage"))
# )

sample_metadata <- rbind(
  sample_metadata1[,c("id_rna","id_acc","id_met","sample","stage","stage_lineage")],
  sample_metadata2[,c("id_rna","id_acc","id_met","sample","stage")] %>% .[,stage_lineage:=stage]
)
