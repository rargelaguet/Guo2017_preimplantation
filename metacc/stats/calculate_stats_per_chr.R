# Define settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/metacc/stats/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
} else {
  stop("Computer not recognised")
}

# Load methylation data and summarise
met_dt <- lapply(opts$annos, function(n) {
  print(n)
  
  file <- sprintf("%s/%s.bed.gz",io$features.dir,n)
  if (file.exists(file)) {
    anno_dt <- fread(file, select=c(1,5)) %>% 
      setnames(c("chr","id")) %>% .[chr!="MT"]
    
    file = sprintf("%s/%s.tsv.gz",io$met_data_parsed,n)
    if (file.exists(file)) {
      fread(file, showProgress = F, header = F,
        select = c("V1"="character","V2"="character","V3"="factor","V4"="integer","V5"="integer","V6"="integer")
        ) %>% setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>%
        merge(anno_dt,by=c("id")) %>%
        .[,.(rate=100*(sum(Nmet)/sum(Ntotal))),by = c("id_met","anno","chr")]
    }
  }
}) %>% rbindlist %>% merge(sample_metadata[,c("id_met","sample")]) %>% .[,id_met:=NULL]


# Load accessibility data and summarise
acc_dt <- lapply(opts$annos, function(n) {
  print(n)
  
  file <- sprintf("%s/%s.bed.gz",io$features.dir,n)
  if (file.exists(file)) {
    anno_dt <- fread(file, select=c(1,5)) %>% 
      setnames(c("chr","id")) %>% .[chr!="MT"]
    
    file = sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n)
    if (file.exists(file)) {
      fread(file, showProgress = F, header = F,
            select = c("V1"="character","V2"="character","V3"="factor","V4"="integer","V5"="integer","V6"="integer")
      ) %>% setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>%
        merge(anno_dt,by=c("id")) %>%
        .[,.(rate=100*(sum(Nmet)/sum(Ntotal))),by = c("id_acc","anno","chr")]
    }
  }
}) %>% rbindlist %>% merge(sample_metadata[,c("id_acc","sample")]) %>% .[,id_acc:=NULL]

# Concatenate
metacc_dt <- rbind(met_dt[,context:="CG"], acc_dt[,context:="GC"])

# Save
fwrite(metacc_dt, paste0(io$outdir,"/metacc_stats_anno_chr.tsv.gz"), sep="\t")
