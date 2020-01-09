load_data <- function(id, anno, met_cells, acc_cells, min.cpg = 1, min.gpc = 5) {
  
  # Load DNA methylation data
  met_dt <- fread(
    file = sprintf("%s/%s.tsv.gz",io$met_data_parsed,anno), 
    showProgress=F, header=F,
    select = c("V1"="factor", "V2"= "factor", "V4"="integer", "V5"="integer", "V6"="integer")
  ) %>% .[V1%in%met_cells & V2%in%id] %>% droplevels %>% setnames(c("id_met","id","Nmet","N","rate"))
  
  
  # Load DNA accessibility data
  acc_dt <- fread(
    file = sprintf("%s/%s.tsv.gz",io$acc_data_parsed,anno), 
    showProgress=F, header=F,
    select = c("V1"="factor", "V2"= "factor", "V4"="integer", "V5"="integer", "V6"="integer")
  ) %>% .[V1%in%acc_cells & V2%in%id] %>% droplevels %>% setnames(c("id_acc","id","Nmet","N","rate"))
  
  # Filter by coverage
  met_dt <- met_dt[N>=min.cpg]
  acc_dt <- acc_dt[N>=min.gpc]
  
  return (list("met"=met_dt, "acc"=acc_dt))
}