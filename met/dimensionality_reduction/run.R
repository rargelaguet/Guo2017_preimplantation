library(MOFA2)
library(data.table)
library(purrr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/met/dimensionality_reduction/load_settings.R")
  reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required=TRUE)
} else if (grepl("ricard",Sys.info()['nodename'])) {
    source("/homes/ricard/Guo2017_preimplantation/met/dimensionality_reduction/load_settings.R")
    # reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required=TRUE)
} else {
  stop("Computer not recognised")
}

###############################
## Load DNA methylation data ##
###############################

print("Loading data...")

met_dt <- lapply(names(opts$annos), function(n)
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n), select=c(1,2,3,5,6)) %>% 
    setnames(c("id_met","id","anno","N","rate")) %>% 
    .[id_met%in%opts$cells] %>%
    .[N>=opts$min.CpGs] %>% .[,N:=NULL]
) %>% rbindlist

################################################
## Merge methylation data and sample metadata ##
################################################

met_dt <- met_dt %>% merge(sample_metadata, by="id_met")

#######################################
## Calculate M value from Beta value ##
#######################################

met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#################
## Filter data ##
#################

# Filter features by coverage
nsamples <- length(unique(met_dt$id_met))
met_dt <- met_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

############################
## Regress out covariates ##
############################

# Global methylation rate
# foo <- met_dt[,.(mean=mean(m)),by=c("id_met","anno")]
# met_dt <- met_dt %>% merge(foo, by=c("id_met","anno")) %>%
#   .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n = opts$nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

###########################
## Prepare data for MOFA ##
###########################

data <- met_dt %>% .[,c("id_met","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>% .[,c("group"):=list("group")]

##########
## Save ##
##########

# file <- paste0(io$outdir,"/data.txt")
# fwrite(data, file, col.names=T, quote=F, sep="\t")
# system(sprintf("pigz -f %s", file))


####################
## Fit MOFA model ##
####################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)

# model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42


#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train the model
model <- run_mofa(MOFAobject, io$outfile)

