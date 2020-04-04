suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(MOFA2))

############################
## Define I/O and options ##
############################

source("/Users/ricard/Guo2017_preimplantation/mofa/load_settings.R")

if (file.exists(paste0(io$basedir,"/mofa/rds/model.rds"))) {
  
  model <- readRDS(paste0(io$basedir,"/mofa/rds/model.rds"))
  
} else {
  
  ################
  ## Load model ##
  ################
  
  file <- paste0(io$basedir,"/mofa/hdf5/model.hdf5")
  model <- load_model(file)
  
  ##################
  ## Rename views ##
  ##################
  
  # opts$rename.views <- c(
  #   "met_prom_2000_2000$" = "Promoter methylation",
  #   "acc_prom_200_200$" = "Promoter accessibility"
  # )
  # views(model) = stringr::str_replace_all(views(model), opts$rename.views)
  # model <- subset_views(model, views=unname(opts$rename.views) )
  
  #########################
  ## Add sample metadata ##
  #########################
  
  cells <- as.character(unname(unlist(MOFA2::samples(model))))
  
  sample_metadata_filt <- sample_metadata %>% setkey(sample) %>% .[cells] %>% 
    .[,group:="group"]
  
  stopifnot(all(cells==sample_metadata_filt$sample))
  samples_metadata(model) <- sample_metadata_filt
  
  ####################
  ## Subset factors ##
  ####################
  
  # r2 <- model@cache$variance_explained$r2_per_factor
  # factors <- sapply(r2, function(x) x[,"RNA expression"]>0.005)
  # model <- subset_factors(model, which(apply(factors,1,sum)>=1))
  # 
  # factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")
  
  ##########
  ## Save ##
  ##########
  
  # saveRDS(model, paste0(io$basedir,"/mofa/rds/model.rds"))
  
}
