suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(MOFA2))

############################
## Define I/O and options ##
############################

source("/Users/ricard/Guo_2017/met/dimensionality_reduction/load_settings.R")

################
## Load model ##
################

file <- paste0(io$basedir,"/mofa/hdf5/model_prom_2000_2000_E5_E6.hdf5")
model <- load_model(file)

##################
## Rename views ##
##################

# opts$rename.views <- c(
#   "met_Promoters" = "Promoter methylation",
#   # "met_Gene bodies" = "Genebody methylation",
#   # "met_E3.5 enhancers" = "E3.5 enhancer methylation",
#   "met_E7.5 enhancers" = "Enhancer methylation",
#   # "met_E10.5 enhancers" = "E10.5 enhancer methylation",
#   
#   "RNA$" = "RNA expression"
# )
# 
# views(model) = stringr::str_replace_all(views(model), opts$rename.views)
# 
# model <- subset_views(model, views=unname(opts$rename.views) )

#########################
## Add sample metadata ##
#########################

cells <- as.character(unname(unlist(MOFA2::samples(model))))

sample_metadata_filt <- sample_metadata %>% setkey(id_met) %>% .[cells] %>%
  .[,sample:=id_met] %>%
  .[,group:="group"]

stopifnot(all(cells==sample_metadata_filt$sample))
samples_metadata(model) <- sample_metadata_filt

####################
## Subset factors ##
####################

# r2 <- model@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# model <- subset_factors(model, which(apply(factors,1,sum)>=1))
# factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")

plot_variance_explained(model)
plot_factor(model, factors = 1:2, color_by = "stage", add_violin = T, dodge=T)
plot_factors(model, factors = 1:5, color_by = "stage")
plot_factors(model, factors = c(1,2), color_by = "stage")
plot_weights(model, view=1, factor = 1, nfeatures = 0, scale = F)

cor(model@samples_metadata$mean, model@expectations$Z$group)
# plot_data_scatter(model, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)


