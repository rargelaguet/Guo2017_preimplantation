source("/Users/ricard/Guo_2017/settings.R")

# 
opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

opts$annos <- c(
  "prom_2000_2000" = "Promoters"
)

# Define window positions and characteristics
opts$positions <- c(
  "prom_2000_2000"="center"
)

opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 100

# Define which cells to  use
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages]

opts$met.cells <- sample_metadata[,id_met]
opts$acc.cells <- sample_metadata[,id_acc]

