if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
} else if (grepl("yoda",Sys.info()['nodename'])) {
  source("/homes/ricard/Guo2017_preimplantation/settings.R")
} else {
  stop("Computer not recognised")
}


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
opts$met.tile <- 100
opts$acc.tile <- 50

# Define which cells to  use
sample_metadata <- sample_metadata %>% 
  .[,c("sample","id_acc","id_met","stage")] %>%
  .[stage%in%opts$stages]

opts$met.cells <- sample_metadata[,id_met] %>% as.character
opts$acc.cells <- sample_metadata[,id_acc] %>% as.character

