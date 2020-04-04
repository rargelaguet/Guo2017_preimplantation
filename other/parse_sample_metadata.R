library(data.table)
library(purrr)

io <- list()
io$basedir <- "/Users/ricard/data/Guo2017_preimplantation"
io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

sample_metadata <- fread(io$metadata) 

# Check for duplicate cells
# Define assay (scMT-seq, scNMT-seq, scRNA-seq)
strsplit(sample_metadata$sample, split="_")
sample_metadata[,embryo:="scRNA-seq"]
sample_metadata[!is.na(id_met) & !is.na(id_acc),assay:="scNMT-seq"]
sample_metadata[!is.na(id_met) & is.na(id_acc),assay:="scMT-seq"]
table(sample_metadata$assay)


fwrite(sample_metadata, file=io$metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)



###

table(sample_metadata[stage=="E6.5",lineage10x_2])

sample_metadata[!is.na(id_acc)] %>% View
