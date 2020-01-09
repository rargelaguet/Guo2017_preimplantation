io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

sample_metadata <- fread(io$metadata) 

bar <- sample_metadata %>% merge(foo,by="sample", all.x = T)
foo <- sample_metadata[grepl("embryo",sample_metadata$sample)]



fwrite(bar, io$metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)


