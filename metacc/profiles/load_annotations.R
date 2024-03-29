
# Load genomic annotations
anno_list <- list()
for (anno in names(opts$annos)) {
  tmp <- fread(
    file = sprintf("%s/%s.bed.gz",io$features.dir,anno),
    colClasses = c("V1"="factor", "V2"="integer", "V3"="integer", "V4"="factor", "V5"="factor", "V6"="factor")
  ) %>% setnames(c("chr","start","end","strand","id","anno"))
  
  # Define central position for the window approach
  if (opts$positions[anno] == "start") {
    tmp <- rbind(tmp[strand=="+",.(chr,start,strand,id,anno)] %>% .[,center:=start] %>% .[,c("start"):=NULL], 
                 tmp[strand=="-",.(chr,end,strand,id,anno)] %>% .[,center:=end] %>% .[,c("end"):=NULL]) 
  }
  if (opts$positions[anno] == "center") {
    stopifnot(all(tmp[,end] > tmp[,start]))
    tmp <- tmp[,.(chr,start,end,strand,id,anno)][,center:=round(end+start)/2][,c("start","end"):=NULL]
  }
  if (opts$positions[anno] == "end") {
    tmp <- rbind(tmp[strand=="+",.(chr,end,strand,id,anno)][,center:=end][,c("end"):=NULL], 
                 tmp[strand=="-",.(chr,start,strand,id,anno)][,center:=start][,c("start"):=NULL])
  }
  anno_list[[anno]] <- tmp %>% .[, c("start","end") := list(center-opts$window_size,center+opts$window_size)]
}

anno_df <- rbindlist(anno_list) %>% 
  .[,chr:=as.factor(sub("chr","",chr))]

integer.cols <- c("start","end","center")
anno_df %>%  .[,(integer.cols):=lapply(.SD, as.integer),.SDcols=(integer.cols)]
  
anno_df %>% setkey(chr,start,end)

rm(anno_list)
