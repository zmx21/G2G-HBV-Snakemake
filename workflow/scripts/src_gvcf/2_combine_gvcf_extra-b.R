## Aim: find all files in vcfs_combine.list that are not stored in any of the 
##      vcfs_combine_split{COUNTER}.list files, and then store it in 
##      vcfs_combine_split31.list


dat_full <- read.table("/scratch/rueger/gvcf_tmp/vcfs_combine.list", stringsAsFactors = FALSE)

dat_list <- NULL
for (i in 1:30){
  if(i < 10) {
    tmp <- read.table(paste0("/scratch/rueger/gvcf_tmp/vcfs_combine_split0", i,".list"), stringsAsFactors = FALSE)
  }else{
    tmp <- read.table(paste0("/scratch/rueger/gvcf_tmp/vcfs_combine_split", i,".list"), stringsAsFactors = FALSE)
  }
  dat_list <- rbind(dat_list, tmp)
}

dat_remain <- dat_full[!(dat_full$V1 %in% dat_list$V1), ]
write.table(dat_remain, paste0("/scratch/rueger/gvcf_tmp/vcfs_combine_split31.list"), col.names = FALSE, row.names = FALSE, quote = FALSE)
