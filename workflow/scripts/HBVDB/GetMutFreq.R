GetAllIDs <- function(clustal_path){
  #Read Alignment
  aln <- seqinr::read.alignment(clustal_path,format = 'clustal')
  names <- aln$nam
  
  #Parse IDs
  IDs <- sapply(names,function(x) strsplit(strsplit(x,split = '\\|')[[1]][3],split = '_')[[1]][1])
  
  #write ID file
  write(IDs,file = paste0(clustal_path,'.id.txt'))
  return(list(aln=aln,IDs=IDs))
}
ParseGenbank <- function(gb_path){
  #Read in gb file
  gb_file <- readLines(gb_path)

  #Split entries of gb file
  split_points <- which(gb_file=='//')
  split_points <- c(1,split_points)
  gb_split <- lapply(2:length(split_points),function(i)gb_file[(split_points[i-1]):(split_points[i])])
  
  seq_df <- data.frame(ID = character(),Country = character(),stringsAsFactors = F)
  #Parse gb file
  for(i in 1:length(gb_split)){
    #Extract entry
    cur_entry <- gb_split[[i]]
    #Extract accession number
    cur_accession <- gsub(x=cur_entry[grep(pattern = 'ACCESSION   ',x=cur_entry)],pattern = 'ACCESSION   ',replacement = '')
    
    country_entry <- grep(pattern = '/country=',x=cur_entry)
    if(length(country_entry) == 0){
      #if country field unavailable, set as NA
      cur_country <- NA
    }else{
      #Extract country field
      cur_country <- strsplit(strsplit(cur_entry[country_entry],split = '\"')[[1]][2],split = '\\"')[[1]][1]
      #Ignore region info
      if(grepl(pattern = ':',x=cur_country)){
        cur_country <- strsplit(cur_country,split = ":")[[1]][1]
      }
    }
    
    seq_df <- rbind(seq_df,data.frame(ID = cur_accession,Country = cur_country))
  }
  return(seq_df)
}
GetResidualOfInterest <- function(aln,positions,metadata,regions){
  residual_info <- vector(mode = 'list',length = length(positions))
  names(residual_info) <- names(positions)
  
  for(i in 1:length(positions)){
    cur_GT <- names(positions)[i]
    cur_aln <- aln[[cur_GT]]$aln$seq
    cur_IDs <- aln[[cur_GT]]$IDs
    cur_metadata <- metadata[[cur_GT]]
    
    df <- data.frame(ID = cur_IDs,stringsAsFactors = F)
    
    for(k in 1:length(positions[[i]])){
      cur_pos <- positions[[i]][k]
      cur_residuals <- sapply(cur_aln,function(x) strsplit(x,split = '')[[1]][cur_pos])
      cur_df <- data.frame(residual = cur_residuals)
      colnames(cur_df) <- paste0('Mut_',k)
      df <- cbind(df,cur_df)
    }
    residual_info[[i]] <- dplyr::left_join(df,cur_metadata,by=c('ID'='ID'))
    residual_info[[i]]$Region <- unlist(sapply(residual_info[[i]]$Country,function(x) {
      if(!is.na(x)){
        return(names(regions)[sapply(regions,function(y) x%in%y)])
      }else{
        return(NA)
      }
    }))
  }
  return(residual_info)
}
regions <- list('East Asia' = c('Japan','China','Taiwan','South Korea','Hong Kong'),
                'South East Asia' = c('Cambodia','Indonesia','Laos','Malaysia','Myanmar','Philippines','Thailand','Viet Nam'),
                'South Asia' = c('Bangladesh','India'),
                'Central Asia' = c('Uzbekistan','Mongolia','Nepal','Tajikistan'),
                'Europe' = c('Belgium','Russia','United Kingdom','Greenland','Spain','Belarus','France','Germany','Italy','Latvia',
                             'Netherlands','Poland','Serbia','Turkey','Sweden','Serbia','Estonia'),
                'Middle East and North Africa' = c('Tunisia','Syria','Pakistan','Lebanon','Iran','Egypt','United Arab Emirates'),
                'Sub-Saharan Africa' = c('Botswana','Cameroon','Cape Verde','Central African Republic','Democratic Republic of the Congo',
                                         'Ethiopia','Gabon','Gambia','Sudan','Guinea','Kenya','Malawi','Nigeria','Rwanda','Senegal','Somalia',
                                         'South Africa','Tanzania','Uganda','Zimbabwe'),
                'North America'= c('USA','Canada'),
                'Central and South America'=c('Bolivia','Brazil','Panama','Argentina','Colombia','Cuba','Haiti','Uruguay','Venezuela'),
                'Polynesia and Micronesia' = c('New Caledonia','Papua New Guinea','Tonga','Fiji'),
                'Australia and New Zealand' = c('Australia','New Zealand'))

genotypes <- c('A','B','C','D')
aln <- lapply(genotypes,function(x) GetAllIDs(paste0('~/G2G-HBV/HBV_DB/MSA/',x,'_LHBs.clu')))
names(aln) <- genotypes
metadata <- lapply(genotypes,function(x) {ParseGenbank(paste0('~/G2G-HBV/HBV_DB/MSA/',x,'_LHBs.gb'))})
names(metadata) <- genotypes


pos_of_interest <- list('A' = c(17,35,51),'B' = c(17,35,51),'C' = c(20,38,54),'D'=c(7,25,41))
residual_info <- GetResidualOfInterest(aln,pos_of_interest,metadata,regions)

