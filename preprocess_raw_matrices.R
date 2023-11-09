

paths <- list("Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190903_IC71\\190903_IC71_run1\\IC71_190903_001.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190903_IC71\\190903_IC71_run2\\IC71_190903_002.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190903_IC71\\190903_IC71_run3\\IC71_190903_003.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190904_IC71\\190904_IC71_run1\\IC71_190904_001.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190904_IC71\\190904_IC71_run2\\IC71_190904_002.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190904_IC71\\190904_IC71_run3\\IC71_190904_003.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190906_IC71\\190906_IC71_run1\\IC71_190906_001.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190906_IC71\\190906_IC71_run2\\IC71_190906_002.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190906_IC71\\190906_IC71_run3\\IC71_190906_003.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190906_IC71\\190906_IC71_run4\\IC71_190906_004.mat2",
              "Y:\\livneh\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_7\\IC71\\190906_IC71\\190906_IC71_run5\\IC71_190906_005.mat2")

preprocess_matrices <- function(path, old=F){
  
  #output_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\temp\\"
  output_path <- "Y:\\livneh\\itayta\\data"
  mt <- readMat(path)
  cellsort <- mt$cellsort
  fmat <- c()
  # 
  # if (dim(cellsort)[1] == 1) {
  #   s <- which(dimnames(mt$cellsort[1,1,1]$timecourse)[[1]] == "dff.axon")  
  # } else {
  #   s <- which(dimnames(mt$cellsort[12,1,1]$timecourse)[[1]] == "dff.axon")  
  # }
  
  s <- which(dimnames(mt$cellsort[,1,1]$timecourse)[[1]] == "dff.axon") 
  
  for (i in 1:dim(cellsort)[3]) {
    # if (dim(cellsort)[1] == 1) {
    #   fmat <- rbind(fmat, as.vector(cellsort[1,1,i]$timecourse[[s]]))
    # } else {
    #   fmat <- rbind(fmat, as.vector(cellsort[12,1,i]$timecourse[[s]]))
    # }
    
    fmat <- rbind(fmat, as.vector(cellsort[,1,i]$timecourse[[s]]))
    print(sprintf("%d",i))
  }
  
  
  rm(mt)
  
  if(old){
    base_str <- str_split(path, "\\\\")
    base_str <- base_str[[1]][len(base_str[[1]])]
    
    split_base_str <- unlist(str_split(base_str, "_"))
    print(split_base_str)
    nrun <- str_split(split_base_str[3], "run")[[1]][2]
    
    name <- sprintf("%s_%s_00%s", split_base_str[2], split_base_str[1], nrun)  
    final_output_path <- sprintf("%s\\%s.R", output_path, str_replace_all(name, "/", "_"))
    
  } else {
    fpath <- paste(output_path, str_split(path, "\\\\")[[1]][11], sep="\\")
    final_path <- paste(strsplit(fpath, ".mat")[[1]][1], ".R", sep="")
    final_output_path <- final_path[[1]]
  }
  
  print(sprintf("Saving to %s", final_output_path))
  print(dim(fmat))
  save(fmat, file=final_output_path)
  
  return(final_output_path)
}

get_path_by_prefix_vector <- function(prefix_vec) {
  return(sprintf("%s\\%s.R", output_path, prefix_vec))
}

# 
# 
# fpaths <- rep(F, times=50)
# 
# for (mice in names(day_list)) {
#   relevant_days <- day_list[[mice]]
# 
#   mice_ind <- grep(mice, paths_all)
#   mice_paths <- paths_all[mice_ind]
#   
#   for (r_day in relevant_days) {
#     fpaths[mice_ind] <- fpaths[mice_ind]  | grepl(r_day, mice_paths)
#   }
# }

