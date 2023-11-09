create_scam_dalshtaim <- function(path, 
                                  h1_chunk, 
                                  h2_chunk, 
                                  window_size=15, 
                                  activity_threshold=.75,
                                  matrix_subpath="reduced_matrices_full_day",
                                  red_mat_name="DalShtaim") {
  
  
  runs <- list.files(path, recursive=F)
  runs <- runs[which(grepl("\\.R", runs))]
  runs <- sapply(str_split(runs, ".R"), function(l){return(l[[1]])})
  
  matrices_list <- lapply(runs, 
                          function(r) {load(sprintf("%s\\%s.R", path, r)); 
                            return(fmat)})
  
  are_na_matrices <- unlist(lapply(matrices_list, function(mt) {all(is.na(mt))}))
  
  if (sum(are_na_matrices) > 0) {
    print("Removing matrix at indexes:")
    print(which(are_na_matrices))
    matrices_list = matrices_list[which(!are_na_matrices)]
  }
  
  are_same_dim_1 <- unlist(lapply(matrices_list, function(mt) {dim(mt)[1]}))
  are_same_dim_2 <- unlist(lapply(matrices_list, function(mt) {dim(mt)[2]}))
  
  if (!all(are_same_dim_1 == are_same_dim_1[1])) {
    print("Incompatible dimensions!")
    return(list())
  }


  avgd_mat_list <- lapply(matrices_list, function(mat) {time_bin_average(mat, window_size)})
  to_keep_list <- list()  
  
  for (ind in list(h1_chunk, h2_chunk)) {
    
    mean_avgd <- lapply(matrices_list[ind], function(mat) {rowMeans(mat)})
    
    outliers  <- lapply(mean_avgd, function(avg) {which(abs(avg) > activity_threshold & !is.nan(avg) & !is.na(avg))})
    outliers <- unique(unlist(outliers))
    
    to_keep <- 1:nrow(matrices_list[[1]])
    to_keep <- to_keep[which(!to_keep %in% outliers)]

    to_keep_list <- append(to_keep_list,
                           list(to_keep))
  }
  
  
  final_ind <- rep(TRUE, times = nrow(matrices_list[[1]]))
  for (to_keep in to_keep_list) {
    final_ind <- final_ind & (1:nrow(matrices_list[[1]]) %in% to_keep )
  }
  
  final_mat_h1 <- do.call(cbind,
                       lapply(avgd_mat_list[h1_chunk], 
                              function(mat) {return(mat[final_ind,])}))
  
  final_mat_h2 <- do.call(cbind,
                          lapply(avgd_mat_list[h2_chunk], 
                                 function(mat) {return(mat[final_ind,])}))
  
  final_mat_h1_zs  <- smooth_ca_trace(final_mat_h1)
  final_mat_h2_zs  <- smooth_ca_trace(final_mat_h2)
  
  
  final_mat_all <- cbind(final_mat_h1_zs, final_mat_h2_zs)
  
  #idx <- sample(1:nrow(final_mat_h1_zs), 1); plot(rollmean(c(final_mat_h1_zs[idx,], final_mat_h2_zs[idx,]), 500), type='l')
  final_mat_all <- final_mat_all[!apply(final_mat_all, 1, function(r) {sum(is.na(r)) > 0}),]
  
  reduced_20 <- reduce_dim(t(final_mat_all), "lem", knn1=0.075, ndim=20)
  reduced_6 <- reduce_dim(reduced_20, "lem", knn1=0.025, ndim=6)
  
  reduced <- reduced_6
  
  dir.create(sprintf("%s\\%s", path, matrix_subpath))
  save(reduced, file=sprintf("%s\\%s\\%s.R", path, matrix_subpath, red_mat_name))
  print(sprintf("Saving mat! %s", "DalShtaim (scam)"))
}
