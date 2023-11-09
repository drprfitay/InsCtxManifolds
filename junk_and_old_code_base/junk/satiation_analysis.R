

get_licking_vec <-  function(path, lick_file_idx, verbose=F)  {
 
  search_path <- sprintf("%s//licking",path)
  licking_files = list.files(search_path)[(grep("IC.*.licking", list.files(search_path)))] 
  
  if(verbose){
    print(licking_files)
    print(lick_file_idx)
  }
  
  regexp_str = "00[%d]{1}"
  
  # Unsafe
  if (lick_file_idx > 9) {
    regexp_str = "01[%d]{1}"
    lick_file_idx = lick_file_idx - 10
  }
  
  requested_file <- licking_files[which(regexpr(sprintf(regexp_str, lick_file_idx), licking_files) != -1)]
  
  if (verbose) {
    print(requested_file)
  }
  
  lick_vec <- readMat(sprintf("%s//%s", search_path, requested_file))
  
  return(as.vector(lick_vec$lick.vec))
}

get_lick_bouts <- function(lick_vec, window_size=15, interval=30, unbinned_bouts=F) {
  
  licks_indices <- which(lick_vec > 0)
  licks_interval <- diff(licks_indices)
  lick_bouts <- licks_indices[which(licks_interval > interval)  + 1]
  
  if (!unbinned_bouts) {
    binned_bouts = unlist(lapply(lick_bouts, function(i) {get_binned_index(i, window_size)}))
  } else {
    binned_bouts <- unlist(lick_bouts)
  }
  
  return(binned_bouts)  
  
}

plot_licking_bouts <- function(licking_vector, unbinned_licking_bouts) {
  
  g <- 
    ggplot()
  
  for (bout_i in 1:len(unbinned_bouts)){
    indices_of_bout <- (unbinned_bouts[bout_i] - 1 * 30):((unbinned_bouts[bout_i] + 3 * 30))
    
    
    lv <- licking_vector[indices_of_bout]
    lv[indices_of_bout > len(licking_vector)] <- 0
    
    plot_df <- data.frame(x=(1:len(lv))[lv != 0],
                          y=rep((len(unbinned_bouts) + 1) - bout_i, times=sum(lv > 0)))
    
    g <- g + geom_point(data=plot_df, aes(x=x, y=y), alpha=.5)
  }
  
  g <- 
    g + theme(line=element_blank(),
              rect=element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.margin = margin(-10, 0, -10, 0, "pt")) +
    ylab("Bout (#)") +
    xlab("Time (seconds")
  
  
  return(g)  
}



satiation_run_indices_all = c(5,4,3, 3)
satiation_paths <- get_satiation_paths()

for (i in 1:len(am)) {
  
  metadata <- am[[i]] 
  path  <-  satiation_paths[i]
  satiation_run_idx <- satiation_run_indices_all[i] 
  licking_vec <- get_licking_vec(path, satiation_run_idx)
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)

  frames_per_mat <- lapply(actual_runs,
                         function(r_path)
                         {
                           load(sprintf("%s\\%s", path, r_path))
                           return(ncol(fmat))
                         })


  binned_frames_per_mat <- unlist(frames_per_mat) / 15
  binned_frames_to_add <- cumsum(binned_frames_per_mat) - binned_frames_per_mat[1]




  unbinned_bouts <- get_lick_bouts(licking_vec, interval=30, unbinned_bouts = T)
  binned_bouts <- get_lick_bouts(licking_vec, interval=30)
  bouts_indices <-  binned_bouts + binned_frames_to_add[satiation_run_idx]
  
  # overlp <- lapply(bouts_indices, function(ind) {(ind - 1):(ind + 4)})
  # bouts_indices <- bouts_indices[!unlist(lapply(2:len(overlp), function(idx) {sum(sapply(overlp[[idx - 1]], function(tp)  {sum(tp %in% overlp[[idx]] > 0)})) > 0}))]
  bouts_mat <- matrix(rep(0, times=5 * len(bouts_indices)), nrow=len(bouts_indices))

  for (idx in 1:len(bouts_indices)) {
    bout <- bouts_indices[idx]
    bout_trial_ind <- (bout - 1):(bout + 3)
    bouts_mat[idx,] <- metadata$cluster_mat$labs[bout_trial_ind]
  }

  consecutive_seconds <- frames_per_mat[[satiation_run_idx]] / 15 / 20
  consecutive_mat <- matrix(rep(0, times=20 * (consecutive_seconds)), nrow=consecutive_seconds)
  
  for (j in 1:nrow(consecutive_mat)) {
    consecutive_mat[j,] <- metadata$cluster_mat$labs[((binned_frames_to_add[satiation_run_idx] + (j  - 1) * 20) + 1):(binned_frames_to_add[satiation_run_idx] + j * 20)]
  }
    

  bouts <- pheatmap(bouts_mat, cluster_rows=F, cluster_cols=F, legend = F, border_col=NA)
  
  pheatmap(consecutive_mat, cluster_rows=F, cluster_cols=F)


  reward_mat <- metadata$trials_mat[metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0,]

  phr <- pheatmap(reward_mat, cluster_rows=F, cluster_cols=F, legend = F, border_col=NA, show_colnames = F, show_rownames = F)
  
  all_mat <- rbind(reward_mat, bouts_mat)

  annot_df <- data.frame(Type=rep(c("Structured", "Unstructured"), times=c(nrow(reward_mat), nrow(bouts_mat))))

  rownames(annot_df) <- 1:nrow(all_mat)
  rownames(all_mat) <- 1:nrow(all_mat)
  pheatmap(all_mat, cluster_rows=F, cluster_cols=F, show_rownames = F, annotation_row = annot_df, border_col=NA)

  reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]
  plot_grid(phr[[4]], plot_structured_licking_bouts(path, reward_trials, frames_per_mat))
  plot_grid(bouts[[4]], plot_licking_bouts(licking_vec, unbinned_bouts))
}


plot_structured_licking_bouts <-  function(path, reward_trials_stim_mat, frames_per_mat) {
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  run_indices <- unique(reward_trials_stim_mat[,ncol(reward_trials_stim_mat)])

  lick_vectors_all <- 
  lapply(run_indices, function(run_ind) {print(run_ind);get_licking_vec(path, run_ind)})
  
  names(lick_vectors_all) <- run_indices
  
  
  g <- ggplot() 
  
  for (reward_trial_idx in 1:nrow(reward_trials_stim_mat)) {
    
    relevant_trial <- reward_trials_stim_mat[reward_trial_idx, ]
    trial_run <- relevant_trial[len(relevant_trial)]
    relevant_lick_vec <- lick_vectors_all[[as.character(trial_run)]]
    frame_in_run <- relevant_trial["Frames"] - frames_to_add[trial_run]
    
    
    lv <- relevant_lick_vec[frame_in_run:(frame_in_run + 300)]
   
    
    plot_df <- data.frame(x=(1:len(lv))[lv != 0],
                          y=rep((nrow(reward_trials_stim_mat) + 1) - reward_trial_idx, times=sum(lv > 0)))
    
    g <- g + geom_point(data=plot_df, aes(x=x, y=y), alpha=.5) 
  }
  
  g <- 
    g + theme(line=element_blank(),
              rect=element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.margin = margin(-10, 0, -10, 0, "pt")) +
    ylab("Bout (#)") +
    xlab("Time (seconds")
  
  
  
  return(g)  
}


get_task_vs_freely_lick_distribution <- function(metadata, path, satiation_run_idx, extra_s=3, bouts_interval=30) {
  
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
                             
  reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]
  run_indices <- unique(reward_trials[,ncol(reward_trials)])
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  
  binned_frames_per_mat <- unlist(frames_per_mat) / 15
  binned_frames_to_add <- cumsum(binned_frames_per_mat) - binned_frames_per_mat[1]
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {get_licking_vec(path, run_ind)})
  
  
  names(lick_vectors_all) <- run_indices
  
  
  task_lick_bouts_mat <- c()
  task_licking_mat <- c()
  neuronal_cluster_task_licks_dist <- list()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))]
    
    
    if( sum(licks_in_trial_unbinned > 0) <= 0 ) {
      print("No licks in trial continue")
      next
    }
    
    lick_in_trial_binned = which(time_bin_average_vec(licks_in_trial_unbinned, 15) > 0)
    
    lick_in_trial_binned <- unique(c(lick_in_trial_binned, max(lick_in_trial_binned): (max(lick_in_trial_binned) + extra_s)))
    
    lick_in_trial_binned <- lick_in_trial_binned[lick_in_trial_binned <= 20]
    
    cluster_labels_of_trial <- metadata$trials_mat[as.character(which(metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0)[reward_trial_idx]),]
    
    neuronal_cluster_task_licks_dist <- append(neuronal_cluster_task_licks_dist, list(cluster_labels_of_trial[lick_in_trial_binned]))
    
    
    
    unbinned_first_lick_in_respone <- which(licks_in_trial_unbinned > 0)[which(which(licks_in_trial_unbinned > 0) > 60)[1]]
    binned_first_lick_in_response <- get_binned_index(unbinned_first_lick_in_respone, 15)
    
    if (binned_first_lick_in_response + 3 > 20) {
      binned_first_lick_in_response = 20 - 3
    }
    
    task_lick_bouts_mat <- rbind(task_lick_bouts_mat,
                            cluster_labels_of_trial[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 8)])
    
    task_licking_mat <- rbind(task_licking_mat,
                              as.numeric(licks_in_trial_unbinned[(unbinned_first_lick_in_respone - 60):(unbinned_first_lick_in_respone + 120)]>0))
    
  }
  
  
  
  
  freely_licking_vec <- get_licking_vec(path, satiation_run_idx)
  freely_licking_indices <- which(freely_licking_vec > 0) + frames_to_add[satiation_run_idx]
  binned_freely_licking <- unique(unlist(lapply(freely_licking_indices, function(i) {get_binned_index(i, 15)})))
  #binned_freely_licking <- unique(unlist(lapply(binned_freely_licking, function(i) {i:(i+extra_s)})))
  
  neuronal_cluster_task_licks_dist <- table(unlist(neuronal_cluster_task_licks_dist))
  neuronal_cluster_freely_licks_dist <- table(metadata$cluster_mat$labs[binned_freely_licking])
  
  
  freely_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  task_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  
  names(freely_licks_cluster) <- c(-1, 1:(len(freely_licks_cluster) - 1))
  names(task_licks_cluster) <- c(-1, 1:(len(task_licks_cluster) - 1))
  
  task_licks_cluster[names(neuronal_cluster_task_licks_dist)] <- neuronal_cluster_task_licks_dist
  freely_licks_cluster[names(neuronal_cluster_freely_licks_dist)] <- neuronal_cluster_freely_licks_dist
  
  freely_binned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval)
  freely_unbinned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval, unbinned_bouts = T)
  freely_bouts_indices <-  freely_binned_bouts + binned_frames_to_add[satiation_run_idx]
  
  freely_lick_bouts_mat <- c()
  freely_licking_mat <- c()
  freely_indices_mat <- c()
  
  for (idx in 1:len(freely_binned_bouts)) {
    bout <- freely_binned_bouts[idx]
    bout_trial_ind <- (bout - 3):(bout + 8)
    freely_lick_bouts_mat <- rbind(freely_lick_bouts_mat,
                              metadata$cluster_mat$labs[bout_trial_ind])
    
    freely_licking_mat <- rbind(freely_licking_mat, 
                                as.numeric(freely_licking_vec[(freely_unbinned_bouts[idx] - 60):(freely_unbinned_bouts[idx] + 120)] > 0))
    
    freely_indices_mat <- rbind(freely_indices_mat,
                            (freely_unbinned_bouts[idx] - 60):(freely_unbinned_bouts[idx] + 120))    
  }

  
  is_done = FALSE
  tmp_freely_indices_mat = freely_indices_mat
  rownames(tmp_freely_indices_mat) <- 1:nrow(freely_indices_mat)
  ri = 1
  
  while (!is_done) {
    
    if (ri == nrow(tmp_freely_indices_mat)) {
      is_done = TRUE
      print("Done!")
      next
    }
    
    overlap = sum(tmp_freely_indices_mat[ri,] %in% tmp_freely_indices_mat[ri+1,])
    if (overlap > 0) {
      print(sprintf("Removing row %d due to %d overlap to row %d", ri+1, overlap, ri))
      tmp_freely_indices_mat <- tmp_freely_indices_mat[-(ri+1),]
    } else {
      ri <- ri +1
    }
  }
  
  #keep = rowSums(apply(freely_indices_mat, 1, function(r) {apply(freely_indices_mat, 1, function(r2) {sum(max(r) %in% r2)})})) <= 1
  keep = as.numeric(rownames(tmp_freely_indices_mat))
  
  freely_licking_mat <- freely_licking_mat[keep,]
  freely_lick_bouts_mat <- freely_lick_bouts_mat[keep,]
  
  return(list(freely=freely_licks_cluster,
              task=task_licks_cluster,
              freely_licking_mat=freely_licking_mat,
              task_licking_mat=task_licking_mat,
              freely_cluster_mat=freely_lick_bouts_mat,
              task_cluster_mat=task_lick_bouts_mat))
}


get_task_vs_agrp_lick_distribution <- function(metadata, path, agrp_indices, extra_s=3) {
  
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  
  reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]
  run_indices <- unique(reward_trials[,ncol(reward_trials)])
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {get_licking_vec(path, run_ind)})
  
  
  names(lick_vectors_all) <- run_indices
  
  
  pre_AGRP <- list()
  post_AGRP <- list()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))]
    
    print(run %in% mainpulation_run_indices)
    
    #print(run)
    if(sum(licks_in_trial_unbinned > 0) <= 0 ) {
      #print("No licks in trial continue")
      next
    }
    
    
    lick_in_trial_binned = which(time_bin_average_vec(licks_in_trial_unbinned, 15) > 0)
    
    lick_in_trial_binned <- unique(c(lick_in_trial_binned, max(lick_in_trial_binned): (max(lick_in_trial_binned) + extra_s)))
    
    lick_in_trial_binned <- lick_in_trial_binned[lick_in_trial_binned <= 20]
    
    
    cluster_labels_of_trial <- metadata$trials_mat[as.character(which(metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0)[reward_trial_idx]),]
    
    
    
    
    if (run %in% agrp_indices) {
      post_AGRP <- append(post_AGRP, list(cluster_labels_of_trial[lick_in_trial_binned]))  
      print(sprintf("Run %d in agrp", run))
    } else {

      pre_AGRP <- append(pre_AGRP, list(cluster_labels_of_trial[lick_in_trial_binned]))
    }
    
  }
  
  
  
  

  
  pre_AGRP_cluster_task_licks_dist <- table(unlist(pre_AGRP))
  post_AGRP_cluster_task_licks_dist <- table(unlist(post_AGRP))
  
  
  
  pre_AGRP_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  post_AGRP_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  
  names(pre_AGRP_licks_cluster) <- c(-1, 1:(len(pre_AGRP_licks_cluster) - 1))
  names(post_AGRP_licks_cluster) <- c(-1, 1:(len(post_AGRP_licks_cluster) - 1))
  
  pre_AGRP_licks_cluster[names(pre_AGRP_cluster_task_licks_dist)] <- pre_AGRP_cluster_task_licks_dist
  post_AGRP_licks_cluster[names(post_AGRP_cluster_task_licks_dist)] <- post_AGRP_cluster_task_licks_dist
  
  
  return(list(pre_AGRP=pre_AGRP_licks_cluster,
              post_AGRP=post_AGRP_licks_cluster))
}

all_results <- list()
satiation_run_indices_all = c(5,4,3, 3)
satiation_run_indices_interval <- c(30, 10, 30, 10)
satiation_paths <- get_satiation_paths()
satiation_metadata <- metadata_sacross_mice_decoding_build_metadata(satiation_paths)


JSD_satiation <- c()
KL_satiation <- c()
cor_satiation <- c() 

write_path <- sprintf("%s\\figure_3\\", figures_base_path)
dir.create(write_path)

write_path <- sprintf("%s\\structure_annotation", write_path)
dir.create(write_path)  
datasets_names = paste("SATIATION_", get_datasets_names(satiation_paths, sep="_"), sep="")

for (i in 1:len(satiation_metadata)) {
  res = get_task_vs_freely_lick_distribution(satiation_metadata[[i]], 
                                             satiation_paths[i], 
                                             satiation_run_indices_all[i], 
                                             bouts_interval = satiation_run_indices_interval[i], 
                                             extra_s = 3)
  
  
  # barplot(res$freely / sum(res$freely))
  # barplot(res$task / sum(res$task), col=adjustcolor("red", alpha=.5), add=T)
  
  KL_satiation <- c(KL_satiation,
                    KL_dist((res$freely / sum(res$freely)),
                            (res$task / sum(res$task))))
  
  JSD_satiation <- c(JSD_satiation,
                      JSD_dist((res$freely / sum(res$freely)),
                               (res$task / sum(res$task))))
  
  
  cor_satiation <- c(cor_satiation, cor(res$freely, res$task))
  
  cmt <- rbind(res$task_cluster_mat, res$freely_cluster_mat)
  lmt <- rbind(res$task_licking_mat, res$freely_licking_mat)
  structured_unstructured <- data.frame(type=rep(c("Structured", "Unstructed"), times=c(nrow(res$task_cluster_mat),
                                                                                  nrow(res$freely_cluster_mat))))
  
  rownames(structured_unstructured) <- 1:nrow(cmt)
  rownames(cmt) <- 1:nrow(cmt)
  colnames(cmt) <- c()
  cluster_color_label <- spec_cg(len(unique(satiation_metadata[[i]]$cluster_mat$labs)))
  names(cluster_color_label) <- c(-1, 1:(len(unique(satiation_metadata[[i]]$cluster_mat$labs)) - 1))
  
  phlmt <- pheatmap(lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
  phcmt <- pheatmap(cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_row = structured_unstructured, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                    col=cluster_color_label )
  plt <- plot_grid(phcmt[[4]], phlmt[[4]])
  print(cor_satiation)
  
  
  for (size_name in names(sizes)) {
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    pdf(sprintf("%s\\%s\\%s\\no_overlap_satiation_lick_mat.pdf",
                write_path,
                datasets_names[i],
                size_name),
        height=sizes[[size_name]][["width"]] * 3,
        width=sizes[[size_name]][["width"]] * 1.5)
    # unit="in",
    # res=1500)
    
    plot(plt)
    dev.off()
  }
}


AGRP_chunks = list(`1`=5:9,
                   `2`=6:10,
                   `3`=5:9,
                   `4`=3:4)

agrp_paths <- get_agrp_paths()
agrp_metadata <- across_mice_decoding_build_metadata(agrp_paths)

KL_agrp <- c()
JSD_agrp <- c()
cor_agrp <- c()

for (i in 1:len(agrp_metadata)) {
  res = get_task_vs_agrp_lick_distribution(agrp_metadata[[i]], agrp_paths[i], AGRP_chunks[[i]], extra_s = 3)
  
  barplot(res$pre_AGRP / sum(res$pre_AGRP))
  barplot(res$post_AGRP / sum(res$post_AGRP), col=adjustcolor("red", alpha=.5), add=T)
  
  KL_agrp <- c(KL_agrp,
                    KL_dist((res$pre_AGRP / sum(res$pre_AGRP)),
                            (res$post_AGRP / sum(res$post_AGRP))))
  
  JSD_agrp <- c(JSD_agrp,
                     JSD_dist((res$pre_AGRP / sum(res$pre_AGRP)),
                              (res$post_AGRP / sum(res$post_AGRP))))
  
  
  cor_agrp <- c(cor_agrp, cor(res$pre_AGRP, res$post_AGRP))
  print(cor_agrp)
}






get_task_lick_distribution <- function(metadata, path, false_alarms=F, CR=F) {
  
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  if (false_alarms) {
    reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & 
                                                (metadata$stim_master_mat[,"Response"] == 5 | metadata$stim_master_mat[,"Response"] == 3),]  
  } else if (CR) {
    reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & 
                                                (metadata$stim_master_mat[,"Response"] == 4 | metadata$stim_master_mat[,"Response"] == 2),]  
  } else {
    reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]  
  }
  
  
  run_indices <- unique(reward_trials[,ncol(reward_trials)])
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  
  binned_frames_per_mat <- unlist(frames_per_mat) / 15
  binned_frames_to_add <- cumsum(binned_frames_per_mat) - binned_frames_per_mat[1]
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {get_licking_vec(path, run_ind)})
  
  
  names(lick_vectors_all) <- run_indices
  
  
  task_lick_bouts_mat <- c()
  task_licking_mat <- c()
  tt = c()
  idxvc <- c()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    licking_indices_of_interest <- trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))
    
    if (max(licking_indices_of_interest) > len(lick_vectors_all[[as.character(run)]])) {
      print("Exceeding max!")
      next
    }
    
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][licking_indices_of_interest]
    
    
    if( sum(licks_in_trial_unbinned > 0) <= 0 ) {
      print(sprintf("%d. No licks in trial continue", reward_trial_idx))
      next
    }
    
    if (false_alarms) {
      cluster_labels_of_trial <- metadata$trials_mat[as.character(which((metadata$annot_df[,1] == 4 | metadata$annot_df[,1] == 5) & metadata$annot_df[,2] == 1)[reward_trial_idx]),]  
    } else if (CR) {
      cluster_labels_of_trial <- metadata$trials_mat[as.character(which((metadata$annot_df[,1] == 4 | metadata$annot_df[,1] == 5) & metadata$annot_df[,2] == 0)[reward_trial_idx]),]  
    } else {
      cluster_labels_of_trial <- metadata$trials_mat[as.character(which(metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0)[reward_trial_idx]),]
    }
    
    
    
    
    if(all((!(which(licks_in_trial_unbinned > 0)) > 60) | (which(licks_in_trial_unbinned > 0) > 120))) {
      print(sprintf("%d. Only anticipatory licks?", reward_trial_idx))
      
      next
    }
    
    unbinned_first_lick_in_respone <- which(licks_in_trial_unbinned > 0)[which(which(licks_in_trial_unbinned > 0) > 60)[1]]
    binned_first_lick_in_response <- get_binned_index(unbinned_first_lick_in_respone, 15)
    
    if (binned_first_lick_in_response + 3 > 20) {
      binned_first_lick_in_response = 20 - 3
    }
    
    task_lick_bouts_mat <- rbind(task_lick_bouts_mat,
                                 cluster_labels_of_trial[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 8)])
    
    task_licking_mat <- rbind(task_licking_mat,
                              as.numeric(licks_in_trial_unbinned[(unbinned_first_lick_in_respone - 60):(unbinned_first_lick_in_respone + 120)]>0))
    
    tt <- c(tt, relevant_trial["TrialType"])
    
    idxvc <- c(idxvc, reward_trial_idx)
    
  }

  return(list(task_licking_mat=task_licking_mat,
              task_cluster_mat=task_lick_bouts_mat,
              tt=tt))
}
sizes = list(big=c(width=2.5))

paths <- get_thirsty_quenched_paths()
metadata_all <- across_mice_decoding_build_metadata(paths)

write_path <- sprintf("%s\\figure_3\\", figures_base_path)
dir.create(write_path)

write_path <- sprintf("%s\\structure_annotation", write_path)
dir.create(write_path)  
datasets_names = paste("TQ_", get_datasets_names(paths, sep="_"), sep="")

cor_res <- list()
for (i in 1:len(metadata_all)) {
res <- get_task_lick_distribution(metadata_all[[i]], paths[i])

cmt <- res$task_cluster_mat
lmt <- res$task_licking_mat

colnames(cmt) <- c()

cluster_color_label <- spec_cg(len(unique(metadata_all[[i]]$cluster_mat$labs)))
names(cluster_color_label) <- c(-1, 1:(len(unique(metadata_all[[i]]$cluster_mat$labs)) - 1))

colnames(cmt) <- rep("", times=ncol(cmt))
colnames(cmt)[4] <- "Onset"
colnames(lmt) <- rep("", times=ncol(lmt))
colnames(lmt)[61] <- "Onset"

#lick_onset_order <- order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, sum), decreasing=T)
lick_onset_order <-  order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)}), decreasing = F)

phlmt <- pheatmap(lmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
phcmt <- pheatmap(cmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                  col=cluster_color_label )
#plt <- 
  plot_grid(phcmt[[4]], phlmt[[4]])

tbld <- table(cmt[,4])
tbld <- tbld[names(tbld) != -1]
first_clust <- as.numeric(names(which.max(tbld)))
second_clust <- as.numeric(names(sort(tbld, decreasing = T)[2]))


cluster_onset <- apply(cmt, 1, function(trial) {which(trial == first_clust)[1]})
second_cluster_onset <-  apply(cmt, 1, function(trial) {which(trial == second_clust)[1]})
lick_onset <- apply(lmt, 1, function(trial) {which(trial == 1)[1]})
num_licks <- apply(lmt[,1:61], 1, sum)

second_clust_lick_onset <- lick_onset[!is.na(second_cluster_onset)]
second_cluster_onset  <- second_cluster_onset[!is.na(second_cluster_onset)]

lick_onset <- lick_onset[!is.na(cluster_onset)]
nl <- num_licks[!is.na(cluster_onset)]
cluster_onset <- cluster_onset[!is.na(cluster_onset)]

cor_res <- append(cor_res,
                list(list(lo=lick_onset,
                          slo=second_clust_lick_onset,
                   co=cluster_onset,
                   sco=second_cluster_onset,
                   cor=cor(lick_onset, cluster_onset),
                   pv=cor.testcor(lick_onset, cluster_onset),
                   scor=cor(second_clust_lick_onset, second_cluster_onset),
                   nl=cor(nl,cluster_onset),
                   spv=cor.test(second_clust_lick_onset, second_cluster_onset))))

for (size_name in names(sizes)) {
  dir.create(sprintf("%s\\%s",
                     write_path,
                     size_name))
  pdf(sprintf("%s\\%s\\%s\\sorted_lick_mat.pdf",
              write_path,
              datasets_names[i],
              size_name),
      height=sizes[[size_name]][["width"]] * 3,
      width=sizes[[size_name]][["width"]] * 1.5)
  # unit="in",
  # res=1500)
  
  plot(plt)
  dev.off()
}


# all_transitions <- unlist(apply(cmt, 1, 
#                                 function(row) {
#                                   rs <- unlist(lapply(2:len(row), 
#                                                       function(clust_i) {if(row[clust_i-1] != row[clust_i]  & !is.na(row[clust_i])) {paste(row[clust_i-1], row[clust_i], sep="_")}else{NA}}
#                                                       ))
#                                   return(rs[!is.na(rs)])
#                                          }))
# 
# 
# fa_res <-  get_task_lick_distribution(metadata_all[[i]], paths[i], false_alarms = T)
# 
# 
# fa_cmt <- fa_res$task_cluster_mat
# 
# 
# fa_all_transitions <- unlist(apply(fa_cmt, 1, 
#                                 function(row) {
#                                   rs <- unlist(lapply(2:len(row), 
#                                                       function(clust_i) {if(row[clust_i-1] != row[clust_i]  & !is.na(row[clust_i])) {paste(row[clust_i-1], row[clust_i], sep="_")}else{NA}}
#                                   ))
#                                   return(rs[!is.na(rs)])
#                                 }))
}


































for (i in 1:len(am)) {
  
  path <- get_hungry_sated_paths()[i]
  metadata <- metadata_HS[[i]] 

  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  
  

  
  
  reward_mat <- metadata$trials_mat[metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0,]
  
  phr <- pheatmap(reward_mat, cluster_rows=F, cluster_cols=F, legend = F, border_col=NA, show_colnames = F, show_rownames = F)
  reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]
  pt <- plot_structured_licking_bouts(path, reward_trials, frames_per_mat)
  
  plot_grid(phr[[4]], pt)
  
}



