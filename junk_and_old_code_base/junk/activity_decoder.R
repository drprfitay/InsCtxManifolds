
activity_across_mice_decoding_all_trials_decode <- function()  {
  #output_path <- sprintf("%s//single_trial_analysis//", output_path_f)
  write_path <- sprintf("%s\\figure_5\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\SFO_across_mice_decoding", write_path)
  dir.create(write_path)  
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm")
  
  sizes = list(big=c(width=2.5,
                     height=2.5),
               big_1=c(width=2.25,
                       height=2.25),
               medium=c(width=2,
                        height=2),
               medium_1=c(width=1.85,
                          height=1.85),
               medium_2=c(width=2,
                          height=2),
               small=c(width=1.75,
                       height=1.755))
  
  SFO_paths <- get_SFO_paths()
  path_list <- c(rep(SFO_paths[3], times=4), rep(SFO_paths[4], times=4))
  # paths_old <- get_thirsty_quenched_paths()
  chunks <- list(c(1,2,4), c(2:4), c(1:3), c(1,3), 2:4,
                 c(1,2,5), c(1,2,4), c(1:3), 2:4, c(2,3,5))
  
  path_list <-  rep(SFO_paths, each=2)
  chunks <- list(c(1,3), 4:5, c(1,4), 5:6, c(1:3), c(6:8), c(1:3), c(7:9))
  
  metadata_SFO <- across_mice_decoding_build_metadata(path_list=path_list, chunk_list=chunks)
  metadata_thirst <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
  
  # mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  # mice_names <- unlist(lapply(1:len(mice_name_indices), 
  #                             function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  # days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  # days <- unlist(lapply(1:len(days_indices),
  #                       function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  # datasets_names <- paste(mice_names, days, sep = " ")
  # 
  # #metadata <- across_mice_decoding_build_metadata()
  # 
  # confusion_matrix <- matrix(rep(1, times=len(path_list) * len(paths_old)),
  #                            nrow=len(path_list))
  # shuffle_confusion_matrix <- matrix(rep(1, times=len(path_list) * len(paths_old)),
  #                                    nrow=len(path_list))
  # 
  
  cos_dist_vec <- function(a,b) {cosine(a,b)[1,1]}
  shuffle <- function(vec) {sample(vec, len(vec))}
  confusion_matrix_list <- list()
  shuffle_matrix_list <- list()
  pvalues_all <- list()
  
  manhattan_dist <- function(a, b){
    dist <- abs(a-b)
    dist <- sum(dist)
    return(dist)
  }
  
  hamming_dist <- function(a,b) {return(sum(a != b))}
  
  Metrics <-  list(cosine=cos_dist_vec)
                   #euc=euc_dist)
  #hamming=hamming_dist,
  #manhattan=manhattan_dist)
  
  Metrics_objective <- list(cosine=which.max,
                            euc=which.min,
                            hamming=which.min,
                            manhattan=which.min)
  
  for (resp_name in Responses) {
    for (metric_name in names(Metrics)) {
      confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
      
      shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
    }
  }
  
  hits_decoding <- list()
  for (metric_name in names(Metrics)) {
    pvalues_all[metric_name] <- c()
    hits_decoding[[metric_name]] <- list()
  }
  
  # cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
  #                                   nrow=len(path_list))
  # shuffle_cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
  #                                           nrow=len(path_list))
  
  window_size= 15
  trial_duration = 20
  
  for (idx_1 in 1:len(metadata_1)) {
    for (idx_2 in 1:len(metadata_2)) {
      
      if (idx_1 == idx_2) {
        next
      }
      
      
      annotated_mice_1 <- metadata_1[[idx_1]]
      annotated_mice_2 <- metadata_2[[idx_2]]
      
      stim_master_1 <- annotated_mice_1$stim_master_mat
      stim_master_2 <- annotated_mice_2$stim_master_mat 
      
      labels_mice_1 <- annotated_mice_1$cluster_mat$labs
      labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
      mask_mice_1 <- annotated_mice_1$mask
      mask_mice_2 <- annotated_mice_2$mask
      
      original_mat_1 <- get_reduced_mat_full_day(get_thirsty_quenched_paths()[idx_1], just_original_mat = T, window_size=15)
      original_mat_2 <- get_reduced_mat_full_day(get_thirsty_quenched_paths()[idx_2], just_original_mat = T, window_size=15)
      
      if (ncol(original_mat_1) > nrow(original_mat_1)) {
        original_mat_1 <- t(original_mat_1)
      }
      
      if (ncol(original_mat_2) > nrow(original_mat_2)) {
        original_mat_2 <- t(original_mat_2)
      }
      mean_activity_1 <- rowMeans(original_mat_1)
      mean_activity_2 <- rowMeans(original_mat_2)
      
      decoding_stats <- list()
      
      for (metric_name in names(Metrics)) {
        decoding_stats[[metric_name]] <- c()
        
      }
      #cosine_decoding_stats <- c()
      
      
      for (shuff_i in c(-1:1)) {
        
        shuffled_stim_mat <- stim_master_1
        #shuffled_stim_mat[,"Response"]
        if (shuff_i != -1) {
          #print("NOT shuff")
          shuffled_stim_mat[,"Frames"] <- shuffle(shuffled_stim_mat[,"Frames"])
        }
        
        trials_mice_1 <- which(shuffled_stim_mat[,"Response"] %in% c(0,1,2,3,4,5))
        trials_labels_mice_1 <- shuffled_stim_mat[trials_mice_1,"Response"]
        trials_mat_mice_1 <- c()
        activity_mat_trials_mice_1  <- c()
        
        for (r_trial in trials_mice_1) {
          trial_binned_index <- as.numeric(get_binned_index(shuffled_stim_mat[r_trial,"Frames"], window_size))
          activity_in_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          activity_in_trial[is.na(activity_in_trial)] <- 0
          activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_in_trial)
        }
        
        #annot_df_mice_1 <- data.frame(Result=paste(annotated_mice_1$annot_df[,1], annotated_mice_1$annot_df[,2]))
        rownames(activity_mat_trials_mice_1) <- 1:nrow(activity_mat_trials_mice_1)
        #rownames(annot_df_mice_1) <- 1:nrow(trials_mat_mice_1)
        
        trials_mice_2 <- which(stim_master_2[,"Response"] %in% c(0,1,2,3,4,5))
        trials_labels_mice_2 <- stim_master_2[trials_mice_2,"Response"]
        trials_mat_mice_2 <- c()
        activity_mat_trials_mice_2  <- c()
        
        for (r_trial in trials_mice_2) {
          trial_binned_index <- as.numeric(get_binned_index(stim_master_2[r_trial,"Frames"], window_size))
          activity_in_trial <- mean_activity_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          activity_in_trial[is.na(activity_in_trial)] <- 0
          activity_mat_trials_mice_2 <- rbind(activity_mat_trials_mice_2, activity_in_trial)
        }
        
        #annot_df_mice_2 <- data.frame(Result=paste(annotated_mice_2$annot_df[,1], annotated_mice_2$annot_df[,2]))
        rownames(activity_mat_trials_mice_1) <- 1:nrow(activity_mat_trials_mice_1)
        #rownames(annot_df_mice_2) <- 1:nrow(trials_mat_mice_2)
        
        
        decoded_vec <- list()
        for (metric_name in names(Metrics)) {
          decoded_vec[[metric_name]] <- c()
        }
        
        for (decoded_trial_idx in 1:nrow(activity_mat_trials_mice_2)) {
          decoded_trial <- activity_mat_trials_mice_2[decoded_trial_idx,]
          
          confidence <-   
            apply(activity_mat_trials_mice_1, 1, 
                  function(trial) {
                    dist <- 
                      lapply(names(Metrics),
                             function(met_name) {
                               return(Metrics[[met_name]](trial, decoded_trial))
                             })
                    return(unlist(dist))
                  })  
          confidence <- t(confidence)
          
          conf_vec <- c()
          
          
          trial_results <- sort(unique(trials_labels_mice_1))
          for (res in trial_results) {
            
            
            conf_vec <- rbind(conf_vec,
                              mean(confidence[trials_labels_mice_1 == res]))
            
          }
          
          colnames(conf_vec) <- names(Metrics)
          
          
          
          for (metric_name in names(Metrics)) {
            decoded_vec[[metric_name]] <- c(decoded_vec[[metric_name]],
                                            trial_results[Metrics_objective[[metric_name]](conf_vec[,metric_name])])
          }
        }
        
        
        
        decoding_statistics <- list()
        for (metric_name in names(Metrics)) { 
          decoding_statistics[[metric_name]] <- c()
        }
        
        
        for (metric_name in names(Metrics)) { 
          for (response_type in c(0,1,2,3,4,5)) {
            response_indices <- trials_labels_mice_2 == response_type
            
            decoding_accuracy = sum(decoded_vec[[metric_name]][response_indices] == trials_labels_mice_2[response_indices]) / 
              len(trials_labels_mice_2[response_indices])
            
            decoding_statistics[[metric_name]] <- 
              c(decoding_statistics[[metric_name]], decoding_accuracy)
            
            
            if (response_type ==0 && shuff_i == -1) {
              #print(list(decoded_vec[[metric_name]][response_indices] == trials_labels_mice_2[response_indices]))
              hits_decoding[[metric_name]] <- append(hits_decoding[[metric_name]],
                                                     list(decoded_vec[[metric_name]][response_indices] == trials_labels_mice_2[response_indices]))
            }
          }
        }
        
        
        
        
        # print(sprintf("Decoding statistics from %d to %d, %f (cosine: %f)", 
        #               idx_1,
        #               idx_2,
        #               decoding_statistics,
        #               cosine_decoding_statistics))
        
        # decoding_stats <- c(decoding_stats,
        #                     decoding_statistics)
        for (metric_name in names(Metrics)) { 
          decoding_stats[[metric_name]] <- rbind(decoding_stats[[metric_name]], 
                                                 decoding_statistics[[metric_name]])
        }
      }
      
      for (metric_name in names(Metrics)) { 
        
        
        colnames(decoding_stats[[metric_name]]) <- Responses
        decoding_statistics <- decoding_stats[[metric_name]][1,]
        shuffled_statistics <- colMeans(decoding_stats[[metric_name]][-1,])
        
        ecdfs <- apply(decoding_stats[[metric_name]][-1,], 2, ecdf)
        pvalues <-  unlist(lapply(1:len(decoding_stats[[metric_name]][1,]), 
                                  function(acc_idx) {
                                    ecdfs[[acc_idx]](decoding_stats[[metric_name]][1,acc_idx])
                                  }))
        
        names(pvalues) <- Responses
        
        pvalues_all[[metric_name]] <- rbind(pvalues_all[[metric_name]],pvalues)
        
        print(sprintf("######################## SHUFFLED DECODING: (Shuffle vs decoder) - %d vs  %d ####", idx_1, idx_2))
        
        for (resp_name in Responses)  {
          
          mat_name <- sprintf("%s_%s", resp_name, metric_name)
          confusion_matrix_list[[mat_name]][idx_1, idx_2] <- decoding_statistics[resp_name]
          shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- shuffled_statistics[resp_name]
          #print(sprintf("%.3f vs %.3f", shuffled_decoding_statistics, decoding_statistics))
          print(sprintf("(%d). %s: %.3f vs %.3f [Pvalue: %f] -> %s", 
                        nrow(pvalues_all[[metric_name]]),
                        resp_name, 
                        shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                        confusion_matrix_list[[mat_name]][idx_1, idx_2],
                        pvalues[[resp_name]],
                        metric_name))
          
          
        }
        
        print("#######################")
        
      }
      # 
      # confusion_matrix[idx_1, idx_2] <- decoding_statistics
      # shuffle_confusion_matrix[idx_1, idx_2] <- shuffled_decoding_statistics
      # cosine_confusion_matrix[idx_1, idx_2] <- cosine_decoding_statistics
      # shuffle_cosine_confusion_matrix[idx_1, idx_2] <- shuffled_cosine
      # miss_confusion_matrix[idx_1, idx_2] <- miss_decoding_statistics
      # hit_confusion_matrix[idx_1, idx_2] <- hit_decoding_statistics
      
      
      # activity_decoded_vec <- c()
      # 
      # 
      # for (decoded_trial_idx in 1:nrow(activity_mat_trials_mice_2)) {
      #   decoded_trial <- activity_mat_trials_mice_2[decoded_trial_idx,]
      #   confidence <- 
      #     apply(activity_mat_trials_mice_1, 1, 
      #           function(trial) {cor(trial, decoded_trial)})  
      #   
      #   
      #   is_hit <- mean(confidence[reward_trials_labels_mice_1 == 0]) 
      #   is_miss <- mean(confidence[reward_trials_labels_mice_1 == 1])
      #   activity_decoded_vec <- c(activity_decoded_vec, ifelse(is_hit > is_miss, 0, 1))
      # }
      # 
      # activity_decoding_statistics <-  sum(activity_decoded_vec == reward_trials_labels_mice_2) / len(reward_trials_labels_mice_2)
      # activity_hit_decoding_statistics <- sum(activity_decoded_vec[reward_trials_labels_mice_2 == 0] == 0) / sum(reward_trials_labels_mice_2 == 0)
      # activity_miss_decoding_statistics <- sum(activity_decoded_vec[reward_trials_labels_mice_2 == 1] == 1) / sum(reward_trials_labels_mice_2 == 1)
      # 
      # print(sprintf("Activity decoding statistics from %d to %d, %f (hits only: %f) (misses only :%f)", 
      #               idx_1,
      #               idx_2,
      #               activity_decoding_statistics,
      #               activity_hit_decoding_statistics,
      #               activity_miss_decoding_statistics))
      # 
      # activity_confusion_matrix[idx_1, idx_2] <- activity_decoding_statistics
      # activity_miss_confusion_matrix[idx_1, idx_2] <- activity_miss_decoding_statistics
      # activity_hit_confusion_matrix[idx_1, idx_2] <- activity_hit_decoding_statistics
      
    }
  }
  
  # rownames(confusion_matrix) <- datasets_names
  # colnames(confusion_matrix) <- datasets_names
  
  dir.create(sprintf("%s\\data\\figure_3\\", base_output_path))
  
  save(file=sprintf("%s\\data\\figure_3\\activity_confusion_matrices.Rda", base_output_path), confusion_matrix_list)
}