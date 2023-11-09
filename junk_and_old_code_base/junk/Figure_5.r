figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"


metadata_14 <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_14_shuffle <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths(), preset="dalshtaimshuff")
metadata_SFO <- across_mice_decoding_build_metadata(path_list=rep(get_SFO_paths(), each=2),
                                                    chunk_list=list(c(1,3), 4:5, c(1,4), 5:6, c(1:3), c(6:8), c(1:3), c(7:9)))

chunks <- list(c(1,3), 4:5, c(1,4), 5:6, c(1:3), c(6:8), c(1:3), c(7:9))

across_mice_decoding_all_trials_decode <- function()  {
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
  
  Metrics <-  list(cosine=cos_dist_vec,
                   euc=euc_dist,
                   hamming=hamming_dist,
                   manhattan=manhattan_dist)
  
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
      
      decoding_stats <- list()
      
      for (metric_name in names(Metrics)) {
        decoding_stats[[metric_name]] <- c()
    
      }
      #cosine_decoding_stats <- c()
      
      
      for (shuff_i in c(-1:200)) {
        
        shuffled_stim_mat <- stim_master_1
        #shuffled_stim_mat[,"Response"]
        if (shuff_i != -1) {
          print("NOT shuff")
          shuffled_stim_mat[,"Frames"] <- shuffle(shuffled_stim_mat[,"Frames"])
        }
        
        trials_mice_1 <- which(shuffled_stim_mat[,"Response"] %in% c(0,1,2,3,4,5))
        trials_labels_mice_1 <- shuffled_stim_mat[trials_mice_1,"Response"]
        trials_mat_mice_1 <- c()
        activity_mat_trials_mice_1  <- c()
        
        for (r_trial in trials_mice_1) {
          trial_binned_index <- as.numeric(get_binned_index(shuffled_stim_mat[r_trial,"Frames"], window_size))
          
          clustered_trial <- labels_mice_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
          }
          
          trials_mat_mice_1 <- rbind(trials_mat_mice_1, clustered_trial)
          
          
          # activity_in_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          # activity_in_trial[is.na(activity_in_trial)] <- 0
          # 
          # activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_in_trial)
          
          
        }
        
        #annot_df_mice_1 <- data.frame(Result=paste(annotated_mice_1$annot_df[,1], annotated_mice_1$annot_df[,2]))
        rownames(trials_mat_mice_1) <- 1:nrow(trials_mat_mice_1)
        #rownames(annot_df_mice_1) <- 1:nrow(trials_mat_mice_1)
        
        
        
        trials_mice_2 <- which(stim_master_2[,"Response"] %in% c(0,1,2,3,4,5))
        trials_labels_mice_2 <- stim_master_2[trials_mice_2,"Response"]
        trials_mat_mice_2 <- c()
        activity_mat_trials_mice_2  <- c()
        
        for (r_trial in trials_mice_2) {
          trial_binned_index <- as.numeric(get_binned_index(stim_master_2[r_trial,"Frames"], window_size))
          
          
          clustered_trial <- labels_mice_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
          }
          
          trials_mat_mice_2 <- rbind(trials_mat_mice_2,clustered_trial)
          
          # 
          # activity_in_trial <- mean_activity_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          # activity_in_trial[is.na(activity_in_trial)] <- 0
          # activity_mat_trials_mice_2 <- rbind(activity_mat_trials_mice_2, activity_in_trial)
        }
        
        #annot_df_mice_2 <- data.frame(Result=paste(annotated_mice_2$annot_df[,1], annotated_mice_2$annot_df[,2]))
        rownames(trials_mat_mice_2) <- 1:nrow(trials_mat_mice_2)
        #rownames(annot_df_mice_2) <- 1:nrow(trials_mat_mice_2)
        
        
        
        translated_mat_mice_1_mice_2 <- c()
        
        
        for (trial_idx in 1:nrow(trials_mat_mice_1)) {
          trial_before <- trials_mat_mice_1[trial_idx,]
          trial_after <- trial_before
          
          for (cluster_idx in 1:len(mask_mice_1)) {
            trial_after[which(trial_before == mask_mice_1[cluster_idx])] <- mask_mice_2[cluster_idx] 
          }
          
          translated_mat_mice_1_mice_2 <- rbind(translated_mat_mice_1_mice_2,
                                                as.numeric(trial_after))
          
        }
        
        
        decoded_vec <- list()
        for (metric_name in names(Metrics)) {
          decoded_vec[[metric_name]] <- c()
        }
        
        for (decoded_trial_idx in 1:nrow(trials_mat_mice_2)) {
          decoded_trial <- trials_mat_mice_2[decoded_trial_idx,]
          
          confidence <-   
            apply(translated_mat_mice_1_mice_2, 1, 
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
                              colMeans(confidence[trials_labels_mice_1 == res,]))
            
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
  
  save(file=sprintf("%s\\confusion_matrices.Rda", write_path), confusion_matrix_list)
  save(file=sprintf("%s\\shuffle_confusion_matrices.Rda", write_path), shuffle_confusion_matrix)
  save(file=sprintf("%s\\pvalues_decoding.Rda", write_path), pvalues_all)
  
  
  behavior_table <- lapply(metadata,
                           function(m) { 
                             error_types <- as.character(c(0,1,2,3,4,5))
                             cont <- rep(0, times=len(error_types))
                             names(cont) <- error_types
                             tb <- table(m$stim_master_mat[,8])
                             cont[names(tb)] <- tb
                             cont <- cont[error_types]
                             return(cont / sum(cont))
                           })
  bmat <- do.call(rbind, behavior_table)
  bmat_reg <- bmat[1:4 * 2 - 1,]
  bmat_shuff <- bmat[1:4 * 2,]
  colnames(bmat_shuff) <- Responses
  colnames(bmat_reg) <- Responses
  melted_reg <- melt(bmat_reg)
  melted_shuff <- melt(bmat_shuff)
  melted_reg <- cbind(melted_reg, rep("Regular", times=24))
  melted_shuff <- cbind(melted_shuff, rep("Shuffle", times=24))
  colnames(melted_reg) <- c("#", "Variable", "Fraction", "Group")
  colnames(melted_shuff) <- c("#", "Variable", "Fraction", "Group")
  rbind(melted_shuff, melted_reg)
  melted_all <- rbind(melted_shuff, melted_reg)
  #gbehav <- 
    ggplot(melted_all) + 
    geom_bar(aes(x=Variable, 
                 fill=Group, 
                 y=Fraction, group=Group), stat="summary", position=position_dodge()) + 
              ylim(0,0.4) + 
      geom_errorbar(aes(x=Variable, y=Fraction, group=Group), 
                    stat="summary", position=position_dodge()) +
        theme_classic() +
        theme(text=element_text(size=14, color="black"), 
              axis.text=element_text(size=14, color="black"),
              axis.ticks = element_line(color="black"),
              plot.title=element_text(size=10),
              legend.position="NA") +
      scale_y_continuous(expand=c(0,0))
    
  for (metric_name in names(Metrics)) {
    plot_pvals = T  
    matrices_names <- names(confusion_matrix_list)[grepl(metric_name, names(confusion_matrix_list))]
    for (trial_response in matrices_names) {
      
      
      
      confusion <- confusion_matrix_list[[trial_response]]
      shuffle  <- shuffle_matrix_list[[trial_response]]
      
      results_df <- data.frame()
      within <- c()
      across <- c()
      within_shuffle <- c()
      across_shuffle <- c()
      
      pvalues_across <- c()
      pvalues_within <- c()
      pvalues_comp_idx <- 1
      
      for (idx1 in 1:len(paths)) {
        for (idx2 in 1:len(paths)) {
          
          if (idx1 == idx2) { 
            next
          }
          
          if (mice_names[idx1] == mice_names[idx2]) {
            print("Within")
            within <- c(within, confusion[idx1, idx2])
            within_shuffle <- c(within_shuffle, shuffle[idx1, idx2])
            results_df <- rbind(results_df, 
                                data.frame(group="Within", decoding=confusion[idx1, idx2]))
            
            results_df <- rbind(results_df, 
                                data.frame(group="Within (shuffle)", decoding=shuffle[idx1, idx2]))
            
            if (plot_pvals) {
              pvalues_within <- rbind(pvalues_within,
                                      pvalues_all[[metric_name]][pvalues_comp_idx,])
            }
          } else {
            print("Across")
            across <- c(across, confusion[idx1, idx2])
            across_shuffle <- c(across_shuffle, shuffle[idx1, idx2])
            results_df <- rbind(results_df, 
                                data.frame(group="Across", decoding=confusion[idx1, idx2]))
            
            results_df <- rbind(results_df, 
                                data.frame(group="Across (shuffle)", decoding=shuffle[idx1, idx2]))
            
            if (plot_pvals) {
              pvalues_across <- rbind(pvalues_across,
                                      pvalues_all[[metric_name]][pvalues_comp_idx,])
            }
          }
          
          pvalues_comp_idx <- pvalues_comp_idx + 1
        }
      }
      
      if (plot_pvals) {
        colnames(pvalues_across) <- Responses 
        colnames(pvalues_within) <- Responses 
        
        pvalues_of_pvalues_across <- lapply(2:ncol(pvalues_across),
                                            function(i) 
                                            {wilcox.test(pvalues_across[,1],
                                                         pvalues_across[,i],
                                                         paired=T, alternative="greater")$p.value})
        
        pvalues_of_pvalues_within <- lapply(2:ncol(pvalues_within),
                                            function(i) 
                                            {wilcox.test(pvalues_within[,1],
                                                         pvalues_within[,i],
                                                         paired=T, alternative="greater")$p.value})
        
        pvalues_of_pvalues_across <- unlist(pvalues_of_pvalues_across) * (ncol(pvalues_across) - 1)
        pvalues_of_pvalues_within <- unlist(pvalues_of_pvalues_within) * (ncol(pvalues_across) - 1)
        
        melted_df_across <- melt(pvalues_across)
        colnames(melted_df_across) <- c("#", "Variable", "Pvalue")
        
        melted_df_within <- melt(pvalues_within)
        colnames(melted_df_within) <- c("#", "Variable", "Pvalue")
        
        gmelted_within <- 
          ggplot(melted_df_within, aes(y=Pvalue, x=Variable)) +
          geom_jitter() +
          #geom_violin(fill="gray50", color=NA) +
          geom_boxplot(width=0.5, 
                       size=1, 
                       color=adjustcolor("gray65", alpha=1), 
                       fill=adjustcolor("gray80", alpha=0.9)) + 
          #stat_summary() + 
          
          #scale_y_continuous(expand=c(0,0)) +
          # geom_bar(stat="summary", width=.45,
          #          fill="gray70") +
          # geom_errorbar(stat="summary", width = .3) +
          theme_classic() +
          theme(text=element_text(size=14, color="black"), 
                axis.text=element_text(size=14, color="black"),
                axis.ticks = element_line(color="black"),
                plot.title=element_text(size=10),
                legend.position="NA")  + 
          ggtitle("Within") + 
          xlab("") +
          ylab("Shuffles defeated (%)")
        
        gmelted_across <- 
          ggplot(melted_df_across, aes(y=Pvalue, x=Variable)) +
          geom_jitter() +
          # geom_violin(fill="gray50", color=NA) +
          # stat_summary() +
          geom_boxplot(width=0.5, 
                       size=1, 
                       color=adjustcolor("gray65", alpha=1), 
                       fill=adjustcolor("gray80", alpha=0.9)) + 
          #scale_y_continuous(expand=c(0,0)) +
          # geom_bar(stat="summary", width=.45,
          #          fill="gray70") +
          # geom_errorbar(stat="summary", width = .3) +
          theme_classic() +
          theme(text=element_text(size=14, color="black"), 
                axis.text=element_text(size=14, color="black"),
                axis.ticks = element_line(color="black"),
                plot.title=element_text(size=10),
                legend.position="NA")  + 
          ggtitle("Across") +
          xlab("") +
          ylab("Shuffles defeated (%)")
      }
      
      
      wilx_across <- wilcox.test(across, across_shuffle, paired=T, alternative="greater")
      wilx_within <- wilcox.test(within, within_shuffle, paired=T, alternative="greater")
      wilx_across_within <- wilcox.test(across, within, paired=F)
      
      ylim_max <- max(results_df$decoding)
      gbox <- 
        ggplot(results_df, aes(x=group,y=decoding)) +
        geom_hline(yintercept=1/6, lty="dashed", size=.5) + 
        geom_jitter(aes(fill=group),
                    position=position_jitterdodge(.5),
                    alpha=0.8,
                    stroke=.5,
                    color=adjustcolor("gray20", alpha=.2),
                    size=1.2) + 
        geom_boxplot(width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9)) + 
        
        theme_classic() +
        xlab("") +
        ylab("Correctly classified (%)") +
        ggtitle(trial_response) + 
        theme(text=element_text(size=14, color="black"),
              plot.title = element_text(size=9, color="black"),
              legend.position="NA",
              legend.title = element_blank(),
              axis.text = element_text(color="black"),
              axis.text.x = element_text(size=10),#element_text(angle = 45, vjust=.6),
              axis.ticks = element_line(color="black")) + 
        
        geom_line(data=data.frame(x=c(1,1,2,2),
                                  y=c(ylim_max + 0.05,
                                      ylim_max + 0.075,
                                      ylim_max + 0.075,
                                      ylim_max + 0.05)),
                  aes(x=x,y=y)) +
        geom_text(x=1.5, y=ylim_max + 0.075, label=signif.num(wilx_across$p.value)) + 
        
        geom_line(data=data.frame(x=c(3,3,4,4),
                                  y=c(ylim_max + 0.05,
                                      ylim_max + 0.075,
                                      ylim_max + 0.075,
                                      ylim_max + 0.05)),
                  aes(x=x,y=y)) +
        geom_text(x=3.5, y=ylim_max + 0.075, label=signif.num(wilx_within$p.value)) +
        
        
        geom_line(data=data.frame(x=c(1,1,3,3),
                                  y=c(ylim_max + 0.125,
                                      ylim_max + 0.15,
                                      ylim_max + 0.15,
                                      ylim_max + 0.125)),
                  aes(x=x,y=y)) +
        geom_text(x=2, y=ylim_max  + 0.15 , label=signif.num(wilx_across_within$p.value)) + 
        ylim(c(0,1.075))
      
      
      for (size_name in names(sizes)) {
        
        dir.create(sprintf("%s\\%s",write_path, size_name))
        dir.create(sprintf("%s\\%s\\%s",write_path, size_name, metric_name))
        
        pdf(sprintf("%s\\%s\\%s\\%s_across_mice_decoding.pdf",
                    write_path,
                    size_name,
                    metric_name,
                    trial_response),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]] * 1.4)
        
        plot(gbox)
        dev.off()
        
        if (plot_pvals) 
        { 
          pdf(sprintf("%s\\%s\\%s\\boxplot_across_pvalues.pdf",
                      write_path,
                      size_name,
                      metric_name),
              height=sizes[[size_name]][["width"]],
              width=sizes[[size_name]][["width"]] * 2)
          
          plot(gmelted_across)
          dev.off()
          
          pdf(sprintf("%s\\%s\\%s\\boxplot_within_pvalues.pdf",
                      write_path,
                      size_name,
                      metric_name),
              height=sizes[[size_name]][["width"]],
              width=sizes[[size_name]][["width"]] * 2)
          
          plot(gmelted_within)
          dev.off()
          
          print("Plotting pvalues!")
          if (size_name == "small")
            plot_pvals =F
        }
      }  
    }
  }
  
  
  colnames(pvalues_all) <- Responses
  melted_df <- melt(pvalues_all)
  colnames(melted_df) <- c("#", "Variable", "Pvalue")
  gmelted <- 
    ggplot(melted_df, aes(y=Pvalue, x=Variable)) +
    geom_violin(fill="gray50", color=NA) +
    stat_summary() + 
    #scale_y_continuous(expand=c(0,0)) +
    # geom_bar(stat="summary", width=.45,
    #          fill="gray70") +
    # geom_errorbar(stat="summary", width = .3) +
    theme_classic() +
    theme(text=element_text(size=14, color="black"), 
          axis.text=element_text(size=14, color="black"),
          axis.ticks = element_line(color="black"),
          plot.title=element_text(size=10),
          legend.position="NA")  + 
    xlab("") +
    ylab("Shuffles defeated (%)")
}


figure_5_structure_similarity_plots <- function() 
{
  
  write_path <- sprintf("%s\\figure_2\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\structure_similarity_plots", write_path)
  dir.create(write_path)  
  
  
  sizes = list(big=c(width=3.5),
               medium=c(width=3),
               medium_3=c(width=2.5),
               medium_2=c(width=2.25),
               medium_1=c(width=2),
               small=c(width=1.75))
  
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  label_datasets_names <- get_datasets_names(paths, sep = " ")
  
  metric = cosine
  
  m1 <- metadata
  m2 <- metadata_old
  
  nclusters=100
  chunk=-1
  nreps=10
  #compare_shuffle = F
  combination_matrix <- matrix(rep(0, times=len(m1) * len(m2)),
                               nrow=len(m1))
  results_all <- list()
  
  for (compare_shuffle in c(F)) {
    
    metric_functions <- list(KL=KL_dist,
                             JSD=JSD_dist,
                             wasserstein=wasserstein_dist,
                             corr=cor,
                             spearman=spearman_cor,
                             euc=euc_dist,
                             cosine=cosine)
    
    
    structure_sim_list <- list()
    
    for (metric_name in names(metric_functions)) {
      
      structure_sim_list[[metric_name]] <- c()
      
    }
    
    for (idx1 in 1:len(m1)) {
      for (idx2 in 1:len(m2)) {
        
        day1_mat <- m1[[idx1]]$red_mat
        day2_mat <- m2[[idx2]]$red_mat
        
        
        
        km_day1 <- kmeans(apply(day1_mat, 2, scale), nclusters, iter.max=300)
        km_day2 <- kmeans(apply(day2_mat, 2, scale), nclusters, iter.max=300)
        
        struct_corr_pair <- list()
        
        for (metric_name in names(metric_functions)) {
          
          struct_corr_pair[[metric_name]] <- c()
          
        }
        
        
        for (i in 1:nreps) {
          
          
          
          dist_matrix_day1 <- as.matrix(dist(km_day1$centers))
          dist_matrix_day2 <- as.matrix(dist(km_day2$centers))
          
          # h1 <- heatmap(dist_matrix_day1, plot=F)
          # h2 <- heatmap(dist_matrix_day2, plot=F)
          
          h1 <- hclust(dist(dist_matrix_day1), method="single")
          h2 <- hclust(dist(dist_matrix_day2), method="single")
          
          
          
          for (metric_name in names(metric_functions)) {
            
            metric <- metric_functions[[metric_name]]
            struct_corr <- metric(c(dist_matrix_day1[h1$order, h1$order]), 
                                  c(dist_matrix_day2[h2$order, h2$order]))
            
            
            struct_corr_pair[[metric_name]] <- c(struct_corr_pair[[metric_name]], 
                                                 struct_corr)  
          }
          
          
        }
        
        for (metric_name in names(metric_functions)) {
          
          
          print(sprintf("METRIC ################## %s", metric_name))
          print(struct_corr_pair[[metric_name]])
          
          structure_sim_list[[metric_name]] <- 
            rbind(structure_sim_list[[metric_name]],
                  struct_corr_pair[[metric_name]])  
          
          combination_matrix[idx1, idx2] <- mean(struct_corr_pair[["JSD"]])
        }
      }
    }
    
    
    name = ifelse(compare_shuffle, "Shuffle", "Regular")
    
    results_all[[name]] <- structure_sim_list
  }
  
  
  alternatives = list(KL="less",
                      JSD="less",
                      wasserstein="less",
                      corr="greater",
                      spearman="greater",
                      euc="less",
                      cosine="greater")
  
  self_comp_values = list(KL=0,
                          JSD=0,
                          wasserstein=0,
                          corr=1,
                          spearman=1,
                          euc=0,
                          cosine=1)
  
  mice_names <- unlist(lapply(str_split(datasets_names, "_"), function(ls) {ls[[1]]}))
  
  within_across_indices <- apply(combination_matrix, 2, function(c) {mice_names[c[1]] == mice_names[c[2]]})
  within_across_names <- rep("Across", times=len(within_across_indices))
  within_across_names[within_across_indices] <- "Within"
  
  regular_labels <- rep("Regular", times=nrow(results_all$Regular[[1]]))
  shuffle_labels <- rep("Shuffle", times=nrow(results_all$Shuffle[[1]]))
  
  regular_within_across_labels <- paste(regular_labels, within_across_names, sep=" - ")
  shuffle_within_across_labels <- paste(shuffle_labels, within_across_names, sep=" - ")
  
  annotation_df <- data.frame(Mice=mice_names)
  rownames(annotation_df) = label_datasets_names
  
  for (metric_name in names(metric_functions)) {
    plot_df <- data.frame(Similarity=rowMeans(results_all$Regular[[metric_name]]),
                          Group=rep(rep(c("Regular", "SFO"), times=4), each=14))
    
    ylim_max = max(plot_df$Similarity)
    ylim_min = min(plot_df$Similarity)
    
    wilx = wilcox.test(plot_df[plot_df$Group == "Regular","Similarity"],
                       plot_df[plot_df$Group != "Regular","Similarity"],
                       paired=T,
                       alternative=alternatives[[metric_name]])

    g <- 
      ggplot(plot_df) +
      geom_jitter(aes(y=Similarity, x=Group, group=Group, fill=Group),
                  position=position_jitterdodge(.5),
                  color="gray20", alpha=0.8) + 
      geom_boxplot(aes(y=Similarity, x=Group, group=Group),
                   width=0.5, 
                   size=1, 
                   color=adjustcolor("gray65", alpha=1), 
                   fill=adjustcolor("gray80", alpha=0.9)) +
      theme_classic() +
      theme(text=element_text(size=14, color="black"), 
            axis.text=element_text(size=14, color="black"),
            axis.ticks = element_line(color="black"),
            plot.title=element_text(size=10),
            legend.position="NA") +
      geom_line(data=data.frame(x=c(1,1,2,2),
                                y=c(ylim_max * 1.005,
                                    ylim_max * 1.007,
                                    ylim_max * 1.007,
                                    ylim_max * 1.005)),
                aes(x=x,y=y)) +
      geom_text(x=1.5, y=ylim_max * 1.008, label=signif.num(wilx$p.value)) +
      ggtitle(metric_name) + 
      xlab("") + 
      ylab("Structure similarity")
    
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\all",
                         write_path))
      dir.create(sprintf("%s\\all\\%s",
                         write_path,
                         size_name))
      pdf(sprintf("%s\\all\\%s\\structure_similarity_%s.pdf",
                  write_path,
                  size_name,
                  metric_name),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]]) 
      
      plot(g)
      dev.off()
    }
    
   
  }
}

figure_5_pairwise_topological_permutation_tests <- function()
{
  write_path <- sprintf("%s\\figure_2\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\topological_similarity_plots", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=3.5),
               medium=c(width=3),
               medium_3=c(width=2.5),
               medium_2=c(width=2.25),
               medium_1=c(width=2),
               small=c(width=1.75))
  

  
  shuffled_pariwise <- matrix(rep(0, times=len(m1) * len(m2)), nrow=len(m2))
  pairwise <- matrix(rep(0, times=len(m1) * len(m2)), nrow=len(m2))
  
  id_perm_test = F
  metric_functions <- list(KL=KL_dist,
                           JSD=JSD_dist,
                           wasserstein=wasserstein_dist)
  
  
    pvalue_matrices_betti_0 <- list()
    pvalue_matrices_betti_1 <- list()
    
    
    for (metric_name in names(metric_functions)) {
      pvalue_matrices_betti_0[[metric_name]] <- matrix(rep(0, times=len(m1) * len(m2)), nrow=len(m1))
      pvalue_matrices_betti_1[[metric_name]] <- matrix(rep(0, times=len(m1) * len(m2)), nrow=len(m1))

    }
  
  remove_main_diagonal=T
  
  

  
  comparision_list <- list(
                           before_SFO=list(m1=metadata_SFO[1:4 * 2 -1],
                                           m2=metadata_14,
                                           diag=F),
                           after_SFO=list(m1=metadata_SFO[1:4 * 2],
                                          m2=metadata_14,
                                          diag=F),
                           within=list(m1=metadata_14,
                                       m2=metadata_14,
                                       diag=T),
                           shuffle=list(m1=metadata_14,
                                       m2=metadata_14_shuffle,
                                       diag=T))
  
  
  
                          
  all_results_all <- list()
  
  for (comp_name in names(comparision_list)) {
    comp <- comparision_list[[comp_name]]
    m1 <- comp$m1
    m2 <- comp$m2
    remove_main_diagonal <- comp$diag
    nreps = 20
    results_all <- data.frame()
  
  for (idx1 in 1:len(m1)) {
    for (idx2 in 1:len(m2))  { 
    
      
      if (remove_main_diagonal) {
        if (idx1 == idx2) {
          next
        }
      }
    

    orig_reg_mat_1 <- m1[[idx1]]$red_mat
    orig_reg_mat_2 <- m2[[idx2]]$red_mat


    reg_mat_1 <- kmeans(orig_reg_mat_1, centers=80, iter.max=500)$centers
    reg_mat_2 <- kmeans(orig_reg_mat_2, centers=80, iter.max=500)$centers
    
    mtall <- rbind(reg_mat_1,
                   reg_mat_2)

    
   
      phomology_similarity <- list()

      permutation_phomology_similarity <- list()

      
      phom_reg_mat_1 <- calculate_homology(reg_mat_1)
      phom_reg_mat_2 <- calculate_homology(reg_mat_2)

      
      # Barcode 0 
      life_death_0_mat_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 0,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 0,2]
      life_death_0_mat_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 0,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 0,2]

      
      # Barcode 1
      life_death_1_mat_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 1,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 1, 2]
      life_death_1_mat_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 1, 3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 1, 2]

      
      
      
      
      for (metric_name in names(metric_functions)) {
        
        phomology_similarity[[metric_name]] <- list()

        
        permutation_phomology_similarity[[metric_name]] <- list()
        permutation_phomology_similarity[[metric_name]][["betti_0"]] <- c()
        permutation_phomology_similarity[[metric_name]][["betti_1"]] <- c()

        
        phomology_similarity[[metric_name]][["betti_0"]] <- metric_functions[[metric_name]](life_death_0_mat_1,
                                                                                            life_death_0_mat_2)
        phomology_similarity[[metric_name]][["betti_1"]] <- metric_functions[[metric_name]](life_death_1_mat_1,
                                                                                            life_death_1_mat_2)
        

        
        print(sprintf("Dataset (%d) vs (%d) b0: %.3f, b1: %.3f [%s]",
                      idx1,
                      idx2,
                      phomology_similarity[[metric_name]][["betti_0"]],
                      phomology_similarity[[metric_name]][["betti_1"]],
                      metric_name))
        

        print("------")
        print("------")
      }
    
    
    for (i in 1:1) {
      
      logical_ind <- 1:nrow(mtall) %in% sample(1:nrow(mtall), nrow(reg_mat_1))
      ind1 <- which(logical_ind)
      ind2 <- which(!logical_ind)
      

        perm_phom_reg_mat_1 <- calculate_homology(mtall[ind1,])
        perm_phom_reg_mat_2 <- calculate_homology(mtall[ind2,])

        # Barcode 0 
        perm_life_death_0_mat_1 <- 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 0,3] - 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 0,2]
        
        perm_life_death_0_mat_2 <- 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 0,3] - 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 0,2]
        

        # Barcode 1
        perm_life_death_1_mat_1 <- 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 1,3] - 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 1,2]
        
        perm_life_death_1_mat_2 <- 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 1,3] - 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 1,2]
        

        
        for (metric_name in names(metric_functions)) {
          
          permutation_phomology_similarity[[metric_name]][["betti_0"]] <- 
            c(permutation_phomology_similarity[[metric_name]][["betti_0"]],
              metric_functions[[metric_name]](perm_life_death_0_mat_1, perm_life_death_0_mat_2))
          
          permutation_phomology_similarity[[metric_name]][["betti_1"]] <- 
            c(permutation_phomology_similarity[[metric_name]][["betti_1"]],
              metric_functions[[metric_name]](perm_life_death_1_mat_1,perm_life_death_1_mat_2))
        }
    }
    
    
    for (metric_name in names(metric_functions)) {
      results_all <- rbind(results_all,
                           data.frame(path1=idx1, 
                                      path2=idx2, 
                                      betti0=phomology_similarity[[metric_name]]$betti_0,
                                      betti1=phomology_similarity[[metric_name]]$betti_1,
                                      isShuffle="Regular",
                                      metric=metric_name))
      
      
      pval_func_betti_0 <- ecdf(permutation_phomology_similarity[[metric_name]]$betti_0)
      pval_func_betti_1 <- ecdf(permutation_phomology_similarity[[metric_name]]$betti_1)

      
      
      print(sprintf("Pvalues b0 (%.3f) b1(%.3f) %s",
                    1 - pval_func_betti_0(phomology_similarity[[metric_name]]$betti_0),
                    1 - pval_func_betti_1(phomology_similarity[[metric_name]]$betti_1),

                    metric_name))
      
      pvalue_matrices_betti_0[[metric_name]][idx1, idx2] <- 1 - pval_func_betti_0(phomology_similarity[[metric_name]]$betti_0)
      pvalue_matrices_betti_1[[metric_name]][idx1, idx2] <- 1 - pval_func_betti_1(phomology_similarity[[metric_name]]$betti_1)

    }
    
    pheatmap(pvalue_matrices_betti_0[["JSD"]], breaks=c(0,0.05,0.1,0.15,1), col=rdylbu_cg(4))
    
    }
    
  }

    for (barcode_name in names(barcodes)){
      barcode =  barcodes[[barcode_name]]
      
      for (metric_name in names(metric_functions)) {
        regular_df <- results_all[results_all$metric == metric_name,]
        
        
        tmp_df <- 
          ddply(regular_df, .(path2), function(path_df) {metric_df = path_df[path_df$metric == metric_name,]; return(mean(metric_df[[barcode]]))})
        
        regular_labels <- rep(comp_name, nrow(tmp_df))
        plot_df <- data.frame(phom=tmp_df[,"V1"],
                              all_group=c(regular_labels))
        
        all_results_all[[paste(barcode_name, metric_name, sep="_")]] <- 
          rbind(all_results_all[[paste(barcode_name, metric_name, sep="_")]],
                plot_df)
      }
    }
  }

  
  for (barcode_name in names(barcodes)){
    barcode =  barcodes[[barcode_name]]
    
    for (metric_name in names(metric_functions)) {
      
      plot_df <- all_results_all[[paste(barcode_name, metric_name, sep="_")]]
      
      plot_df$all_group <- factor(plot_df$all_group,
                                  levels=c("within",
                                           "before_SFO",
                                           "after_SFO",
                                           "shuffle"))
      
      #plot_df <- plot_df[plot_df$all_group != "shuffle",]
      ylim_max = max(plot_df$phom)
      ylim_min = min(plot_df$phom)

      wilx_before_after <-
             t.test(plot_df[plot_df$all_group == "before_SFO",]$phom,
                         plot_df[plot_df$all_group == "after_SFO",]$phom,
                         paired=F,
                         alternative="less")

      wilx_within_after <-
             t.test(plot_df[plot_df$all_group == "within",]$phom,
                         plot_df[plot_df$all_group == "after_SFO",]$phom,
                         paired=F,
                         alternative="less")

      wilx_within_before <-
             t.test(plot_df[plot_df$all_group == "within",]$phom,
                         plot_df[plot_df$all_group == "before_SFO",]$phom,
                         paired=F,
                         alternative="less")


      wilx_within_shuffle <- 
        t.test(plot_df[plot_df$all_group == "within",]$phom,
               plot_df[plot_df$all_group == "shuffle",]$phom,
               paired=F,
               alternative="less")
      
      g <- 
        ggplot(plot_df) +
        geom_jitter(aes(y=phom, x=all_group, group=all_group, fill=all_group),
                    position=position_jitterdodge(.5),
                    color="gray20", alpha=0.8) + 
        geom_boxplot(aes(y=phom, x=all_group, group=all_group),
                     width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9)) +
        theme_classic() +
        theme(text=element_text(size=14, color="black"), 
              axis.text=element_text(size=14, color="black"),
              axis.ticks = element_line(color="black"),
              plot.title=element_text(size=10),
              legend.position="NA") +
        geom_line(data=data.frame(x=c(1,1,2,2),
                                  y=c(ylim_max * .75,
                                      ylim_max * .85,
                                      ylim_max * .85,
                                      ylim_max * .75)),
                  aes(x=x,y=y)) +
       geom_text(x=1.5, y=ylim_max * .85, label=signif.num(wilx_within_before$p.value * 3)) +

          geom_line(data=data.frame(x=c(1,1,3,3),
                                    y=c(ylim_max * .85,
                                        ylim_max * .95,
                                        ylim_max * .95,
                                        ylim_max * .85)),
                    aes(x=x,y=y)) +
          geom_text(x=2, y=ylim_max * .95, label=signif.num(wilx_within_after$p.value * 3)) +
          geom_line(data=data.frame(x=c(2,2,3,3),
                                    y=c(ylim_max * .8,
                                        ylim_max * .9,
                                        ylim_max * .9,
                                        ylim_max * .8)),
                    aes(x=x,y=y)) +
          geom_text(x=2.5, y=ylim_max * .9, label=signif.num(wilx_before_after$p.value * 3)) +
       #  
        geom_line(data=data.frame(x=c(1,1,4,4),
                                  y=c(ylim_max * 1.05,
                                      ylim_max * 1.075,
                                      ylim_max * 1.075,
                                      ylim_max * 1.05)),
                  aes(x=x,y=y)) +
        geom_text(x=2, y=ylim_max * 1.075, label=signif.num(wilx_within_shuffle$p.value * 1)) +
        ggtitle(metric_name) + 
        xlab("") + 
        ylab(barcodes_labs[[barcode_name]])
      
      
      for (size_name in names(sizes)) {
        
        dir.create(sprintf("%s\\all",
                           write_path))
        dir.create(sprintf("%s\\all\\%s",
                           write_path,
                           size_name))
        pdf(sprintf("%s\\all\\%s\\shuff_real_topological_similarity_%s_%s.pdf",
                    write_path,
                    size_name,
                    barcode_name,
                    metric_name),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]]) 
        
        plot(g)
        dev.off()
      }
      # 
      # plot_df_wthin_across <- 
      #   data.frame(Similarity=c(regular_df[[barcode]],shuffle_df[[barcode]]),
      #              Group=factor(c(regular_within_across_labels, shuffle_within_across_labels),
      #                           levels=c(unique(regular_within_across_labels)[2],
      #                                    unique(shuffle_within_across_labels)[2],
      #                                    unique(regular_within_across_labels)[1],
      #                                    unique(shuffle_within_across_labels)[1])))
      # 
      # #ylim_max = max(plot_df$Similarity)
      # 
      # wilx_within = wilcox.test(regular_df[[barcode]][within_across_indices],
      #                           shuffle_df[[barcode]][within_across_indices],
      #                           paired=T,
      #                           alternative="less")
      # 
      # wilx_across = wilcox.test(regular_df[[barcode]][!within_across_indices],
      #                           shuffle_df[[barcode]][!within_across_indices],
      #                           paired=T,
      #                           alternative="less")
      # 
      # wilx_across_within = wilcox.test(regular_df[[barcode]][!within_across_indices],
      #                                  regular_df[[barcode]][within_across_indices],
      #                                  paired=F)
      # g_wa <- 
      #   ggplot(plot_df_wthin_across) +
      #   geom_jitter(aes(y=Similarity, x=Group, group=Group, fill=Group),
      #               position=position_jitterdodge(.5),
      #               color="gray20", alpha=0.8) + 
      #   geom_boxplot(aes(y=Similarity, x=Group, group=Group),
      #                width=0.5, 
      #                size=1, 
      #                color=adjustcolor("gray65", alpha=1), 
      #                fill=adjustcolor("gray80", alpha=0.9)) +
      #   theme_classic() +
      #   theme(text=element_text(size=14, color="black"), 
      #         axis.ticks = element_line(color="black"),
      #         axis.text=element_text(size=14, color="black"),
      #         plot.title=element_text(size=10),
      #         legend.position="NA") +
      #   geom_line(data=data.frame(x=c(1,1,2,2),
      #                             y=c(ylim_max * 1.005,
      #                                 ylim_max * 1.007,
      #                                 ylim_max * 1.007,
      #                                 ylim_max * 1.005)),
      #             aes(x=x,y=y)) +
      #   geom_text(x=1.5, y=ylim_max * 1.008, label=signif.num(wilx_across$p.value)) +
      #   
      #   geom_line(data=data.frame(x=c(3,3,4,4),
      #                             y=c(ylim_max * 1.005,
      #                                 ylim_max * 1.007,
      #                                 ylim_max * 1.007,
      #                                 ylim_max * 1.005)),
      #             aes(x=x,y=y)) +
      #   geom_text(x=3.5, y=ylim_max * 1.008, label=signif.num(wilx_within$p.value)) +
      #   ggtitle(metric_name) + 
      #   
      #   geom_line(data=data.frame(x=c(1,1,3,3),
      #                             y=c(ylim_max * 1.010,
      #                                 ylim_max * 1.012,
      #                                 ylim_max * 1.012,
      #                                 ylim_max * 1.010)),
      #             aes(x=x,y=y)) +
      #   geom_text(x=2, y=ylim_max * 1.013, label=signif.num(wilx_across_within$p.value)) +
      #   ggtitle(metric_name) + 
      #   xlab("") + 
      #   ylab(barcodes_labs[[barcode_name]])
      # 
      # 
      # for (size_name in names(sizes)) {
      #   
      #   dir.create(sprintf("%s\\across_within",
      #                      write_path))
      #   dir.create(sprintf("%s\\across_within\\%s",
      #                      write_path,
      #                      size_name))
      #   pdf(sprintf("%s\\across_within\\%s\\topological_similarity_%s_%s.pdf",
      #               write_path,
      #               size_name,
      #               barcode_name,
      #               metric_name),
      #       height=sizes[[size_name]][["width"]],
      #       width=sizes[[size_name]][["width"]] * 1.65) 
      #   
      #   plot(g_wa)
      #   dev.off()
      # }
      
    }
  }
}


across_mice_decoding_all_trials_decode <- function(paths_group_1, chunks_group_1,
                                                   paths_group_2, chunks_group_2,
                                                   nshuffles,
                                                   remove_diagonal=F)  {
 
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm")
  

  
  # path_list <- c(rep(SFO_paths[3], times=4), rep(SFO_paths[4], times=4))
  # paths_old <- get_thirsty_quenched_paths()
  # chunks <- list(c(1,4), c(1,3), c(1:3), c(6,8), 
  #                c(1,5), c(1,4), c(1:3), c(7:9))
  
  # path_list <-  rep(SFO_paths, each=2)
  # chunks <- list(c(1,3), 4:5, c(1,4), 5:6, c(1:3), c(6:8), c(1:3), c(7:9))
  
  
  metadata_group1 <- across_mice_decoding_build_metadata(path_list = paths_group_1,
                                                         chunk_list = chunks_group_1)
  metadata_group2 <- across_mice_decoding_build_metadata(path_list = paths_group_2,
                                                         chunk_list = chunks_group_2)
  
  # mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  # mice_names <- unlist(lapply(1:len(mice_name_indices), 
  #                             function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  # days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  # days <- unlist(lapply(1:len(days_indices),
  #                       function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  # datasets_names <- paste(mice_names, days, sep = " ")
  
  #metadata <- across_mice_decoding_build_metadata()
  
  # confusion_matrix <- matrix(rep(1, times=len(path_list) * len(paths_old)),
  #                            nrow=len(path_list))
  # shuffle_confusion_matrix <- matrix(rep(1, times=len(path_list) * len(paths_old)),
  #                                    nrow=len(path_list))
  
  
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
  
  Metrics <-  list(cosine=cos_dist_vec,
                   hamming=hamming_dist,
                   manhattan=manhattan_dist)
  
  Metrics_objective <- list(cosine=which.max,

                            hamming=which.min,
                            manhattan=which.min)
  
  for (resp_name in Responses) {
    for (metric_name in names(Metrics)) {
      confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
        matrix(rep(1, times=len(metadata_group1) * len(metadata_group2)),nrow=len(metadata_group1))
      
      shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
        matrix(rep(1, times=len(metadata_group1) * len(metadata_group2)),nrow=len(metadata_group1))
    }
  }
  
  for (metric_name in names(Metrics)) {
    pvalues_all[metric_name] <- c()
  }
  

  window_size= 15
  trial_duration = 20
  
  for (idx_1 in 1:len(metadata_group1)) {
    for (idx_2 in 1:len(metadata_group2)) {
      
      
      if (remove_diagonal) {
        if(idx_1 == idx_2) {
          print(sprintf("Skipping %d-%d!", idx_1, idx_2))
          next
        }
      }
      annotated_mice_1 <- metadata_group1[[idx_1]]
      
      annotated_mice_2 <- metadata_group2[[idx_2]]
      

      stim_master_1 <- annotated_mice_1$stim_master_mat
      stim_master_2 <- annotated_mice_2$stim_master_mat 
      
      labels_mice_1 <- annotated_mice_1$cluster_mat$labs
      labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
      mask_mice_1 <- annotated_mice_1$mask
      mask_mice_2 <- annotated_mice_2$mask
      
      decoding_stats <- list()
      
      for (metric_name in names(Metrics)) {
        decoding_stats[[metric_name]] <- c()
      }

      
      
      for (shuff_i in c(-1:nshuffles)) {
        
        shuffled_stim_mat <- stim_master_1
        shuffled_stim_mat[,"Response"]
        if (shuff_i != -1) {
          #print("NOT shuff")
          shuffled_stim_mat[,1] <- shuffle(shuffled_stim_mat[,1])
        }
        
        trials_mice_1 <- which(shuffled_stim_mat[,8] %in% c(0,1,2,3,4,5))
        trials_labels_mice_1 <- shuffled_stim_mat[trials_mice_1,8]
        trials_mat_mice_1 <- c()
        activity_mat_trials_mice_1  <- c()
        
        for (r_trial in trials_mice_1) {
          trial_binned_index <- as.numeric(get_binned_index(shuffled_stim_mat[r_trial,1], window_size))
          
          clustered_trial <- labels_mice_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
          }
          
          trials_mat_mice_1 <- rbind(trials_mat_mice_1, clustered_trial)
          
        
          
          
        }

        rownames(trials_mat_mice_1) <- 1:nrow(trials_mat_mice_1)

        
        
        
        trials_mice_2 <- which(stim_master_2[,8] %in% c(0,1,2,3,4,5))
        trials_labels_mice_2 <- stim_master_2[trials_mice_2,8]
        trials_mat_mice_2 <- c()
        activity_mat_trials_mice_2  <- c()
        
        for (r_trial in trials_mice_2) {
          trial_binned_index <- as.numeric(get_binned_index(stim_master_2[r_trial,1], window_size))
          
          
          clustered_trial <- labels_mice_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
          }
          
          trials_mat_mice_2 <- rbind(trials_mat_mice_2,clustered_trial)
          

        }

        rownames(trials_mat_mice_2) <- 1:nrow(trials_mat_mice_2)

        
        
        
        translated_mat_mice_1_mice_2 <- c()
        
        
        for (trial_idx in 1:nrow(trials_mat_mice_1)) {
          trial_before <- trials_mat_mice_1[trial_idx,]
          trial_after <- trial_before
          
          for (cluster_idx in 1:len(mask_mice_1)) {
            trial_after[which(trial_after == mask_mice_1[cluster_idx])] <- mask_mice_2[cluster_idx] 
          }
          
          translated_mat_mice_1_mice_2 <- rbind(translated_mat_mice_1_mice_2,
                                                as.numeric(trial_after))
          
        }
        
        
        decoded_vec <- list()
        for (metric_name in names(Metrics)) {
          decoded_vec[[metric_name]] <- c()
        }
        
        for (decoded_trial_idx in 1:nrow(trials_mat_mice_2)) {
          decoded_trial <- trials_mat_mice_2[decoded_trial_idx,]
          
          confidence <-   
            apply(translated_mat_mice_1_mice_2, 1, 
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
                              colMeans(confidence[trials_labels_mice_1 == res,]))
            
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
          }
        }
        

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

      
    }
  }
  


  
  return(list(pvalues_all=pvalues_all,
              shuffles=shuffle_matrix_list,
              decodede=confusion_matrix_list))

}


plot_decoding <- function(decoding_df, external_label=c(), ext="", compare_label=F, average_over_decoded=T) {
  write_path <- sprintf("%s\\figure_5\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\across_mice_decoding", write_path)
  dir.create(write_path)  
  
  if (average_over_decoded) {
    write_path <- sprintf("%s\\%s", write_path, "averaged_over_decoded")
    dir.create(write_path)
  }
  
  write_path <- sprintf("%s\\%s", write_path, ext)
  dir.create(write_path)
  
  
  
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
  
  
  decoded_trials <- names(decoding_df$decodede)
  decoded_datasets_n <- ncol(decoding_df$decodede$Hit_cosine)

  
  for (decoded_trial_type in decoded_trials) {
    
    if (average_over_decoded) {
      if (len(external_label) > 0) {
        regular <- c()
        shuffles <- c()
        
        for (lbl in unique(external_label)) {
          print(which(external_label == lbl))
          regular <-  c(regular,
                        colMeans(decoding_df$decodede[[decoded_trial_type]][which(external_label == lbl), ]))
          shuffles <-  c(shuffles,
                         colMeans(decoding_df$shuffles[[decoded_trial_type]][which(external_label == lbl), ]))
        }
      } else {
        regular <-  colMeans(decoding_df$decodede[[decoded_trial_type]])
        shuffles <-  colMeans(decoding_df$shuffles[[decoded_trial_type]])
      }
    } else {
      regular <-  c(decoding_df$decodede[[decoded_trial_type]])
      shuffles <-  c(decoding_df$shuffles[[decoded_trial_type]])
    }
    
    # plot_df <- data.frame(decoded=c(regular, shuffles),
                          group=rep(c("Regular", "Shuffle"), each=len(shuffles)))
    
    if (len(external_label) > 0) {
       if (!average_over_decoded) {
         ext_label = rep(external_label, times=decoded_datasets_n)
       } else {
         ext_label <- rep(unique(external_label), each=decoded_datasets_n)
       }
      
      
       ext_label <- c(ext_label, ext_label)
       plot_df$group <- paste(plot_df$group, ext_label, sep=" - ")
    }
    
    gcomp <- 
      ggplot(plot_df, aes(x=group, y=decoded)) +
      geom_point() +
      geom_boxplot(width=0.5, 
                   size=1, 
                   color=adjustcolor("gray65", alpha=1), 
                   fill=adjustcolor("gray80", alpha=0.9),
                   outlier.shape=NA) + 
      ylim(0,1) + 
      theme_classic() +
      xlab("") +
      ylab("Correctly classified (%)") +
      ggtitle(decoded_trial_type) + 
      theme(text=element_text(size=14, color="black"),
            plot.title = element_text(size=9, color="black"),
            legend.position="NA",
            legend.title = element_blank(),
            axis.text = element_text(color="black"),
            axis.text.x = element_text(size=10),#element_text(angle = 45, vjust=.6),
            axis.ticks = element_line(color="black")) 
    
    
    if (compare_label) { 
      
      gcomp_with_lines <- 
        ggplot(plot_df, aes(x=group, y=decoded)) +
        geom_hline(yintercept=1/6, lty="dashed", size=.5)
      
      if (average_over_decoded) {
      comp_df <- cbind(regular[1:decoded_datasets_n],
                       regular[(decoded_datasets_n + 1):(decoded_datasets_n * 2)])
      
      shuffle_comp_df <- cbind(shuffles[1:decoded_datasets_n],
                               shuffles[(decoded_datasets_n + 1):(decoded_datasets_n * 2)])
      } else {
        all_labels_comp <- ext_label[1:(len(ext_label) / 2)]
        comp_df <- cbind(regular[all_labels_comp == unique(external_label)[1]],
                         regular[all_labels_comp == unique(external_label)[2]])
        

      }
      for (idx in 1:nrow(comp_df)) {
        gcomp_with_lines <- 
          gcomp_with_lines + 
          geom_line(data=data.frame(x=c(paste("Regular", unique(external_label)[1], sep=" - "),
                                        paste("Regular", unique(external_label)[2], sep=" - ")),
                                    y=c(comp_df[idx,1], comp_df[idx,2])), 
                    aes(x=x,y=y,group=1),
                    alpha=.3) #+
          # geom_line(data=data.frame(x=c(paste("Shuffle", unique(external_label)[1], sep=" - "),
          #                               paste("Shuffle", unique(external_label)[2], sep=" - ")),
          #                           y=c(shuffle_comp_df[idx,1], shuffle_comp_df[idx,2])), 
          #           aes(x=x,y=y,group=1),
          #           alpha=.3)
      }
      
      wilk <- wilcox.test(comp_df[,1], comp_df[,2], paired=T, correct = F)

      gcomp_lines_f <- 
        gcomp_with_lines +       
        geom_point() +
        geom_boxplot(width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9),
                     outlier.shape=NA) + 
        ylim(0,1) + 
        theme_classic() +
        xlab("") +
        ylab("Correctly classified (%)") +
        ggtitle(decoded_trial_type) + 
        theme(text=element_text(size=14, color="black"),
              plot.title = element_text(size=9, color="black"),
              legend.position="NA",
              legend.title = element_blank(),
              axis.text = element_text(color="black"),
              axis.text.x = element_text(size=10),#element_text(angle = 45, vjust=.6),
              axis.ticks = element_line(color="black")) +
        geom_text(data=data.frame(x=1.5, y=0.8, label=signif.num(wilk$p.value)),
                  aes(x=x,y=y,label=label))
      
    } else {
    
    gcomp_with_lines <- 
      ggplot(plot_df, aes(x=group, y=decoded)) +
      geom_hline(yintercept=1/6, lty="dashed", size=.5)
    
    
    for (idx in 1:len(regular)) {
      gcomp_with_lines <- 
        gcomp_with_lines + 
        geom_line(data=data.frame(x=c("Regular", "Shuffle"),
                                  y=c(regular[idx], shuffles[idx])), 
                  aes(x=x,y=y,group=1),
                  alpha=.3)
    }
    
    wilk <- wilcox.test(regular, shuffles, paired=T, correct = T)
    gcomp_lines_f <- 
      gcomp_with_lines +       
      geom_point() +
      geom_boxplot(width=0.5, 
                   size=1, 
                   color=adjustcolor("gray65", alpha=1), 
                   fill=adjustcolor("gray80", alpha=0.9),
                   outlier.shape=NA) + 
      ylim(0,1) + 
      theme_classic() +
      xlab("") +
      ylab("Correctly classified (%)") +
      ggtitle(decoded_trial_type) + 
      theme(text=element_text(size=14, color="black"),
            plot.title = element_text(size=9, color="black"),
            legend.position="NA",
            legend.title = element_blank(),
            axis.text = element_text(color="black"),
            axis.text.x = element_text(size=10),#element_text(angle = 45, vjust=.6),
            axis.ticks = element_line(color="black"))  + 
      geom_text(data=data.frame(x=1.5, y=0.8, label=signif.num(wilk$p.value)),
                aes(x=x,y=y,label=label))
    }
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\%s",write_path, size_name))
      
      pdf(sprintf("%s\\%s\\lines_%s.pdf",
                  write_path,
                  size_name,
                  decoded_trial_type),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]] * 1.4)
      
      plot(gcomp_lines_f)
      dev.off()
      
      pdf(sprintf("%s\\%s\\box_%s.pdf",
                  write_path,
                  size_name,
                  decoded_trial_type),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]] * 1.4)
      
      plot(gcomp)
      dev.off()
      
      
    }  
  }
 
  if (len(external_label > 0)) {
    all_external_labels  <- unique(external_label)
    
    indices <- list()
    
    
    for (lbl in all_external_labels) {
      which_lbls <- which(lbl == external_label)
      indices[[lbl]] <- unlist(lapply(which_lbls, function(i) {(decoded_datasets_n * (i - 1) + 1):(decoded_datasets_n * i )}))
    }
  } else {
    indices <- list(all=1:nrow(decoding_df$pvalues_all$cosine))
  }
  
  for (ind_name in names(indices)) {
    ind_to_use <- indices[[ind_name]]
    print("Using:")
    print(ind_to_use)
    print(ind_name)
    
    for (metric_type in names(decoding_df$pvalues_all)) {
      
      gm <- 
      ggplot(melt(decoding_df$pvalues_all[[metric_type]][ind_to_use,]),
             aes(x=Var2,y=value))  + 
        geom_jitter(alpha=.2, size=.5) +
        geom_violin(fill="gray50", color=NA) +
        stat_summary() + 
        theme_classic() +
        theme(text=element_text(size=14, color="black"), 
              axis.text=element_text(size=14, color="black"),
              axis.ticks = element_line(color="black"),
              plot.title=element_text(size=10),
              legend.position="NA")  + 
        xlab("") + 
        ylab("% defeated shuffles") +
        ggtitle(metric_type)
      
      
      for (size_name in names(sizes)) {
        
        dir.create(sprintf("%s\\%s",write_path, size_name))
        dir.create(sprintf("%s\\%s\\%s",write_path, size_name, metric_type))
        
        pdf(sprintf("%s\\%s\\%s\\%sdecoding_pvalues.pdf",
                    write_path,
                    size_name,
                    metric_type,
                    ifelse(ind_name == "all", "", paste(ind_name, "_"))),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]] * 1.4)
        
        plot(gm)
        dev.off()
  
    
      }  
    }
  }
}
# 




