
figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"




spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))


across_mice_decoding_all_trials_decode <- function()  {
  #output_path <- sprintf("%s//single_trial_analysis//", output_path_f)
  write_path <- sprintf("%s\\figure_4\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\across_mice_decoding", write_path)
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
  
  
  path_list <- get_thirsty_quenched_paths()
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  metadata <- across_mice_decoding_build_metadata()
  
  confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                             nrow=len(path_list))
  shuffle_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                                     nrow=len(path_list))
  
  
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
        matrix(rep(1, times=len(path_list) ** 2),nrow=len(path_list))
    
    shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
      matrix(rep(1, times=len(path_list) ** 2),nrow=len(path_list))
    }
  }
  
  for (metric_name in names(Metrics)) {
    pvalues_all[metric_name] <- c()
  }
  
  # cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
  #                                   nrow=len(path_list))
  # shuffle_cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
  #                                           nrow=len(path_list))

  activity_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                                      nrow=len(path_list))
  
  window_size= 15
  trial_duration = 20
  
  for (idx_1 in 1:len(path_list)) {
    for (idx_2 in 1:len(path_list)) {
      if (idx_1 == idx_2) {
        next
      }
      
      # p1 <- path_list[idx_1]
      # p2 <- path_list[idx_2]
      # 
      # orig_mat_1 <- get_reduced_mat_full_day(p1, window_size=15, just_original_mat = T)
      # orig_mat_2 <- get_reduced_mat_full_day(p2, window_size=15, just_original_mat = T)
      # 
      # if (nrow(orig_mat_2) > ncol(orig_mat_2)) {
      #   orig_mat_2 <- t(orig_mat_2)
      # }
      # 
      # if (nrow(orig_mat_1) > ncol(orig_mat_1)) {
      #   orig_mat_1 <- t(orig_mat_1)
      # }
      # 
      # mean_activity_1 <- colMeans(orig_mat_1)
      # mean_activity_2 <- colMeans(orig_mat_2)
      # 
      
      
      # annotated_struct_path_1 <-  sprintf("%s\\%s\\structure_annotation_all_trials_%d.Rda", 
      #                                     output_path,
      #                                     str_replace(datasets_names[idx_1]," ", "_"),
      #                                     7)
      # 
      # load(annotated_struct_path_1)
      annotated_mice_1 <- metadata[[idx_1]]
      
      # annotated_struct_path_2 <-  sprintf("%s\\%s\\structure_annotation_all_trials_%d.Rda", 
      #                                     output_path,
      #                                     str_replace(datasets_names[idx_2]," ", "_"),
      #                                     7)
      # 
      # load(annotated_struct_path_2)
      annotated_mice_2 <- metadata[[idx_2]]
      
      #stim_master_1 <- get_stim_mat(p1,just_mat = T, window_size = window_size, chunk=-1)
      #stim_master_2 <- get_stim_mat(p2,just_mat = T, window_size = window_size, chunk=-1)
      stim_master_1 <- annotated_mice_1$stim_master_mat
      stim_master_2 <- annotated_mice_2$stim_master_mat 
      
      labels_mice_1 <- annotated_mice_1$cluster_mat$labs
      labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
      mask_mice_1 <- annotated_mice_1$prob_annot$EntireTrials3$all
      mask_mice_2 <- annotated_mice_2$prob_annot$EntireTrials3$all
      
      decoding_stats <- list()
      
      for (metric_name in names(Metrics)) {
        decoding_stats[[metric_name]] <- c()
      }
      #cosine_decoding_stats <- c()
      
      
      for (shuff_i in c(-1:200)) {
        
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
          
          
          # activity_in_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          # activity_in_trial[is.na(activity_in_trial)] <- 0
          # 
          # activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_in_trial)
          
          
        }
        
        #annot_df_mice_1 <- data.frame(Result=paste(annotated_mice_1$annot_df[,1], annotated_mice_1$annot_df[,2]))
        rownames(trials_mat_mice_1) <- 1:nrow(trials_mat_mice_1)
        #rownames(annot_df_mice_1) <- 1:nrow(trials_mat_mice_1)
        
        
        
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


across_mice_decoding_build_metadata <- function(path_list=get_thirsty_quenched_paths(), chunk_list=list(), control=F) {
  
  #path_list <- get_thirsty_quenched_paths()
  results <- list()
  
  for (idx in 1:len(path_list)) {
    
    if (len(chunk_list) > 0) {
      chunk = chunk_list[[idx]]
    } else {
      chunk = -1
    }
    
    p <- path_list[[idx]]
    all_mat <- get_mat_with_preset(p, "dalshtaim", chunk=chunk)
    if (control) {
      stim_master_mat <- get_stim_mat_control(p,just_mat = T, window_size = 15)
    } else {
      stim_master_mat <- get_stim_mat(p,just_mat = T, window_size = 15, chunk=chunk)
    }
    clust_mat_all <- get_clusters_mat_kmeans(p, chunk=chunk)
    
    num_of_clusters <- clust_mat_all$nc
    cluster_labels <- clust_mat_all$labs
    
    relevant_trials_mat <- stim_master_mat[stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    num_of_trials <- nrow(relevant_trials_mat)
    
    chronological_mat <- matrix(rep(-1, times=num_of_trials * 20), nrow=num_of_trials)
    prob_matrix <- matrix(rep(0, times=(num_of_clusters + 1) ** 2), nrow=(num_of_clusters + 1))
    
    colnames(prob_matrix) <- c(-1, 1:num_of_clusters)
    rownames(prob_matrix) <- c(-1, 1:num_of_clusters)
    
    colnames(chronological_mat) <- 1:20
    rownames(chronological_mat) <- 1:num_of_trials
    
    colnames(chronological_mat) <- sprintf("t =  %.1f", 1:20 / 2)
    
    trials <- 1:num_of_trials
    
    
    for (trial in trials) {
      trial_ind <- get_binned_index(relevant_trials_mat[trial,1], 15)
      trial_ind <- (trial_ind):(trial_ind + 19)
      
      clusters_in_trial <- cluster_labels[trial_ind]
      
      if (sum(is.na(clusters_in_trial)) > 0) {
        clusters_in_trial[is.na(clusters_in_trial)] <- 
          clusters_in_trial[!is.na(clusters_in_trial)][len(clusters_in_trial[!is.na(clusters_in_trial)])]
      }
      
      
      
      chronological_mat[trial,] <- clusters_in_trial
    }
    
    #title <- sprintf("%s - %d Clusters", datasets_names[[idx]], num_of_clusters)
    
    
    for (trial_idx in 1:len(trials)) {
      for (clust_idx in 1:(ncol(chronological_mat)  - 1)) {
        first_clust <- as.character(chronological_mat[trial_idx, clust_idx])
        second_clust <- as.character(chronological_mat[trial_idx, clust_idx + 1])
        #print(sprintf("From %s to %s", first_clust, second_clust))
        prob_matrix[first_clust, second_clust] <- prob_matrix[first_clust, second_clust] + 1
        
      }
    }
    
    prb <- prob_matrix / sum(prob_matrix)
    
    
    
    hc <- hclust(dist(prb), method="ward.D2")
    print(rownames(prb)[hc$order])
    

    annot_df <- data.frame(trialtype=relevant_trials_mat[,"TrialType"], response=relevant_trials_mat[,"Response"] %% 2)
    rownames(annot_df) <- 1:num_of_trials
    
    

    
    
    
    
    final_list <- list()
    final_list$cluster_mat <- clust_mat_all
    final_list$mask <- rownames(prb)[hc$order]
    final_list$trials_mat <- chronological_mat
    final_list$annot_df <- annot_df
    final_list$probability_mat <- prb
    final_list$stim_master_mat <- stim_master_mat
    final_list$red_mat <- all_mat
    
    results[[idx]] <- final_list
    

  }
  
  return(results)
}










# 
# #decoded_vec <- c()
# cosine_decoded_vec <- c()
# 
# for (decoded_trial_idx in 1:nrow(trials_mat_mice_2)) {
#   decoded_trial <- trials_mat_mice_2[decoded_trial_idx,]
#   
#   cosine_confidence <-   
#     apply(translated_mat_mice_1_mice_2, 1, 
#           function(trial) {cos_dist_vec(trial, decoded_trial)})  
#   # confidence <- 
#   #   apply(translated_mat_mice_1_mice_2, 1, 
#   #         function(trial) {euc_dist(trial, decoded_trial)})  
#   
#   
#   #conf_vec <- c()
#   cosine_conf_vec <- c()
#   trial_results <- sort(unique(trials_labels_mice_1))
#   for (res in trial_results) {
#     # conf_vec <- c(conf_vec,
#     #               mean(confidence[trials_labels_mice_1 == res]))
#     
#     cosine_conf_vec <- c(cosine_conf_vec,
#                          mean(cosine_confidence[trials_labels_mice_1 == res]))
#     
#   }
#   
#   
#   # is_hit <- mean(confidence[reward_trials_labels_mice_1 == 0]) 
#   # is_miss <- mean(confidence[reward_trials_labels_mice_1 == 1])
#   cosine_decoded_vec <- c(cosine_decoded_vec, trial_results[which.max(cosine_conf_vec)])
#   #decoded_vec <- c(decoded_vec, trial_results[which.min(conf_vec)])
# }
# 
# 
# #decoding_statistics <-  sum(decoded_vec == trials_labels_mice_2) / len(trials_labels_mice_2)
# 
# cosine_decoding_statistics <- c()
# for (response_type in c(0,1,2,3,4,5)) {
#   response_indices <- trials_labels_mice_2 == response_type
#   cosine_decoding_statistics <- 
#     c(cosine_decoding_statistics,
#       sum(cosine_decoded_vec[response_indices] == trials_labels_mice_2[response_indices]) / len(trials_labels_mice_2[response_indices]))
# }
# 
# 
# 
# 
# # print(sprintf("Decoding statistics from %d to %d, %f (cosine: %f)", 
# #               idx_1,
# #               idx_2,
# #               decoding_statistics,
# #               cosine_decoding_statistics))
# 
# # decoding_stats <- c(decoding_stats,
# #                     decoding_statistics)
# 
# cosine_decoding_stats <- rbind(cosine_decoding_stats, cosine_decoding_statistics)
# }
# 
# 
# #print(table(trials_labels_mice_2[decoded_vec == trials_labels_mice_2]))
# 
# #decoding_statistics <- decoding_stats[1]
# colnames(cosine_decoding_stats) <- Responses
# cosine_decoding_statistics <- cosine_decoding_stats[1,]
# #shuffled_decoding_statistics <- mean(decoding_stats[-1])
# shuffled_cosine <- colMeans(cosine_decoding_stats[-1,])
# 
# ecdfs <- apply(cosine_decoding_stats[-1,], 1, ecdf)
# pvalues <-  unlist(lapply(1:len(cosine_decoding_stats[1,]), function(acc_idx) {ecdfs[[acc_idx]](cosine_decoding_stats[1,acc_idx])}))
# names(pvalues) <- Responses
# pvalues_all <- rbind(pvalues_all,
#                      pvalues)
# 
# print(sprintf("######################## SHUFFLED DECODING: (Shuffle vs decoder) - %d vs  %d ####", idx_1, idx_2))
# for (resp_name in Responses)  {
#   
#   confusion_matrix_list[[resp_name]][idx_1, idx_2] <- cosine_decoding_statistics[resp_name]
#   shuffle_matrix_list[[resp_name]][idx_1, idx_2] <- shuffled_cosine[resp_name]
#   #print(sprintf("%.3f vs %.3f", shuffled_decoding_statistics, decoding_statistics))
#   print(sprintf("(%d). %s: %.3f vs %.3f [Pvalue: %f]", 
#                 nrow(pvalues_all),
#                 resp_name, 
#                 shuffle_matrix_list[[resp_name]][idx_1, idx_2],
#                 confusion_matrix_list[[resp_name]][idx_1, idx_2],
#                 pvalues[[resp_name]]))
#   
#   
# }
# print("#######################")