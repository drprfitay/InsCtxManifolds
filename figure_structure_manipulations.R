
figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\figures_structure_manipulations\\", base_output_path) 

dir.create(base_output_path)
dir.create(output_path)

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_HS <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
metadata_control <- across_mice_decoding_build_metadata(get_control_paths(), control=T)
metadata_SFO <- across_mice_decoding_build_metadata(get_SFO_paths())
metadata_AGRP <- across_mice_decoding_build_metadata(get_agrp_paths())
metadata_hypertonicsaline <- across_mice_decoding_build_metadata(get_hypertonic_saline_paths())

SFO_chunks = list(`1`=4:5,
                  `2`=5:6,
                  `3`=6:8,
                  `4`=7:9)


AGRP_chunks = list(`1`=5:9,
                   `2`=6:10,
                   `3`=5:9,
                   `4`=3:4)

hypertonic_chunks <- list(`1`=4:5,
                          `2`=4,
                          `3`=3:4)


across_mice_decoding_all_trials_decode_manipulations <- function()
{
  
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm")
  cos_dist_vec <- function(a,b) {cosine(a,b)[1,1]}
  shuffle <- function(vec) {sample(vec, len(vec))}
  pre_manipulation_confusion_matrix_list <- list()
  pre_manipulation_shuffle_matrix_list <- list()
  post_manipulation_confusion_matrix_list <- list()
  post_manipulation_shuffle_matrix_list <- list()
  pre_pvalues_all <- list()
  post_pvalues_all <- list()
  
  manhattan_dist <- function(a, b){
    dist <- abs(a-b)
    dist <- sum(dist)
    return(dist)
  }
  
  cosine_qg <- function(a,b) {a <- paste(a, collapse="");
  b <- paste(b, collapse="");
  stringdist(a,b, "cosine", q=1)}
  hamming_dist <- function(a,b) {return(sum(a != b))}
  
  
  lv_dist <- function(a,b) {a <- paste(a, collapse="");
  b <- paste(b, collapse="");
  stringdist(a,b, "lv")}
  
  Metrics <-  list(cosine=cos_dist_vec)
  #jaccard=jaccard)
  #lv=lv_dist)
  #hamming=hamming_dist,
  #manhattan=manhattan_dist)
  
  manipulation_chunks <- SFO_chunks
  Metrics_objective <- list(cosine=which.max,
                            euc=which.min,
                            lv=which.min,
                            jaccard=which.max,
                            cosineqg=which.max,
                            hamming=which.min,
                            manhattan=which.min)
  
  for (resp_name in Responses) {
    for (metric_name in names(Metrics)) {
      pre_manipulation_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
      
      pre_manipulation_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
      
      post_manipulation_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
      
      post_manipulation_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
        matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
    }
  }
  
  window_size= 15
  trial_duration = 20
  
  for (idx_1 in 1:len(metadata_1)) {
    for (idx_2 in 1:len(metadata_2)) {
      
      # if (idx_1 == idx_2) {
      #   next
      # }
      # 
      
      annotated_mice_1 <- metadata_1[[idx_1]]
      annotated_mice_2 <- metadata_2[[idx_2]]
      
      stim_master_1 <- annotated_mice_1$stim_master_mat
      stim_master_2 <- annotated_mice_2$stim_master_mat 
      
      labels_mice_1 <- annotated_mice_1$cluster_mat$labs
      labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
      mask_mice_1 <- annotated_mice_1$mask
      mask_mice_2 <- annotated_mice_2$mask
      
      pre_decoding_stats <- list()
      post_decoding_stats <- list()
      
      for (metric_name in names(Metrics)) {
        pre_decoding_stats[[metric_name]] <- c()
        post_decoding_stats[[metric_name]] <- c()
        
      }

      trials_mice_2 <- which(stim_master_2[,"Response"] %in% c(0,1,2,3,4,5))
      trials_labels_mice_2 <- stim_master_2[trials_mice_2,"Response"]
      pre_manipulation_trials <- !stim_master_2[trials_mice_2,ncol(stim_master_2)] %in% manipulation_chunks[[idx_2]]
      post_manipulation_trials <- stim_master_2[trials_mice_2,ncol(stim_master_2)] %in% manipulation_chunks[[idx_2]]
      
      trials_mat_mice_2 <- c()

      
      for (r_trial in trials_mice_2) {
        trial_binned_index <- as.numeric(get_binned_index(stim_master_2[r_trial,"Frames"], window_size))
        clustered_trial <- labels_mice_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
        
        if (sum(is.na(clustered_trial)) > 0 ) {
          clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
        }
        
        trials_mat_mice_2 <- rbind(trials_mat_mice_2,clustered_trial)
      }
      
      rownames(trials_mat_mice_2) <- 1:nrow(trials_mat_mice_2)
      
      
      for (shuff_i in c(-1:200)) {
        
        shuffled_stim_mat <- stim_master_1
        #shuffled_stim_mat[,"Response"]
        if (shuff_i != -1) {
          #print("NOT shuff")
          if (shuff_i %% 50 == 0) {
            print(shuff_i)
          }
          
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
                              mean(confidence[trials_labels_mice_1 == res]))
            
          }
          
          colnames(conf_vec) <- names(Metrics)
          
          for (metric_name in names(Metrics)) {
            decoded_vec[[metric_name]] <- c(decoded_vec[[metric_name]],
                                            trial_results[Metrics_objective[[metric_name]](conf_vec[,metric_name])])
          }
        }
        
        pre_decoding_statistics <- list()
        post_decoding_statistics <- list()
        
        for (metric_name in names(Metrics)) { 
          pre_decoding_statistics[[metric_name]] <- c()
          post_decoding_statistics[[metric_name]] <- c()
        }
        
        
        for (metric_name in names(Metrics)) { 
          for (response_type in c(0)) {
            response_indices <- trials_labels_mice_2 == response_type
            
            pre_indices <- response_indices & pre_manipulation_trials
            post_indices <- response_indices & post_manipulation_trials
            
            pre_decoding_accuracy = sum(decoded_vec[[metric_name]][pre_indices] == trials_labels_mice_2[pre_indices]) / len(trials_labels_mice_2[pre_indices])
            post_decoding_accuracy = sum(decoded_vec[[metric_name]][post_indices] == trials_labels_mice_2[post_indices]) / len(trials_labels_mice_2[post_indices])
            
            pre_decoding_statistics[[metric_name]] <-  c(pre_decoding_statistics[[metric_name]], pre_decoding_accuracy)
            post_decoding_statistics[[metric_name]] <-  c(post_decoding_statistics[[metric_name]], post_decoding_accuracy)
          }
        }
        
        
        for (metric_name in names(Metrics)) { 
          pre_decoding_stats[[metric_name]] <- rbind(pre_decoding_stats[[metric_name]], pre_decoding_statistics[[metric_name]])
          post_decoding_stats[[metric_name]] <- rbind(post_decoding_stats[[metric_name]], post_decoding_statistics[[metric_name]])
        }
      }
      
      for (metric_name in names(Metrics)) { 
        
        
        colnames(pre_decoding_stats[[metric_name]]) <- Responses[1]
        colnames(post_decoding_stats[[metric_name]]) <- Responses[1]
        
        pre_decoding_statistics <- pre_decoding_stats[[metric_name]][1,]
        pre_shuffled_statistics <- mean(pre_decoding_stats[[metric_name]][-1,])
        post_decoding_statistics <- post_decoding_stats[[metric_name]][1,]
        post_shuffled_statistics <- mean(post_decoding_stats[[metric_name]][-1,])
        
        pre_ecdfs <- list(ecdf(pre_decoding_stats[[metric_name]][-1,]))#apply(pre_decoding_stats[[metric_name]][-1,], 2, ecdf)
        post_ecdfs <- list(ecdf(pre_decoding_stats[[metric_name]][-1,]))#apply(post_decoding_stats[[metric_name]][-1,], 2, ecdf)
        
        pre_pvalues <-  unlist(lapply(1:len(pre_decoding_stats[[metric_name]][1,]), 
                                  function(acc_idx) {
                                    pre_ecdfs[[acc_idx]](pre_decoding_stats[[metric_name]][1,acc_idx])
                                  }))
        
        post_pvalues <- unlist(lapply(1:len(post_decoding_stats[[metric_name]][1,]), 
                                      function(acc_idx) {
                                        post_ecdfs[[acc_idx]](post_decoding_stats[[metric_name]][1,acc_idx])
                                      }))
        
        names(pre_pvalues) <- Responses[1]
        names(post_pvalues) <- Responses[1]
        
        pre_pvalues_all[[metric_name]] <- rbind(pre_pvalues_all[[metric_name]], pre_pvalues)
        post_pvalues_all[[metric_name]] <- rbind(post_pvalues_all[[metric_name]], post_pvalues)
        
        print(sprintf("######################## SHUFFLED DECODING: (Shuffle vs decoder) - %d vs  %d ####", idx_1, idx_2))
        
        for (resp_name in Responses[1])  {
          
          mat_name <- sprintf("%s_%s", resp_name, metric_name)
          
          pre_manipulation_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- pre_decoding_statistics[resp_name]
          pre_manipulation_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- pre_shuffled_statistics[resp_name]
          post_manipulation_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- post_decoding_statistics[resp_name]
          post_manipulation_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- post_shuffled_statistics[resp_name]
          
          #print(sprintf("%.3f vs %.3f", shuffled_decoding_statistics, decoding_statistics))
          print(sprintf("(%d). %s: Pr[(S) %.3f vs %.3f], Po[(S) %.3f vs %.3f]  [Pvalue: %f - %f] -> %s", 
                        nrow(pre_pvalues_all[[metric_name]]),
                        resp_name, 
                        pre_manipulation_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                        pre_manipulation_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                        post_manipulation_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                        post_manipulation_confusion_matrix_list[[mat_name]][idx_1, idx_2],                        
                        pre_pvalues[[resp_name]],
                        post_pvalues[[resp_name]],
                        metric_name))
          
          
        }
        
        print("#######################")
        
      }
    }
  }

  dir.create(sprintf("%s\\data\\figure_5\\", base_output_path))
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_confusion_matrices.Rda", base_output_path), pre_manipulation_confusion_matrix_list)
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_shuffle_confusion_matrices.Rda", base_output_path), pre_manipulation_shuffle_matrix_list)
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_pvalues_decoding.Rda", base_output_path), pre_pvalues_all)
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_confusion_matrices.Rda", base_output_path), post_manipulation_confusion_matrix_list)
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_shuffle_confusion_matrices.Rda", base_output_path), post_manipulation_shuffle_matrix_list)
  save(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_pvalues_decoding.Rda", base_output_path), post_pvalues_all)
}


plot_AGRP_SFO_hypertonic_trial_dynamics <- function()
{
  
  manipulation_paths <- c(get_SFO_paths(),
                          get_agrp_paths(),
                          get_hypertonic_saline_paths())
  
  manipulation_prefix <- rep(c("SFO", "AGRP", "SALINE"), 
                             times=c(len(get_SFO_paths()),
                                     len(get_agrp_paths()),
                                     len(get_hypertonic_saline_paths())))
  
  datasets_names <- get_datasets_names(manipulation_paths, sep="_")
  datasets_names <- paste(manipulation_prefix, datasets_names, sep="_")
  
  metadata <- list()
  metadata <- append(metadata, metadata_SFO)
  metadata <- append(metadata, metadata_AGRP)
  metadata <- append(metadata, metadata_hypertonicsaline)
  
  base_figure_path <- sprintf("%s\\AGRP_SFO_hypertonic_trial_dynamics", output_path)
  dir.create(base_figure_path)
  
  for (idx in 1:len(metadata)) { 
    
    
    work_mat <- metadata[[idx]]$trials_mat
    work_annot <- metadata[[idx]]$annot_df
    
    non_pavlov_trials <- work_annot[,1] != 1
    
    work_mat <- work_mat[non_pavlov_trials,]
    work_annot <- work_annot[non_pavlov_trials,]
    
    rownames(work_mat) <- 1:nrow(work_mat)
    rownames(work_annot) <- 1:nrow(work_mat)
    
    reward_trials <- work_annot[,1] == 3
    sorted_trials <-  order(paste(work_annot[,1], work_annot[,2]))
    cluster_color_label <- spec_cg(len(unique(metadata[[idx]]$cluster_mat$labs)))
    
    
    trials_mat_reward_only <- pheatmap(work_mat[reward_trials,], 
                                       cluster_rows=F, 
                                       cluster_cols=F, 
                                       col=cluster_color_label, 
                                       show_rownames = F,
                                       show_colnames = F,
                                       annotation_row = work_annot[reward_trials,])
    
    annoCol<-list(trialtype=c(`1`="blue", `3`="red", `4`="gray70", `5`="grey60"), response=c(`0`="white",`1`="black"))
    trials_mat_clustered <- pheatmap(work_mat,
                                     cluster_rows=T, 
                                     cluster_cols=F, 
                                     show_rownames = F,
                                     show_colnames = F,
                                     col=cluster_color_label, 
                                     annotation_row = work_annot, 
                                     annotation_colors =  annoCol)
    
    
    trials_mat_sorted <- pheatmap(work_mat[sorted_trials,], 
                                  cluster_rows=F, 
                                  cluster_cols=F, 
                                  show_rownames = F,
                                  show_colnames = F,
                                  col=cluster_color_label, 
                                  annotation_row = work_annot[sorted_trials,])
    
    trials_mat_all <- pheatmap(metadata[[idx]]$trials_mat, 
                               cluster_rows=F, 
                               cluster_cols=F, 
                               show_rownames = F,
                               show_colnames = F,
                               col=cluster_color_label, 
                               annotation_row = work_annot)
    
    
    
    prob_trans <- 
      pheatmap(metadata[[idx]]$probability_mat, cluster_rows=F, cluster_cols=F,
               breaks=seq(0,.025, length.out=100),
               color = rev(colorspace::heat_hcl(100)))
    
    
    prob_trans_sorted <- 
      pheatmap(metadata[[idx]]$probability_mat[metadata[[idx]]$mask,
                                               metadata[[idx]]$mask], cluster_rows=F, cluster_cols=F,
               breaks=seq(0,.025, length.out=100),
               color = rev(colorspace::heat_hcl(100)))
    
    
    write_path <- sprintf("%s\\%s\\",
                          base_figure_path,
                          datasets_names[idx])
    dir.create(write_path)
    
    for (size_name in names(a4_sizes)) {
      
      size <- a4_sizes[[size_name]]
      dir.create(sprintf("%s\\%s",
                         write_path,
                         size_name))
      
      pdf(sprintf("%s\\%s\\trials_mat_all.pdf",
                  write_path,
                  size_name),
          height=size * 3,
          width=size * 1.5)
      # unit="in",
      # res=1500)
      
      plot(trials_mat_all[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\trials_mat_reward.pdf",
                  write_path,
                  size_name),
          height=size * 3,
          width=size * 1.5)
      # unit="in",
      # res=1500)
      
      plot(trials_mat_reward_only[[4]])
      dev.off()
      
      
      pdf(sprintf("%s\\%s\\trials_mat_sorted.pdf",
                  write_path,
                  size_name),
          height=size * 3,
          width=size * 1.5)
      # unit="in",
      # res=1500)
      
      plot(trials_mat_sorted[[4]])
      dev.off()
      
      
      pdf(sprintf("%s\\%s\\trials_mat_clustered.pdf",
                  write_path,
                  size_name),
          height=size * 3,
          width=size * 1.5)
      # unit="in",
      # res=1500)
      
      plot(trials_mat_clustered[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\transition_probability.pdf",
                  write_path,
                  size_name),
          height=size,
          width=size)
      # unit="in",
      # res=1500)
      
      plot(prob_trans[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\transition_probability_sorted.pdf",
                  write_path,
                  size_name),
          height=size,
          width=size)
      # unit="in",
      # res=1500)
      
      plot(prob_trans_sorted[[4]])
      dev.off()
      
      
      
    }
  }
}

plot_AGRP_SFO_hypertonic_trial_similarity <- function()
{
  
  
  base_figure_path <- sprintf("%s\\AGRP_SFO_hypertonic_trial_similarity", output_path)
  dir.create(base_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  
  
  
  metadata_list <- list(SFO=metadata_SFO,
                        AGRP=metadata_AGRP,
                        Hypertonic=metadata_hypertonicsaline)
  
  
  similarity_df <- data.frame()
  
  similarity_dist <- list(SFO=list(),
                          AGRP=list(),
                          Hypertonic=list())
  
  
  first_similarity_dist <- list(SFO=list(),
                                AGRP=list(),
                                Hypertonic=list())
  
  second_similarity_dist <- list(SFO=list(),
                                 AGRP=list(),
                                 Hypertonic=list())
  
  
  for(metadata_name in names(metadata_list)){
    
    half_similarity <- c()
    metadata <- metadata_list[[metadata_name]]
    
    print("------------")
    print("------------")
    print("------------")
    print("------------")
    
    
    hits_dec_all <- c()
    
    for (i in 1:len(metadata)) {
      
      stim_master_mat <- metadata[[i]]$stim_master_mat
      relevant_trials <- which(stim_master_mat[,"TrialType"] %in% (c(1,3,4,5)))
      trials_mat <- metadata[[i]]$trials_mat
      reward_trials <- which(metadata[[i]]$annot_df[,1] == 3 & metadata[[i]]$annot_df[,2] == 0)
      
      
      
      if (metadata_name == "SFO") {
        first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],15] %in% SFO_chunks[[as.character(i)]]]
        second_half_trials <- reward_trials[stim_master_mat[relevant_trials[reward_trials],15] %in% SFO_chunks[[as.character(i)]]]
        print("Using this SFO")
        print(SFO_chunks[[as.character(i)]])
      } else if (metadata_name == "AGRP") {
        first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],9] %in% AGRP_chunks[[as.character(i)]]]
        second_half_trials <- reward_trials[stim_master_mat[relevant_trials[reward_trials],9] %in% AGRP_chunks[[as.character(i)]]]
        print("Using this AGRP")
        print(AGRP_chunks[[as.character(i)]])
      } else {
        first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],15] %in% hypertonic_chunks[[as.character(i)]]]
        second_half_trials <- reward_trials[stim_master_mat[relevant_trials[reward_trials],15] %in% hypertonic_chunks[[as.character(i)]]]
        print("Using this hypertonic")
        print(SFO_chunks[[as.character(i)]])
      }
      
      
      first_half_mat <- trials_mat[first_half_trials,]
      second_half_mat <- trials_mat[second_half_trials,]
      
      
      similarity_matrix <- cosine(t(rbind(first_half_mat, second_half_mat)))
      
      first_to_second_sim <- similarity_matrix[1:len(first_half_trials),(len(first_half_trials)):nrow(similarity_matrix)]
      second_to_first_sim <- similarity_matrix[(len(first_half_trials)):nrow(similarity_matrix),1:len(first_half_trials)]
      
      for (main_diag_i in 1:nrow(similarity_matrix)) {similarity_matrix[main_diag_i,main_diag_i] <- NA}
      
      mean_sim <- median(c(c(first_to_second_sim), c(second_to_first_sim)), na.rm=T)
      
      similarity_dist[[metadata_name]] = append(similarity_dist[[metadata_name]], list(similarity_matrix))
      first_similarity_dist[[metadata_name]] = append(first_similarity_dist[[metadata_name]], list(similarity_matrix[1:len(first_half_trials),
                                                                                                                     1:len(first_half_trials)]))
      second_similarity_dist[[metadata_name]] = append(second_similarity_dist[[metadata_name]], list(similarity_matrix[1:len(first_half_trials), 
                                                                                                                       (len(first_half_trials) + 1):nrow(similarity_matrix)]))
      
      print(sprintf("all %f %f %d", median(similarity_matrix, na.rm=T), mean_sim, i))
      
      half_similarity <- c(half_similarity,
                           mean_sim)
      
    }
    
    similarity_df <- rbind(similarity_df,
                           data.frame(similarity=half_similarity,
                                      x=rep(metadata_name, times=len(half_similarity))))
  }
  
  
  fdf <- 
    lapply(names(similarity_dist),
           function(nm) {
             
             pooled_trial_sim_dist <- 
               unlist(lapply(similarity_dist[[nm]], c))
             
             h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=50), plot=F) 

             return(data.frame(Frac=c(0,h$counts / sum(h$counts)), Breaks=h$breaks, Group=rep(nm, times=len(h$breaks))))
           })
  
  hist_df <- do.call(rbind, fdf)
  
  first_half_fdf <- lapply(names(first_similarity_dist),
                           function(nm) {
                             
                             pooled_trial_sim_dist <- 
                               unlist(lapply(first_similarity_dist[[nm]], c))
                             
                             h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=20), plot=F) 
                             
                             return(data.frame(Frac=c(0,h$counts / sum(h$counts)), Breaks=h$breaks, Group=rep(nm, times=len(h$breaks))))
                           })
  
  first_hist_df <- do.call(rbind, first_half_fdf)
  
  
  second_half_fdf <- lapply(names(second_similarity_dist),
                            function(nm) {
                              
                              pooled_trial_sim_dist <- 
                                unlist(lapply(second_similarity_dist[[nm]], c))
                              
                              h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=20), plot=F) 
                              
                              return(data.frame(Frac=c(0,h$counts  / sum(h$counts)), Breaks=h$breaks, Group=rep(nm, times=len(h$breaks))))
                            })
  
  second_hist_df <- do.call(rbind, second_half_fdf)
  
  
  manip_names <- c("SFO",
                   "AGRP",
                   "Hypertonic")
  
  all_vals <- 
  lapply(manip_names, function(manipulation_name) {c(unlist(lapply(second_similarity_dist[[manipulation_name]], mean, na.rm=T)), 
                                                     unlist(lapply(first_similarity_dist[[manipulation_name]], mean, na.rm=T)))})
  
  all_vals <- unlist(all_vals)
  
  statistics_df <- data.frame()
  for (manipulation_name in manip_names) {
  
    pre_to_post <- second_hist_df[second_hist_df$Group == manipulation_name,]
    pre_to_pre <- first_hist_df[first_hist_df$Group == manipulation_name,]
    pre_to_pre$Group <- "Pre vs Pre"
    pre_to_post$Group <- "Pre vs Post"
    
    comp_df <- rbind(pre_to_pre, pre_to_post)
    gsim <- 
    ggplot(comp_df, aes(x=Breaks, y=Frac)) +
      geom_line(aes(color=Group),  size=.85, alpha=.85) + 
      #geom_ribbon(aes(fill=group),  color=NA, alpha=.25, stat="summary") +
      big_text_base_plot_theme + 
      scale_y_continuous(expand=c(0,0)) + 
      xlab("Cosine similarity")  +
      ylab("Fraction")
    
    
    tmp_df <- 
      cbind(unlist(lapply(second_similarity_dist[[manipulation_name]], mean, na.rm=T)), 
            unlist(lapply(first_similarity_dist[[manipulation_name]], mean, na.rm=T)))
    colnames(tmp_df) <- c("Pre vs Pre", "Pre vs Post")
    
    melted <- melt(tmp_df)
    colnames(melted) <- c("#", "Group", "Sim")
    
    
    ## Statistics
    
    wilc <- t.test(tmp_df[,1], tmp_df[,2], correct = F, paired=T, alternative="two.sided")

    
    stat_df <-  
      data.frame(Comp=manipulation_name,
                 Statistics=wilc$statistic,
                 Pvalue=wilc$p.value,
                 Signif=signif.num(wilc$p.value),
                 Method=wilc$method,
                 Alternative=wilc$alternative,
                 Mean_prepre=mean(tmp_df[,1]),
                 Mean_prepost=mean(tmp_df[,2]),
                 Median_prepre=median(tmp_df[,1]),
                 Median_prepost=median(tmp_df[,2]),
                 SD_prepre=sd(tmp_df[,1]),
                 SD_prepost=sd(tmp_df[,2]),
                 SEM_prepre=sem(tmp_df[,1]),
                 SEM_prepost=sem(tmp_df[,2]),
                 N_prepre=len(tmp_df[,1]),
                 N_prepost=len(tmp_df[,2]))
    
    statistics_df <- rbind(statistics_df,
                           stat_df)
    
    gbox <- 
      ggplot(melted) +
      geom_boxplot(aes(x=Group,y=Sim), fill=NA, width=.5, outlier.shape=NA) 

    
    
    for (row_idx in 1:nrow(tmp_df)) {
      tmp_df[row_idx,]
      row_df <- data.frame(x=c(1,2),
                           y=tmp_df[row_idx,])  
      
      gbox <- gbox + geom_line(data=row_df, aes(x=x, y=y))
    }
    
    gboxf <- 
      gbox + 
      geom_point(aes(x=Group,y=Sim), size=1.25, position=position_dodge(1)) +
      theme_light() + 
      big_text_base_plot_theme + 
      xlab("") +
      ylim(min(all_vals),max(all_vals)) +
      ylab("Cosine similarity")
    
    
    combined_plot <- 
      plot_grid(gsim, gboxf, rel_widths=c(2,1), nrow=1)
      
   
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_trial_similarity.pdf",base_figure_path, size_name, manipulation_name), height=size, width=size)
      plot(gsim)
      dev.off()
      
      pdf(sprintf("%s\\%s_%s_trial_similarity_boxplot_comp.pdf",base_figure_path, size_name, manipulation_name), height=size, width=size)
      plot(gboxf)
      dev.off()
      
      pdf(sprintf("%s\\%s_%s_combined_trial_similarity_boxplot_comp.pdf",base_figure_path, size_name, manipulation_name), height=size/1.5, width=size)
      plot(combined_plot)
      dev.off()
    }
    
  }
  
  
  write.csv(file=sprintf("%s//ttest_statistics.csv",base_figure_path), statistics_df)
  
  
}

plot_AGRP_SFO_hypertonic_decoding <- function()
{
  
  
  base_figure_path <- sprintf("%s\\AGRP_SFO_hypertonic_decoding", output_path)
  dir.create(base_figure_path)

  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_pre_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_pre_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_pre_pvalues_decoding.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_post_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_post_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\hunger_AGRP_post_pvalues_decoding.Rda", base_output_path), verbose=T)
  
  hunger_AGRP_pre = pre_manipulation_confusion_matrix_list
  hunger_AGRP_post =  post_manipulation_confusion_matrix_list
  hunger_AGRP_pre_shuffle = pre_manipulation_shuffle_matrix_list
  hunger_AGRP_post_shuffle = post_manipulation_shuffle_matrix_list
  
  
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_pre_pvalues_decoding.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_SFO_post_pvalues_decoding.Rda", base_output_path), verbose=T)
  
  thirst_SFO_pre = pre_manipulation_confusion_matrix_list
  thirst_SFO_post =  post_manipulation_confusion_matrix_list
  thirst_SFO_pre_shuffle = pre_manipulation_shuffle_matrix_list
  thirst_SFO_post_shuffle = post_manipulation_shuffle_matrix_list
  
  
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_pre_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_pre_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_pre_pvalues_decoding.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_post_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_post_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_5\\thirst_hypertonic_post_pvalues_decoding.Rda", base_output_path), verbose=T)

  
  thirst_hypertonic_pre = pre_manipulation_confusion_matrix_list
  thirst_hypertonic_post =  post_manipulation_confusion_matrix_list
  thirst_hypertonic_pre_shuffle = pre_manipulation_shuffle_matrix_list
  thirst_hypertonic_post_shuffle = post_manipulation_shuffle_matrix_list
  
  
  
  decoding_results = list(AGRP=list(pre=hunger_AGRP_pre,
                                    post=hunger_AGRP_post,
                                    pre_shuffle=hunger_AGRP_pre_shuffle,
                                    post_shuffle=hunger_AGRP_post_shuffle),
                          SFO=list(pre=thirst_SFO_pre,
                                   post=thirst_SFO_post,
                                   pre_shuffle=thirst_SFO_pre_shuffle,
                                   post_shuffle=thirst_SFO_post_shuffle),
                          hypertonic=list(pre=thirst_hypertonic_pre,
                                          post=thirst_hypertonic_post,
                                          pre_shuffle=thirst_hypertonic_pre_shuffle,
                                          post_shuffle=thirst_hypertonic_post_shuffle))
  
  
  statistics_df <- data.frame()
  values_df <- data.frame()
  
  for (manipulation_name in names(decoding_results)) {
    
    manip_decoding <- decoding_results[[manipulation_name]]
    
    
    pre <- rowMeans(manip_decoding$pre$Hit_cosine)
    post <- rowMeans(manip_decoding$post$Hit_cosine)
    pre_shuffle <- rowMeans(manip_decoding$pre_shuffle$Hit_cosine)
    post_shuffle <- rowMeans(manip_decoding$post_shuffle$Hit_cosine)
    
    pre_df <- data.frame(Decoding=pre, Group=rep("pre", times=len(pre)))
    post_df <- data.frame(Decoding=post, Group=rep("post", times=len(post)))
    pre_shuffle_df <- data.frame(Decoding=pre_shuffle, Group=rep("pre_shuffle", times=len(pre_shuffle)))
    post_shuffle_df <- data.frame(Decoding=post_shuffle, Group=rep("post_shuffle", times=len(post_shuffle)))
    
    df_final <- rbind(pre_df,
                      post_df,
                      pre_shuffle_df,
                      post_shuffle_df)
    
    plot_groups <- c("pre", "post", "pre_shuffle", "post_shuffle")
    gbox <- 
      ggplot(df_final, 
             aes(x=factor(Group, levels=plot_groups), 
                 y=Decoding), 
             group=Group) +
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA)
    
    
    for (row_idx in 1:nrow(pre_df)) {

      row_df <- data.frame(x=c(1,2), y=c(pre_df[row_idx,"Decoding"], post_df[row_idx,"Decoding"]))
      shuff_row_df <- data.frame(x=c(3,4), y=c(pre_shuffle_df[row_idx,"Decoding"], post_shuffle_df[row_idx,"Decoding"]))
      gbox <- gbox + geom_line(data=row_df, aes(x=x, y=y))
      gbox <- gbox + geom_line(data=shuff_row_df, aes(x=x, y=y))
    }
    
    gboxf <- 
      gbox +
      geom_point() + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Decoding accuracy (%)") + 
      ylim(0,1) +
      xlab("")  +
      theme(plot.title = element_text(size=9))
  
    
    dt <- dunn.test::dunn.test(df_final$Decoding, df_final$Group, "bonferroni")
    
    stat_df <- as.data.frame(dt)
    stat_df$SignifCorrected <- signif.num(stat_df$P.adjusted)
    stat_df$Signif <- signif.num(stat_df$P.adjusted)
    stat_df$Group <- manipulation_name
    
    wilc <- wilcox.test(pre_df[,1], post_df[,1], correct = F, paired=T, alternative="greater")
    wilc_shuff <- wilcox.test(pre_shuffle_df[,1], 
                         post_shuffle_df[,1], 
                         correct = F, paired=T, alternative="two.sided")
    
    
    wilc_pre <- wilcox.test(pre_df[,1], pre_shuffle_df[,1], correct = F, paired=T, alternative="greater")
    wilc_post <- wilcox.test(post_df[,1], post_shuffle_df[,1], correct = F, paired=T, alternative="greater")
  
    
    
    val_df <-  
      data.frame(Comp=manipulation_name,
                 Pvalue=wilc$p.value,
                 Method=wilc$method,
                 Statistic=wilc$statistic,
                 Alternative=wilc$alternative,
                 DatName=wilc$data.name,
                 Pvalue_shuffle=wilc_shuff$p.value,
                 Method_shuffle=wilc_shuff$method,
                 Statistic_shuffle=wilc_shuff$statistic,
                 Alternative_shuffle=wilc_shuff$alternative,
                 DatName_shuffle=wilc_shuff$data.name,
                 
                 Pvalue_pre=wilc_pre$p.value,
                 Method_pre=wilc_pre$method,
                 Statistic_pre=wilc_pre$statistic,
                 Alternative_pre=wilc_pre$alternative,
                 DatName_pre=wilc_pre$data.name,
                 
                 Pvalue_post=wilc_post$p.value,
                 Method_post=wilc_post$method,
                 Statistic_post=wilc_post$statistic,
                 Alternative_post=wilc_post$alternative,
                 DatName_post=wilc_post$data.name,
                 
                 Mean_pre=mean(pre_df[,1]),
                 Mean_post=mean(post_df[,1]),
                 Median_pre=median(pre_df[,1]),
                 Median_post=median(post_df[,1]),
                 SD_pre=sd(pre_df[,1]),
                 SD_post=sd(post_df[,1]),
                 SEM_pre=sem(pre_df[,1]),
                 SEM_post=sem(post_df[,1]),
                 N_pre=len(pre_df[,1]),
                 N_post=len(post_df[,1]),
                 Mean_pre_shuffle=mean(pre_shuffle_df[,1]),
                 Mean_post_shuffle=mean(post_shuffle_df[,1]),
                 Median_pre_shuffle=median(pre_shuffle_df[,1]),
                 Median_post_shuffle=median(post_shuffle_df[,1]),
                 SD_pre_shuffle=sd(pre_shuffle_df[,1]),
                 SD_post_shuffle=sd(post_shuffle_df[,1]),
                 SEM_pre_shuffle=sem(pre_shuffle_df[,1]),
                 SEM_post_shuffle=sem(post_shuffle_df[,1]),
                 N_pre_shuffle=len(pre_shuffle_df[,1]),
                 N_post_shuffle=len(post_shuffle_df[,1]))
    
    statistics_df <- rbind(statistics_df, stat_df)
    values_df <- rbind(values_df, val_df)

    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_box.pdf",base_figure_path, size_name, manipulation_name), height=size, width=size)
      plot(gboxf)
      dev.off()

    }
    }
    

  values_df$signif = signif.num(values_df$Pvalue)
  values_df$signif_shuff = signif.num(values_df$Pvalue_shuffle)
  values_df$signif_pre =signif.num(values_df$Pvalue_pre)
  values_df$signif_post =signif.num(values_df$Pvalue_post)
  
  values_df$corrected_signif = signif.num(values_df$Pvalue * 4)
  values_df$corrected_signif_shuff = signif.num(values_df$Pvalue_shuffle * 4)
  values_df$corrected_signif_pre =signif.num(values_df$Pvalue_pre * 4)
  values_df$corrected_signif_post =signif.num(values_df$Pvalue_post * 4)
  
  
  write.csv(file=sprintf("%s//statistics.csv",base_figure_path), statistics_df)
  write.csv(file=sprintf("%s//values.csv",base_figure_path), t(values_df))
  
  
}
