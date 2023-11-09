figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
metadata <- across_mice_decoding_build_metadata()

spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
bupu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "BuPu")))
ylgn_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "YlGn")))
blues_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "Blues")))
rdylbu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "RdYlBu")))


combined_max <- function(vec, combs = c(1:2)) {
  vec <- c(vec)
  ret <- max(unlist(lapply(combs,
                              function(nc) {
                                return(max(colSums(combn(vec, nc))))
                              })))
  
  return(ret)
}
cross_tabulate <- function(vec1, vec2) {
  
  if (len(unique(vec2)) >= len(unique(vec1))) {
    tmp = vec1
    vec1 <- vec2
    vec2 <- tmp
  }
  
  df <- data.frame(var1=vec1,var2=vec2)
  df_f <- df %>%
    group_by(var1, var2) %>%
    tally() %>%
    spread(var1, n)  
  df_f <- as.data.frame(df_f)
  df_f[is.na(df_f)] <- 0
  colnames(df_f) <- df_f[,1]
  df_f <- df_f[,-1]
  
  
  return(df_f)
}
calculate_MI_from_cont <- function(tabulated) {

  joint_prob <-  tabulated / sum(tabulated) 
  
  X <- colSums(joint_prob) / sum(joint_prob)
  Y <- rowSums(joint_prob) / sum(joint_prob)
  
  
  
  sum_all <- c()
  
  for (i in 1:nrow(joint_prob)) {
    for (j in 1:ncol(joint_prob)) {
      sum_all <- c(sum_all, 
                   joint_prob[i,j] * log2(((joint_prob[i,j])/(Y[i] * X[j])) + 10 ^ -30))
    }
  }
   
  return(sum(sum_all))
}


spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))


figure_3_colored_dimensionality_dataset <- function() {
  
  write_path <- sprintf("%s\\figure_3\\", figures_base_path)
  dir.create(write_path)
  
  subsample=F
  
  if (subsample) {
    write_path <- sprintf("%s\\structure_colored_subsampled", write_path)
    dir.create(write_path)  
  } else {
    write_path <- sprintf("%s\\structure_colored_new", write_path)
    dir.create(write_path)  
  }
  
  sizes = list(big=c(width=2.5),
               big_1=c(width=2.2),
               medium=c(width=2),
               medium_2=c(width=1.5),
               medium_1=c(width=1.25),
               small=c(width=1))
  
  structure_sizes <-  list(#medium_1=c(width=1,height=1),
                                  small=c(width=.75, height=.75))
  
  paths <- rep(get_thirsty_quenched_paths(), each=1)
  
  #datasets_names <- paste(rep(c("","SFO_"), times=4), get_datasets_names(paths, sep = "_"), sep="")
  datasets_names <- paste(rep("FULL_", times=4), get_datasets_names(paths, sep="_"), sep="")
  Responses = c("0"="Hit", "1"="Miss", "2"="CR (N)", "3"="FA (N)", "4"="CR (A)", "5"="FA (A)",  "15"="Pupil", "16"="Consumed water")
  frame_rate=30
  binning_window_size=15
  #max_windows = c(0:5)
  max_windows <- c(3)
  begin_offset <- (c(6))
  plot_structure = F
  pupil_q <- .15
  thirsty_q <- .15
  cont_breaks <- 15
  n_subsamples = 100
  
  ylabs <- list(occ_no_main="Occupancy (%)",
                occ_main="Occupancy (%)",
                var_prob="Probability (%)",
                var_prob_main="Probability (%)",
                entropy_no_main="Entropy",
                entropy_main="Entropy",
                MI_main="Mutual information",
                MI_no_main="Mutual information")
  
  all_results_dfs <- list()
  #metadata <- metadata_SFO
  metadata <- SFO_paths
  plot_structure = T
  
  for (offset in rev(begin_offset)) {
  for (max_window_size in rev(max_windows)){
    
      
      # if (offset >= max_window_size) {
      #   next
      # }
        
      #all_windows <- seq(1, (max_window_size - offset), by=.5)
      entropy_window_size = max_window_size - offset
      occ_main_all <- c()
      occ_all <- c()
      occ_of_var_main_all <- c()
      occ_of_var_all <- c()
      # entropy_of_var_main_all <- c()
      # entropy_of_var_all <- c()
      entropy_main_all <- c()
      entropy_all <- c()
      MI_all <- c()
      MI_main_all <- c()
      
      for (path_idx in 1:len(metadata)) {
        
        if (path_idx == 4 || path_idx == 14 || path_idx == 7) {
           next
        }
         
        dataset_nm <- datasets_names[[path_idx]]
        reduced_mat <- metadata[[path_idx]]$red_mat
        work_stim_mat <- metadata[[path_idx]]$stim_master_mat
        
        base <- rep(0, times=nrow(reduced_mat))
        base_color <- rep(adjustcolor("gray70", alpha=0.3), times=nrow(reduced_mat))
        
        cluster_map <- metadata[[path_idx]]$cluster_mat
        cluster_labels <- cluster_map$labs
        no_main_indices <- cluster_labels != -1
        cluster_labels_no_main <- cluster_labels[no_main_indices]
        all_clusters <- sort(unique(cluster_labels))
        cluster_prob <- table(cluster_labels_no_main) / len(cluster_labels_no_main)
        cluster_prob_with_main <- table(cluster_labels) / len(cluster_labels)
        
        pupil_vec_all <- get_pupil_files(paths[[path_idx]], window_size=15)
        
        drink_ind <- work_stim_mat[!is.nan(work_stim_mat[,6]),"Frames"]
        drink_ind <- unique(get_binned_index(drink_ind, binning_window_size))
        consumed_water_vec <- base
        consumed_water_vec[drink_ind] <- 1
        consumed_water_vec <- cumsum(consumed_water_vec)

        ongoing_trials <- which(metadata[[path_idx]]$stim_master_mat[,7] %in% c(4,5) & metadata[[path_idx]]$stim_master_mat[,8]%%2 == 0)
        ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[path_idx]]$stim_master_mat)]
        ongoing_trials_frames = metadata[[path_idx]]$stim_master_mat[ongoing_trials + 1,1]
        ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - max_window_size * 2):(get_binned_index(i, 15)-1)})))
        ongoing_ind_logic <- 1:nrow(reduced_mat) %in% ongoing_ind
        
        indices_all <- list()
        
        for(resp in as.numeric(names(Responses))) {
          resp_name <- Responses[as.character(resp)]
          
          indices_all[[resp_name]] <- base

          if (resp < 15) {

            trials <- work_stim_mat[work_stim_mat[,"Response"] == resp,"Frames"]
            trial_indices <- lapply(trials,
                                    function(trial_start) 
                                    {(trial_start + offset * frame_rate):(trial_start + (offset + entropy_window_size) * frame_rate - 1)})
            
            
            trial_indices_all <- lapply(trial_indices,
                                        function(trial_ind) {
                                          return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, binning_window_size)})))
                                        })
            
            trial_indices_all <- unique(unlist(trial_indices_all))
            trial_indices_all <- trial_indices_all[trial_indices_all <= len(base)] 
            
            indices_all[[resp_name]][trial_indices_all] <- 1
          }
        }
        
        min_timepoints <- min(unlist(lapply(indices_all, sum)))
        
        max_occ <- c() # Fraction of cluster that is explained by event
        max_occ_with_main <- c()
        
        max_occ_of_var <- c() # Fraction of event that is explained by cluster
        max_occ_of_var_with_main <- c()
        
        entropy <- c() # Distribution of event in different clusters
        entropy_with_main <- c()
        
        MI <- c()
        MI_with_main <- c()
        
        
        for(resp in as.numeric(names(Responses))) {
          resp_name <- Responses[as.character(resp)]
          print(resp_name)
          # if (subsample) {
          # 
          #   subsample_max_occ <- c()
          #   subsample_max_occ_with_main <- c()
          #   subsample_max_occ_of_var <- c()
          #   subsample_max_occ_of_var_with_main <- c()
          #   subsample_entropy <- c()
          #   subsample_entropy_with_main <- c()
          #   # subsample_entropy_of_var <- c()
          #   # subsample_entropy_of_var_with_main <- c()
          #   subsample_MI <- c()
          #   subsample_MI_with_main <- c()
          # 
          #   original_indices <- which(indices_all[[resp_name]] == 1)
          #   subsample_ind <- sort(sample(original_indices, min_timepoints))
          #   print(sprintf("%f - %d", mean(subsample_ind), len(subsample_ind)))
          #   
          #   for (subsample_rep in 1:n_subsamples) {
          # 
          #     subsample_ind <- sort(sample(original_indices, min_timepoints))
          #     #print(sprintf("%f - %d", mean(subsample_ind), len(subsample_ind)))
          #     variable_indices <- base
          # 
          #     variable_indices[subsample_ind] <- 1
          #     variable_indices_no_main <- variable_indices[no_main_indices]
          #     var_contingencies <- table(cluster_labels[variable_indices == 1] )
          # 
          #     all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
          #     all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
          # 
          #     names(all_var_contingencies) <- all_clusters[-1]
          #     names(all_var_contingencies_with_main) <- all_clusters
          # 
          # 
          #     all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
          #     all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies
          # 
          #     var_probability <- all_var_contingencies / sum(all_var_contingencies)
          #     var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)
          # 
          #     clust_occ_with_main <- all_var_contingencies_with_main / table(cluster_labels)
          #     clust_occ <- all_var_contingencies / table(cluster_labels_no_main)
          #     #clust_occ_with_main <- clust_occ_with_main / (sum(clust_occ_with_main)  + 10^-30)
          #     #clust_occ <- clust_occ / (sum(clust_occ)   + 10^-30)
          #     
          #     
          #     tabulated_main <- cross_tabulate(cluster_labels,
          #                                      variable_indices)
          #     
          #     tabulated_no_main <- cross_tabulate(cluster_labels_no_main,
          #                                         variable_indices_no_main)
          #     
          #     subsample_MI <- c(subsample_MI, calculate_MI_from_cont(tabulated_no_main))
          #     subsample_MI_with_main <- c(subsample_MI_with_main, calculate_MI_from_cont(tabulated_main))
          # 
          #     subsample_max_occ <- c(subsample_max_occ,
          #                           max(clust_occ))
          #     subsample_max_occ_with_main <- c(subsample_max_occ_with_main,
          #                                      max(clust_occ_with_main))
          #     subsample_max_occ_of_var <- c(subsample_max_occ_of_var,
          #                                    max(var_probability))
          #     subsample_max_occ_of_var_with_main <- c(subsample_max_occ_of_var_with_main,
          #                                              max(var_probability_with_main))
          #     subsample_entropy <- c(subsample_entropy,
          #                            calculate_entropy(var_probability))
          #     subsample_entropy_with_main <- c(subsample_entropy_with_main,
          #                                     calculate_entropy(var_probability_with_main))
          #     # subsample_entropy_of_var <- c(subsample_entropy_of_var,
          #     #                               calculate_entropy(clust_occ))
          #     # subsample_entropy_of_var_with_main <- c(subsample_entropy_of_var_with_main,
          #     #                                         calculate_entropy(clust_occ_with_main))
          #   }
          # 
          #   print(sprintf("%s - %f - %f - %f - %f",
          #                 resp_name,
          #                 mean(subsample_max_occ),
          #                 median(subsample_max_occ),
          #                 min(subsample_max_occ),
          #                 max(subsample_max_occ)))
          #   
          #   
          #   print(sprintf("%s - %f - %f - %f - %f",
          #                 resp_name,
          #                 mean(subsample_MI),
          #                 median(subsample_MI),
          #                 min(subsample_MI),
          #                 max(subsample_MI)))
          #   max_occ <- c(max_occ, mean(subsample_max_occ))
          #   max_occ_with_main <- c(max_occ_with_main, mean(subsample_max_occ_with_main))
          #   max_occ_of_var <- c(max_occ_of_var, mean(subsample_max_occ_of_var))
          #   max_occ_of_var_with_main <- c(max_occ_of_var_with_main, mean(subsample_max_occ_of_var_with_main))
          #   entropy <- c(entropy, mean(subsample_entropy))
          #   entropy_with_main <- c(entropy_with_main, mean(subsample_entropy_with_main))
          #   # entropy_of_var <- c(entropy_of_var, mean(subsample_entropy_of_var))
          #   # entropy_of_var_with_main <- c(entropy_of_var_with_main, mean(subsample_entropy_of_var_with_main))
          #   
          #   MI <- c(MI, mean(subsample_MI))
          #   MI_with_main <- c(MI_with_main, mean(subsample_MI_with_main))
          # 
          # } else {
          #   
            if (resp >= 15) {
              
              if (resp_name == "Pupil") {
                vec_to_discretize <- pupil_vec_all$smoothed
              } else {
                vec_to_discretize <- consumed_water_vec
              }
              
              
              
              discretized_vec <- 
                as.numeric(cut(vec_to_discretize[ongoing_ind], 
                               breaks = seq(min(vec_to_discretize[ongoing_ind] - 1), 
                                            max(vec_to_discretize[ongoing_ind]), 
                                            length.out=cont_breaks + 1)))
              
              discretized_vec_no_main <- 
                as.numeric(cut(vec_to_discretize[no_main_indices & ongoing_ind_logic], 
                               breaks = seq(min(vec_to_discretize[no_main_indices & ongoing_ind_logic] - 1), 
                                            max(vec_to_discretize[no_main_indices & ongoing_ind_logic]), 
                                            length.out=cont_breaks + 1)))
              
              
              working_cluster_labels <- rep(0, times=len(all_clusters) - 1)
              working_cluster_labels_with_main <- rep(0, times=len(all_clusters))
              names(working_cluster_labels) <- all_clusters[-1]
              names(working_cluster_labels_with_main) <- all_clusters
              
              tmp <- table(cluster_labels[ongoing_ind_logic & no_main_indices])
              working_cluster_labels[names(tmp)] <- tmp
              
              tmp <- table(cluster_labels[ongoing_ind])
              working_cluster_labels_with_main[names(tmp)] <- tmp
              

              parameterization_variables <- 
              lapply(1:cont_breaks,
                     function(b) {

                       var_contingencies <- table(cluster_labels[ongoing_ind][discretized_vec == b])
                       
                       all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
                       all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
                       

                       
                       names(all_var_contingencies) <- all_clusters[-1]
                       names(all_var_contingencies_with_main) <- all_clusters

                       
                       
                       all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
                       all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies
                       
                       clust_occ_with_main <- all_var_contingencies_with_main / working_cluster_labels_with_main
                       clust_occ <- all_var_contingencies / working_cluster_labels
                       var_probability <- all_var_contingencies / sum(all_var_contingencies)
                       var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)
                       
                       return(list(clust_occ_with_main=clust_occ_with_main,
                                   clust_occ=clust_occ,
                                   var_probability=var_probability,
                                   var_probability_with_main=var_probability_with_main))
                     })
              
              
              max_occ <- c(max_occ, 
                           max(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$clust_occ})), 2, median, na.rm=T)))
              
              max_occ_with_main <- c(max_occ_with_main, 
                                     max(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$clust_occ_with_main})), 2, median, na.rm=T)))
              
              max_occ_of_var <- c(max_occ_of_var, 
                                  max(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$var_probability})), 2, median, na.rm=T)))
              
              max_occ_of_var_with_main <- c(max_occ_of_var_with_main, 
                                            max(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$var_probability_with_main})), 2, median, na.rm=T)))
              
              
              
              tabulated_main <- cross_tabulate(discretized_vec,
                                               cluster_labels[ongoing_ind])
              
              tabulated_no_main <- cross_tabulate(discretized_vec_no_main,
                                                  cluster_labels[ongoing_ind_logic & no_main_indices])
              
              MI <- c(MI,
                      calculate_MI_from_cont(tabulated_no_main))
              
              MI_with_main <- c(MI_with_main,
                                calculate_MI_from_cont(tabulated_main))
              
              
            } else {
              
            variable_indices <- indices_all[[resp_name]]
            variable_indices_no_main <- variable_indices[no_main_indices]
            var_contingencies <- table(cluster_labels[variable_indices == 1] )

            all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
            all_var_contingencies_with_main <- rep(0, times=len(all_clusters))

            names(all_var_contingencies) <- all_clusters[-1]
            names(all_var_contingencies_with_main) <- all_clusters


            all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
            all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies

            var_probability <- all_var_contingencies / sum(all_var_contingencies)
            var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)

            clust_occ_with_main <- all_var_contingencies_with_main / table(cluster_labels)
            clust_occ <- all_var_contingencies / table(cluster_labels_no_main)
            #clust_occ_with_main <- clust_occ_with_main / (sum(clust_occ_with_main)  + 10^-30)
            #clust_occ <- clust_occ / (sum(clust_occ)   + 10^-30)

            max_occ <- c(max_occ, max(clust_occ))
            max_occ_with_main <- c(max_occ_with_main, max(clust_occ_with_main))

            max_occ_of_var <- c(max_occ_of_var, max(var_probability))
            max_occ_of_var_with_main <- c(max_occ_of_var_with_main, max(var_probability_with_main))

            entropy <- c(entropy,calculate_entropy(var_probability))
            entropy_with_main <- c(entropy_with_main, calculate_entropy(var_probability_with_main))


            # entropy_of_var <- c(entropy_of_var, calculate_entropy(clust_occ))
            # entropy_of_var_with_main <- c(entropy_of_var_with_main, calculate_entropy(clust_occ_with_main))

            
            tabulated_main <- cross_tabulate(cluster_labels,
                                             variable_indices)
            
            tabulated_no_main <- cross_tabulate(cluster_labels_no_main,
                                                variable_indices_no_main)
            
            MI <- c(MI,
                    calculate_MI_from_cont(tabulated_no_main))
            
            MI_with_main <- c(MI_with_main,
                              calculate_MI_from_cont(tabulated_main))
            }
          #}
          

        
        }
        
        occ_all <- rbind(occ_all,max_occ)
        occ_main_all <- rbind(occ_main_all,max_occ_with_main)
        
        occ_of_var_all <- rbind(occ_of_var_all,max_occ_of_var)
        occ_of_var_main_all <- rbind(occ_of_var_main_all,max_occ_of_var_with_main)

        entropy_all <- rbind(entropy_all, entropy)
        entropy_main_all <- rbind(entropy_main_all, entropy_with_main)
      
        
        MI_main_all <- rbind(MI_main_all, MI_with_main)
        MI_all <- rbind(MI_all, MI)
        # print(max_occ)
        # print(max_occ_of_var)
        # print(MI_with_main)
        print(MI_main_all)
      }
      
      
      all_results_dfs[[sprintf("%d_%d", max_window_size, offset)]] <- 
                    list(occ_no_main=occ_all,
                         occ_main=occ_main_all,
                         var_prob=occ_of_var_all,
                         var_prob_main=occ_of_var_main_all,
                         entropy_no_main=entropy_all,
                         entropy_main=entropy_main_all,
                         MI_main=MI_all,
                         MI_no_main=MI_main_all)
      
      
      # colnames(entropy_all) <- Responses
      # colnames(entropy_main_all) <- Responses
      # colnames(occ_all) <- Responses
      # colnames(occ_main_all) <- Responses
      # colnames(occ_of_var_all) <- Responses
      # colnames(occ_of_var_main_all)<- Responses
      # colnames(MI_all) <- Responses
      # colnames(MI_main_all) <- Responses
      # 
      #  melted_occ <- melt(occ_all)
      # melted_occ_main <- melt(occ_main_all)
      # melted_var_prob <- melt(occ_of_var_all)
      # melted_var_prob_main <- melt(occ_of_var_main_all)
      # melted_entropy <- melt(entropy_all)
      # melted_entropy_main <- melt(entropy_main_all)
      # melted_MI <- melt(MI_all)
      # melted_MI_main <- melt(MI_main_all)
      # 
      # melted_dfs <- list(occ_no_main=melted_occ,
      #                    occ_main=melted_occ_main,
      #                    var_prob=melted_var_prob,
      #                    var_prob_main=melted_var_prob_main,
      #                    entropy_no_main=melted_entropy,
      #                    entropy_main=melted_entropy_main,
      #                    MI_main=melted_MI,
      #                    MI_no_main=melted_MI_main)
      
  
      
      
      # for (melted_df_name in names(melted_dfs)) {
      #   
      #   melted_df <- melted_dfs[[melted_df_name]]
      #   colnames(melted_df) <- c("#", "Group", "Occ")
      #   
      #   gmelted <- 
      #     ggplot(melted_df, aes(y=Occ, x=Group)) +
      #     scale_y_continuous(expand=c(0,0)) +
      #     geom_bar(stat="summary", width=.45,
      #              fill="gray70") +
      #     geom_errorbar(stat="summary", width = .3) +
      #     geom_jitter(position=position_jitterdodge(.5), aes(fill=Group),
      #                 size=1.5) + 
      #     theme_classic() +
      #     theme(text=element_text(size=14, color="black"), 
      #           axis.text=element_text(size=14, color="black"),
      #           axis.ticks = element_line(color="black"),
      #           plot.title=element_text(size=10),
      #           legend.position="NA")  + 
      #     xlab("") +
      #     ylab(ylabs[[melted_df_name]])
      #   
      #   
      #   for (size_name in names(sizes)) {
      #     dir.create(sprintf("%s\\%s\\",
      #                        write_path,
      #                        size_name))
      #     dir.create(sprintf("%s\\%s\\%d_%d",
      #                        write_path,
      #                        size_name,
      #                        max_window_size,
      #                        offset))
      #     
      #     
      #     
      #     pdf(sprintf("%s\\%s\\%d_%d\\comp_%s.pdf",
      #                 write_path,
      #                 size_name,
      #                 max_window_size,
      #                 offset,
      #                 melted_df_name),
      #         height=sizes[[size_name]][["width"]],
      #         width=sizes[[size_name]][["width"]]) 
      #     
      #     plot(gmelted)
      #     dev.off()
      #     
      #     
      #   }
      # }
    }
  }
  
  # for (max_window_size in max_windows) {
  #   
  #   all_confs <- names(all_results_dfs)
  #   #all_confs <- all_confs[11:len(all_confs)]
  #   df_names <- names(all_results_dfs[[1]])
  # 
  # 
  #   relevant_confs <- all_confs[unlist(lapply(str_split(all_confs, "_"), 
  #                                             function(str) {
  #                                               str[[1]] == max_window_size #& str[[1]] > str[[2]]
  #                                               }))]
  #   relevant_confs <- sort(relevant_confs)
  #   
  #   
  #   relevant_dfs <- all_results_dfs[relevant_confs]
  #   
  # 
  #   for (df_name in df_names) {
  #     
  #     conf_df <- data.frame()
  #     
  # 
  #     for (conf_name in names(relevant_dfs)) {
  #       offset_to_use <- as.numeric(str_split(conf_name, "_")[[1]])[[2]]
  #       
  #       work_df <-relevant_dfs[[conf_name]][[df_name]]
  # 
  #       mean_df <- apply(work_df, 2, mean, na.rm=T)
  #       sd_df <- apply(work_df, 2, sem)
  #       
  #       fdf <- cbind(mean_df, sd_df)
  #       fdf <- as.data.frame(fdf)
  #       fdf <- cbind(fdf, Responses)
  #       fdf <- cbind(fdf, rep(offset_to_use, times=len(Responses)))
  #       colnames(fdf) <- c("Mean", "Sd", "Var", "Window")
  #       
  #       conf_df <- rbind(conf_df,
  #                        fdf)
  #       
  #     }
  #     
  #     gtime_all <- 
  #       ggplot(conf_df, aes(x=Window, y=Mean, group=Var)) + 
  #       geom_line(aes(color=Var), size=.8) + 
  #       geom_point(aes(color=Var), size=1.25) + 
  #       #geom_errorbar(aes(color=Var, ymin=Mean - Sd, ymax=Mean + Sd), width=.2, alpha=1) +
  #       geom_ribbon(aes(fill=Var, ymin=Mean - Sd, ymax=Mean + Sd), color=NA, alpha=.2) +
  #       theme_classic() +
  #       theme(text=element_text(size=14, color="black"), 
  #             axis.text=element_text(size=14, color="black"),
  #             axis.ticks = element_line(color="black"),
  #             plot.title=element_text(size=10),
  #             legend.position="NA")  + 
  #       xlab("Offset size (s)") +
  #       ylab(ylabs[[df_name]]) +
  #       ggtitle(sprintf("Offset: %d %s", offset, ifelse(grepl("main", df_name), "(Main)", ""))) + 
  #       geom_vline(xintercept=0, size=.6, color="gray60", linetype="dashed")
  #       
  #     
  #     panels <- list()
  #     ymax <- max(conf_df$Mean + conf_df$Sd)
  #     ymin <- min(conf_df$Mean - conf_df$Sd)
  #     
  #     variables <- unique(conf_df$Var)
  #     
  #     for (var_name in variables) {
  #       
  #       work_conf_df <- conf_df[conf_df$Var == var_name,]
  #       gtime_panel <- 
  #           ggplot(work_conf_df, aes(x=Window, y=Mean, group=Var)) + 
  #           geom_line(aes(color=Var), size=.8) + 
  #           geom_point(aes(color=Var), size=1.25) + 
  #           #geom_errorbar(aes(color=Var, ymin=Mean - Sd, ymax=Mean + Sd), width=.2, alpha=1) +
  #           geom_ribbon(aes(fill=Var, ymin=Mean - Sd, ymax=Mean + Sd), color=NA, alpha=.2) +
  #           theme_classic() +
  #           theme(text=element_text(size=14, color="black"), 
  #                 axis.text=element_text(size=14, color="black"),
  #                 axis.ticks = element_line(color="black"),
  #                 plot.title=element_text(size=10),
  #                 legend.position="NA")  + 
  #           xlab("Offset (s)") +
  #           ylab(ylabs[[df_name]]) +
  #           #ggtitle(sprintf("%s %s", var_name, ifelse(grepl("main", df_name), "(Main)", ""))) +
  #           ggtitle(var_name) +
  #           ylim(ymin, ymax) + 
  #           geom_vline(xintercept=0, size=.6, color="gray60", linetype="dashed")
  #       
  #       panels <- append(panels, list(gtime_panel))
  #       
  #     }
  #     
  #     panels$nrow <- 2
  #     panels_plot <- do.call(plot_grid, panels)
  #     
  #     for (size_name in names(sizes)) {
  #       dir.create(sprintf("%s\\%s\\",
  #                          write_path,
  #                          size_name))
  #       
  #       dir.create(sprintf("%s\\%s\\window_%d",
  #                          write_path,
  #                          size_name,
  #                          max_window_size))
  #       
  #       
  #       
  #       pdf(sprintf("%s\\%s\\window_%d\\comp_by_time_%s.pdf",
  #                   write_path,
  #                   size_name,
  #                   max_window_size,
  #                   df_name),
  #           height=sizes[[size_name]][["width"]],
  #           width=sizes[[size_name]][["width"]]) 
  #       
  #       plot(gtime_all)
  #       dev.off()
  #       
  #       
  #       pdf(sprintf("%s\\%s\\window_%d\\panels_comp_by_time_%s.pdf",
  #                   write_path,
  #                   size_name,
  #                   max_window_size,
  #                   df_name),
  #           height=sizes[[size_name]][["width"]] * 1.5,
  #           width=sizes[[size_name]][["width"]] * 3) 
  #       
  #       plot(panels_plot)
  #       dev.off()
  #     }
  #   }
  # }
  
  save(file=sprintf("%s//results_df.Rda", write_path), all_results_dfs)
}


figure_3_entropy_maps <- function() {
  write_path <- sprintf("%s\\figure_3\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\entropy", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=2.5,
                     height=2.5),
               big_1=c(width=2.25,
                       height=2.25),
               medium=c(width=2,
                        height=2),
               medium_2=c(width=1.75,
                          height=1.75),
               medium_1=c(width=1.5,
                          height=1.5),
               small=c(width=1.25,
                       height=1.25))
  paths <- get_thirsty_quenched_paths()
  
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  calculate_entropy <- function(pmf) {- sum(pmf * log2(pmf + 10 ** -30))}
  

  
  trial_slices <- list(c(begin=0, end=1), c(begin=0, end=2), c(begin=0, end=3), c(begin=0, end=4), c(begin=0, end=5), c(begin=0, end=6), c(begin=0, end=7), c(begin=0, end=8), c(begin=0, end=9),
                       c(begin=1, end=2), c(begin=1, end=3), c(begin=1, end=4), c(begin=1, end=5), c(begin=1, end=6), c(begin=1, end=7), c(begin=1, end=8), c(begin=1, end=9),
                       c(begin=2, end=3), c(begin=2, end=4), c(begin=2, end=5), c(begin=2, end=6), c(begin=2, end=7), c(begin=2, end=8), c(begin=2, end=9),
                       c(begin=3, end=4), c(begin=3, end=5), c(begin=3, end=6), c(begin=3, end=7), c(begin=3, end=8), c(begin=3, end=9))

  for (trial_slice in trial_slices) {
    
    trial_begin  = trial_slice["begin"]
    trial_end  = trial_slice["end"]
    
    entropy_df <- c()
    posterior_entropy_df <- c()
    entropy_df_with_main <- c()
    posterior_entropy_df_with_main <- c()
    
    for (i in 1:len(paths)) {
      
      if (i == 4 || i == 14) {
        next
      }
      
      dataset_nm <- datasets_names[[i]]
      reduced_mat <- get_mat_with_preset(paths[[i]], "dalshtaim")  
      cluster_map <- get_clusters_mat_kmeans(paths[[i]])
      variables <- get_behavior_variables(paths[[i]], begin=trial_begin, end=trial_end)
      cluster_labels <- cluster_map$labs
      no_main_indices <- cluster_labels != -1
      cluster_labels_no_main <- cluster_labels[no_main_indices]
      all_clusters <- sort(unique(cluster_labels))
      
      cluster_prob <- table(cluster_labels_no_main) / len(cluster_labels_no_main)
      cluster_prob_with_main <- table(cluster_labels) / len(cluster_labels)
      
      # Disjoint pupil vector into two vectors
      pupil_vec <- variables$variables$Pupil
      low_pupil_vec <- as.numeric(pupil_vec == 1)
      high_pupil_vec <- as.numeric(pupil_vec == 2)
      
      #names(variables$variables)[which(names(variables$variables) == "Pupil")] <- "LowPupil"
      variables$variables$LowPupil <- low_pupil_vec
      variables$variables$HighPupil <- high_pupil_vec
      
      entropy_of_variables <- c()
      entropy_of_variables_with_main <- c()
      posterior_entropy_of_variables <- c()
      posterior_entropy_of_variables_with_main <- c()
      
      
      all_var_names <- names(variables$variables)
      all_var_names <- all_var_names[all_var_names != "Pupil"]
      
      var_prob_plot_lists <- list(no_main=list(),
                                  main=list(),
                                  posterior_no_main=list(),
                                  posterior_main=list())
      for (var_name in all_var_names) {
        # Var was disjointed
        if (var_name == "Pupil") {
          next
        }
        
        variable_indices <- variables$variables[[var_name]]
        variable_indices_no_main <- variable_indices[no_main_indices]
        var_contingencies <- table(cluster_labels[variable_indices == 1] )
        all_var_prob_pmf <- table(variable_indices) 
        all_var_prob_no_main_pmf <- table(variable_indices_no_main)
        all_var_prob   <- all_var_prob_pmf[2] / sum(all_var_prob_pmf)
        all_var_prob_no_main   <- all_var_prob_no_main_pmf[2] / sum(all_var_prob_no_main_pmf)
      
        
        
        # SANITY CHECK
        assert(all(table(cluster_labels_no_main[variable_indices_no_main == 1]) == 
                     var_contingencies[-1]))
        
        all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
        all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
        
        names(all_var_contingencies) <- all_clusters[-1]
        names(all_var_contingencies_with_main) <- all_clusters
        
        likelihood_prob <- lapply(all_clusters[-1], function(c) {table(variable_indices_no_main[which(cluster_labels_no_main == c)])})
        likelihood <- unlist(lapply(likelihood_prob, 
                                     function(v) {
                                       if (all(names(v) == "0")) {return(0)} 
                                       else if (all(names(v) == "1")) {return(1)}
                                       else {return(v[2] / sum(v))}
                                     }))
        names(likelihood) <- all_clusters[-1]
        
        likelihood_prob_with_main <- lapply(all_clusters, function(c) {table(variable_indices[which(cluster_labels == c)])})
        likelihood_with_main <- unlist(lapply(likelihood_prob_with_main, 
                                              function(v) {
                                                if (all(names(v) == "0")) {return(0)} 
                                                else if (all(names(v) == "1")) {return(1)}
                                                else {return(v[2] / sum(v))}
                                              }))
        
        for (clust in all_clusters) {
          
          # In case there werent any frames of variables in that perticular cluster
          if(!as.character(clust) %in% names(var_contingencies)) {
            next
          }
          
          if (clust == -1) {
            all_var_contingencies_with_main[as.character(clust)]  <- 
              var_contingencies[as.character(clust)]
          } else {
            all_var_contingencies_with_main[as.character(clust)]  <- 
              var_contingencies[as.character(clust)]
            
            all_var_contingencies[as.character(clust)]  <- 
              var_contingencies[as.character(clust)]
          }
        }
        
        var_probability <- all_var_contingencies / sum(all_var_contingencies)
        var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)  
        posterior_probability <- (likelihood * cluster_prob) / all_var_prob_no_main
        posterior_probability_with_main <- (likelihood_with_main * cluster_prob_with_main) / all_var_prob
        
        
        names(posterior_probability) <- all_clusters[-1]
        names(posterior_probability_with_main) <- all_clusters
        
        entropy_of_variables <- c(entropy_of_variables,
                                  calculate_entropy(var_probability))
        
        entropy_of_variables_with_main <- c(entropy_of_variables_with_main,
                                            calculate_entropy(var_probability_with_main))
        
        posterior_entropy_of_variables <- c(posterior_entropy_of_variables,
                                            calculate_entropy(posterior_probability))
        
        posterior_entropy_of_variables_with_main <- c(posterior_entropy_of_variables_with_main,
                                                      calculate_entropy(posterior_probability_with_main))
        
        
        prob_list <- list(no_main=list(prob_vec=var_probability,
                                         clust=as.numeric(names(var_probability))),
                          main=list(prob_vec=var_probability_with_main,
                                    clust=as.numeric(names(var_probability_with_main))),
                          posterior_no_main=list(prob_vec=posterior_probability,
                                                 clust=as.numeric(names(var_probability))),
                          posterior_main=list(prob_vec=posterior_probability_with_main,
                                              clust=as.numeric(names(posterior_probability_with_main))))
        
        for (prob_name in names(prob_list)) {
          
          
          if(!var_name %in% all_var_names[c(1,2,3,4,5,6,11,12)]) {
            next
          }
          
          
          #tmp_prob_v <- 
          #names(tmp_prob_v) <- c()
          df_to_use <- data.frame(prob=as.numeric(prob_list[[prob_name]]$prob_vec),
                                  x=prob_list[[prob_name]]$clust)
          colnames(df_to_use) <- c("prob", "x")

          
          gvarprob <- 
            ggplot(df_to_use, aes(x=factor(x), y=prob),) +
            geom_bar(stat="summary", width=.45, fill="gray70") + 
              scale_y_continuous(expand=c(0,0))  +
              xlab("") +
              ylab("Probability") +
              theme(text=element_text(size=14, color="black"),
                    legend.position="NA",
                    axis.text=element_text(size=14, color="black"),
                    axis.ticks=element_line(color="black"))  +
              theme_classic() +
              ggtitle(var_name)
            
            
          var_prob_plot_lists[[prob_name]]  <- append(var_prob_plot_lists[[prob_name]],
                                                      list(gvarprob))
                  
        }
      }
      
      names(entropy_of_variables) <- all_var_names
      names(entropy_of_variables_with_main) <- all_var_names
      names(posterior_entropy_of_variables) <- all_var_names
      names(posterior_entropy_of_variables_with_main) <- all_var_names
      
      entropy_df <- rbind(entropy_df,
                          entropy_of_variables)
      
      entropy_df_with_main <- rbind(entropy_df_with_main,
                                    entropy_of_variables_with_main)
      
      posterior_entropy_df <- rbind(posterior_entropy_df,
                                    posterior_entropy_of_variables)
      
      posterior_entropy_df_with_main <- rbind(posterior_entropy_df_with_main,
                                              posterior_entropy_of_variables_with_main)
      
      
      
      for (plot_conf in names(var_prob_plot_lists)) {
        plt <- var_prob_plot_lists[[plot_conf]]
        plt$nrow <- 2
        gf <- do.call(plot_grid, plt)
        for (size_name in names(sizes)) {
          
          dir.create(sprintf("%s\\%s",
                             write_path,
                             size_name))
          dir.create(sprintf("%s\\%s\\%s\\",
                             write_path,
                             size_name,
                             plot_conf))
          
          dir.create(sprintf("%s\\%s\\%s\\var_probabilities",
                             write_path,
                             size_name,
                             plot_conf))
          
          dir.create(sprintf("%s\\%s\\%s\\var_probabilities\\%s",
                             write_path,
                             size_name,
                             plot_conf,
                             dataset_nm))
          
          
          pdf(sprintf("%s\\%s\\%s\\var_probabilities\\%s\\var_prob_begin%d_end%d.pdf",
                      write_path,
                      size_name,
                      plot_conf,
                      dataset_nm,
                      trial_begin,
                      trial_end),
              height=sizes[[size_name]][["width"]],
              width=sizes[[size_name]][["width"]] * 2) 
          
          plot(gf)
          dev.off()
          
        }
      }
    } 
    
    
    
    entropy_df_relevant  <- entropy_df[, c(1,2,3,4,5,6,11,12)]
    entropy_df_with_main_relevant  <- entropy_df_with_main[, c(1,2,3,4,5,6,11,12)]
    posterior_entropy_df_relevant  <- posterior_entropy_df[, c(1,2,3,4,5,6,11,12)]
    posterior_entropy_df_with_main_relevant  <- posterior_entropy_df_with_main[, c(1,2,3,4,5,6,11,12)]
    
    
    # Bonferroni correction
    pvalues_smaller_hit <- 
      unlist(lapply(2:(ncol(entropy_df_relevant)), 
                    function(i) {
                      wilcox.test(entropy_df_relevant[,1], 
                                  entropy_df_relevant[,i], 
                                  alternative = "less", 
                                  paired=T)$p.value * (ncol(entropy_df_relevant) - 1)
                    }))
    
    
    # Bonferroni correction
    pvalues_smaller_hit_with_main <- 
      unlist(lapply(2:(ncol(entropy_df_with_main_relevant)), 
                    function(i) {
                      wilcox.test(entropy_df_with_main_relevant[,1], 
                                  entropy_df_with_main_relevant[,i], 
                                  alternative = "greater", 
                                  paired=T)$p.value * (ncol(entropy_df_with_main_relevant) - 1)
                    }))
    
    
    pvalues_posterior_smaller_hit <- 
      unlist(lapply(2:(ncol(posterior_entropy_df_relevant)), 
                    function(i) {
                      wilcox.test(posterior_entropy_df_relevant[,1], 
                                  posterior_entropy_df_relevant[,i], 
                                  alternative = "less", 
                                  paired=T)$p.value * (ncol(posterior_entropy_df_relevant) - 1)
                    }))
    
    
    # Bonferroni correction
    pvalues_posterior_smaller_hit_with_main <- 
      unlist(lapply(2:(ncol(posterior_entropy_df_with_main_relevant)), 
                    function(i) {
                      wilcox.test(posterior_entropy_df_with_main_relevant[,1], 
                                  posterior_entropy_df_with_main_relevant[,i], 
                                  alternative = "greater", 
                                  paired=T)$p.value * (ncol(posterior_entropy_df_with_main_relevant) - 1)
                    }))
    
      
    pvalues_smaller_hit[pvalues_smaller_hit > 1] <- 1
    pvalues_smaller_hit_with_main[pvalues_smaller_hit_with_main > 1] <- 1
    pvalues_posterior_smaller_hit_with_main[pvalues_posterior_smaller_hit_with_main > 1] <- 1
    pvalues_posterior_smaller_hit[pvalues_posterior_smaller_hit > 1] <- 1
      
    
    
    entropy_dfs_confs <- 
    list(no_main=list(df=entropy_df_relevant, pvals=pvalues_smaller_hit),
         main=list(df=entropy_df_with_main_relevant, pvals=pvalues_smaller_hit_with_main),
         posterior_no_main=list(df=posterior_entropy_df_relevant, pvals=pvalues_posterior_smaller_hit),
         posterior_main=list(df=posterior_entropy_df_with_main_relevant, pvals=pvalues_posterior_smaller_hit_with_main))
    
    
    for (plot_conf in names(entropy_dfs_confs)) {
      
      conf <- entropy_dfs_confs[[plot_conf]]
      df_to_use <- conf$df
      pvals <- conf$pvals
      melted_df <- melt(df_to_use)[,c(2:3)]
      colnames(melted_df) <- c("Group", "Entropy")
      
      gbox <- 
        ggplot(melted_df, aes(x=Group,y=Entropy)) +
        
        geom_boxplot(width=0.5, 
                     size=.5, 
                     #aes(fill=Group, col=Group)) + 
                     color=adjustcolor("gray65", alpha=1),   
                     fill=adjustcolor("gray80", alpha=0.9)) +
        
        geom_jitter(position=position_jitterdodge(.5), aes(color=Group),
                    size=1.5) + 
        #scale_color_brewer(palette="Spectral") + 
        scale_color_manual(values=adjustcolor(spec_cg(ncol(df_to_use)), alpha=.75)) + 
        scale_fill_manual(values=adjustcolor(spec_cg(ncol(df_to_use)), alpha=.5)) + 
        theme_classic() +
        xlab("") +
        #ylab("H(Experimental Variable)") +
        ylab("Entropy") + 
        theme(text=element_text(size=14, color="black"),
              #legend.position="NA",
              legend.title = element_blank(),
              axis.text = element_text(color="black"),
              axis.text.x = element_blank(),#element_text(angle = 45, vjust=.6),
              axis.ticks = element_line(color="black")) 
      
      level_jump = 0.135
      max_e <- max(melted_df$Entropy) + level_jump - .05
      
      
      gboxf <- gbox
      
      for (pv_idx in 1:len(pvals)) {
        pv <- pvals[pv_idx]
        line_level <- max_e + (level_jump * (pv_idx - 1))
        
        line_df <- data.frame(x=c(1, 1, pv_idx + 1, pv_idx + 1),
                              y=c(line_level - (level_jump - .05), line_level - (level_jump - .1), 
                                  line_level -  (level_jump - .1), line_level - (level_jump - .05)))
        gboxf <- gboxf + geom_line(data=line_df, aes(x=x,y=y)) +
          geom_text(label = signif.num(pvals[pv_idx]),
                    x = (1 + pv_idx + 1) / 2,
                    y = line_level -  (level_jump - .125))
        
      }
      
      
      for (size_name in names(sizes)) {
        
# 
#         dir.create(sprintf("%s\\%s",
#                            write_path,
#                            size_name))
#         dir.create(sprintf("%s\\%s\\%s\\",
#                            write_path,
#                            size_name,
#                            plot_conf))
        
        dir.create(sprintf("%s\\%s\\%s\\entropy_comp",
                           write_path,
                           size_name,
                           plot_conf))
        
        pdf(sprintf("%s\\%s\\%s\\entropy_comp\\entropy_begin%d_end%d.pdf",
                    write_path,
                    size_name,
                    plot_conf,
                    trial_begin,
                    trial_end),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]] * 2.2) 
        
        plot(gboxf)
        dev.off()
        
        pdf(sprintf("%s\\%s\\%s\\entropy_comp\\no_legend_entropy_begin%d_end%d.pdf",
                    write_path,
                    size_name,
                    plot_conf,
                    trial_begin,
                    trial_end),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]]) 
        
        plot(gboxf + theme(legend.position="NA"))
        dev.off()
      }
    }
  }
}

figure_3_continuous_entropy_maps <- function() {
  write_path <- sprintf("%s\\figure_3\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\entropy_by_time", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=3.5,
                     height=3.5),
               big_1=c(width=3.25,
                       height=3.25),
               medium=c(width=2.75,
                        height=2.75),
               medium_2=c(width=2.5,
                          height=2.5),
               medium_1=c(width=2,
                          height=2),
               small=c(width=1.75,
                       height=1.755))
  paths <- get_thirsty_quenched_paths()
  
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  calculate_entropy <- function(pmf) {- sum(pmf * log2(pmf + 10 ** -30))}
  
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm",  "15"="HighPupil", "16"="LowPupil")
  
  
  max_windows = c(8,9,10)
  begin_offset <- c(0,1,2,3,4)
  
  for (max_window_size in max_windows){
    for (offset in begin_offset) {
    
      all_windows <- seq(1, (max_window_size - offset), by=.5)
      
      entropy_list <- list()
      posterior_entropy_list <- list()
      entropy_wth_main_list <- list()
      posterior_entropy_with_main_list <- list()
      
      
      for(resp in as.numeric(names(Responses))) {
        resp_name <- Responses[as.character(resp)]
        
        entropy_list[[resp_name]] <- c()
        posterior_entropy_list[[resp_name]] <- c()
        entropy_wth_main_list[[resp_name]] <- c()
        posterior_entropy_with_main_list[[resp_name]] <- c() 
      }
      
      for (i in 1:len(paths)) {
        
        if (i == 4 || i == 14) {
          next
        }
        
        dataset_nm <- datasets_names[[i]]
        reduced_mat <- get_mat_with_preset(paths[[i]], "dalshtaim")  
        
        base <- rep(0, times=nrow(reduced_mat))
        
        cluster_map <- get_clusters_mat_kmeans(paths[[i]])
        variables <- get_behavior_variables(paths[[i]], begin=0, end=5)
        cluster_labels <- cluster_map$labs
        no_main_indices <- cluster_labels != -1
        cluster_labels_no_main <- cluster_labels[no_main_indices]
        all_clusters <- sort(unique(cluster_labels))
        
        cluster_prob <- table(cluster_labels_no_main) / len(cluster_labels_no_main)
        cluster_prob_with_main <- table(cluster_labels) / len(cluster_labels)
        
        # Disjoint pupil vector into two vectors
        pupil_vec <- variables$variables$Pupil
        low_pupil_vec <- as.numeric(pupil_vec == 1)
        high_pupil_vec <- as.numeric(pupil_vec == 2)
        
        #names(variables$variables)[which(names(variables$variables) == "Pupil")] <- "LowPupil"
        variables$variables$LowPupil <- low_pupil_vec
        variables$variables$HighPupil <- high_pupil_vec
        
        
        
        all_var_names <- names(variables$variables)
        all_var_names <- all_var_names[all_var_names != "Pupil"]
        
        var_prob_plot_lists <- list(no_main=list(),
                                    main=list(),
                                    posterior_no_main=list(),
                                    posterior_main=list())
        
        stim_master_mat <- variables$stim_mater_mat
        
        
        
        
        
        for(resp in as.numeric(names(Responses))) {
          resp_name <- Responses[as.character(resp)]
          # Var was disjointed
          entropy_of_variables <- c()
          entropy_of_variables_with_main <- c()
          posterior_entropy_of_variables <- c()
          posterior_entropy_of_variables_with_main <- c()
          
          
          for (entropy_window_size in all_windows) {
            
            if (grepl("Pupil", resp_name)) {
              
              trials <- stim_master_mat[stim_master_mat[,"Response"] %in% c(2:4),"Frames"]
              
              
              trial_indices <- lapply(trials,
                                      function(trial_start) 
                                      {(trial_start + max_window_size * frame_rate):(trial_start + (max_window_size - entropy_window_size) * frame_rate + 1)})
              
              
              trial_indices_all <- lapply(trial_indices,
                                          function(trial_ind) {
                                            return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, window_size)})))
                                          })
              
              trial_indices_all <- unique(unlist(trial_indices_all))
              trial_indices_all <- trial_indices_all[trial_indices_all <= len(base)] 
              
              variable_indices <- base
              
              if (resp_name == "HighPupil") {
                variable_indices[trial_indices_all] <- high_pupil_vec[trial_indices_all]
              } else {
                variable_indices[trial_indices_all] <- low_pupil_vec[trial_indices_all]
              }
              
            } else {
              
              trials <- stim_master_mat[stim_master_mat[,"Response"] == resp,"Frames"]
              
              
              trial_indices <- lapply(trials,
                                      function(trial_start) 
                                      {(trial_start + offset * frame_rate):(trial_start + (offset + entropy_window_size) * frame_rate - 1)})
              
              
              trial_indices_all <- lapply(trial_indices,
                                          function(trial_ind) {
                                            return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, window_size)})))
                                          })
              
              trial_indices_all <- unique(unlist(trial_indices_all))
              
              variable_indices <- base
              variable_indices[trial_indices_all] <- 1
            }
            
            
            
            
            
            variable_indices_no_main <- variable_indices[no_main_indices]
            var_contingencies <- table(cluster_labels[variable_indices == 1] )
            all_var_prob_pmf <- table(variable_indices) 
            all_var_prob_no_main_pmf <- table(variable_indices_no_main)
            all_var_prob   <- all_var_prob_pmf[2] / sum(all_var_prob_pmf)
            all_var_prob_no_main   <- all_var_prob_no_main_pmf[2] / sum(all_var_prob_no_main_pmf)
            
            
            print(sprintf("(%s - %s - %f), Running from %f to %f (Frames sanity: %d - %d - %f - %f)",
                          resp_name,
                          as.character(resp),
                          entropy_window_size,
                          offset * frame_rate,
                          (offset + entropy_window_size) * frame_rate,
                          sum(variable_indices),
                          sum(variable_indices_no_main),
                          all_var_prob,
                          all_var_prob_no_main))
            
            
            # SANITY CHECK
            assert(all(table(cluster_labels_no_main[variable_indices_no_main == 1]) == 
                         var_contingencies[-1]))
            
            all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
            all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
            
            names(all_var_contingencies) <- all_clusters[-1]
            names(all_var_contingencies_with_main) <- all_clusters
            
            likelihood_prob <- lapply(all_clusters[-1], function(c) {table(variable_indices_no_main[which(cluster_labels_no_main == c)])})
            likelihood <- unlist(lapply(likelihood_prob, 
                                        function(v) {
                                          if (all(names(v) == "0")) {return(0)} 
                                          else if (all(names(v) == "1")) {return(1)}
                                          else {return(v[2] / sum(v))}
                                        }))
            names(likelihood) <- all_clusters[-1]
            
            likelihood_prob_with_main <- lapply(all_clusters, function(c) {table(variable_indices[which(cluster_labels == c)])})
            likelihood_with_main <- unlist(lapply(likelihood_prob_with_main, 
                                                  function(v) {
                                                    if (all(names(v) == "0")) {return(0)} 
                                                    else if (all(names(v) == "1")) {return(1)}
                                                    else {return(v[2] / sum(v))}
                                                  }))
            
            for (clust in all_clusters) {
              
              # In case there werent any frames of variables in that perticular cluster
              if(!as.character(clust) %in% names(var_contingencies)) {
                next
              }
              
              if (clust == -1) {
                all_var_contingencies_with_main[as.character(clust)]  <- 
                  var_contingencies[as.character(clust)]
              } else {
                all_var_contingencies_with_main[as.character(clust)]  <- 
                  var_contingencies[as.character(clust)]
                
                all_var_contingencies[as.character(clust)]  <- 
                  var_contingencies[as.character(clust)]
              }
            }
            
            var_probability <- all_var_contingencies / sum(all_var_contingencies)
            var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)  
            posterior_probability <- (likelihood * cluster_prob) / all_var_prob_no_main
            posterior_probability_with_main <- (likelihood_with_main * cluster_prob_with_main) / all_var_prob
            
            
            names(posterior_probability) <- all_clusters[-1]
            names(posterior_probability_with_main) <- all_clusters
            
            entropy_of_variables <- c(entropy_of_variables,
                                      calculate_entropy(var_probability))
            
            entropy_of_variables_with_main <- c(entropy_of_variables_with_main,
                                                calculate_entropy(var_probability_with_main))
            
            posterior_entropy_of_variables <- c(posterior_entropy_of_variables,
                                                calculate_entropy(posterior_probability))
            
            posterior_entropy_of_variables_with_main <- c(posterior_entropy_of_variables_with_main,
                                                          calculate_entropy(posterior_probability_with_main))
          }
          
          entropy_list[[resp_name]] <- rbind(entropy_list[[resp_name]],
                                             entropy_of_variables)
          posterior_entropy_list[[resp_name]] <- rbind(posterior_entropy_list[[resp_name]],
                                                       posterior_entropy_of_variables)
          entropy_wth_main_list[[resp_name]] <- rbind(entropy_wth_main_list[[resp_name]],
                                                      entropy_of_variables_with_main)
          posterior_entropy_with_main_list[[resp_name]] <- rbind(posterior_entropy_with_main_list[[resp_name]],
                                                                 posterior_entropy_of_variables_with_main)
          
        }
      } 
      
      
      list_of_dfs <- list(entropy_no_main=entropy_list,
                          posterior_no_main=posterior_entropy_list,
                          entropy_with_main=entropy_wth_main_list,
                          posterior_with_main=posterior_entropy_with_main_list)
      
      
      class_titles <- list(entropy_no_main="Entropy (no main cluster)",
                           posterior_no_main="Posterior entropy (no main cluster)",
                           entropy_with_main="Entropy (all)",
                           posterior_with_main="Posterior entropy (all)")
      for (class_name in names(list_of_dfs)) {
        
        df_list <- list_of_dfs[[class_name]]
        
        
        mean_sd <- 
          lapply(df_list, 
                 function(df)
                 {
                   
                   rbind(apply(df, 2, mean, na.rm=T),
                         apply(df, 2, sem))
                   
                 })
        
        mean_sd_window <- 
          lapply(mean_sd,
                 function(mean_sd_df) {
                   t(rbind(mean_sd_df,
                           all_windows))
                 })
        
        df_all <- do.call(rbind, mean_sd_window)
        df_all <- as.data.frame(df_all)
        df_all <- cbind(df_all, rep(names(df_list), each=ncol(mean_sd[[1]])))
        
        colnames(df_all) <- c("mean", "sd", "window", "group")
        
        gf <- 
          ggplot(df_all) + 
          geom_ribbon(aes(x=window, ymin=mean-sd, ymax=mean+sd, fill=group), color="NA", alpha=.2) + 
          geom_line(aes(x=window, y=mean, color=group), size=1) +
          #scale_color_manual(values=adjustcolor(spec_cg(ncol(df_to_use)), alpha=.75)) + 
          #scale_fill_manual(values=adjustcolor(spec_cg(ncol(df_to_use)), alpha=.5)) + 
          scale_fill_brewer(palette="Dark2") + 
          scale_color_brewer(palette="Dark2") + 
          theme_classic() +
          xlab("Window size (sec)") +
          #ylab("H(Experimental Variable)") +
          ylab("Entropy") + 
          theme(text=element_text(size=14, color="black"),
                #legend.position="NA",
                legend.title = element_blank(),
                axis.text = element_text(color="black"),
                #axis.text.x = element_blank(),#element_text(angle = 45, vjust=.6),
                axis.ticks = element_line(color="black")) +
          ggtitle(class_titles[[class_name]])
        
        
        for (size_name in names(sizes)) {
          
          
          dir.create(sprintf("%s\\%s",write_path, size_name))
          
          
          pdf(sprintf("%s\\%s\\%s_begin_offset%d_window%d.pdf",
                      write_path,
                      size_name,
                      class_name,
                      offset,
                      max_window_size),
              height=sizes[[size_name]][["width"]],
              width=sizes[[size_name]][["width"]] * 1.5) 
          
          plot(gf)
          dev.off()
        }  
      }
    }
  }
}

figure_3_continuous_entropy_shuffles <- function() {
  write_path <- sprintf("%s\\figure_3\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\shuffle_distributions", write_path)
  dir.create(write_path)  
  
  write_path <- sprintf("%s\\all_regular_new\\", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=3.5,
                     height=3.5),
               big_1=c(width=3.25,
                       height=3.25),
               medium=c(width=2.75,
                        height=2.75),
               medium_2=c(width=2.5,
                          height=2.5),
               medium_1=c(width=2,
                          height=2),
               small=c(width=1.75,
                       height=1.755))
  paths <- get_thirsty_quenched_paths()
  
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  calculate_entropy <- function(pmf) {- sum(pmf * log2(pmf + 10 ** -30))}
  
  
  Responses = c("0"="Hit", "1"="Miss", "2"="CR (N)", "3"="FA (N)", "4"="CR (A)", "5"="FA (A)",  "15"="Pupil", "16"="Consumed water")
  frame_rate=30
  binning_window_size=15
  #max_windows = c(0:5)
  #max_windows <- c(4,3,2)
  #begin_offset <- c(0,6,7)
  max_windows <- c(3)
  begin_offset <- c(6)
  plot_structure = F
  thirsty_q <- .15
  max_fun <- combined_max
  
  cont_breaks = 7
  n_shuffles = 150
  n_subsamples = 100
  
  ylabs <- list(occ_no_main="Occupancy (%)",
                occ_main="Occupancy (%)",
                var_prob="Probability (%)",
                var_prob_main="Probability (%)",
                entropy_no_main="Entropy",
                entropy_main="Entropy",
                MI_main="Mutual information",
                MI_no_main="Mutual information")
  
  all_results <- list()
  all_pval_results <- list()
  
  
    for (max_window_size in max_windows){
      for (offset in rev(begin_offset)) {
      
      #all_windows <- seq(1, (max_window_size - offset), by=.5)
      entropy_window_size = max_window_size - offset

      all_shuffle_matrices <- list()
      all_pvals <- list()
      
      for (path_idx in c(1:14)) {
        
        # if (path_idx == 4 || path_idx == 14) {
        #   next
        # }
        
        dataset_nm <- datasets_names[[path_idx]]
        reduced_mat <- metadata[[path_idx]]$red_mat
        stim_master_mat <- metadata[[path_idx]]$stim_master_mat
        
        base <- rep(0, times=nrow(reduced_mat))
        base_color <- rep(adjustcolor("gray70", alpha=0.3), times=nrow(reduced_mat))
        
        cluster_map <- metadata[[path_idx]]$cluster_mat
        cluster_labels <- cluster_map$labs
        no_main_indices <- cluster_labels != -1
        cluster_labels_no_main <- cluster_labels[no_main_indices]
        all_clusters <- sort(unique(cluster_labels))
        cluster_prob <- table(cluster_labels_no_main) / len(cluster_labels_no_main)
        cluster_prob_with_main <- table(cluster_labels) / len(cluster_labels)
        
        pupil_vec_all <- get_pupil_files(paths[[path_idx]], window_size=15, load_all = T)
        
        drink_ind <- stim_master_mat[!is.nan(stim_master_mat[,6]),"Frames"]
        drink_ind <- unique(get_binned_index(drink_ind, binning_window_size))
        consumed_water_vec <- base
        consumed_water_vec[drink_ind] <- 1
        consumed_water_vec <- cumsum(consumed_water_vec)
        
        
        ongoing_trials <- which(metadata[[path_idx]]$stim_master_mat[,7] %in% c(4,5) & metadata[[path_idx]]$stim_master_mat[,8]%%2 == 0)
        ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[path_idx]]$stim_master_mat)]
        
        if (offset == 6 || offset == 7) {
          ongoing_trials_frames = metadata[[path_idx]]$stim_master_mat[ongoing_trials + 1,1]
          ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - max_window_size * 2):(get_binned_index(i, 15)-1)})))
        } else {
          print ("USING OFFSET <6")     
          ongoing_trials_frames <- lapply(metadata[[path_idx]]$stim_master_mat[ongoing_trials,1],
                                  function(trial_start) 
                                  {(trial_start + offset * frame_rate):(trial_start + (offset + entropy_window_size) * frame_rate - 1)})
          
          
          ongoing_ind <- lapply(ongoing_trials_frames,
                                      function(trial_ind) {
                                        return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, binning_window_size)})))
                                      })
          
          
          ongoing_ind <- unique(unlist(ongoing_ind))
          ongoing_ind <- ongoing_ind[ongoing_ind <= len(base)] 
        }
        ongoing_ind_logic <- 1:nrow(reduced_mat) %in% ongoing_ind
        ongoing_ind_no_main <- ongoing_ind_logic & no_main_indices
        
        #all_windows <- seq(1, (max_window_size - offset), by=.5)
        
        occ_main_all <- c()
        occ_all <- c()
        occ_of_var_main_all <- c()
        occ_of_var_all <- c()

        MI_all <- c()
        MI_main_all <- c()
        
        
        for (shuff_i in c(-1, 1:n_shuffles)) {
          iteration_occ_main_all <- c()
          iteration_occ_all <- c()
          iteration_occ_of_var_main_all <- c()
          iteration_occ_of_var_all <- c()

          iteration_MI_all <- c()
          iteration_MI_main_all <- c()
          
          shuffled_stim_mat <- stim_master_mat
          shuffled_stim_mat[,"Frames"] <- sample(stim_master_mat[,"Frames"], nrow(stim_master_mat))
          
          if (shuff_i == -1) {
            print("Original mat!")
            work_stim_mat <- stim_master_mat
          } else {
            if (shuff_i %% 50 == 0) {
              print(sprintf("shuffle!!! %d", shuff_i))
            }
            work_stim_mat <- shuffled_stim_mat
            #block_shuffle_ind <- block_shuffle(nrow(reduced_mat), 20)
          }
          
          
          if (offset == 6 || offset == 7) {
            shuffle_ongoing_trials_frames = work_stim_mat[ongoing_trials + 1,1]
            shuffle_ongoing_ind <- unique(unlist(lapply(shuffle_ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - max_window_size * 2):(get_binned_index(i, 15)-1)})))
          } else {
            #print ("USING OFFSET <6")     
            shuffle_ongoing_trials_frames <- lapply(work_stim_mat[ongoing_trials,1],
                                            function(trial_start) 
                                            {(trial_start + offset * frame_rate):(trial_start + (offset + entropy_window_size) * frame_rate - 1)})
            
            
            shuffle_ongoing_ind <- lapply(ongoing_trials_frames,
                                  function(trial_ind) {
                                    return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, binning_window_size)})))
                                  })
            
            
            shuffle_ongoing_ind <- unique(unlist(shuffle_ongoing_ind))
            shuffle_ongoing_ind <- shuffle_ongoing_ind[shuffle_ongoing_ind <= len(base)] 
          }

          shuffle_ongoing_ind_logic <- 1:nrow(reduced_mat) %in% shuffle_ongoing_ind
          #shuffle_ongoing_ind_no_main <- sample(1:nrow(reduced_mat),  sum(ongoing_ind_logic & no_main_indices))
          shuffle_ongoing_ind_no_main <- shuffle_ongoing_ind_logic & no_main_indices

          ind_to_use <- 1:nrow(reduced_mat)
          ind_to_use_no_main <- ind_to_use[no_main_indices]
          
          if (shuff_i == -1) {
            ind_to_use_ongoing <- ongoing_ind
            ind_to_use_ongoing_no_main <- ongoing_ind_no_main
          } else {
            ind_to_use_ongoing <- shuffle_ongoing_ind
            ind_to_use_ongoing_no_main <- shuffle_ongoing_ind[shuffle_ongoing_ind %in% which(no_main_indices)]#shuffle_ongoing_ind_no_main
          }
          
          
          indices_all <- list()
          
          for(resp in as.numeric(names(Responses))) {
            resp_name <- Responses[as.character(resp)]
            
            indices_all[[resp_name]] <- base
            
            #if (grepl("Pupil", resp_name)) {
            if (resp < 15){
              
              trials <- work_stim_mat[work_stim_mat[,"Response"] == resp,"Frames"]
              
              
              trial_indices <- lapply(trials,
                                      function(trial_start) 
                                      {(trial_start + offset * frame_rate):(trial_start + (offset + entropy_window_size) * frame_rate - 1)})
              
              
              trial_indices_all <- lapply(trial_indices,
                                          function(trial_ind) {
                                            return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, binning_window_size)})))
                                          })
              
              trial_indices_all <- unique(unlist(trial_indices_all))
              trial_indices_all <- trial_indices_all[trial_indices_all <= len(base)] 
              
              indices_all[[resp_name]][trial_indices_all] <- 1
            }
            
            
          }
          
          # min_timepoints <- min(unlist(lapply(indices_all, sum)))
          # max_occ <- c() # Fraction of cluster that is explained by event
          # max_occ_with_main <- c()
          # max_occ_of_var <- c() # Fraction of event that is explained by cluster
          # max_occ_of_var_with_main <- c()
          # 
          # MI <- c()
          # MI_with_main <- c()
          
          
          for(resp in as.numeric(names(Responses))) {
            
            resp_name <- Responses[as.character(resp)]
            
            if (resp >= 15) {
              
              if (resp_name == "Pupil") {
                vec_to_discretize <- pupil_vec_all$smoothed
              } else {
                vec_to_discretize <- consumed_water_vec
              }
              
              
              
              
              discretized_vec <- 
                as.numeric(cut(vec_to_discretize[ind_to_use_ongoing], 
                               breaks = seq(min(vec_to_discretize[ind_to_use_ongoing] - 1), 
                                            max(vec_to_discretize[ind_to_use_ongoing]), 
                                            length.out=cont_breaks + 1)))
              
              discretized_vec_no_main <- 
                as.numeric(cut(vec_to_discretize[ind_to_use_ongoing_no_main], 
                               breaks = seq(min(vec_to_discretize[ind_to_use_ongoing_no_main] - 1), 
                                            max(vec_to_discretize[ind_to_use_ongoing_no_main]), 
                                            length.out=cont_breaks + 1)))
              
              
              working_cluster_labels <- rep(0, times=len(all_clusters) - 1)
              working_cluster_labels_with_main <- rep(0, times=len(all_clusters))
              names(working_cluster_labels) <- all_clusters[-1]
              names(working_cluster_labels_with_main) <- all_clusters
              
              tmp <- table(cluster_labels[ind_to_use_ongoing_no_main])
              working_cluster_labels[names(tmp)] <- tmp
              
              tmp <- table(cluster_labels[ind_to_use_ongoing])
              working_cluster_labels_with_main[names(tmp)] <- tmp
              
              
              parameterization_variables <- 
                lapply(1:cont_breaks,
                       function(b) {
                         var_contingencies <- table(cluster_labels[ind_to_use_ongoing][discretized_vec == b])
                         
                         all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
                         all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
                         
                         names(all_var_contingencies) <- all_clusters[-1]
                         names(all_var_contingencies_with_main) <- all_clusters
                         
                         all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
                         all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies
                         
                         clust_occ_with_main <- all_var_contingencies_with_main / table(cluster_labels)
                         clust_occ <- all_var_contingencies / table(cluster_labels_no_main)
                         var_probability <- all_var_contingencies / sum(all_var_contingencies)
                         var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)
                         
                         return(list(raw_main=all_var_contingencies_with_main,
                                     raw=all_var_contingencies,
                                     clust_occ_with_main=clust_occ_with_main,
                                     clust_occ=clust_occ,
                                     var_probability=var_probability,
                                     var_probability_with_main=var_probability_with_main))
                       })
              
              
              
              #iteration_occ_all <-
              #   c(iteration_occ_all, 
              #     max_fun(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$clust_occ})), 2, median, na.rm=T)))
              # 
              # iteration_occ_main_all <- 
              #   c(iteration_occ_main_all, 
              #     max_fun(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$clust_occ_with_main})), 2, median, na.rm=T)))
              # 
              # iteration_occ_of_var_all <- 
              #     c(iteration_occ_of_var_all, 
              #       max_fun(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$var_probability})), 2, median, na.rm=T)))
              # 
              # iteration_occ_of_var_main_all <- 
              #     c(iteration_occ_of_var_main_all, 
              #       max_fun(apply(do.call(rbind,lapply(parameterization_variables, function(vc) {vc$var_probability_with_main})), 2, median, na.rm=T)))
              
              raw_main_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$raw_main}));
              raw_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$raw}));
              prob_var_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$var_probability}));
              prob_var_main_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$var_probability_with_main}));
              clust_occ_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$clust_occ}));
              clust_occ_main_table <- do.call(rbind, lapply(parameterization_variables, function(vc) {vc$clust_occ_with_main}));
              
              # prob_var_table <- prob_var_table/ sum(prob_var_table)
              # prob_var_main_table <- prob_var_main_table/ sum(prob_var_main_table)
              # clust_occ_table <- clust_occ_table/ sum(clust_occ_table)
              # clust_occ_main_table <- clust_occ_main_table/ sum(clust_occ_main_table)
              # 
              # 
              # iteration_occ_all <- c(iteration_occ_all, max_fun(clust_occ_table))
              # iteration_occ_main_all <- c(iteration_occ_main_all, max_fun(clust_occ_main_table))
              # iteration_occ_of_var_all <- c(iteration_occ_of_var_all, max_fun(prob_var_table))
              # iteration_occ_of_var_main_all <- c(iteration_occ_of_var_main_all, max_fun(prob_var_main_table))
              
              
              #iteration_occ_all <- c(iteration_occ_all, mean(apply(clust_occ_table, 2, function(v) { max_fun(v)})))
              #iteration_occ_main_all <- c(iteration_occ_main_all, mean(apply(clust_occ_main_table, 2, function(v) { max_fun(v)})))
              iteration_occ_all <- c(iteration_occ_all, mean(sort(apply(clust_occ_table, 2, max), decreasing = T)[1:2]))
              iteration_occ_main_all <- c(iteration_occ_main_all, mean(sort(apply(clust_occ_main_table, 2, max), decreasing = T)[1:2]))
              iteration_occ_of_var_all <- c(iteration_occ_of_var_all, mean(apply(prob_var_table, 2, function(v) { max_fun(v)})))
              iteration_occ_of_var_main_all <- c(iteration_occ_of_var_main_all, mean(apply(prob_var_main_table, 2, function(v) { max_fun(v)})))
              
              tabulated_main <- cross_tabulate(discretized_vec,
                                               cluster_labels[sort(ind_to_use_ongoing)])
              
              tabulated_no_main <- cross_tabulate(discretized_vec_no_main,
                                                  cluster_labels[sort(ind_to_use_ongoing_no_main)])
              
              iteration_MI_all <- 
                c(iteration_MI_all,
                  calculate_MI_from_cont(tabulated_no_main))
              
              iteration_MI_main_all <- 
                c(iteration_MI_main_all,
                  calculate_MI_from_cont(tabulated_main))
              
              
            } else {
              
              variable_indices <- indices_all[[resp_name]][ind_to_use]
              variable_indices_no_main <- variable_indices[ind_to_use_no_main]
              var_contingencies <- table(cluster_labels[variable_indices == 1] )
              
              all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
              all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
              
              names(all_var_contingencies) <- all_clusters[-1]
              names(all_var_contingencies_with_main) <- all_clusters
              
              
              all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
              all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies
              
              var_probability <- all_var_contingencies / sum(all_var_contingencies)
              var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)
              
              clust_occ_with_main <- all_var_contingencies_with_main / table(cluster_labels)
              clust_occ <- all_var_contingencies / table(cluster_labels_no_main)
              
          
              tabulated_main <- cross_tabulate(cluster_labels,
                                               variable_indices)
              
              tabulated_no_main <- cross_tabulate(cluster_labels_no_main,
                                                  variable_indices_no_main)
              
              #iteration_occ_all <- c(iteration_occ_all, max_fun(clust_occ))
              #iteration_occ_main_all <- c(iteration_occ_main_all, max_fun(clust_occ_with_main))
              iteration_occ_all <- c(iteration_occ_all, mean(sort(clust_occ, decreasing=T)[1:2]))
              iteration_occ_main_all <- c(iteration_occ_main_all, mean(sort(clust_occ_with_main, decreasing=T)[1:2]))
              iteration_occ_of_var_all <- c(iteration_occ_of_var_all, max_fun(var_probability))
              iteration_occ_of_var_main_all <- c(iteration_occ_of_var_main_all, max_fun(var_probability_with_main))
              iteration_MI_all <- c(iteration_MI_all, calculate_MI_from_cont(tabulated_no_main))
              iteration_MI_main_all <- c(iteration_MI_main_all, calculate_MI_from_cont(tabulated_main))
            }
              # 
              # 
              # variable_indices <- indices_all[[resp_name]]
              # variable_indices_no_main <- variable_indices[no_main_indices]
              # var_contingencies <- table(cluster_labels[variable_indices == 1] )
              # 
              # all_var_contingencies <- rep(0, times=len(all_clusters) - 1)
              # all_var_contingencies_with_main <- rep(0, times=len(all_clusters))
              # 
              # names(all_var_contingencies) <- all_clusters[-1]
              # names(all_var_contingencies_with_main) <- all_clusters
              # 
              # 
              # all_var_contingencies[names(var_contingencies[-1])] <- var_contingencies[-1]
              # all_var_contingencies_with_main[names(var_contingencies)] <- var_contingencies
              # 
              # var_probability <- all_var_contingencies / sum(all_var_contingencies)
              # var_probability_with_main <- all_var_contingencies_with_main / sum(all_var_contingencies_with_main)
              # 
              # clust_occ_with_main <- all_var_contingencies_with_main / table(cluster_labels)
              # clust_occ <- all_var_contingencies / table(cluster_labels_no_main)
              # 
              # 
              # 
              # iteration_occ_all <- c(iteration_occ_all, max(clust_occ))
              # iteration_occ_main_all <- c(iteration_occ_main_all, max(clust_occ_with_main))
              # 
              # iteration_occ_of_var_all <- c(iteration_occ_of_var_all, max(var_probability))
              # iteration_occ_of_var_main_all <- c(iteration_occ_of_var_main_all, max(var_probability_with_main))
              # 
              # iteration_entropy_all <- c(iteration_entropy_all, calculate_entropy(var_probability))
              # iteration_entropy_main_all <- c(iteration_entropy_main_all, calculate_entropy(var_probability_with_main))
              # 
              # tabulated_main <- cross_tabulate(cluster_labels,variable_indices)
              # tabulated_no_main <- cross_tabulate(cluster_labels_no_main, variable_indices_no_main)
              # 
              # iteration_MI_all <- c(iteration_MI_all, calculate_MI_from_cont(tabulated_no_main))
              # iteration_MI_main_all <- c(iteration_MI_main_all, calculate_MI_from_cont(tabulated_main))

          }
          
          if (shuff_i == -1) {
            mdf <- rbind(c(iteration_MI_main_all),
                         c(iteration_occ_main_all),
                         c(iteration_occ_of_var_main_all),
                         c(iteration_MI_all),
                         c(iteration_occ_all),
                         c(iteration_occ_of_var_all))
            
            rownames(mdf) <- c("MI_main",
                               "clust occ_main",
                               "occ_of_var_main",
                               "MI",
                               "clust occ",
                               "occ_of_var")
            
            colnames(mdf) <- Responses
            
            print("####################")
            print(mdf)
            print("####################")
          }

          occ_main_all <- rbind(occ_main_all, iteration_occ_main_all)
          occ_all <- rbind(occ_all, iteration_occ_all)
          occ_of_var_main_all <- rbind(occ_of_var_main_all, iteration_occ_of_var_main_all)
          occ_of_var_all <- rbind(occ_of_var_all, iteration_occ_of_var_all)
          MI_all <- rbind(MI_all, iteration_MI_all)
          MI_main_all <- rbind(MI_main_all, iteration_MI_main_all)
        }
        
        
        shuffle_results <- list(shuffle_occ_main_all=occ_main_all,
                                shuffle_occ_all=occ_all,
                                shuffle_occ_of_var_main_all=occ_of_var_main_all,
                                shuffle_occ_of_var_all=occ_of_var_all,
                                 shuffle_MI_all=MI_all,
                                shuffle_MI_main_all=MI_main_all)
        
        for (shuff_mat_name in names(shuffle_results)) {
          shuff_mat <- shuffle_results[[shuff_mat_name]]
          non_shuffle <- shuff_mat[1,]
          shuff_only <- shuff_mat[-1,]
          print(dim(shuff_only))
          pvals <- unlist(lapply(1:ncol(shuff_mat), function(ci) {if(is.na(non_shuffle[ci])) {return(NA)}; return(ecdf(shuff_mat[,ci])(non_shuffle[ci]))}))
          print(pvals)
          all_pvals[[shuff_mat_name]] <- rbind(all_pvals[[shuff_mat_name]],
                                               pvals)
        }

        all_shuffle_matrices[[as.character(path_idx)]] <- shuffle_results
      } 
      
      all_results[[sprintf("w%d_o%d",max_window_size, offset)]] <- all_shuffle_matrices
      
      all_pval_results[[sprintf("w%d_o%d",max_window_size, offset)]] <- all_pvals
    
        
        
    }
  }
  
  save(file=sprintf("%s/figure_3_parameterization_comp.Rda", figures_base_path), all_results)
  save(file=sprintf("%s/figure_3_parameterization_pval.Rda", figures_base_path), all_pval_results)
  # 
  ylabs <- list(shuffle_occ_main_all="Fraction of variable (%)",
                shuffle_occ_all="Fraction of variable (%)",
                shuffle_occ_of_var_main_all="Fraction of cluster (%)",
                shuffle_occ_of_var_all="Fraction of cluster (%)",
                shuffle_entropy_main_all="Entropy",
                shuffle_entropy_all="Entropy",
                shuffle_MI_all="Mutual information",
                shuffle_MI_main_all="Mutual information")
  
  
  save(file=sprintf("%s\\pvalue_distributions.Rda", write_path), all_results)
  
  all_confs <- names(all_results)
  df_names <- names(all_results[[1]][[1]])
  
  for (conf_name in all_confs) {
    
    
    
    
    
    relevant_dfs <- all_results[[conf_name]]
    
    
    for (df_name in df_names) {
      
      conf_df <- data.frame()

      work_df <- c()
      shuffle_work_df <- c()
        
      for (dataset_idx in 1:len(relevant_dfs)) {
          

          work_df <- rbind(work_df, relevant_dfs[[dataset_idx]][[df_name]][1,])
          shuffle_work_df <- rbind(shuffle_work_df, colMeans(relevant_dfs[[dataset_idx]][[df_name]][-1,]))
      }

        
        
      mean_df <- apply(work_df, 2, mean, na.rm=T)
      sd_df <- apply(work_df, 2, sem)
      fdf <- cbind(mean_df, sd_df)
      fdf <- as.data.frame(fdf)
      fdf <- cbind(fdf, Responses)
      colnames(fdf) <- c("Mean", "Sd", "Var")
        
      shuffle_mean_df <- apply(shuffle_work_df, 2, mean, na.rm=T)
      shuffle_sd_df <- apply(shuffle_work_df, 2, sem)
      shuffle_fdf <- cbind(shuffle_mean_df, shuffle_sd_df)
      shuffle_fdf <- as.data.frame(shuffle_fdf)
      shuffle_fdf <- cbind(shuffle_fdf, Responses)
      colnames(shuffle_fdf) <- c("Mean", "Sd", "Var")
        
      df_all <-  rbind(fdf, shuffle_fdf)
    
  
      variables <- Responses
      colnames(work_df) <- variables
      colnames(shuffle_work_df) <- variables
      melt_reg <-  melt(work_df)
      melt_shuff <-  melt(shuffle_work_df)
      melt_shuff$group = rep("Shuffle", times=nrow(melt_shuff))
      melt_reg$group = rep("Regular", times=nrow(melt_reg))
      
      
      all_df <- rbind(melt_reg[,-1])#, melt_shuff[,-1])
      
      top_var_melt_reg <-  melt(work_df[,c("Hit", "Pupil", "Consumed water")])
      top_var_melt_shuff <-  melt(shuffle_work_df[,c("Hit", "Pupil", "Consumed water")])
      top_var_melt_shuff$group = rep("Shuffle", times=nrow(top_var_melt_shuff))
      top_var_melt_reg$group = rep("Regular", times=nrow(top_var_melt_reg))
      top_var_all_df <- rbind(top_var_melt_reg[,-1], top_var_melt_shuff[,-1])
      
      colnames(all_df) <- c("Var", "Val", "Group")
      colnames(top_var_all_df) <- c("Var", "Val", "Group")
      # (width=0.5, 
      #   size=.5, 
      #   #aes(fill=Group, col=Group)) + 
      #   color),   
      #   fill=)
      
      wilcox_statistics <- 
      lapply(1:3, function(i) {
        var  <- c("Hit", "Pupil", "Consumed water")[i]
        wilk <- 
        wilcox.test(work_df[,var],
                    shuffle_work_df[,var],
                    alternative="greater",
                    signed=T,
                    correct = F)
        
        stat_df <- 
        data.frame(w=wilk$statistic,
                   pval=wilk$p.value,
                   adj_pval=wilk$p.value * 3,
                   alternative=wilk$alternative,
                   method=wilk$method,
                   signif=signif.num(wilk$p.value),
                   adj_signif=signif.num(wilk$p.value * 3),
                   var=var)
        return(stat_df)
      })
      
      wilcox_statistics <- do.call(rbind, wilcox_statistics)
      
      gbox_top_var <- 
        ggplot(top_var_all_df, aes(x=Var, y=Val)) + 
        geom_jitter(aes(fill=Group, color=Group, group=Group), position=position_dodge(.95), fill="black", color="black") +
        geom_boxplot(aes(fill=Group, color=Group), width=.5, size=.25, position=position_dodge(.95), outlier.shape = NA) +
        scale_fill_manual(breaks=c("Regular", "Shuffle"),
                          values=c(adjustcolor(c("gray80", "#ff7d7d"), alpha=0.9))) + 
        scale_color_manual(breaks=c("Regular", "Shuffle"),
                           values=c(adjustcolor(c("gray55", "#e81c1c"), alpha=1))) + 
        theme_classic() +
        theme(text=element_text(size=14, color="black"), 
              axis.text=element_text(size=14, color="black"),
              axis.ticks = element_line(color="black"),
              plot.title=element_text(size=10),
              legend.position="NA")  + 
        ggtitle(sprintf("%s %s", conf_name, df_name)) + 
        xlab("") +
        #ylab(ylabs[[conf_name]]) +
        ylab(ylabs[[df_name]])
      
      
      gbox <- 
        ggplot(all_df, aes(x=Var, y=Val)) + 
              geom_jitter(aes(fill=Group, color=Group, group=Group), position=position_dodge(.95), fill="black", color="black") +
              geom_boxplot(aes(fill=Group, color=Group), width=.5, size=.25, position=position_dodge(.95), outlier.shape = NA) +
              scale_fill_manual(breaks=c("Regular", "Shuffle"),
                                values=c(adjustcolor(c("gray80", "#ff7d7d"), alpha=0.9))) + 
              scale_color_manual(breaks=c("Regular", "Shuffle"),
                          values=c(adjustcolor(c("gray55", "#e81c1c"), alpha=1))) + 
              theme_classic() +
              theme(text=element_text(size=14, color="black"), 
              axis.text=element_text(size=14, color="black"),
              axis.ticks = element_line(color="black"),
              plot.title=element_text(size=10),
              legend.position="NA")  + 
              ggtitle(sprintf("%s %s", conf_name, df_name)) + 
              xlab("") +
              #ylab(ylabs[[conf_name]]) +
              ylab(ylabs[[df_name]])
      
      
      
      
      statistics <- 
        as.data.frame(dunn.test::dunn.test(all_df[,2], paste(all_df[,1], all_df[,3]), 
                                           method="bonferroni"))
      
      statistics <- cbind(statistics, signif.num(statistics$P.adjusted))
        
      #         
      # 
      # gtime_all <- 
      #   ggplot(conf_df, aes(x=Window, y=Mean, group=Var)) + 
      #   geom_line(aes(color=Var), size=.8) + 
      #   geom_point(aes(color=Var), size=1.25) + 
      #   #geom_errorbar(aes(color=Var, ymin=Mean - Sd, ymax=Mean + Sd), width=.2, alpha=1) +
      #   geom_ribbon(aes(fill=Var, ymin=Mean - Sd, ymax=Mean + Sd), color=NA, alpha=.2) +
      #   theme_classic() +
      #   theme(text=element_text(size=14, color="black"), 
      #         axis.text=element_text(size=14, color="black"),
      #         axis.ticks = element_line(color="black"),
      #         plot.title=element_text(size=10),
      #         legend.position="NA")  + 
      #   xlab("Offset size (s)") +
      #   ylab(ylabs[[df_name]]) +
      #   ggtitle(sprintf("Offset: %d %s", offset, ifelse(grepl("main", df_name), "(Main)", ""))) + 
      #   geom_vline(xintercept=0, size=.6, color="gray60", linetype="dashed")
      # 
      # 
      # 
      # gtime_all_shuffle <- 
      #   ggplot(shuffle_conf_df, aes(x=Window, y=Mean, group=Var)) + 
      #   geom_line(aes(color=Var), size=.8) + 
      #   geom_point(aes(color=Var), size=1.25) + 
      #   #geom_errorbar(aes(color=Var, ymin=Mean - Sd, ymax=Mean + Sd), width=.2, alpha=1) +
      #   geom_ribbon(aes(fill=Var, ymin=Mean - Sd, ymax=Mean + Sd), color=NA, alpha=.2) +
      #   theme_classic() +
      #   theme(text=element_text(size=14, color="black"), 
      #         axis.text=element_text(size=14, color="black"),
      #         axis.ticks = element_line(color="black"),
      #         plot.title=element_text(size=10),
      #         legend.position="NA")  + 
      #   xlab("Offset size (s)") +
      #   ylab(ylabs[[df_name]]) +
      #   ggtitle(sprintf("Offset: %d %s", offset, ifelse(grepl("main", df_name), "(Main)", ""))) + 
      #   geom_vline(xintercept=0, size=.6, color="gray60", linetype="dashed")
      
      # panels <- list()
      # conf_ymax <- max(conf_df$Mean + conf_df$Sd)
      # conf_ymin <- min(conf_df$Mean - conf_df$Sd)
      # shuffle_conf_ymax <- max(shuffle_conf_df$Mean + shuffle_conf_df$Sd)
      # shuffle_conf_ymin <- min(shuffle_conf_df$Mean - shuffle_conf_df$Sd)
      # 
      # ymin <- min(shuffle_conf_ymin, conf_ymin)
      # ymax <- max(shuffle_conf_ymax, conf_ymax)
      # 
      # 
      # 
      # for (var_name in variables) {
      #   
      #   work_conf_df <- conf_df[conf_df$Var == var_name,]
      #   shuffle_work_conf_df <-  shuffle_conf_df[shuffle_conf_df$Var == var_name,]
      #   
      #   work_df_f <- rbind(work_conf_df,shuffle_work_conf_df)
      #   work_df_f$Group = c(rep("Regular", times=nrow(work_conf_df)), rep("Shuffle", times=nrow(shuffle_work_conf_df)))
      #   
      #   
      #   gtime_panel <- 
      #     ggplot(work_df_f, aes(x=Window, y=Mean, group=interaction(Var, Group))) + 
      #     geom_line(aes(color=interaction(Var, Group)), size=.8) + 
      #     geom_point(aes(color=interaction(Var, Group)), size=1.25) + 
      #     #geom_errorbar(aes(color=Var, ymin=Mean - Sd, ymax=Mean + Sd), width=.2, alpha=1) +
      #     geom_ribbon(aes(fill=interaction(Var, Group), ymin=Mean - Sd, ymax=Mean + Sd), color=NA, alpha=.2) +
      #     theme_classic() +
      #     theme(text=element_text(size=14, color="black"), 
      #           axis.text=element_text(size=14, color="black"),
      #           axis.ticks = element_line(color="black"),
      #           plot.title=element_text(size=10),
      #           legend.position="NA")  + 
      #     xlab("Offset (s)") +
      #     ylab(ylabs[[df_name]]) +
      #     #ggtitle(sprintf("%s %s", var_name, ifelse(grepl("main", df_name), "(Main)", ""))) +
      #     ggtitle(var_name) +
      #     ylim(ymin, ymax)
      #     geom_vline(xintercept=0, size=.6, color="gray60", linetype="dashed")
      #   
      #   panels <- append(panels, list(gtime_panel))
      #   
      # }
      # 
      # panels$nrow <- 2
      # panels_plot <- do.call(plot_grid, panels)
      
      for (size_name in names(sizes)) {
        dir.create(sprintf("%s\\%s\\",
                           write_path,
                           size_name))
        
        dir.create(sprintf("%s\\%s\\%s",
                           write_path,
                           size_name,
                           conf_name))
        
        
        
        pdf(sprintf("%s\\%s\\%s\\%s_comparision.pdf",
                    write_path,
                    size_name,
                    conf_name,
                    df_name),
            height=sizes[[size_name]][["width"]],
            width=2 * sizes[[size_name]][["width"]]) 
        
        plot(gbox)
        dev.off()
        
        pdf(sprintf("%s\\%s\\%s\\%s_comparision_top_var.pdf",
                    write_path,
                    size_name,
                    conf_name,
                    df_name),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]]) 
        
        plot(gbox_top_var)
        dev.off()
        
        write.csv(file=sprintf("%s\\%s\\%s\\%s_comparision.csv",
                               write_path,
                               size_name,
                               conf_name,
                               df_name),
                  statistics)
        
        
        write.csv(file=sprintf("%s\\%s\\%s\\%s_comparision_top_var.csv",
                               write_path,
                               size_name,
                               conf_name,
                               df_name),
                  wilcox_statistics)
        # 
        # pdf(sprintf("%s\\%s\\window_%d\\panels_comp_by_time_%s.pdf",
        #             write_path,
        #             size_name,
        #             max_window_size,
        #             df_name),
        #     height=sizes[[size_name]][["width"]] * 1.5,
        #     width=sizes[[size_name]][["width"]] * 3) 
        # 
        # plot(panels_plot)
        # dev.off()
        # 
        # pdf(sprintf("%s\\%s\\window_%d\\boxplot_%s.pdf",
        #             write_path,
        #             size_name,
        #             max_window_size,
        #             df_name),
        #     height=sizes[[size_name]][["width"]],
        #     width=sizes[[size_name]][["width"]]) 
        # 
        # plot(gbox)
        # dev.off()
      }
    }
  }
  
  # for (conf_name in names(all_results)) {
  #   pval_dfs <- all_results[[conf_name]]
  #   
  #   for (pval_df_name in names(pval_dfs)) {
  #     df <- pval_dfs[[pval_df_name]] 
  #     
  #     
  #     colnames(df)  <-  Responses
  #     
  #     melted_df <- melt(df)
  #     colnames(melted_df) <- c("#", "Variable", "Pvalue")
  #     # geom_violin <- 
  #     #   ggplot(melted_df,aes(x=Variable, y=Pvalue)) + 
  #     #   geom_violin() +
  #     #   geom_dotplot(stackdir="center", binaxis='y', binwidth=.02) + 
  #     #   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #     #                geom="pointrange", color="red") +
  #     #   scale_y_continuous(expand=c(0,0))  +
  #     #   xlab("Variable") +
  #     #   ylab("Shuffles defeated (%)") +
  #     #   theme(text=element_text(size=14, color="black"),
  #     #         legend.position="NA",
  #     #         axis.text=element_text(size=14, color="black"),
  #     #         axis.ticks=element_line(color="black"))  +
  #     #   theme_classic() +
  #     #   ggtitle(pval_df_name)
  #     # 
  #     # gbox <- 
  #     #   ggplot(melted_df, aes(x=Variable,y=Pvalue, group=Variable)) +
  #     #   geom_boxplot(width=0.5, 
  #     #                size=1, 
  #     #                color=adjustcolor("gray65", alpha=1), 
  #     #                fill=adjustcolor("gray80", alpha=0.9)) + 
  #     #   geom_jitter(aes(fill=Variable),
  #     #               position=position_jitterdodge(.5),
  #     #               alpha=0.8,
  #     #               stroke=.5,
  #     #               color=adjustcolor("gray20", alpha=.2),shape=21,
  #     #               size=1.2) + 
  #     #   
  #     #   #scale_color_brewer(palette="Spectral") + 
  #     #   #scale_color_brewer(palette="Spectral") + 
  #     #   #scale_color_manual(values=adjustcolor(spec_cg(ncol(df)), alpha=.75)) + 
  #     #   scale_fill_manual(values=adjustcolor(spec_cg(ncol(df)), alpha=.5)) + 
  #     #   
  #     #   theme_classic() +
  #     #   xlab("Variable") +
  #     #   ylab("Shuffles defeated (%)") +
  #     #   theme(text=element_text(size=14, color="black"),
  #     #         #legend.position="NA",
  #     #         legend.title = element_blank(),
  #     #         axis.text = element_text(color="black"),
  #     #         axis.text.x = element_blank(),#element_text(angle = 45, vjust=.6),
  #     #         axis.ticks = element_line(color="black")) 
  #     
  #     gmelted <- 
  #       ggplot(melted_df, aes(y=Pvalue, x=Variable)) +
  #       scale_y_continuous(expand=c(0,0)) +
  #       geom_bar(stat="summary", width=.45,
  #                fill="gray70") +
  #       geom_errorbar(stat="summary", width = .3) +
  #       geom_jitter(position=position_jitterdodge(.5), aes(fill=Variable),
  #                   size=1.5) + 
  #       theme_classic() +
  #       theme(text=element_text(size=14, color="black"), 
  #             axis.text=element_text(size=14, color="black"),
  #             axis.ticks = element_line(color="black"),
  #             plot.title=element_text(size=10),
  #             legend.position="NA")  + 
  #       xlab("") +
  #       ylab("Shuffles defeated (%)")
  #     
  #     for (size_name in names(sizes)) {
  #       
  #       dir.create(sprintf("%s\\%s",write_path, size_name))
  #       dir.create(sprintf("%s\\%s\\%s",write_path, size_name, conf_name))
  #       
  #       
  #       pdf(sprintf("%s\\%s\\%s\\shuffles_distribution_%s.pdf",
  #                   write_path,
  #                   size_name,
  #                   conf_name,
  #                   pval_df_name),
  #           height=sizes[[size_name]][["width"]],
  #           width=sizes[[size_name]][["width"]])
  #       
  #       plot(gmelted)
  #       dev.off()
  #       
  #       # pdf(sprintf("%s\\%s\\shuffles-distribution_%s_%s_box.pdf",
  #       #             write_path,
  #       #             size_name,
  #       #             conf_name,
  #       #             pval_df_name),
  #       #     height=sizes[[size_name]][["width"]],
  #       #     width=sizes[[size_name]][["width"]] * 3)
  #       # 
  #       # plot(gbox)
  #       # dev.off()
  #       # 
  #       # 
  #       # pdf(sprintf("%s\\%s\\shuffle_distribution_%s_%s_bar.pdf",
  #       #             write_path,
  #       #             size_name,
  #       #             conf_name,
  #       #             pval_df_name),
  #       #     height=sizes[[size_name]][["width"]],
  #       #     width=sizes[[size_name]][["width"]] * 3)
  #       # 
  #       # plot(gbox)
  #       # dev.off()
  #     }  
  #   }
  # }    
}

figure_3_metadata_annotation <- function(metadata, metadata_names) {
  
  write_path <- sprintf("%s\\figure_3\\", figures_base_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\structure_annotation", write_path)
  dir.create(write_path)  
  
  
  sizes = list(big=c(width=2.5))
  # big_1=c(width=2.2),
  # medium=c(width=2),
  # medium_2=c(width=1.5),
  # medium_1=c(width=1.25),
  # small=c(width=1))
  
  
  structure_sizes <-  list(small=c(width=.75, height=.75))
  pupil = F
  thirst_q = F
  #metadata_names <- get_datasets_names(control_paths, sep="_", control=T)
  write_path_all <- write_path
  for (idx in 1:len(metadata)) { 
    
    if (pupil && thirsty_q) {
      if (idx == 4 || idx == 14 || idx == 7 ) {
        next
      }
    }
    write_path <- write_path_all
    write_path <- sprintf("%s\\%s", write_path, metadata_names[[idx]])
    dir.create(write_path)
    
    cluster_color_label <- spec_cg(len(unique(metadata[[idx]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[idx]]$cluster_mat$labs)) - 1))
    
    
    base_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color[unlist(lapply(metadata[[idx]]$stim_master_mat[!is.nan(metadata[[idx]]$stim_master_mat[,"Reward"]),"Reward"], function(i) {get_binned_index(i,15)}))] <- 1
    cumsum_color <- cumsum(cumsum_color)
    
    reward_col = base_color
    reward_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                               metadata[[idx]]$stim_master_mat[,"Response"] == 0,"Reward"], function(i) {get_binned_index(i,15)}))] <- 1
    
    
    
    
    cumsum_color_f = adjustcolor(rev(rdylbu_cg(max(cumsum_color+1))[cumsum_color + 1]), alpha=.3)
    reward_col[reward_col == 0] <- "gray70"
    
    reward_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(reward_col))
    reward_col_f[reward_col == 1] <- "red"
    #reward_col_f <- adjustcolor(reward_col, alpha=.3)
    
    reward_trials <- metadata[[idx]]$annot_df[,1] == 3
    sorted_trials <-  order(paste(metadata[[idx]]$annot_df[,1], metadata[[idx]]$annot_df[,2]))
    
    ongoing_trials <- which(metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(4,5) & metadata[[idx]]$stim_master_mat[,"Response"]%%2 == 0)
    ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[idx]]$stim_master_mat)]
    
    ongoing_trials_frames = metadata[[idx]]$stim_master_mat[ongoing_trials + 1,1]
    
    ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - 6):(get_binned_index(i, 15)-1)})))
    
    ongoing_col <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    ongoing_col[ongoing_ind] <-  adjustcolor(c(rdylbu_cg(max(cumsum_color+1))[cumsum_color + 1]), alpha=.7)[ongoing_ind]
    
    
    if (pupil) {
      pupil_vec_all <- get_pupil_files(paths[[idx]], window_size=15)
      pupil_col <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
      pupil_col[ongoing_ind] <-  adjustcolor(c(spec_cg(8)[as.numeric(cut(pupil_vec_all$smoothed, breaks=8))]), alpha=.7)[ongoing_ind]
    }
    
    clusters_col_all=cluster_color_label[as.character(metadata[[idx]]$cluster_mat$labs)]
    clusters_col_ongoing <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    clusters_col_ongoing[ongoing_ind] <- clusters_col_all[ongoing_ind]
    
    trials_mat_reward_only <- pheatmap(metadata[[idx]]$trials_mat[reward_trials,], 
                                       cluster_rows=F, 
                                       cluster_cols=F, 
                                       col=cluster_color_label, 
                                       annotation_row = metadata[[idx]]$annot_df[reward_trials,])
    
    annoCol<-list(trialtype=c(`1`="blue", `3`="red", `4`="gray70", `5`="grey60"), response=c(`0`="white",`1`="black"))
    trials_mat_clustered <- pheatmap(metadata[[idx]]$trials_mat,
                                     cluster_rows=T, 
                                     cluster_cols=F, 
                                     col=cluster_color_label, 
                                     annotation_row = metadata[[idx]]$annot_df, 
                                     annotation_colors =  annoCol)
    
    
    trials_mat_sorted <- pheatmap(metadata[[idx]]$trials_mat[sorted_trials,], 
                                  cluster_rows=F, 
                                  cluster_cols=F, 
                                  col=cluster_color_label, 
                                  annotation_row = metadata[[idx]]$annot_df[sorted_trials,])
    
    trials_mat_all <- pheatmap(metadata[[idx]]$trials_mat, 
                               cluster_rows=F, 
                               cluster_cols=F, 
                               col=cluster_color_label, 
                               annotation_row = metadata[[idx]]$annot_df)
    
    
    
    
    
    stim_mat_for_trials_mat <-  metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    
    original_ongoing_trials <- which(metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(4,5) & metadata[[idx]]$stim_master_mat[,"Response"] %% 2 == 0)
    original_ongoing_trials <- original_ongoing_trials[original_ongoing_trials < nrow(metadata[[idx]]$stim_master_mat)]
    
    ongoing_trials <- which(stim_mat_for_trials_mat[,"TrialType"] %in% c(4,5) & stim_mat_for_trials_mat[,"Response"]%%2 == 0)
    ongoing_trials <- ongoing_trials[ongoing_trials < nrow(stim_mat_for_trials_mat)]
    
    
    
    reward_succeeded_by_ongoing <- 
      ongoing_trials[metadata[[idx]]$stim_master_mat[original_ongoing_trials - 1,"TrialType"] == 3 & 
                       metadata[[idx]]$stim_master_mat[original_ongoing_trials - 1,"Response"] == 0] - 1
    
    extended_all <- metadata[[idx]]$extended_trials_mat
    
    
    just_reward_heatmap <- 
      pheatmap(extended_all[reward_succeeded_by_ongoing,], cluster_rows=F, cluster_cols=F, border_color = NA, col=cluster_color_label)
    extended_heatmap <- 
      pheatmap(metadata[[idx]]$trials_mat[ongoing_trials,], cluster_rows=F, cluster_cols=F, border_color = NA, col=cluster_color_label)
    
    
    for (size_name in names(sizes)) {
      dir.create(sprintf("%s\\%s",
                         write_path,
                         size_name))
      pdf(sprintf("%s\\%s\\trials_mat_all.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 1.5)
          # unit="in",
          # res=1500)
      
      plot(trials_mat_all[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\trials_mat_reward.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 1.5)
          # unit="in",
          # res=1500)
      
      plot(trials_mat_reward_only[[4]])
      dev.off()
      
      
      pdf(sprintf("%s\\%s\\trials_mat_sorted.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 1.5)
          # unit="in",
          # res=1500)
      
      plot(trials_mat_sorted[[4]])
      dev.off()
      
      
      pdf(sprintf("%s\\%s\\trials_mat_clustered.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 1.5)
          # unit="in",
          # res=1500)
      
      plot(trials_mat_clustered[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\extended_ITI_heatmap.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 3)
          # unit="in",
          # res=1500)
      
      plot(just_reward_heatmap[[4]])
      dev.off()
      
      pdf(sprintf("%s\\%s\\extended_all_reward.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]] * 3,
          width=sizes[[size_name]][["width"]] * 3)
          # unit="in",
          # res=1500)
      
      plot(extended_heatmap[[4]])
      dev.off()
    }
    
    if (pupil) {
      color_list <- list(pupil=pupil_col, ongoing=ongoing_col, reward=reward_col_f, cumsum=cumsum_color_f, clusters=clusters_col_all,
                         ongoing_by_cluster=clusters_col_ongoing)
    } else {
      color_list <- list(ongoing=ongoing_col, reward=reward_col_f, cumsum=cumsum_color_f, clusters=clusters_col_all,
                         ongoing_by_cluster=clusters_col_ongoing)  
    }
    
    
    for (col_name in names(color_list)) {
      for (main_cluster in c(T,F)) {
        col_to_use <- color_list[[col_name]]
        reduced_mat <- metadata[[idx]]$red_mat

        if (main_cluster) {
          main_cluster_ind <- metadata[[idx]]$cluster_mat$labs == -1
          reduced_mat <- reduced_mat[main_cluster_ind,]
          col_to_use <- col_to_use[main_cluster_ind]
          print(sprintf("Using main cluster for color :%s", col_name))
        }

        p_lem2_list <- pair_plots_colored(reduced_mat,
                                          col_to_use,
                                          return_grid = F,
                                          pt_size=.25,
                                          stroke_alpha = 0.0001)

        p_lem2_list_f <- lapply(1:len(p_lem2_list),
                                function(i) {return(p_lem2_list[[i]] +
                                                      theme(#line=element_blank(),
                                                        #rect=element_rect(color="white"),
                                                        plot.margin=margin(t = 2, r = 2, b = 0, l = 0, unit = "pt"),
                                                        panel.border=element_blank(),
                                                        plot.background = element_rect("white"),
                                                        axis.title=element_blank()))})

        p_lem2_5_no_annot <- p_lem2_list_f

        p_lem2_5_no_annot$nrow <- 5
        pf_5 <- do.call(arrangeGrob, p_lem2_5_no_annot)



        for (size_name in names(structure_sizes)) {


          png(sprintf("%s\\%s%s.png",
                      write_path,
                      ifelse(main_cluster,"Main_", ""),
                      col_name),
              height=structure_sizes[[size_name]][["height"]] * 5,
              width=structure_sizes[[size_name]][["width"]] * 3,
              unit="in",
              res=1500)

          plot(pf_5)
          dev.off()
        }
      }
    }
  }
}
