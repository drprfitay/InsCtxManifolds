

figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\figures_structure_parameterization\\", base_output_path) 

dir.create(base_output_path)
dir.create(output_path)

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_HS <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
metadata_control <- across_mice_decoding_build_metadata(get_control_paths(), control=T)

Responses = c("0"="Hit", "1"="Miss", "2"="CR (N)", "3"="FA (N)", "4"="CR (A)", "5"="FA (A)",  "15"="Pupil", "16"="Consumed water")

manifold_parameterization <- function()
{
   
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  calculate_entropy <- function(pmf) {- sum(pmf * log2(pmf + 10 ** -30))}
  
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
  
  metadata <- metadata_TT
  
  for (max_window_size in max_windows){
    for (offset in rev(begin_offset)) {
      
      entropy_window_size = max_window_size - offset
      all_shuffle_matrices <- list()
      all_pvals <- list()
      
      for (path_idx in c(1:len(paths))) {
        
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
  
  dir.create(sprintf("%s\\data\\figure_2\\", base_output_path))
  
  for (conf_name in names(all_results)) {
    param_comp <- all_results[[conf_name]]
    param_pvalue <- all_pval_results[[conf_name]]
    save(file=sprintf("%s\\data\\figure_2\\%s_TT_parameterization_comp.Rda", base_output_path, conf_name), param_comp)
    save(file=sprintf("%s\\data\\figure_2\\%s_TT_parameterization_pval.Rda", base_output_path, conf_name), param_pvalue)    
  }
}

plot_TT_structure_parameterization_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\TT_structure_parameterization_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
  
  sizes = list(medium=c(width=1,
                        height=1))
               # small=c(width=.75,
               #         height=.75))
  
  metadata <- metadata_TT
  additive_secs_before = 2
  additive_secs_after = 5
  
  for (idx in c(1:len(paths))) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[idx]])
    
    dir.create(figure_path)
    
    
    cluster_color_label <- spec_cg(len(unique(metadata[[idx]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[idx]]$cluster_mat$labs)) - 1))
    
    
    base_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color[unlist(lapply(metadata[[idx]]$stim_master_mat[!is.nan(metadata[[idx]]$stim_master_mat[,"Reward"]),"Reward"], function(i) {get_binned_index(i,15)}))] <- 1
    cumsum_color <- cumsum(cumsum_color)
    
    
    reward_col = base_color
    reward_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 0,"Frames"], 
                                                             function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    miss_col = base_color
    miss_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                           metadata[[idx]]$stim_master_mat[,"Response"] == 1,"Frames"], 
                                                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    CR_A_col = base_color
    CR_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                           metadata[[idx]]$stim_master_mat[,"Response"] == 4, "Frames"], 
                                                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_A_col = base_color
    FA_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                           metadata[[idx]]$stim_master_mat[,"Response"] == 5, "Frames"],
                                                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    CR_N_col = base_color
    CR_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                           metadata[[idx]]$stim_master_mat[,"Response"] == 2, "Frames"],
                                                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_N_col = base_color
    FA_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                           metadata[[idx]]$stim_master_mat[,"Response"] == 3, "Frames"],
                                                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    
    
    reward_col[reward_col == 0] <- "gray70"
    reward_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(reward_col))
    reward_col_f[reward_col == 1] <- "red"
    
    miss_col[miss_col == 0] <- "gray70"
    miss_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(miss_col))
    miss_col_f[miss_col == 1] <- "red"
    
    CR_A_col[CR_A_col == 0] <- "gray70"
    CR_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_A_col))
    CR_A_col_f[CR_A_col == 1] <- "red"
    
    FA_A_col[FA_A_col == 0] <- "gray70"
    FA_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_A_col))
    FA_A_col_f[FA_A_col == 1] <- "red"
    
    CR_N_col[CR_N_col == 0] <- "gray70"
    CR_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_N_col))
    CR_N_col_f[CR_N_col == 1] <- "red"
    
    FA_N_col[FA_N_col == 0] <- "gray70"
    FA_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_N_col))
    FA_N_col_f[FA_N_col == 1] <- "red"
    
    
    
    
    cumsum_color_f = adjustcolor((viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.3)
    #reward_col_f <- adjustcolor(reward_col, alpha=.3)
    
    reward_trials <- metadata[[idx]]$annot_df[,1] == 3
    sorted_trials <-  order(paste(metadata[[idx]]$annot_df[,1], metadata[[idx]]$annot_df[,2]))
    
    ongoing_trials <- which(metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(4,5) & metadata[[idx]]$stim_master_mat[,"Response"]%%2 == 0)
    ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[idx]]$stim_master_mat)]
    
    ongoing_trials_frames = metadata[[idx]]$stim_master_mat[ongoing_trials + 1,1]
    
    ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - 6):(get_binned_index(i, 15)-1)})))
    
    ongoing_col <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    ongoing_col[ongoing_ind] <-  adjustcolor(c(viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.7)[ongoing_ind]
    
    
    pupil_vec_all <- get_pupil_files(paths[[idx]], window_size=15, load_all = T)
    plot(pupil_vec_all$smoothed)
    pupil_col_ongoing <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    pupil_col_ongoing[ongoing_ind] <-  adjustcolor(c(viridis(7)[as.numeric(cut(pupil_vec_all$smoothed, breaks=7))]), alpha=.7)[ongoing_ind]
    pupil_col_all <-  adjustcolor(c(viridis(7)[as.numeric(cut(pupil_vec_all$smoothed, breaks=7))]), alpha=.7)
    
    clusters_col_all=cluster_color_label[as.character(metadata[[idx]]$cluster_mat$labs)]
    clusters_col_ongoing <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    clusters_col_ongoing[ongoing_ind] <- clusters_col_all[ongoing_ind]
    
    reduced_mat <- metadata[[idx]]$red_mat
    
    color_list <- list(#pupil=pupil_col_all,
                       #pupil_col_ongoing=pupil_col_ongoing)
                       #ongoing=ongoing_col, 
                       reward=reward_col_f,
                       miss=miss_col_f,
                       CR_A=CR_A_col_f,
                       FA_A=FA_A_col_f,
                       CR_N=CR_N_col_f,
                       FA_N=FA_N_col_f)
                       # cumsum=cumsum_color_f,
                       # clusters=clusters_col_all,
                       # ongoing_by_cluster=clusters_col_ongoing)
    
    
    for (col_name in names(color_list)) {
      col_to_use <- color_list[[col_name]]
      
      
      
      p_lem2_list <- pair_plots_colored(reduced_mat, 
                                        col_to_use,
                                        return_grid = F,
                                        pt_size=.25, 
                                        stroke_alpha = 0.0001,
                                        just_annotation = F)  
      
      
      p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                              col_to_use,
                                              return_grid = F,
                                              pt_size=.55, 
                                              stroke_alpha = 0.0001,
                                              just_annotation = T)  
      
      p_lem2_list <- lapply(1:len(p_lem2_list), 
                            function(i) {return(p_lem2_list[[i]] + 
                                                  theme(line=element_blank(),
                                                        #rect=element_rect(color="white"),
                                                        plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                        panel.border=element_rect(color="white"),
                                                        plot.background = element_rect("white"),
                                                        axis.title=element_text(color="white")))})
      
      
      
      
      for (size_name in names(sizes)) {
        
        
        
        p_lem2_no_annot <- p_lem2_list
        p_lem2_no_annot$nrow <- 3
        pf_structure <- do.call(plot_grid, p_lem2_no_annot)
        
        
        p_lem2_annot_only <- p_lem2_list_annot
        p_lem2_annot_only$nrow <- 3
        pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
        
        
        
        png(sprintf("%s\\%s_%s_lem2_example_no_annot.png",
                    figure_path,
                    size_name,
                    col_name),
            height=sizes[[size_name]][["height"]] * 3,
            width=sizes[[size_name]][["width"]] * 5,
            unit="in",
            res=1500)
        
        plot(pf_structure)
        dev.off()
        
        # Just once
        if (col_name == "pupil") {
          
          pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                      figure_path,
                      size_name),
              height=sizes[[size_name]][["height"]] * 3,
              width=sizes[[size_name]][["width"]] * 5)
          
          plot(pf_annot_only)
          dev.off()
        }
        
      }
    }
  }
}

plot_HS_structure_parameterization_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\HS_structure_parameterization_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_hungry_sated_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
  
  sizes = list(medium=c(width=1,
                        height=1))
  # small=c(width=.75,
  #         height=.75))
  
  metadata <- metadata_HS
  additive_secs_before = 2
  additive_secs_after = 5
  
  for (idx in c(1:len(paths))) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[idx]])
    
    dir.create(figure_path)
    
    
    cluster_color_label <- spec_cg(len(unique(metadata[[idx]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[idx]]$cluster_mat$labs)) - 1))
    
    base_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color[unlist(lapply(metadata[[idx]]$stim_master_mat[!is.nan(metadata[[idx]]$stim_master_mat[,"Reward"]),"Reward"], function(i) {get_binned_index(i,15)}))] <- 1
    cumsum_color <- cumsum(cumsum_color)
    
    
    reward_col = base_color
    reward_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                               metadata[[idx]]$stim_master_mat[,"Response"] == 0,"Frames"], 
                             function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    miss_col = base_color
    miss_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 1,"Frames"], 
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    CR_A_col = base_color
    CR_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 4, "Frames"], 
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_A_col = base_color
    FA_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 5, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    CR_N_col = base_color
    CR_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 2, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_N_col = base_color
    FA_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 3, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    
    reward_col[reward_col == 0] <- "gray70"
    reward_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(reward_col))
    reward_col_f[reward_col == 1] <- "red"
    
    miss_col[miss_col == 0] <- "gray70"
    miss_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(miss_col))
    miss_col_f[miss_col == 1] <- "red"
    
    CR_A_col[CR_A_col == 0] <- "gray70"
    CR_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_A_col))
    CR_A_col_f[CR_A_col == 1] <- "red"
    
    FA_A_col[FA_A_col == 0] <- "gray70"
    FA_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_A_col))
    FA_A_col_f[FA_A_col == 1] <- "red"
    
    CR_N_col[CR_N_col == 0] <- "gray70"
    CR_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_N_col))
    CR_N_col_f[CR_N_col == 1] <- "red"
    
    FA_N_col[FA_N_col == 0] <- "gray70"
    FA_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_N_col))
    FA_N_col_f[FA_N_col == 1] <- "red"
    
    
    
    
    cumsum_color_f = adjustcolor((viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.3)
    #reward_col_f <- adjustcolor(reward_col, alpha=.3)
    
    reward_trials <- metadata[[idx]]$annot_df[,1] == 3
    sorted_trials <-  order(paste(metadata[[idx]]$annot_df[,1], metadata[[idx]]$annot_df[,2]))
    
    ongoing_trials <- which(metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(4,5) & metadata[[idx]]$stim_master_mat[,"Response"]%%2 == 0)
    ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[idx]]$stim_master_mat)]
    
    ongoing_trials_frames = metadata[[idx]]$stim_master_mat[ongoing_trials + 1,1]
    
    ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - 6):(get_binned_index(i, 15)-1)})))
    
    ongoing_col <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    ongoing_col[ongoing_ind] <-  adjustcolor(c(viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.7)[ongoing_ind]
    
    clusters_col_all=cluster_color_label[as.character(metadata[[idx]]$cluster_mat$labs)]
    clusters_col_ongoing <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    clusters_col_ongoing[ongoing_ind] <- clusters_col_all[ongoing_ind]
    
    reduced_mat <- metadata[[idx]]$red_mat
    
    color_list <- list(ongoing=ongoing_col, 
                       reward=reward_col_f,
                       miss=miss_col_f,
                       CR_A=CR_A_col_f,
                       FA_A=FA_A_col_f,
                       CR_N=CR_N_col_f,
                       FA_N=FA_N_col_f,
                       cumsum=cumsum_color_f,
                       clusters=clusters_col_all,
                       ongoing_by_cluster=clusters_col_ongoing)
    
    
    for (col_name in names(color_list)) {
      col_to_use <- color_list[[col_name]]
      
      
      
      p_lem2_list <- pair_plots_colored(reduced_mat, 
                                        col_to_use,
                                        return_grid = F,
                                        pt_size=.25, 
                                        stroke_alpha = 0.0001,
                                        just_annotation = F)  
      
      
      p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                              col_to_use,
                                              return_grid = F,
                                              pt_size=.55, 
                                              stroke_alpha = 0.0001,
                                              just_annotation = T)  
      
      p_lem2_list <- lapply(1:len(p_lem2_list), 
                            function(i) {return(p_lem2_list[[i]] + 
                                                  theme(line=element_blank(),
                                                        #rect=element_rect(color="white"),
                                                        plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                        panel.border=element_rect(color="white"),
                                                        plot.background = element_rect("white"),
                                                        axis.title=element_text(color="white")))})
      
      
      
      
      for (size_name in names(sizes)) {
        
        
        
        p_lem2_no_annot <- p_lem2_list
        p_lem2_no_annot$nrow <- 3
        pf_structure <- do.call(plot_grid, p_lem2_no_annot)
        
        
        p_lem2_annot_only <- p_lem2_list_annot
        p_lem2_annot_only$nrow <- 3
        pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
        
        
        
        png(sprintf("%s\\%s_%s_lem2_example_no_annot.png",
                    figure_path,
                    size_name,
                    col_name),
            height=sizes[[size_name]][["height"]] * 3,
            width=sizes[[size_name]][["width"]] * 5,
            unit="in",
            res=1500)
        
        plot(pf_structure)
        dev.off()
        
        # Just once
        if (col_name == "pupil") {
          
          pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                      figure_path,
                      size_name),
              height=sizes[[size_name]][["height"]] * 3,
              width=sizes[[size_name]][["width"]] * 5)
          
          plot(pf_annot_only)
          dev.off()
        }
        
      }
    }
  }
}

plot_control_structure_parameterization_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\control_structure_parameterization_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_control_paths()
  datasets_names <- get_datasets_names(paths, sep="_", control=T)
  
  sizes = list(medium=c(width=1,
                        height=1))
  # small=c(width=.75,
  #         height=.75))
  
  metadata <- metadata_control
  additive_secs_before = 2
  additive_secs_after = 5
  
  for (idx in c(1:len(paths))) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[idx]])
    
    dir.create(figure_path)
    
    
    cluster_color_label <- spec_cg(len(unique(metadata[[idx]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[idx]]$cluster_mat$labs)) - 1))
    
    base_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color <- rep(0, times=nrow(metadata[[idx]]$red_mat))
    cumsum_color[unlist(lapply(metadata[[idx]]$stim_master_mat[!is.nan(metadata[[idx]]$stim_master_mat[,"Reward"]),"Reward"], function(i) {get_binned_index(i,15)}))] <- 1
    cumsum_color <- cumsum(cumsum_color)
    
    
    reward_col = base_color
    reward_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                               metadata[[idx]]$stim_master_mat[,"Response"] == 0,"Frames"], 
                             function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    miss_col = base_color
    miss_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 3 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 1,"Frames"], 
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    CR_A_col = base_color
    CR_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 4, "Frames"], 
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_A_col = base_color
    FA_A_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 4 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 5, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    CR_N_col = base_color
    CR_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 2, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    FA_N_col = base_color
    FA_N_col[unlist(lapply(metadata[[idx]]$stim_master_mat[metadata[[idx]]$stim_master_mat[,"TrialType"] == 5 & 
                                                             metadata[[idx]]$stim_master_mat[,"Response"] == 3, "Frames"],
                           function(i) {idx <- get_binned_index(i,15); return((idx+additive_secs_before):(idx+additive_secs_after))}))] <- 1
    
    
    reward_col[reward_col == 0] <- "gray70"
    reward_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(reward_col))
    reward_col_f[reward_col == 1] <- "red"
    
    miss_col[miss_col == 0] <- "gray70"
    miss_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(miss_col))
    miss_col_f[miss_col == 1] <- "red"
    
    CR_A_col[CR_A_col == 0] <- "gray70"
    CR_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_A_col))
    CR_A_col_f[CR_A_col == 1] <- "red"
    
    FA_A_col[FA_A_col == 0] <- "gray70"
    FA_A_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_A_col))
    FA_A_col_f[FA_A_col == 1] <- "red"
    
    CR_N_col[CR_N_col == 0] <- "gray70"
    CR_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(CR_N_col))
    CR_N_col_f[CR_N_col == 1] <- "red"
    
    FA_N_col[FA_N_col == 0] <- "gray70"
    FA_N_col_f <- rep(adjustcolor("gray70", alpha=.2), times=len(FA_N_col))
    FA_N_col_f[FA_N_col == 1] <- "red"
    
    
    
    
    cumsum_color_f = adjustcolor((viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.3)
    #reward_col_f <- adjustcolor(reward_col, alpha=.3)
    
    reward_trials <- metadata[[idx]]$annot_df[,1] == 3
    sorted_trials <-  order(paste(metadata[[idx]]$annot_df[,1], metadata[[idx]]$annot_df[,2]))
    
    ongoing_trials <- which(metadata[[idx]]$stim_master_mat[,"TrialType"] %in% c(4,5) & metadata[[idx]]$stim_master_mat[,"Response"]%%2 == 0)
    ongoing_trials <- ongoing_trials[ongoing_trials < nrow(metadata[[idx]]$stim_master_mat)]
    
    ongoing_trials_frames = metadata[[idx]]$stim_master_mat[ongoing_trials + 1,1]
    
    ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - 6):(get_binned_index(i, 15)-1)})))
    
    ongoing_col <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    ongoing_col[ongoing_ind] <-  adjustcolor(c(viridis(max(cumsum_color+1))[cumsum_color + 1]), alpha=.7)[ongoing_ind]
    
    clusters_col_all=cluster_color_label[as.character(metadata[[idx]]$cluster_mat$labs)]
    clusters_col_ongoing <- rep(adjustcolor("gray70", alpha=.05), times=len(base_color))
    clusters_col_ongoing[ongoing_ind] <- clusters_col_all[ongoing_ind]
    
    reduced_mat <- metadata[[idx]]$red_mat
    
    color_list <- list(ongoing=ongoing_col, 
                       reward=reward_col_f,
                       miss=miss_col_f,
                       CR_A=CR_A_col_f,
                       FA_A=FA_A_col_f,
                       CR_N=CR_N_col_f,
                       FA_N=FA_N_col_f,
                       cumsum=cumsum_color_f,
                       clusters=clusters_col_all,
                       ongoing_by_cluster=clusters_col_ongoing)
    
    
    for (col_name in names(color_list)) {
      col_to_use <- color_list[[col_name]]
      
      
      
      p_lem2_list <- pair_plots_colored(reduced_mat, 
                                        col_to_use,
                                        return_grid = F,
                                        pt_size=.25, 
                                        stroke_alpha = 0.0001,
                                        just_annotation = F)  
      
      
      p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                              col_to_use,
                                              return_grid = F,
                                              pt_size=.55, 
                                              stroke_alpha = 0.0001,
                                              just_annotation = T)  
      
      p_lem2_list <- lapply(1:len(p_lem2_list), 
                            function(i) {return(p_lem2_list[[i]] + 
                                                  theme(line=element_blank(),
                                                        #rect=element_rect(color="white"),
                                                        plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                        panel.border=element_rect(color="white"),
                                                        plot.background = element_rect("white"),
                                                        axis.title=element_text(color="white")))})
      
      
      
      
      for (size_name in names(sizes)) {
        
        
        
        p_lem2_no_annot <- p_lem2_list
        p_lem2_no_annot$nrow <- 3
        pf_structure <- do.call(plot_grid, p_lem2_no_annot)
        
        
        p_lem2_annot_only <- p_lem2_list_annot
        p_lem2_annot_only$nrow <- 3
        pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
        
        
        
        png(sprintf("%s\\%s_%s_lem2_example_no_annot.png",
                    figure_path,
                    size_name,
                    col_name),
            height=sizes[[size_name]][["height"]] * 3,
            width=sizes[[size_name]][["width"]] * 5,
            unit="in",
            res=1500)
        
        plot(pf_structure)
        dev.off()
        
        # Just once
        if (col_name == "pupil") {
          
          pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                      figure_path,
                      size_name),
              height=sizes[[size_name]][["height"]] * 3,
              width=sizes[[size_name]][["width"]] * 5)
          
          plot(pf_annot_only)
          dev.off()
        }
        
      }
    }
  }
}

plot_TT_parameterization_boxplots <- function()
{
  base_figure_path <- sprintf("%s\\TT_parameterization_boxplots", output_path)
  dir.create(base_figure_path)

  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path)
  
  output_name = "w3_o6_TT_parameterization_comp.Rda"
  load(sprintf("%s\\data\\figure_2\\%s", base_output_path, output_name), verbose=T)
  param_boxplots <- param_comp
  
  add_lines_to_plot <- function(plot_df, plot, logit=F) {
    
    new_plot <- plot
    for (var in unique(plot_df$Var)) {
      xval <- which(levels(plot_df$Var) == var)
      var_plot_df <- plot_df[plot_df$Var == var,]
      
      for (point_idx in unique(var_plot_df[["#"]])) {
        point_comp <- var_plot_df[var_plot_df[["#"]] == point_idx,] 
        
        data_value <- point_comp[point_comp$Group == "Data","MI"]
        shuff_value <- point_comp[point_comp$Group == "Shuffle","MI"]
        
        if (logit) {
          data_value <- log2(data_value)
          shuff_value <- log2(shuff_value)
        }
        
        new_plot <- 
          new_plot + 
          geom_line(data=data.frame(x=c(xval-.25,
                                        xval+.25), 
                                    y=c(data_value,
                                        shuff_value)),
                    aes(x=x,y=y))
      }
    }
    
    return(new_plot)
  }
  
  
  MI_all <- do.call(rbind,lapply(param_boxplots, function(df) {df$shuffle_MI_all[1,]}))
  MI_main_all <- do.call(rbind,lapply(param_boxplots, function(df) {df$shuffle_MI_main_all[1,]}))
  shuffle_MI_all <- do.call(rbind,lapply(param_boxplots, function(df) {colMeans(df$shuffle_MI_all[-1,])}))
  shuffle_MI_main_all <- do.call(rbind,lapply(param_boxplots, function(df) {colMeans(df$shuffle_MI_main_all[-1,])}))
  
  colnames(MI_all) <- Responses
  colnames(MI_main_all) <- Responses
  colnames(shuffle_MI_all) <- Responses
  colnames(shuffle_MI_main_all) <- Responses
  
  melted_MI_all <- melt(MI_all)
  melted_MI_main_all <- melt(MI_main_all)
  melted_shuffle_MI_all <- melt(shuffle_MI_all)
  melted_shuffle_MI_main_all <- melt(shuffle_MI_main_all)
  
  colnames(melted_MI_all) <- c("#", "Var", "MI")
  colnames(melted_MI_main_all) <- c("#", "Var", "MI")
  colnames(melted_shuffle_MI_all) <- c("#", "Var", "MI")
  colnames(melted_shuffle_MI_main_all) <- c("#", "Var", "MI")
  
  
  melted_MI_all$Var <- factor(melted_MI_all$Var, levels=colnames(MI_all)[order(colMeans(MI_all), decreasing = T)])
  melted_MI_main_all$Var <- factor(melted_MI_main_all$Var, levels=colnames(MI_main_all)[order(colMeans(MI_main_all), decreasing = T)])
  melted_shuffle_MI_all$Var <- factor(melted_shuffle_MI_all$Var, colnames(MI_all)[order(colMeans(MI_all), decreasing = T)])
  melted_shuffle_MI_main_all$Var <- factor(melted_shuffle_MI_main_all$Var, colnames(MI_main_all)[order(colMeans(MI_main_all), decreasing = T)])

  
  gMI_all <- 
    ggplot(melted_MI_all, aes(x=Var, y=(MI))) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.5)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape = NA)+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  gMI_all_log2 <- 
    ggplot(melted_MI_all, aes(x=Var, y=log2(MI))) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.5)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape = NA)+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("log2(Mutual Information)") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  gMI_main_all <- 
    ggplot(melted_MI_main_all, aes(x=Var, y=(MI))) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.5)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape = NA)+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  gMI_main_all_log2 <- 
    ggplot(melted_MI_main_all, aes(x=Var, y=log2(MI))) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.5)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape = NA)+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("log2(Mutual Information)") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  
  MI_main_all_data_vs_shuffle <- rbind(melted_MI_main_all,melted_shuffle_MI_main_all)
  MI_main_all_data_vs_shuffle$Group <- rep(c("Data", "Shuffle"), times=c(nrow(melted_MI_main_all),nrow(melted_shuffle_MI_main_all)))
  MI_main_all_data_vs_shuffle <- MI_main_all_data_vs_shuffle[MI_main_all_data_vs_shuffle$Var %in% levels(melted_MI_main_all$Var)[1:3],]
  
  MI_all_data_vs_shuffle <- rbind(melted_MI_all,melted_shuffle_MI_all)
  MI_all_data_vs_shuffle$Group <- rep(c("Data", "Shuffle"), times=c(nrow(melted_MI_all),nrow(melted_shuffle_MI_all)))
  MI_all_data_vs_shuffle <- MI_all_data_vs_shuffle[MI_all_data_vs_shuffle$Var %in% levels(melted_MI_all$Var)[1:3],]
  
  
  
  gMI_main_all_data_vs_shuff <- 
    ggplot(MI_main_all_data_vs_shuffle, aes(x=Var, y=MI)) +
    geom_boxplot(outlier.shape = NA, position=position_dodge(1), aes(color=Group))
  
  
  gMI_main_all_data_vs_shuff <- add_lines_to_plot(MI_main_all_data_vs_shuffle, gMI_main_all_data_vs_shuff)
  gMI_main_all_data_vs_shuff <- 
    gMI_main_all_data_vs_shuff + 
    geom_point(size=1.75, position=position_dodge(1), aes(color=Group))+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  
  
  gMI_main_all_data_vs_shuff_log2 <- 
    ggplot(MI_main_all_data_vs_shuffle, aes(x=Var, y=log2(MI))) +
    geom_boxplot(outlier.shape = NA, position=position_dodge(1), aes(color=Group))
  
  
  gMI_main_all_data_vs_shuff_log2 <- add_lines_to_plot(MI_main_all_data_vs_shuffle, gMI_main_all_data_vs_shuff_log2, logit = T)
  gMI_main_all_data_vs_shuff_log2 <- 
    gMI_main_all_data_vs_shuff_log2 + 
    geom_point(size=1.75, position=position_dodge(1), aes(color=Group))+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  
  
  
  gMI_all_data_vs_shuff <- 
    ggplot(MI_all_data_vs_shuffle, aes(x=Var, y=(MI))) +
    geom_boxplot(outlier.shape = NA, position=position_dodge(1), aes(color=Group))
  
  
  gMI_all_data_vs_shuff <- add_lines_to_plot(MI_all_data_vs_shuffle, gMI_all_data_vs_shuff)
  gMI_all_data_vs_shuff <- 
    gMI_all_data_vs_shuff + 
    geom_point(size=1.75, position=position_dodge(1), aes(color=Group))+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
  
  
  gMI_all_data_vs_shuff_log2 <- 
    ggplot(MI_all_data_vs_shuffle, aes(x=Var, y=log2(MI))) +
    geom_boxplot(outlier.shape = NA, position=position_dodge(1), aes(color=Group))
  
  
  gMI_all_data_vs_shuff_log2 <- add_lines_to_plot(MI_all_data_vs_shuffle, gMI_all_data_vs_shuff_log2, logit = T)
  gMI_all_data_vs_shuff_log2 <- 
    gMI_all_data_vs_shuff_log2 + 
    geom_point(size=1.75, position=position_dodge(1), aes(color=Group))+
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Mutual Information") + 
    xlab("")  +
    theme(plot.title = element_text(size=9))
  
    
    ### Shuffle statistics
    pupil_s <- melted_shuffle_MI_all[melted_shuffle_MI_all$Var == "Pupil",]
    pupil_d <- melted_MI_all[melted_MI_all$Var == "Pupil",]
    pupil_main_s <- melted_shuffle_MI_main_all[melted_shuffle_MI_main_all$Var == "Pupil",]
    pupil_main_d <- melted_MI_main_all[melted_MI_main_all$Var == "Pupil",]    
    hit_s <- melted_shuffle_MI_all[melted_shuffle_MI_all$Var == "Hit",]
    hit_d <- melted_MI_all[melted_MI_all$Var == "Hit",]
    hit_main_s <- melted_shuffle_MI_main_all[melted_shuffle_MI_main_all$Var == "Hit",]
    hit_main_d <- melted_MI_main_all[melted_MI_main_all$Var == "Hit",]
    consumed_water_s <- melted_shuffle_MI_main_all[melted_shuffle_MI_all$Var == "Consumed water",]
    consumed_water_d <- melted_MI_all[melted_MI_all$Var == "Consumed water",]
    consumed_water_main_s <- melted_shuffle_MI_main_all[melted_shuffle_MI_main_all$Var == "Consumed water",]
    consumed_water_main_d <- melted_MI_main_all[melted_MI_main_all$Var == "Consumed water",]
    
    
    shuffle_statistics_df <- data.frame()
      
    comparisions <- list(list(pupil_s, pupil_d, F),
                         list(pupil_main_s, pupil_main_d, T),
                         list(hit_s, hit_d, F),
                         list(hit_main_s, hit_main_d, T),
                         list(consumed_water_s,consumed_water_d, F),
                         list(consumed_water_main_s, consumed_water_main_d, T))
    
    
    for (comp in comparisions) {
      df_shuffle <- comp[[1]]
      df_data <- comp[[2]]
      with_main <- comp[[3]]
      
      wilc <- wilcox.test(df_data$MI, df_shuffle$MI, paired=T, correct = F, alternative="greater")
      
      comp_df <-  data.frame(method=wilc$method,
                             statistic=wilc$statistic,
                             pvalue=wilc$p.value,
                             alternative=wilc$alternative,
                             mean_shuffle=mean(df_shuffle$MI),
                             mean_data=mean(df_data$MI),
                             sem_shuffle=sem(df_shuffle$MI),
                             sem_data=sem(df_data$MI),
                             sd_shuffle=sd(df_shuffle$MI),
                             sd_data=sd(df_data$MI),
                             N_shuffle=len(df_shuffle$MI),
                             N_data=len(df_data$MI),
                             var_shuffle=df_shuffle$Var[1],
                             var_data=df_data$Var[1],
                             signif=signif.num(wilc$p.value),
                             corrected_pval=wilc$p.value * 3,
                             corrected_signif = signif.num(wilc$p.value * 3),
                             with_main=with_main)
      
      shuffle_statistics_df <- rbind(shuffle_statistics_df,
                                     comp_df)
                             
    }
    
    
    ### Raw MI comp
    
    dt_no_main <- dunn.test::dunn.test(melted_MI_all$MI, melted_MI_all$Var, "bonferroni")
    dt_main <- dunn.test::dunn.test(melted_MI_main_all$MI, melted_MI_main_all$Var, "bonferroni")
    dt_no_main <- as.data.frame(dt_no_main)
    dt_main <- as.data.frame(dt_main)
    
    dt_no_main$signif.adjusted <- signif.num(dt_no_main$P.adjusted)
    dt_no_main$signif <- signif.num(dt_no_main$P)
    dt_main$signif.adjusted <- signif.num(dt_main$P.adjusted)
    dt_main$signif <- signif.num(dt_main$P)
    
    
    
    l1 <-  list()
    l2 <- list()
    
    for (j in c(1:8)) {
      with_main_res_df <- do.call(rbind,lapply(c(1:8), 
                                     function(i) {
                                       wilc <- wilcox.test(MI_main_all[,j], MI_main_all[,i], alternative="greater", correct = F, paired=T);
                                       df <- data.frame(comp=sprintf("%d-%d", j,i), 
                                                        pval=wilc$p.value,
                                                        corrected=wilc$p.value * 7,
                                                        statistic=wilc$statistic, 
                                                        alternative=wilc$alternative, 
                                                        method=wilc$method, 
                                                        mean_1=mean(MI_main_all[,j]),
                                                        mean_2=mean(MI_main_all[,i]),
                                                        sem_1=sem(MI_main_all[,j]),
                                                        sem_2=sem(MI_main_all[,i]),
                                                        N_1=len(MI_main_all[,j]),
                                                        N_2=len(MI_main_all[,i]),
                                                        sd_1=sd(MI_main_all[,j]),
                                                        sd_2=sd(MI_main_all[,i]));
                                       return(df)}))
      
      l1 <- append(l1, list(with_main_res_df))
      
      res_df <- do.call(rbind,lapply(c(1:8), 
                                      function(i) {
                                        wilc <- wilcox.test(MI_all[,j], MI_all[,i], alternative="greater", correct = F, paired=T);
                                          df <- data.frame(comp=sprintf("%d-%d", j,i), 
                                                           pval=wilc$p.value, 
                                                           corrected=wilc$p.value * 7,
                                                           statistic=wilc$statistic, 
                                                           alternative=wilc$alternative, 
                                                           method=wilc$method, 
                                                           mean_1=mean(MI_all[,j]),
                                                           mean_2=mean(MI_all[,i]),
                                                           sem_1=sem(MI_all[,j]),
                                                           sem_2=sem(MI_all[,i]),
                                                           N_1=len(MI_all[,j]),
                                                           N_2=len(MI_all[,i]),
                                                           sd_1=sd(MI_all[,j]),
                                                           sd_2=sd(MI_all[,i]));
                                          return(df)}))
      
      l2 <- append(l2, list(res_df))      
    }
    
    
    MI_main_all_wilcox_df <- as.data.frame(do.call(rbind, l1))
    MI_all_wilcox_df <- as.data.frame(do.call(rbind, l2))
    
    
    MI_all_wilcox_df$Comp <- unlist(lapply(MI_all_wilcox_df$comp, function(comp) {split <- str_split(comp, "-"); sprintf("%s - %s", Responses[[as.numeric(split[[1]][[1]])]], Responses[[as.numeric(split[[1]][[2]])]])}))
    MI_main_all_wilcox_df$Comp <- unlist(lapply(MI_main_all_wilcox_df$comp, function(comp) {split <- str_split(comp, "-"); sprintf("%s - %s", Responses[[as.numeric(split[[1]][[1]])]], Responses[[as.numeric(split[[1]][[2]])]])}))
    
    MI_all_wilcox_df$signif <- signif.num(MI_all_wilcox_df$pval)
    MI_main_all_wilcox_df$signif <- signif.num(MI_main_all_wilcox_df$pval)
    MI_all_wilcox_df$corrected_signif <- signif.num(MI_all_wilcox_df$corrected)
    MI_main_all_wilcox_df$corrected_signif <- signif.num(MI_main_all_wilcox_df$corrected)
    
    
    write.csv(file=sprintf("%s\\MI_wilcox_all.csv", statistics_figure_path),
              MI_all_wilcox_df)
    
    write.csv(file=sprintf("%s\\MI_main_wilcox_all.csv", statistics_figure_path),
              MI_main_all_wilcox_df)    
    
    write.csv(file=sprintf("%s\\shuffle_statistics.csv", statistics_figure_path),
              shuffle_statistics_df)    
    
    write.csv(file=sprintf("%s\\dunn_test_MI.csv", statistics_figure_path),
              as.data.frame(dt_no_main))    
    
    write.csv(file=sprintf("%s\\dunn_test_MI_main.csv", statistics_figure_path),
              as.data.frame(dt_main))    
    
    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      
      pdf(sprintf("%s\\%s_MI_all_data_vs_shuff.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_all_data_vs_shuff)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_all_data_vs_shuff_log2.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_all_data_vs_shuff_log2)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_main_all_data_vs_shuff.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_main_all_data_vs_shuff)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_main_all_data_vs_shuff_log2.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_main_all_data_vs_shuff_log2)
      dev.off()      

      pdf(sprintf("%s\\%s_MI_main_all.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_main_all)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_main_all_log2.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_main_all_log2)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_all.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_all)
      dev.off()
      
      pdf(sprintf("%s\\%s_MI_all_log2.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gMI_all_log2)
      dev.off()
    }
}

plot_cluster_cumulative_rewards <- function()
{
  
  
  base_figure_path <- sprintf("%s\\cluster_cumulative_rewards", output_path)
  dir.create(base_figure_path)
  
  cl1_all <- c()
  cl2_all <- c()
  reward_cl1_all <- c()
  reward_cl2_all <- c()
  
  cmsm_all <- c()
  cmsm_f_all <- c()
  
  for (idx in 1:len(metadata_TT)) {
    m_obj  <- metadata_TT[[idx]]
    
    
    used_trials <- m_obj$stim_master_mat[,"TrialType"] %in% c(1,3,4,5)
    stim_master_mat <- m_obj$stim_master_mat
    cumsum_all <- rep(0, times=nrow(stim_master_mat))
    drink_ind <- which(!is.nan(stim_master_mat[,6]))
    cumsum_all[drink_ind] <- 1
    
    cumsum_f <- cumsum(cumsum_all)
    
    used_cumsum <- cumsum_f[used_trials] / max(cumsum_f[used_trials])
    
    cmsm_all <- rbind(cmsm_all,
                      used_cumsum[floor(seq(1, len(used_cumsum), length.out=200))])
    
    cmsm_f_all <- rbind(cmsm_f_all,
                        cumsum_f[floor(seq(1, len(cumsum_f), length.out=200))] / max(cumsum_f[floor(seq(1, len(cumsum_f), length.out=200))]))
    
    zero_q <- 1
    twenty_q <- which(used_cumsum >= .2)[1]
    forty_q <- which(used_cumsum >=  .4)[1]
    sixty_q <- which(used_cumsum >= .6)[1]
    eighty_q<- which(used_cumsum >= .8)[1]
    oneh_q <- sum(used_trials)
    
    
    thirst_f <- m_obj$trials_mat[zero_q:twenty_q,15:20]
    second_f <- m_obj$trials_mat[(twenty_q + 1):forty_q,15:20]
    third_f <-  m_obj$trials_mat[(forty_q + 1):sixty_q,15:20]
    fourth_f <- m_obj$trials_mat[(sixty_q + 1):eighty_q,15:20]
    fifth_f <- m_obj$trials_mat[(eighty_q + 1):oneh_q,15:20]
    
    
    tabled_thirst_f <- rep(0, times=8)
    names(tabled_thirst_f) <- c(-1,1:7)
    tabled_thirst_f[names(table(thirst_f))] <- table(thirst_f)
    tabled_second_f <- rep(0, times=8)
    names(tabled_second_f) <- c(-1,1:7)
    tabled_second_f[names(table(second_f))] <- table(second_f)
    tabled_third_f <-  rep(0, times=8)
    names(tabled_third_f) <- c(-1,1:7)
    tabled_third_f[names(table(third_f))] <- table(third_f)
    tabled_fourth_f <- rep(0, times=8)
    names(tabled_fourth_f) <- c(-1,1:7)
    tabled_fourth_f[names(table(fourth_f))] <- table(fourth_f)
    tabled_fifth_f <-  rep(0, times=8)
    names(tabled_fifth_f) <- c(-1,1:7)
    tabled_fifth_f[names(table(fifth_f))] <- table(fifth_f)
    
    tmp <- rbind(tabled_thirst_f,
                 tabled_second_f,
                 tabled_third_f,
                 tabled_fourth_f,
                 tabled_fifth_f)
    
    
    reward_thirst_f <- m_obj$trials_mat[zero_q:twenty_q, 8:13]
    reward_second_f <- m_obj$trials_mat[(twenty_q + 1):forty_q, 8:13]
    reward_third_f <-  m_obj$trials_mat[(forty_q + 1):sixty_q, 8:13]
    reward_fourth_f <- m_obj$trials_mat[(sixty_q + 1):eighty_q, 8:13]
    reward_fifth_f <- m_obj$trials_mat[(eighty_q + 1):oneh_q, 8:13]
    
    
    reward_tabled_thirst_f <- rep(0, times=8)
    names(reward_tabled_thirst_f) <- c(-1,1:7)
    reward_tabled_thirst_f[names(table(reward_thirst_f))] <- table(reward_thirst_f)
    
    reward_tabled_second_f <- rep(0, times=8)
    names(reward_tabled_second_f) <- c(-1,1:7)
    reward_tabled_second_f[names(table(reward_second_f))] <- table(reward_second_f)
    
    reward_tabled_third_f <-  rep(0, times=8)
    names(reward_tabled_third_f) <- c(-1,1:7)
    reward_tabled_third_f[names(table(reward_third_f))] <- table(reward_third_f)
    
    reward_tabled_fourth_f <- rep(0, times=8)
    names(reward_tabled_fourth_f) <- c(-1,1:7)
    reward_tabled_fourth_f[names(table(reward_fourth_f))] <- table(reward_fourth_f)
    
    reward_tabled_fifth_f <-  rep(0, times=8)
    names(reward_tabled_fifth_f) <- c(-1,1:7)
    reward_tabled_fifth_f[names(table(reward_fifth_f))] <- table(reward_fifth_f)
    
    reward_tmp <- rbind(reward_tabled_thirst_f,
                        reward_tabled_second_f,
                        reward_tabled_third_f,
                        reward_tabled_fourth_f,
                        reward_tabled_fifth_f)
    
    
    tabled_reward <- table(m_obj$trials_mat[m_obj$annot_df[,1] == 3 & m_obj$annot_df[,2] == 0,8:13])
    reward_cl1 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[1]
    reward_cl2 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[2]
    cl1 <- names(sort(tabled_fifth_f, decreasing=T)[names(sort(tabled_fifth_f, decreasing=T)) != -1])[1]
    cl2 <- names(sort(tabled_fifth_f, decreasing=T)[names(sort(tabled_fifth_f, decreasing=T)) != -1])[2]
    
    contingency_list <- list(tabled_thirst_f,
                             tabled_second_f,
                             tabled_third_f,
                             tabled_fourth_f,
                             tabled_fifth_f)
    
    reward_contingency_list <- list(reward_tabled_thirst_f,
                                    reward_tabled_second_f,
                                    reward_tabled_third_f,
                                    reward_tabled_fourth_f,
                                    reward_tabled_fifth_f)
    
    cl1_vec <- c()
    cl2_vec <- c()
    
    for (tabled in contingency_list) {
      
      cl1_vec <- c(cl1_vec, tabled[cl1])
      cl2_vec <- c(cl2_vec, tabled[cl2])
      
    }
    
    reward_cl1_vec <- c()
    reward_cl2_vec <- c()
    
    for (tabled in reward_contingency_list) {
      
      reward_cl1_vec <- c(reward_cl1_vec, tabled[reward_cl1])
      reward_cl2_vec <- c(reward_cl2_vec, tabled[reward_cl2])
      
    }
    
    
    cl1_all <- rbind(cl1_all, cl1_vec)
    cl2_all <- rbind(cl2_all, cl2_vec)
    
    reward_cl1_all <- rbind(reward_cl1_all, reward_cl1_vec)
    reward_cl2_all <- rbind(reward_cl2_all, reward_cl2_vec)
  }
  
  cl1_all <- t(apply(cl1_all, 1, function(r) {r/sum(r)}))
  cl2_all <- t(apply(cl2_all, 1, function(r) {r/sum(r)}))
  reward_cl1_all <- t(apply(reward_cl1_all, 1, function(r) {r/sum(r)}))
  reward_cl2_all <- t(apply(reward_cl2_all, 1, function(r) {r/sum(r)}))
  
  colnames(cl1_all) <- seq(0.2,1,by=.2)
  colnames(cl2_all) <- seq(0.2,1,by=.2)
  colnames(reward_cl1_all) <- seq(0.2,1,by=.2)
  colnames(reward_cl2_all) <- seq(0.2,1,by=.2)
  
  mdf_cl1 <- as.data.frame(melt(cl1_all))
  mdf_cl2 <- as.data.frame(melt(cl2_all))
  mdf_reward_cl1 <- as.data.frame(melt(reward_cl1_all))
  mdf_reward_cl2 <- as.data.frame(melt(reward_cl2_all))
  
  colnames(mdf_cl1) <- c("#", "Thirst", "Fraction")
  colnames(mdf_cl2) <- c("#", "Thirst", "Fraction")
  colnames(mdf_reward_cl1) <- c("#", "Thirst", "Fraction")
  colnames(mdf_reward_cl2) <- c("#", "Thirst", "Fraction")
  
  
  mdf_cl1$Group <- c("ITI cluster 1")
  mdf_cl2$Group <- c("ITI cluster 2")
  mdf_reward_cl1$Group <- c("Reward cluster 1")
  mdf_reward_cl2$Group <- c("Reward cluster 2")
  
  df_final <- rbind(mdf_cl1,
                    mdf_cl2,
                    mdf_reward_cl1,
                    mdf_reward_cl2)
  
  
  
  gclust_legend <- 
    ggplot(df_final) + 
    geom_vline(xintercept=.2, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.4, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.6, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.8, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=1, linetype="dashed", size=.5, col="gray50") +
    geom_line(aes(x=Thirst, y=Fraction, color=Group), stat="summary") +
    geom_ribbon(aes(x=Thirst, y=Fraction, fill=Group), stat="summary", alpha=.5, color=NA) + 
    big_text_base_plot_theme_wl +
    theme(legend.position="top") + 
    ylab("Cluster occupancy")
  
  gclust <- 
    ggplot(df_final) + 
    geom_vline(xintercept=.2, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.4, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.6, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=.8, linetype="dashed", size=.5, col="gray50") +
    geom_vline(xintercept=1, linetype="dashed", size=.5, col="gray50") +
    geom_line(aes(x=Thirst, y=Fraction, color=Group), stat="summary") +
    geom_ribbon(aes(x=Thirst, y=Fraction, fill=Group), stat="summary", alpha=.5, color=NA) + 
    big_text_base_plot_theme +
    ylab("Cluster occupancy")
  
  
  colnames(cmsm_all) <- seq(0,1, length.out=200)
  cmsm_df <- as.data.frame(melt(cmsm_all))
  colnames(cmsm_df) <- c("#", "Duration", "Thirst")
  
  gbehav <- 
    ggplot(cmsm_df) + 
    geom_hline(yintercept=.2, linetype="dashed", size=.5, col="gray50") +
    geom_hline(yintercept=.4, linetype="dashed", size=.5, col="gray50") +
    geom_hline(yintercept=.6, linetype="dashed", size=.5, col="gray50") +
    geom_hline(yintercept=.8, linetype="dashed", size=.5, col="gray50") +
    geom_hline(yintercept=1, linetype="dashed", size=.5, col="gray50") +
    geom_line(aes(x=Duration, y=Thirst), stat="summary") +
    geom_ribbon(aes(x=Duration, y=Thirst), stat="summary", alpha=.5, color=NA) + 
    big_text_base_plot_theme_wl +
    theme(legend.position="top") + 
    xlab("Experimental duration")
  
  gf <- plot_grid(gbehav, gclust)
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    
    pdf(sprintf("%s\\%s_cluster_occupancy_exp_progression.pdf",base_figure_path, size_name), height=size, width=2 * size)
    plot(gf)
    dev.off()

    pdf(sprintf("%s\\%s_cluster_occupancy_exp_progression_legend.pdf",base_figure_path, size_name), height=size, width=size)
    plot(gclust_legend)
    dev.off()
    
  }
}
