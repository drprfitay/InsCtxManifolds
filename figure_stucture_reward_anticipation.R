
figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\figures_structure_reward_anticipation\\", base_output_path) 

dir.create(base_output_path)
dir.create(output_path)

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_HS <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
metadata_satiation <- across_mice_decoding_build_metadata(get_satiation_paths())

metadata_HS_interval_thresholds_sd <- c(.35, .35, .35, .35, .1, .35, .35, .35, .35, .35)
metadata_TT_interval_thresholds_sd <- rep(.1, 14)# c(.35, .35, .35, .35, .35, .35, .35, .35, .35, .35, 35, .35, .35, .35)

satiation_run_indices_all = c(5, 4, 3, 3, 3)
satiation_run_indices_interval <- c(30, 16, 30, 16,30)
satiation_paths <- get_satiation_paths()

get_task_lick_distribution <- function(metadata, path, false_alarms=F, CR=F, verbose=F)
{
  
  
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
  unaligned_task_licking_mat <- c()
  unaligned_task_cluster_mat <- c()
  tt = c()
  idxvc <- c()
  bvec <- c()
  ubvec <- c()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    licking_indices_of_interest <- trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))
    
    if (max(licking_indices_of_interest) > len(lick_vectors_all[[as.character(run)]])) {
      if (verbose) {
        print("Exceeding max!")
      }
      next
    }
    
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][licking_indices_of_interest]
    
    
    if( sum(licks_in_trial_unbinned > 0) <= 0 ) {
      if (verbose) {
        print(sprintf("%d. No licks in trial continue", reward_trial_idx))
      }
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
      if (verbose) {
        print(sprintf("%d. Only anticipatory licks?", reward_trial_idx))
      }
      next
    }
    
    unbinned_first_lick_in_respone <- which(licks_in_trial_unbinned > 0)[which(which(licks_in_trial_unbinned > 0) > 60)[1]]
    binned_first_lick_in_response <- get_binned_index(unbinned_first_lick_in_respone, 15)
    
    if (binned_first_lick_in_response + 3 > 20) {
      binned_first_lick_in_response = 20 - 3
    }
    
    task_lick_bouts_mat <- rbind(task_lick_bouts_mat,
                                 cluster_labels_of_trial[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 12)])
    
    task_licking_mat <- rbind(task_licking_mat,
                              as.numeric(licks_in_trial_unbinned[(unbinned_first_lick_in_respone - 60):(unbinned_first_lick_in_respone + 180)]>0))
    
    
    unaligned_task_licking_mat <- rbind(unaligned_task_licking_mat,
                                        as.numeric(licks_in_trial_unbinned > 0))
    
    unaligned_task_cluster_mat <- rbind(unaligned_task_cluster_mat,
                                        cluster_labels_of_trial)
    
    tt <- c(tt, relevant_trial["TrialType"])
    
    idxvc <- c(idxvc, reward_trial_idx)
    bvec <- c(bvec, binned_first_lick_in_response)
    ubvec <- c(ubvec, unbinned_first_lick_in_respone)
    
  }
  
  return(list(task_licking_mat=task_licking_mat,
              task_cluster_mat=task_lick_bouts_mat,
              unaligned_task_licking_mat=unaligned_task_licking_mat,
              unaligned_task_cluster_mat=unaligned_task_cluster_mat,
              tt=tt,
              trrial_indices=idxvc,
              binned=bvec,
              unbinned=ubvec))
}

plot_satiation_free_vs_dynamics_examples <- function()
{
  base_figure_path <- sprintf("%s\\new_bout_satiation_free_vs_dynamics_examples", output_path)
  dir.create(base_figure_path)
  
  
  pooled_task_similarity <- list()
  pooled_task_freely_similarity <- list()
  pooled_freely_freely_similarity <- list()
  
  datasets_names = paste("SATIATION_", get_datasets_names(satiation_paths, sep="_"), sep="")
  
  for (i in 1:len(metadata_satiation)) {
    res = get_task_vs_freely_lick_distribution(metadata_satiation[[i]], 
                                               satiation_paths[i], 
                                               satiation_run_indices_all[i], 
                                               bouts_interval = satiation_run_indices_interval[i])
    
    
    
    cmt <- rbind(res$task_cluster_mat, res$freely_cluster_mat)
    lmt <- rbind(res$task_licking_mat, res$freely_licking_mat)
    structured_unstructured <- data.frame(type=rep(c("Structured", "Unstructed"), times=c(nrow(res$task_cluster_mat),
                                                                                          nrow(res$freely_cluster_mat))))
    
    if (nrow(res$task_cluster_mat) > nrow(res$freely_cluster_mat)) {
      equal_cmt <- rbind(res$task_cluster_mat[1:nrow(res$freely_cluster_mat),], res$freely_cluster_mat)
      equal_lmt <- rbind(res$task_licking_mat[1:nrow(res$freely_licking_mat),], res$freely_licking_mat)
      bottom_equal_cmt <- rbind(res$task_cluster_mat[(nrow(res$task_cluster_mat) - nrow(res$freely_cluster_mat) + 1):nrow(res$task_cluster_mat),], res$freely_cluster_mat)
      bottom_equal_lmt <- rbind(res$task_licking_mat[(nrow(res$task_licking_mat) - nrow(res$freely_licking_mat) + 1):nrow(res$task_licking_mat),], res$freely_licking_mat)
      equal_structured_unstructured <- data.frame(type=rep(c("Structured", "Unstructed"), times=c(nrow(res$freely_cluster_mat),
                                                                                            nrow(res$freely_cluster_mat))))      
    } else {
      equal_cmt <- rbind(res$task_cluster_mat, res$freely_cluster_mat[1:nrow(res$task_cluster_mat),])
      equal_lmt <- rbind(res$task_licking_mat, res$freely_licking_mat[1:nrow(res$task_licking_mat),])
      bottom_equal_cmt <- rbind(res$task_cluster_mat, res$freely_cluster_mat[(nrow(res$freely_licking_mat) - nrow(res$task_cluster_mat) + 1):nrow(res$freely_licking_mat),])
      bottom_equal_lmt <- rbind(res$task_licking_mat, res$freely_licking_mat[(nrow(res$freely_licking_mat) - nrow(res$task_licking_mat) + 1):nrow(res$freely_licking_mat),])
      equal_structured_unstructured <- data.frame(type=rep(c("Structured", "Unstructed"), times=c(nrow(res$task_cluster_mat),
                                                                                            nrow(res$task_cluster_mat))))
    }

    
    

    
    rownames(structured_unstructured) <- 1:nrow(cmt)
    rownames(cmt) <- 1:nrow(cmt)
    colnames(cmt) <- rep("", times=ncol(cmt))
    colnames(lmt) <- rep("", times=ncol(lmt))
    colnames(cmt)[4] <- "O"
    colnames(lmt)[61] <- "O"
    
    rownames(equal_structured_unstructured) <- 1:nrow(equal_cmt)
    rownames(equal_cmt) <- 1:nrow(equal_cmt)
    colnames(equal_cmt) <- rep("", times=ncol(equal_cmt))
    colnames(equal_lmt) <- rep("", times=ncol(equal_lmt))
    colnames(equal_cmt)[4] <- "O"
    colnames(equal_lmt)[61] <- "O"
    
    rownames(equal_structured_unstructured) <- 1:nrow(bottom_equal_cmt)
    rownames(bottom_equal_cmt) <- 1:nrow(bottom_equal_cmt)
    colnames(bottom_equal_cmt) <- rep("", times=ncol(bottom_equal_cmt))
    colnames(bottom_equal_lmt) <- rep("", times=ncol(bottom_equal_lmt))
    colnames(bottom_equal_cmt)[4] <- "O"
    colnames(bottom_equal_lmt)[61] <- "O"
    
    
    cluster_color_label <- spec_cg(len(unique(metadata_satiation[[i]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata_satiation[[i]]$cluster_mat$labs)) - 1))
    
    phlmt <- pheatmap(lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    phcmt <- pheatmap(cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_row = structured_unstructured, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                      col=cluster_color_label)
    
    pheqlmt <- pheatmap(equal_lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    pheqcmt <- pheatmap(equal_cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_row = equal_structured_unstructured, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                      col=cluster_color_label )
    
    ph_bottom_eqlmt <- pheatmap(bottom_equal_lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    ph_bottom_eqcmt <- pheatmap(bottom_equal_cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_row = equal_structured_unstructured, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                        col=cluster_color_label )
    
    plt <- plot_grid(phcmt[[4]], phlmt[[4]])
    plt_eq <- plot_grid(pheqcmt[[4]], pheqlmt[[4]])
    plt_bottom_eq <- plot_grid(ph_bottom_eqcmt[[4]], ph_bottom_eqlmt[[4]])
    
    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_unequal_mat.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[i]),
          height=size,
          width=size)
      
      plot(plt)
      
      dev.off()
      
      pdf(sprintf("%s\\%s_%s_equal_mat.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[i]),
          height=size,
          width=size)
      
      plot(plt_eq)
      
      dev.off()

      pdf(sprintf("%s\\%s_%s_bottom_equal_mat.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[i]),
          height=size,
          width=size)
      
      plot(plt_bottom_eq)
      
      dev.off()
      
      
      
    }
    
    
    sim_mat <- cosine(t(cmt))
    
    
    task_sim_mat <- sim_mat[1:nrow(res$task_cluster_mat),
                            1:nrow(res$task_cluster_mat)]
    
    free_sim_mat <- sim_mat[(nrow(res$task_cluster_mat) + 1):nrow(cmt), 1:nrow(res$task_cluster_mat)]
    
    free_free_sim_mat <- sim_mat[(nrow(res$task_cluster_mat) + 1):nrow(cmt), (nrow(res$task_cluster_mat) + 1):nrow(cmt)]
    
    pooled_task_similarity <- append(pooled_task_similarity,
                                     list(task_sim_mat[lower.tri(task_sim_mat)]))
    
    pooled_task_freely_similarity <- append(pooled_task_freely_similarity,
                                            list(c(free_sim_mat)))
    
    pooled_freely_freely_similarity <- append(pooled_freely_freely_similarity, list(free_free_sim_mat[lower.tri(free_free_sim_mat)]))
  }
  
  
  pooled_smilarity_list <- list(task=unlist(pooled_task_similarity),
                                freely=unlist(pooled_task_freely_similarity),
                                free_free=unlist(pooled_freely_freely_similarity))

  N_breaks = 20
  hist_df <- 
    lapply(names(pooled_smilarity_list),
           function(lst_name) {
             lst <- pooled_smilarity_list[[lst_name]]
             
             h <- hist(lst, plot=F, seq(min(unlist(pooled_smilarity_list), na.rm=T),
                                        max(unlist(pooled_smilarity_list), na.rm=T),
                                        length.out=(N_breaks + 1)))
             
             df <- 
               data.frame(breaks=h$breaks,
                          freq=c(0, h$counts/sum(h$counts)),
                          group=rep(lst_name, times=(N_breaks + 1)))
             
             return(df)
           })
  
  all_df <- do.call(rbind,hist_df)
  
  
  gsim <- 
    ggplot(all_df, aes(x=breaks, y=freq)) + 
    #geom_bar(aes(fill=group), color=adjustcolor("gray30", alpha=.85), size=.01, alpha=.2, width=.01,  position=position_nudge(), stat="summary") + 
    geom_line(aes(color=group),  size=.85, alpha=.85, stat="summary") + 
    #geom_ribbon(aes(fill=group),  color=NA, alpha=.25, stat="summary") +
    big_text_base_plot_theme + 
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Cosine similarity")  +
    ylab("Fraction")
    
    
    tmp_df <- 
      cbind(unlist(lapply(pooled_task_similarity, mean, na.rm=T)), 
            unlist(lapply(pooled_task_freely_similarity, mean, na.rm=T)))
    colnames(tmp_df) <- c("Task vs Task", "Task vs Freely")
    
    melted <- melt(tmp_df)
    colnames(melted) <- c("#", "Group", "Sim")
    
    
    ## Statistics
    
    wilc <- 
      t.test(tmp_df[,1], tmp_df[,2], correct = F, paired=T, alternative="greater")
    
    
    stat_df <-  
      data.frame(Statistics=wilc$statistic,
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
      ylab("Cosine similarity")
    
    combined_plot <- 
      plot_grid(gsim, gboxf, rel_widths=c(2,1), nrow=1)
    
    
    for (size_name in names(a4_sizes)) {
      dir.create(sprintf("%s\\dist_and_box\\", base_figure_path))
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\dist_and_box\\%s_free_free_trial_similarity.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gsim)
      dev.off()
      
      pdf(sprintf("%s\\dist_and_box\\%s_trial_similarity_boxplot_comp.pdf",base_figure_path, size_name), height=size, width=size)
      plot(gboxf)
      dev.off()
      
      pdf(sprintf("%s\\dist_and_box\\%s_combined_trial_similarity_boxplot_comp.pdf",base_figure_path, size_name), height=size/1.5, width=size)
      plot(combined_plot)
      dev.off()
    }
    
    
    write.csv(file=sprintf("%s//dist_and_box//ttest_statistics.csv",base_figure_path), stat_df)
}

plot_HS_lick_anticipation_examples <- function()
{
  base_figure_path <- sprintf("%s\\new_alignment_HS_lick_anticipation_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_hungry_sated_paths()
  metadata <- metadata_HS
  datasets_names = get_datasets_names(paths, sep="_")
  
  for (i in 1:len(metadata)) {
      res <- get_task_lick_distribution(metadata[[i]], paths[i])
      
      cmt <- res$task_cluster_mat
      lmt <- res$task_licking_mat
      
      unaligned_lmt <- res$unaligned_task_licking_mat
      unaligned_cmt <- res$unaligned_task_cluster_mat
      
      colnames(cmt) <- c()
      
      cluster_color_label <- spec_cg(len(unique(metadata[[i]]$cluster_mat$labs)))
      names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[i]]$cluster_mat$labs)) - 1))
      
      colnames(cmt) <- rep("", times=ncol(cmt))
   #   colnames(cmt)[4] <- "Onset"
      colnames(lmt) <- rep("", times=ncol(lmt))
  #    colnames(lmt)[61] <- "Onset"
      
      colnames(unaligned_cmt) <- rep("", times=ncol(unaligned_cmt))
      #      colnames(unaligned_cmt)[4] <- "Onset"
      colnames(unaligned_lmt) <- rep("", times=ncol(unaligned_lmt))
      #     colnames(unaligned_lmt)[61] <- "Onset"
      
      const=metadata_HS_interval_thresholds_sd[i]; 
      min_frame_interval <- ceiling(1/(mean(unaligned_lmt[1:15,]) + const*sd(unaligned_lmt[1:15,])))
      
      unaligned_licks_per_trial <- lapply(1:nrow(unaligned_lmt), function(ri) {which(unaligned_lmt[ri,] == 1)})
      unaligned_licks_interval_per_trial <- lapply(unaligned_licks_per_trial, diff)
      unaligned_licks_interval_per_trial_under_mnimimum <- lapply(unaligned_licks_interval_per_trial, 
                                                                  function(intervals) {as.numeric(intervals <= min_frame_interval)})
      
      first_lick_index_per_trial <- lapply(unaligned_licks_interval_per_trial_under_mnimimum,
                                            function(under_min_intevals) {which(rollsum(under_min_intevals, 2) == 2)[1]})
      
      first_lick_index_per_trial <- unlist(first_lick_index_per_trial)
      first_lick_per_trial <- lapply(1:nrow(unaligned_lmt), 
                                     function(ri) {
                                       
                                       fl_idx = first_lick_index_per_trial[ri]
                                       
                                         if (!is.na(fl_idx)) {
                                          return(which(unaligned_lmt[ri,] == 1)[fl_idx])
                                         } else {
                                           return(which(unaligned_lmt[ri,] == 1)[which(unaligned_lmt[ri,] == 1) > 60][1])
                                         }
                                       
                                       })
      
      first_lick_per_trial <- unlist(first_lick_per_trial)
      
      

      #lick_onset_order <- order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, sum), decreasing=T)
      # lick_onset_order <-  order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)}), decreasing = F)
      lick_onset_order <- order(first_lick_per_trial, decreasing=F)
      
      phlmt <- pheatmap(lmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
      phcmt <- pheatmap(cmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                        col=cluster_color_label )
      
      
      
      phlmt_unsorted <- pheatmap(lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
      phcmt_unsorted <- pheatmap(cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                        col=cluster_color_label)
      
      
      phlmt_unsorted_all <- pheatmap(unaligned_lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black", `2`="black"), border_col=NA, show_rownames=F,)
      phcmt_unsorted_all <- pheatmap(unaligned_cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                                 col=cluster_color_label)
      
      
      phlmt_sorted_all <- pheatmap(unaligned_lmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black", `2`="black"), border_col=NA, show_rownames=F,)
      phcmt_sorted_all <- pheatmap(unaligned_cmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                                   col=cluster_color_label)
      
      plt <- 
        plot_grid(phcmt[[4]], phlmt[[4]])
      
      plt_unsorted <- 
        plot_grid(phcmt_unsorted[[4]], phlmt_unsorted[[4]])
      
      plt_unsorted_all <- 
        plot_grid(phcmt_unsorted_all[[4]], phlmt_unsorted_all[[4]])
      
      plt_sorted_all <- 
        plot_grid(phcmt_sorted_all[[4]], phlmt_sorted_all[[4]])
      
    
    for (size_name in names(a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[i]]))
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_sorted.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt)
      
      dev.off()
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_unsorted.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_unsorted)
      
      dev.off()
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_unsorted_unaligned.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_unsorted_all)
      
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s_anticipation_sorted_unaligned.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_sorted_all)
      
      dev.off()
    }
  }
}

plot_TQ_lick_anticipation_examples <- function()
{
  base_figure_path <- sprintf("%s\\TQ_lick_anticipation_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  metadata <- metadata_TT
  datasets_names = get_datasets_names(paths, sep="_")
  
  for (i in 1:len(metadata)) {
    res <- get_task_lick_distribution(metadata[[i]], paths[i])
    
    cmt <- res$task_cluster_mat
    lmt <- res$task_licking_mat
    
    unaligned_lmt <- res$unaligned_task_licking_mat
    unaligned_cmt <- res$unaligned_task_cluster_mat
    
    colnames(cmt) <- c()
    
    cluster_color_label <- spec_cg(len(unique(metadata[[i]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[i]]$cluster_mat$labs)) - 1))
    
    colnames(cmt) <- rep("", times=ncol(cmt))
    #   colnames(cmt)[4] <- "Onset"
    colnames(lmt) <- rep("", times=ncol(lmt))
    #    colnames(lmt)[61] <- "Onset"
    
    
    colnames(unaligned_cmt) <- rep("", times=ncol(unaligned_cmt))
    #      colnames(unaligned_cmt)[4] <- "Onset"
    colnames(unaligned_lmt) <- rep("", times=ncol(unaligned_lmt))
    #     colnames(unaligned_lmt)[61] <- "Onset"
    
    const=.35; 
    min_frame_interval <- ceiling(1/(mean(unaligned_lmt[1:15,]) + const*sd(unaligned_lmt[1:15,])))
    
    unaligned_licks_per_trial <- lapply(1:nrow(unaligned_lmt), function(ri) {which(unaligned_lmt[ri,] == 1)})
    unaligned_licks_interval_per_trial <- lapply(unaligned_licks_per_trial, diff)
    unaligned_licks_interval_per_trial_under_mnimimum <- lapply(unaligned_licks_interval_per_trial, 
                                                                function(intervals) {as.numeric(intervals <= min_frame_interval)})
    
    first_lick_index_per_trial <- lapply(unaligned_licks_interval_per_trial_under_mnimimum,
                                         function(under_min_intevals) {which(rollsum(under_min_intevals, 2) == 2)[1]})
    
    first_lick_index_per_trial <- unlist(first_lick_index_per_trial)
    first_lick_per_trial <- lapply(1:nrow(unaligned_lmt), 
                                   function(ri) {
                                     
                                     fl_idx = first_lick_index_per_trial[ri]
                                     
                                     if (!is.na(fl_idx)) {
                                       return(which(unaligned_lmt[ri,] == 1)[fl_idx])
                                     } else {
                                       return(which(unaligned_lmt[ri,] == 1)[which(unaligned_lmt[ri,] == 1) > 60][1])
                                     }
                                     
                                   })
    
    first_lick_per_trial <- unlist(first_lick_per_trial)

    #lick_onset_order <- order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, sum), decreasing=T)
    # lick_onset_order <-  order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)}), decreasing = F)
    lick_onset_order <- order(first_lick_per_trial, decreasing=F)
    
    
    
    #lick_onset_order <- order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, sum), decreasing=T)
    #lick_onset_order <-  order(apply(t(apply(lmt[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)}), decreasing = F)
    
    phlmt <- pheatmap(lmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    phcmt <- pheatmap(cmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                      col=cluster_color_label )
    
    
    
    phlmt_unsorted <- pheatmap(lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    phcmt_unsorted <- pheatmap(cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                               col=cluster_color_label)
    
    
    phlmt_unsorted_all <- pheatmap(unaligned_lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black", `2`="black"), border_col=NA, show_rownames=F,)
    phcmt_unsorted_all <- pheatmap(unaligned_cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                                   col=cluster_color_label)
    
    
    phlmt_sorted_all <- pheatmap(unaligned_lmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black", `2`="black"), border_col=NA, show_rownames=F,)
    phcmt_sorted_all <- pheatmap(unaligned_cmt[lick_onset_order,], cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                                   col=cluster_color_label)
    
    
    plt <- 
      plot_grid(phcmt[[4]], phlmt[[4]])
    
    plt_unsorted <- 
      plot_grid(phcmt_unsorted[[4]], phlmt_unsorted[[4]])
    
    plt_unsorted_all <- 
      plot_grid(phcmt_unsorted_all[[4]], phlmt_unsorted_all[[4]])
    
    plt_sorted_all <- 
      plot_grid(phcmt_sorted_all[[4]], phlmt_sorted_all[[4]])
    
    for (size_name in names(a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[i]]))
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_sorted.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt)
      
      dev.off()
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_unsorted.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_unsorted)
      
      dev.off()
      
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_unsorted_unaligned.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_unsorted_all)
      
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s_anticipation_sorted_unaligned.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(plt_sorted_all)
      
      dev.off()
    }
  }
}

plot_TQ_lick_onset_cluster_onset_anticipation <- function()
{
  base_figure_path <- sprintf("%s\\TQ_lick_onset_cluster_onset_anticipation", output_path)
  dir.create(base_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  metadata <- metadata_TT
  datasets_names = get_datasets_names(paths, sep="_")
  
  all_res <- list()
  
  all_cor <- c()
  all_scor <- c()
  all_tcor <- c()
  
  all_licks <- c()
  all_licks_2nd <- c()
  all_licks_3rd <- c()
  
  all_onsets <- c()
  all_onsets_2nd <- c()
  all_onsets_3rd <- c()
  
  for (i in 1:len(metadata)) {
    res <- get_task_lick_distribution(metadata[[i]], paths[i])
    relevant_trials <- metadata[[i]]$stim_master_mat[metadata[[i]]$stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    reward_trials_to_use <- which(metadata[[i]]$annot_df[,1] == 3 & metadata[[i]]$annot_df[,2] == 0)[res$trrial_indices]
    
    mat <- metadata[[i]]$trials_mat[reward_trials_to_use,]
    
    
    cmt <- res$task_cluster_mat
    lmt <- res$task_licking_mat
    
    unaligned_lmt <- res$unaligned_task_licking_mat
    unaligned_cmt <- res$unaligned_task_cluster_mat
    
    tbld <- table(res$task_cluster_mat[,5:8])
    tbld <- tbld[names(tbld) != -1]
    first_clust <- as.numeric(names(which.max(tbld)))
    second_clust <- as.numeric(names(sort(tbld, decreasing = T)[2]))
    third_clust <- as.numeric(names(sort(tbld, decreasing = T)[3]))
    
    
    cluster_onset <- apply(mat, 1, function(trial) {which(trial == first_clust)[1]})
    second_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == second_clust)[1]})
    thirst_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == third_clust)[1]})
    
    # lcks <- apply(t(apply(res$task_licking_mat[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)})
    # lcks <- res$unbinned - 60 + lcks
    
    const=metadata_TT_interval_thresholds_sd[i];
    min_frame_interval <- ceiling(1/(cmean(unaligned_lmt[1:15,]) + const*sd(unaligned_lmt[1:15,])))

    unaligned_licks_per_trial <- lapply(1:nrow(unaligned_lmt), function(ri) {which(unaligned_lmt[ri,] == 1)})
    unaligned_licks_interval_per_trial <- lapply(unaligned_licks_per_trial, diff)
    unaligned_licks_interval_per_trial_under_mnimimum <- lapply(unaligned_licks_interval_per_trial,
                                                                function(intervals) {as.numeric(intervals <= min_frame_interval)})

    first_lick_index_per_trial <- lapply(unaligned_licks_interval_per_trial_under_mnimimum,
                                         function(under_min_intevals) {which(rollsum(under_min_intevals, 2) == 2)[1]})

    first_lick_index_per_trial <- unlist(first_lick_index_per_trial)
    first_lick_per_trial <- lapply(1:nrow(unaligned_lmt),
                                   function(ri) {

                                     fl_idx = first_lick_index_per_trial[ri]

                                     if (!is.na(fl_idx)) {
                                       return(which(unaligned_lmt[ri,] == 1)[fl_idx])
                                     } else {
                                       return(which(unaligned_lmt[ri,] == 1)[which(unaligned_lmt[ri,] == 1) > 60][1])
                                     }

                                   })

    first_lick_per_trial <- unlist(first_lick_per_trial)




    lcks <- first_lick_per_trial

    cr <- cor(lcks[!is.na(cluster_onset)], 
              cluster_onset[!is.na(cluster_onset)])
    
    scr <- cor(lcks[!is.na(second_cluster_onset)], 
               second_cluster_onset[!is.na(second_cluster_onset)])
    
    tcr <- cor(lcks[!is.na(thirst_cluster_onset)], 
               thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    
    
    all_res <- append(all_res,
                      list(list(licks=lcks[!is.na(cluster_onset)],
                                onsets=cluster_onset[!is.na(cluster_onset)],
                                licks_2nd=lcks[!is.na(second_cluster_onset)],
                                onsets_2nd=second_cluster_onset[!is.na(second_cluster_onset)],
                                licks_3rd=lcks[!is.na(thirst_cluster_onset)],
                                onsets_3rd=thirst_cluster_onset[!is.na(thirst_cluster_onset)])))
    
    all_onsets <- c(all_onsets, cluster_onset[!is.na(cluster_onset)])
    all_onsets_2nd <- c(all_onsets_2nd, second_cluster_onset[!is.na(second_cluster_onset)])
    all_onsets_3rd <- c(all_onsets_3rd, thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    all_licks <- c(all_licks, lcks[!is.na(cluster_onset)])
    all_licks_2nd <- c(all_licks_2nd, lcks[!is.na(second_cluster_onset)])
    all_licks_3rd <- c(all_licks_3rd, lcks[!is.na(thirst_cluster_onset)])
    
    
    all_cor <- c(all_cor, cr)
    all_scor <- c(all_scor, scr)
    all_tcor <- c(all_tcor, tcr)
    
    print(sprintf("1st cluster %.3f, 2nd cluster %.3f, 3rd cluster %.3f", cr, scr, tcr))
  }
  
  vls <- apply(cbind(all_cor, all_scor), 1, which.max)
  combined_licks <- lapply(1:len(vls), function(idx) {if(vls[idx] == 1) {all_res[[idx]]$licks} else if (vls[idx] == 2) {all_res[[idx]]$licks_2nd}})# else {all_res[[idx]]$licks_3rd}})
  combined_onsets <- lapply(1:len(vls), function(idx) {if(vls[idx] == 1) {all_res[[idx]]$onsets} else if (vls[idx] == 2) {all_res[[idx]]$onsets_2nd}})# else {all_res[[idx]]$onsets_3rd}})
  
  
  all_p <- list()
  
  for (di in 1:len(combined_licks)) {
    as <- cbind(unlist(combined_licks[[di]]), unlist(combined_onsets[[di]]))
    
    asd <- data.frame(as)
    colnames(asd) <- c("Licks", "Clusters")
    
  
    lmodel <- lm(formula=Licks~Clusters, data=asd)
    sm <- summary(lmodel)
    print(sm$r.squared)
    gt <- 
    ggplot(asd, aes(x=Clusters / 2, y=Licks / 30)) +
      geom_point() + 
      geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
      ggtitle(sprintf("%.4f - %.4f", sqrt(sm$r.squared), sm$coefficients[2,4])) +
      big_text_base_plot_theme + 
      xlab("Cluster onset (0.5s bin)") + 
      ylab("Lick onset (Frame)") +
      theme(plot.title = element_text(size=9))
    
    for (size_name in names(small_a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[di]]))
      
      size = small_a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_correlation.pdf",
                  base_figure_path,
                  datasets_names[di],
                  size_name),
          height=size,
          width=size)
      
      plot(gt)
      
      dev.off()
    }
  }
}

plot_HS_lick_onset_cluster_onset_anticipation <- function()
{
  base_figure_path <- sprintf("%s\\HS_lick_onset_cluster_onset_anticipation", output_path)
  dir.create(base_figure_path)
  
  paths <- get_hungry_sated_paths()
  metadata <- metadata_HS
  datasets_names = get_datasets_names(paths, sep="_")
  
  all_res <- list()
  
  all_cor <- c()
  all_scor <- c()
  all_tcor <- c()
  
  all_licks <- c()
  all_licks_2nd <- c()
  all_licks_3rd <- c()
  
  all_onsets <- c()
  all_onsets_2nd <- c()
  all_onsets_3rd <- c()
  
  for (i in 1:len(metadata)) {
    res <- get_task_lick_distribution(metadata[[i]], paths[i])
    relevant_trials <- metadata[[i]]$stim_master_mat[metadata[[i]]$stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    reward_trials_to_use <- which(metadata[[i]]$annot_df[,1] == 3 & metadata[[i]]$annot_df[,2] == 0)[res$trrial_indices]
    
    mat <- metadata[[i]]$trials_mat[reward_trials_to_use,]
    
    
    cmt <- res$task_cluster_mat
    lmt <- res$task_licking_mat
    
    unaligned_lmt <- res$unaligned_task_licking_mat
    unaligned_cmt <- res$unaligned_task_cluster_mat
    
    tbld <- table(res$task_cluster_mat[,5:8])
    tbld <- tbld[names(tbld) != -1]
    first_clust <- as.numeric(names(which.max(tbld)))
    second_clust <- as.numeric(names(sort(tbld, decreasing = T)[2]))
    third_clust <- as.numeric(names(sort(tbld, decreasing = T)[3]))
    
    
    cluster_onset <- apply(mat, 1, function(trial) {which(trial == first_clust)[1]})
    second_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == second_clust)[1]})
    thirst_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == third_clust)[1]})
    
    # lcks <- apply(t(apply(res$task_licking_mat[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)})
    # lcks <- res$unbinned - 60 + lcks
    
    const=metadata_HS_interval_thresholds_sd[i];
    min_frame_interval <- ceiling(1/(mean(unaligned_lmt[1:15,]) + const*sd(unaligned_lmt[1:15,])))
    
    unaligned_licks_per_trial <- lapply(1:nrow(unaligned_lmt), function(ri) {which(unaligned_lmt[ri,] == 1)})
    unaligned_licks_interval_per_trial <- lapply(unaligned_licks_per_trial, diff)
    unaligned_licks_interval_per_trial_under_mnimimum <- lapply(unaligned_licks_interval_per_trial,
                                                                function(intervals) {as.numeric(intervals <= min_frame_interval)})
    
    first_lick_index_per_trial <- lapply(unaligned_licks_interval_per_trial_under_mnimimum,
                                         function(under_min_intevals) {which(rollsum(under_min_intevals, 2) == 2)[1]})
    
    first_lick_index_per_trial <- unlist(first_lick_index_per_trial)
    first_lick_per_trial <- lapply(1:nrow(unaligned_lmt),
                                   function(ri) {
                                     
                                     fl_idx = first_lick_index_per_trial[ri]
                                     
                                     if (!is.na(fl_idx)) {
                                       return(which(unaligned_lmt[ri,] == 1)[fl_idx])
                                     } else {
                                       return(which(unaligned_lmt[ri,] == 1)[which(unaligned_lmt[ri,] == 1) > 60][1])
                                     }
                                     
                                   })
    
    first_lick_per_trial <- unlist(first_lick_per_trial)
    
    
    
    
    lcks <- first_lick_per_trial
    
    cr <- cor(lcks[!is.na(cluster_onset)], 
              cluster_onset[!is.na(cluster_onset)])
    
    scr <- cor(lcks[!is.na(second_cluster_onset)], 
               second_cluster_onset[!is.na(second_cluster_onset)])
    
    tcr <- cor(lcks[!is.na(thirst_cluster_onset)], 
               thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    
    
    all_res <- append(all_res,
                      list(list(licks=lcks[!is.na(cluster_onset)],
                                onsets=cluster_onset[!is.na(cluster_onset)],
                                licks_2nd=lcks[!is.na(second_cluster_onset)],
                                onsets_2nd=second_cluster_onset[!is.na(second_cluster_onset)],
                                licks_3rd=lcks[!is.na(thirst_cluster_onset)],
                                onsets_3rd=thirst_cluster_onset[!is.na(thirst_cluster_onset)])))
    
    all_onsets <- c(all_onsets, cluster_onset[!is.na(cluster_onset)])
    all_onsets_2nd <- c(all_onsets_2nd, second_cluster_onset[!is.na(second_cluster_onset)])
    all_onsets_3rd <- c(all_onsets_3rd, thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    all_licks <- c(all_licks, lcks[!is.na(cluster_onset)])
    all_licks_2nd <- c(all_licks_2nd, lcks[!is.na(second_cluster_onset)])
    all_licks_3rd <- c(all_licks_3rd, lcks[!is.na(thirst_cluster_onset)])
    
    
    all_cor <- c(all_cor, cr)
    all_scor <- c(all_scor, scr)
    all_tcor <- c(all_tcor, tcr)
    
    print(sprintf("1st cluster %.3f, 2nd cluster %.3f, 3rd cluster %.3f", cr, scr, tcr))
  }
  
  vls <- apply(cbind(all_cor, all_scor), 1, which.max)
  combined_licks <- lapply(1:len(vls), function(idx) {if(idx==9 || vls[idx] == 1) {all_res[[idx]]$licks} else if (vls[idx] == 2) {all_res[[idx]]$licks_2nd}})# else {all_res[[idx]]$licks_3rd}})
  combined_onsets <- lapply(1:len(vls), function(idx) {if(idx==9 ||vls[idx] == 1) {all_res[[idx]]$onsets} else if (vls[idx] == 2) {all_res[[idx]]$onsets_2nd}})# else {all_res[[idx]]$onsets_3rd}})
  
  
  all_p <- list()
  
  for (di in 1:len(combined_licks)) {
    as <- cbind(unlist(combined_licks[[di]]), unlist(combined_onsets[[di]]))
    
    asd <- data.frame(as)
    colnames(asd) <- c("Licks", "Clusters")
    
    
    lmodel <- lm(formula=Licks~Clusters, data=asd)
    sm <- summary(lmodel)
    print(sm$r.squared)
    gt <- 
      ggplot(asd, aes(x=Clusters / 2, y=Licks / 30)) +
      geom_point() + 
      geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
      ggtitle(sprintf("%.4f - %.4f", sqrt(sm$r.squared), sm$coefficients[2,4])) +
      big_text_base_plot_theme + 
      xlab("Cluster onset (0.5s bin)") + 
      ylab("Lick onset (Frame)") +
      theme(plot.title = element_text(size=9))
    
    for (size_name in names(small_a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[di]]))
      
      size = small_a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_correlation.pdf",
                  base_figure_path,
                  datasets_names[di],
                  size_name),
          height=size,
          width=size)
      
      plot(gt)
      
      dev.off()
    }
  }
}

plot_anticipatory_statistics_all <- function()
{
  base_figure_path <- sprintf("%s\\anticipatory_statistics_all", output_path)
  dir.create(base_figure_path)
  
  paths <- c(get_hungry_sated_paths(), get_thirsty_quenched_paths())
  metadata <- append(metadata_HS, metadata_TT)
  datasets_names = get_datasets_names(paths, sep="_")
  
  all_res <- list()
  
  all_cor <- c()
  all_scor <- c()
  all_tcor <- c()
  
  all_licks <- c()
  all_licks_2nd <- c()
  all_licks_3rd <- c()
  
  all_onsets <- c()
  all_onsets_2nd <- c()
  all_onsets_3rd <- c()
  
  cor_statistics  <- data.frame()
  
  for (i in 1:len(metadata)) {
    res <- get_task_lick_distribution(metadata[[i]], paths[i])
    relevant_trials <- metadata[[i]]$stim_master_mat[metadata[[i]]$stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    reward_trials_to_use <- which(metadata[[i]]$annot_df[,1] == 3 & metadata[[i]]$annot_df[,2] == 0)[res$trrial_indices]
    
    mat <- metadata[[i]]$trials_mat[reward_trials_to_use,]
    
    reward_stim_mat <- metadata[[i]]$stim_master_mat[which(metadata[[i]]$stim_master_mat[,"TrialType"] %in% c(1,3,4,5))[reward_trials_to_use],]
    
    
    cmt <- res$task_cluster_mat
    lmt <- res$task_licking_mat
    
    unaligned_lmt <- res$unaligned_task_licking_mat
    unaligned_cmt <- res$unaligned_task_cluster_mat
    

    
    
    tbld <- table(res$task_cluster_mat[,5:8])
    tbld <- tbld[names(tbld) != -1]
    first_clust <- as.numeric(names(which.max(tbld)))
    second_clust <- as.numeric(names(sort(tbld, decreasing = T)[2]))
    third_clust <- as.numeric(names(sort(tbld, decreasing = T)[3]))
    
    
    cluster_onset <- apply(mat, 1, function(trial) {which(trial == first_clust)[1]})
    second_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == second_clust)[1]})
    thirst_cluster_onset <-  apply(mat, 1, function(trial) {which(trial == third_clust)[1]})
    

    # lcks <- apply(t(apply(res$task_licking_mat[,1:61], 1, function(mr) {rollmean(mr, 10)})), 1, function(r) {rs <- which(r > .1)[1]; ifelse(is.na(rs), 61, rs)})
    # lcks <- res$unbinned - 60 + lcks
    
    const=c(metadata_HS_interval_thresholds_sd,metadata_TT_interval_thresholds_sd)[i];
    min_frame_interval <- ceiling(1/(mean(unaligned_lmt[1:15,]) + const*sd(unaligned_lmt[1:15,])))
    
    unaligned_licks_per_trial <- lapply(1:nrow(unaligned_lmt), function(ri) {which(unaligned_lmt[ri,] == 1)})
    unaligned_licks_interval_per_trial <- lapply(unaligned_licks_per_trial, diff)
    unaligned_licks_interval_per_trial_under_mnimimum <- lapply(unaligned_licks_interval_per_trial,
                                                                function(intervals) {as.numeric(intervals <= min_frame_interval)})
    
    first_lick_index_per_trial <- lapply(unaligned_licks_interval_per_trial_under_mnimimum,
                                         function(under_min_intevals) {which(rollsum(under_min_intevals, 2) == 2)[1]})
    
    first_lick_index_per_trial <- unlist(first_lick_index_per_trial)
    first_lick_per_trial <- lapply(1:nrow(unaligned_lmt),
                                   function(ri) {
                                     
                                     fl_idx = first_lick_index_per_trial[ri]
                                     
                                     if (!is.na(fl_idx)) {
                                       return(which(unaligned_lmt[ri,] == 1)[fl_idx])
                                     } else {
                                       return(which(unaligned_lmt[ri,] == 1)[which(unaligned_lmt[ri,] == 1) > 60][1])
                                     }
                                     
                                   })
    
    first_lick_per_trial <- unlist(first_lick_per_trial)
    
    
    
    
    lcks <- first_lick_per_trial
    
    cr <- cor(lcks[!is.na(cluster_onset)], cluster_onset[!is.na(cluster_onset)])
    scr <- cor(lcks[!is.na(second_cluster_onset)], second_cluster_onset[!is.na(second_cluster_onset)])
    tcr <- cor(lcks[!is.na(thirst_cluster_onset)], thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    cr_test <- cor.test(lcks[!is.na(cluster_onset)], cluster_onset[!is.na(cluster_onset)])
    scr_test <- cor.test(lcks[!is.na(second_cluster_onset)], second_cluster_onset[!is.na(second_cluster_onset)])
    #tcr_test <- cor.test(lcks[!is.na(thirst_cluster_onset)], thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    tmp_unaligned_cmt <- unaligned_cmt
    
    cluster_color_label <- spec_cg(len(unique(metadata[[i]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[i]]$cluster_mat$labs)) - 1))
    
    for(c_idx in 1:len(second_cluster_onset)) {
      tmp_unaligned_cmt[c_idx,] <- rep(0,20);
      c=cluster_onset[c_idx];
      if(is.na(c)){next}; 
      tmp_unaligned_cmt[c_idx, c] <- 50; 
      tmp_unaligned_cmt[c_idx, res$binned[c_idx]] <- -50; 
      tmp_unaligned_cmt[c_idx, ceiling(lcks[c_idx] / 15)] <- -100;
    }
    ph1 <- ph(unaligned_cmt[order(first_lick_per_trial),], show_rownames=F, border_col=NA, legend=F, show_colnames=F, color=cluster_color_label)
    ph2 <- ph(tmp_unaligned_cmt[order(first_lick_per_trial),], show_rownames=F, border_col=NA, legend=F, col=c("50"="red","-50"="blue","-100"="yellow", "0"="green"), show_colnames=F)
    
    combined_plot <- plot_grid(ph1[[4]], ph2[[4]], nrow=1)
    
    for (size_name in names(small_a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[i]]))
      
      size = small_a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_lick_anticipation_onset.pdf",
                  base_figure_path,
                  datasets_names[i],
                  size_name),
          height=size,
          width=size)
      
      plot(combined_plot)
      
      dev.off()
    }
    
    
    all_res <- append(all_res,
                      list(list(licks=lcks[!is.na(cluster_onset)],
                                onsets=cluster_onset[!is.na(cluster_onset)],
                                licks_2nd=lcks[!is.na(second_cluster_onset)],
                                onsets_2nd=second_cluster_onset[!is.na(second_cluster_onset)],
                                licks_3rd=lcks[!is.na(thirst_cluster_onset)],
                                onsets_3rd=thirst_cluster_onset[!is.na(thirst_cluster_onset)],
                                reward_onsets=res$unbinned[!is.na(cluster_onset)],
                                reward_onsets_2nd=res$unbinned[!is.na(second_cluster_onset)],
                                reward_onsets_3rd=res$unbinned[!is.na(thirst_cluster_onset)],
                                cor_test=cr_test,
                                cor_test_2nd=scr_test)))
    
    all_reward_onset <- c(all_reward_onset, (reward_stim_mat[,"Reward"] - reward_stim_mat[,"Frames"])/15)
    all_onsets <- c(all_onsets, cluster_onset[!is.na(cluster_onset)])
    all_onsets_2nd <- c(all_onsets_2nd, second_cluster_onset[!is.na(second_cluster_onset)])
    all_onsets_3rd <- c(all_onsets_3rd, thirst_cluster_onset[!is.na(thirst_cluster_onset)])
    
    all_licks <- c(all_licks, lcks[!is.na(cluster_onset)])
    all_licks_2nd <- c(all_licks_2nd, lcks[!is.na(second_cluster_onset)])
    all_licks_3rd <- c(all_licks_3rd, lcks[!is.na(thirst_cluster_onset)])
    
    
    all_cor <- c(all_cor, cr)
    all_scor <- c(all_scor, scr)
    all_tcor <- c(all_tcor, tcr)
    
    print(sprintf("1st cluster %.3f, 2nd cluster %.3f, 3rd cluster %.3f", cr, scr, tcr))
  }
  
  cor_mat_all <- cbind(all_cor, all_scor)[1:21,]
  max_cor <- apply(cor_mat_all, 1, max)
  vls <- apply(cor_mat_all, 1, which.max)
  combined_licks <- lapply(1:len(vls), function(idx) {if(idx==9 || vls[idx] == 1) {all_res[[idx]]$licks} else if (vls[idx] == 2) {all_res[[idx]]$licks_2nd}})# else {all_res[[idx]]$licks_3rd}})
  combined_onsets <- lapply(1:len(vls), function(idx) {if(idx==9 ||vls[idx] == 1) {all_res[[idx]]$onsets} else if (vls[idx] == 2) {all_res[[idx]]$onsets_2nd}})# else {all_res[[idx]]$onsets_3rd}})
  combined_reward_onsets <- lapply(1:len(vls), function(idx) {if(idx==9 ||vls[idx] == 1) {all_res[[idx]]$reward_onsets} else if (vls[idx] == 2) {all_res[[idx]]$reward_onsets_2nd}})# else {all_res[[idx]]$onsets_3rd}})
  combined_correlation_tests <- lapply(1:len(vls), function(idx) {if(idx==9 ||vls[idx] == 1) {all_res[[idx]]$cor_test} else if (vls[idx] == 2) {all_res[[idx]]$cor_test_2nd}})
   

  hist_df_hunger <- data.frame(x=(unlist(combined_onsets[1:10]) / 2 -unlist(combined_licks[1:10]) / 30))
  hist_df_thirst <- data.frame(x=(unlist(combined_onsets[11:21]) / 2 -unlist(combined_licks[11:21]) / 30))
  hist_df <- data.frame(x=(unlist(combined_onsets) / 2 -unlist(combined_licks) / 30))
  cor_df <- data.frame(y=max_cor, 
                       x=unlist(lapply(combined_correlation_tests, function(crt){crt$p.value})),
                       group=rep(c("Hungry", "Thirsty"), c(10,11)))
  
  
  ghist_hunger <- 
    ggplot(hist_df_hunger, aes(x=x)) + 
    geom_histogram(bins=50, color="gray30", fill="royalblue4", alpha=.35, size=.25) +
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Delta (s)") +
    xlim(-4,11) +
    ylab("Frequency") + 
    geom_vline(xintercept=0, linetype="dashed", size=.25) +
    big_text_base_plot_theme
  
  ghist_thirst <- 
    ggplot(hist_df_thirst, aes(x=x)) + 
    geom_histogram(bins=50, color="gray30", fill="royalblue4", alpha=.35, size=.25) +
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Delta (s)") +
    xlim(-4,11) +
    ylab("Frequency") + 
    geom_vline(xintercept=0, linetype="dashed", size=.25) +
    big_text_base_plot_theme
  
  ghist <- 
    ggplot(hist_df, aes(x=x)) + 
    geom_histogram(bins=50, color="gray30", fill="royalblue4", alpha=.35, size=.25) +
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Delta (s)") +
    ylab("Frequency") + 
    geom_vline(xintercept=0, linetype="dashed", size=.25) +
    big_text_base_plot_theme
  
  gbox <- 
    ggplot(cor_df, aes(x=group, y=y)) +
      geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75) + 
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA)  + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylim(0,.8) +
      ylab("Correlation (r)") +
      xlab("")
      
    gpoint_scatter <- 
      ggplot(cor_df, aes(x=x, y=y, color=group)) +
        geom_point() + 
        theme(legend.position="top") +
        big_text_base_plot_theme +
        ylab("Correlation (r)") +
        xlab("P-value") +
        geom_vline(xintercept=0.05, linetype="dashed", size=.25) +
        xlim(0,.125) + 
        ylim(0,.8)
      
  
  box_scatter_p <- 
    plot_grid(gbox, gpoint_scatter, rel_widths = c(3,2), nrow=1)
  
  
  correlations_stat_df <- lapply(1:len(combined_correlation_tests),
                                 function(crt_obj_idx) {
                                   crt_obj <- combined_correlation_tests[[crt_obj_idx]]
                                   res_df <- data.frame(data_name=datasets_names[crt_obj_idx],
                                                        cor_value=max_cor[crt_obj_idx],
                                                        `statistic`=crt_obj$statistic,
                                                        `parameter`=crt_obj$parameter,
                                                        `p.value`=crt_obj$p.value,
                                                        `method`=crt_obj$method,
                                                        `alternative`=crt_obj$alternative,
                                                        `data.name`=crt_obj$data.name,
                                                        `signif.code`=signif.num(crt_obj$p.value),
                                                        N_licks = len(combined_licks[[crt_obj_idx]]),
                                                        N_clusters = len(combined_onsets[[crt_obj_idx]]),
                                                        mean_licks = mean(combined_licks[[crt_obj_idx]]),
                                                        mean_clusters = mean(combined_onsets[[crt_obj_idx]]),
                                                        median_licks = median(combined_licks[[crt_obj_idx]]),
                                                        median_clusters = median(combined_onsets[[crt_obj_idx]]),
                                                        sd_licks = sd(combined_licks[[crt_obj_idx]]),
                                                        sd_clusters = sd(combined_onsets[[crt_obj_idx]]),
                                                        sem_licks = sem(combined_licks[[crt_obj_idx]]),
                                                        sem_clusters = sem(combined_onsets[[crt_obj_idx]]))
                                 })
  
  
  
  values_df <- data.frame(mean_cor_hunger=mean(cor_df[1:10,1]),
                          mean_cor_thirst=mean(cor_df[11:21,1]),
                          median_cor_hunger=mean(cor_df[1:10,1]),
                          median_cor_thirst=mean(cor_df[11:21,1]),
                          N_cor_hunger=len(cor_df[1:10,1]),
                          N_cor_thirst=len(cor_df[11:21,1]),
                          sem_cor_hunger=sem(cor_df[1:10,1]),
                          sem_cor_thirst=sem(cor_df[11:21,1]),
                          sd_cor_hunger=sd(cor_df[1:10,1]),
                          sd_cor_thirst=sd(cor_df[11:21,1]),
                          mean_delta_hunger=mean(hist_df_hunger[,1]),
                          mean_delta_thirst=mean(hist_df_thirst[,1]),
                          median_delta_hunger=mean(hist_df_hunger[,1]),
                          median_delta_thirst=mean(hist_df_thirst[,1]),
                          N_delta_hunger=len(hist_df_hunger[,1]),
                          N_delta_thirst=len(hist_df_thirst[,1]),
                          sem_delta_hunger=sem(hist_df_hunger[,1]),
                          sem_delta_thirst=sem(hist_df_thirst[,1]),
                          sd_delta_hunger=sd(hist_df_hunger[,1]),
                          sd_delta_thirst=sd(hist_df_thirst[,1]),
                          mean_delta_all=mean(hist_df[,1]),
                          median_delta_all=mean(hist_df[,1]),
                          N_delta_all=len(hist_df[,1]),
                          sem_delta_all=sem(hist_df[,1]),
                          sd_delta_all=sd(hist_df[,1]))
  
  correlations_stat_df <- do.call(rbind,correlations_stat_df)
  
  values_df <- t(values_df)

  for (size_name in names(a4_sizes)) {
    
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_delta_hist.pdf",base_figure_path, size_name), height=size, width=size)
    plot(ghist)
    dev.off()
    
    pdf(sprintf("%s\\%s_delta_hist_hunger.pdf",base_figure_path, size_name), height=size, width=size)
    plot(ghist_hunger)
    dev.off()
    
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_delta_hist_thirst.pdf",base_figure_path, size_name), height=size, width=size)
    plot(ghist_thirst)
    dev.off()
    
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_cor_boxplots.pdf",base_figure_path, size_name), height=size, width=size)
    plot(gbox)
    dev.off()
    
    pdf(sprintf("%s\\%s_cor_scatter_pvalue.pdf",base_figure_path, size_name), height=size, width=size)
    plot(gpoint_scatter)
    dev.off()
    
    pdf(sprintf("%s\\%s_combined_box_pval_scatter.pdf",base_figure_path, size_name), height=size, width=size * (5/3))
    plot(box_scatter_p)
    dev.off()
  }
  
  write.csv(file=sprintf("%s//cor_statistics.csv",base_figure_path), correlations_stat_df)
  write.csv(file=sprintf("%s//values.csv",base_figure_path), values_df)
  
  
  
  for (di in 1:len(combined_licks)) {
    as <- cbind(unlist(combined_licks[[di]]), unlist(combined_onsets[[di]]))
    
    asd <- data.frame(as)
    colnames(asd) <- c("Licks", "Clusters")
    
    
    lmodel <- lm(formula=Licks~Clusters, data=asd)
    sm <- summary(lmodel)
    print(sm$r.squared)
    gt <- 
      ggplot(asd, aes(x=Clusters / 2, y=Licks / 30)) +
      geom_point() + 
      geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
      ggtitle(sprintf("%.4f - %.4f", sqrt(sm$r.squared), sm$coefficients[2,4])) +
      big_text_base_plot_theme + 
      xlab("Cluster onset (0.5s bin)") + 
      ylab("Lick onset (Frame)") +
      theme(plot.title = element_text(size=9))
    
    for (size_name in names(small_a4_sizes)) {
      
      dir.create(sprintf("%s\\%s", base_figure_path, datasets_names[[di]]))
      
      size = small_a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_anticipation_correlation.pdf",
                  base_figure_path,
                  datasets_names[di],
                  size_name),
          height=size,
          width=size)
      
      plot(gt)
      
      dev.off()
    }
  }
}

get_task_lick_distribution_all <- function(m_obj, path)
{
  
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  trials <- m_obj$stim_master_mat[m_obj$stim_master_mat[,"TrialType"] %in% c(3,4,5),]  
  run_indices <- unique(trials[,ncol(trials)])
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  
  binned_frames_per_mat <- unlist(frames_per_mat) / 15
  binned_frames_to_add <- cumsum(binned_frames_per_mat) - binned_frames_per_mat[1]
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {get_licking_vec(path, run_ind)})
  
  
  names(lick_vectors_all) <- run_indices
  
  
  lick_mat <- c()
  cluster_mat <- c()
  aligned_cluster_mat <- c()
  aligned_lick_mat <- c()
  idxvc <- c()
  tt <- c()

  annot_df <- m_obj$annot_df[m_obj$annot_df[,1] %in% c(3,4,5),]
  cluster_trials <- m_obj$trials_mat[m_obj$annot_df[,1] %in% c(3,4,5),]
  
  for (trial_idx in 1:nrow(trials)) {
    
    relevant_trial <- trials[trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    licking_indices_of_interest <- trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))
    
    if (max(licking_indices_of_interest) > len(lick_vectors_all[[as.character(run)]])) {
      if (verbose) {
        print("Exceeding max!")
      }
      next
    }
    
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][licking_indices_of_interest]
    
    
    if( sum(licks_in_trial_unbinned > 0) <= 0 ) {
      if (verbose) {
        print(sprintf("%d. No licks in trial continue", trial_idx))
      }
      next
    }
    
    
    cluster_labels_of_trial <- cluster_trials[trial_idx,]  

    
    
    
    
    if(all((!(which(licks_in_trial_unbinned > 0)) > 60) | (which(licks_in_trial_unbinned > 0) > 120))) {
      if (verbose) {
        print(sprintf("%d. Only anticipatory licks?", trial_idx))
      }
      next
    }
    
    unbinned_first_lick_in_respone <- which(licks_in_trial_unbinned > 0)[which(which(licks_in_trial_unbinned > 0) > 60)[1]]
    binned_first_lick_in_response <- get_binned_index(unbinned_first_lick_in_respone, 15)
    
    if (binned_first_lick_in_response + 3 > 20) {
      binned_first_lick_in_response = 20 - 3
    }
    
    cluster_mat <- rbind(cluster_mat,cluster_labels_of_trial)
    lick_mat <- rbind(lick_mat,as.numeric(licks_in_trial_unbinned))
    
    aligned_cluster_mat <- rbind(aligned_cluster_mat,
                                 cluster_labels_of_trial[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 12)])
    
    aligned_lick_mat <- rbind(aligned_lick_mat,
                              as.numeric(licks_in_trial_unbinned[(unbinned_first_lick_in_respone - 59):(unbinned_first_lick_in_respone + 180)]>0))
    
    tt <- c(tt, relevant_trial["Response"])
    idxvc <- c(idxvc, trial_idx)

  }
  
  
  
  
  return(list(cluster_mat=cluster_mat,
              lick_mat=lick_mat,
              aligned_cluster_mat=aligned_cluster_mat,
              aligned_lick_mat=aligned_lick_mat,
              tt=tt,
              trrial_indices=idxvc))
}

plot_trials_to_hit_similarity_dist <- function()
{
  
  
  base_figure_path <- sprintf("%s\\trials_to_hit_similarity_dist", output_path)
  dir.create(base_figure_path)
  
  fa_pooled_sim <- c()
  cr_pooled_sim <- c()
  miss_pooled_sim <- c()
  for (m_idx in 1:len(metadata_TT)) {
    
    m_obj <- metadata_TT[[m_idx]]
    
    sim_mt <- cosine(t(m_obj$trials_mat))
    
    trials_annotation <- m_obj$annot_df
    
    hit_sim_mt <- sim_mt[trials_annotation[,1] == 3 & trials_annotation[,2] == 0,]
    
    fa_to_hit_sim_mt <- hit_sim_mt[,trials_annotation[,1] %in% c(4:5) & trials_annotation[,2] == 1]
    cr_to_hit_sim_mt <- hit_sim_mt[,trials_annotation[,1] %in% c(4:5) & trials_annotation[,2] == 0]
    miss_to_hit_sim_mt <- hit_sim_mt[,trials_annotation[,1] %in% c(3) & trials_annotation[,2] == 1]
    
    fa_pooled_sim <- c(fa_pooled_sim, c(fa_to_hit_sim_mt))
    cr_pooled_sim <- c(cr_pooled_sim, c(cr_to_hit_sim_mt))
    miss_pooled_sim <- c(miss_pooled_sim, c(miss_to_hit_sim_mt))
    
  }
  
  
  
  pooled_similarity <- list(resp_na_hit=fa_pooled_sim,
                            nresp_all_hit=c(cr_pooled_sim,miss_pooled_sim))
                            #miss_hit=miss_pooled_sim)
  N_breaks <- 20
  hist_df <- 
    lapply(names(pooled_similarity),
           function(lst_name) {
             lst <- pooled_similarity[[lst_name]]
             
             h <- hist(lst, plot=F, seq(min(unlist(pooled_similarity), na.rm=T),
                                        max(unlist(pooled_similarity), na.rm=T),
                                        length.out=(N_breaks + 1)))
             
             df <- 
               data.frame(breaks=h$breaks,
                          freq=c(0, h$counts/sum(h$counts)),
                          group=rep(lst_name, times=(N_breaks + 1)))
             
             return(df)
           })
  
  all_df <- do.call(rbind,hist_df)
  gtrial_sim <- 
  ggplot(all_df, aes(x=breaks, y=freq)) + 
    #geom_bar(aes(fill=group), color=adjustcolor("gray30", alpha=.85), size=.01, alpha=.2, width=.01,  position=position_nudge(), stat="summary") + 
    geom_line(aes(color=group),  size=.85, alpha=.85, stat="summary") + 
    #geom_ribbon(aes(fill=group),  color=NA, alpha=.25, stat="summary") +
    big_text_base_plot_theme + 
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Cosine similarity")  +
    ylab("Fraction") 
    #theme(legend.position="top")
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_trial_sim.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gtrial_sim)
    
    dev.off()
    
  }
}

plot_FA_to_hit_similarity <- function()
{
  
  
  base_figure_path <- sprintf("%s\\FA_to_hit_similarity", output_path)
  dir.create(base_figure_path)
  
  clust_sim_all <- c()
  lick_sim_all <- c()
  cor_all <- c()
  pval_all <- c()
  N_trials_all <- c() 
  paths = get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
  for (m_idx in (1:len(metadata_TT))[1:11]) {
    m_obj <- metadata_TT[[m_idx]]
    res = get_task_lick_distribution_all(m_obj ,paths[m_idx])
    
    cluster_mat <- res$aligned_cluster_mat
    lick_mat <- res$aligned_lick_mat
    tt <- res$tt
    sim_mt <- cosine(t(cluster_mat))
    lick_sim_mt <- 1/(as.matrix(dist(t(apply(lick_mat, 1, function(r) {rollmean(r, 100)})))) + 10^-30)#, method="manhattan"))
    #lick_sim_mt <- cosine((apply(lick_mat, 1, function(r) {rollmean(r, 100)})))#, method="manhattan"))
    #lick_sim_mt <- cor(t(lick_mat))
    
    fa_trials <- which(tt == 3| tt == 5)
    fa_a_trials <- which(tt == 3)
    fa_n_trials <- which(tt == 3)
    
  
    
    hit_trials <-  which(tt==0)
    
    lick_sim <- apply(lick_sim_mt[hit_trials, fa_trials], 2, median)
    clst_sim <- apply(sim_mt[hit_trials, fa_trials], 2, median)
    
    
    #lck_order <- order(first_lick_per_trial[c(fa_trials)])
    lck_order <- (order(lick_sim, decreasing=T))
    annot_df <- data.frame(Type=rep(c("Hit", "Fa"), times=c(len(hit_trials), len(fa_trials))))
    clst_vis <- cluster_mat[c(hit_trials, fa_trials[lck_order]),]
    #lick_vis <- t(apply(lick_mat, 1, function(r) {rollmean(r, 100)}))[c(hit_trials, fa_trials[order(lck_order)]),]
    lick_vis <- lick_mat[c(hit_trials, fa_trials[order(lck_order)]),]
    
    

    rownames(annot_df) <- 1:nrow(annot_df)
    rownames(clst_vis) <- 1:nrow(annot_df)
    rownames(lick_vis) <- 1:nrow(annot_df)
    
    
    cluster_color_label <- spec_cg(len(unique(m_obj$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:((len(unique(m_obj$cluster_mat$labs))) - 1))

    phlmt <- pheatmap(lick_vis,
                      cluster_rows=F,
                      cluster_cols=F,
                      legend=F,
                      col=c(`0`="white",`1`="black"),
                      border_col=NA,
                      annotation_row=annot_df,
                      annotation_names_row=F,
                      show_rownames=F,
                      show_colnames=F,
                      annotation_legend=F)
    phcmt <- pheatmap(clst_vis,
                      cluster_rows=F,
                      cluster_cols=F,
                      legend=F,
                      border_co=NA,
                      annotation_row=annot_df,
                      annotation_names_row=F,
                      annotation_legend=F,
                      show_rownames=F,
                      show_colnames=F,
                      col=cluster_color_label)

    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_cluster_mat.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[m_idx]),
          height=size,
          width=size)
      
      plot(phcmt[[4]])
      
      dev.off()
      
      pdf(sprintf("%s\\%s_%s_lick_mat.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[m_idx]),
          height=size,
          width=size)
      
      plot(phlmt[[4]])
      
      dev.off()
    }

    
    lick_a_sim <- colMeans(lick_sim_mt[hit_trials, fa_trials])
    clst_a_sim <- colMeans(sim_mt[hit_trials, fa_trials])
    lick_n_sim <- colMeans(lick_sim_mt[hit_trials, fa_a_trials])
    clst_n_sim <- colMeans(sim_mt[hit_trials, fa_n_trials])
    
    crt <- cor.test(clst_sim, lick_sim)
    crt_a <- cor(clst_a_sim, lick_a_sim)
    crt_n <- cor(clst_n_sim, lick_n_sim)
    print(sprintf("%.3f A: %.3f, N: %.3f [%.4f]", crt$estimate, crt_a, crt_n, crt$p.value))
    
    #crt <- cor.test(clst_sim[clst_sim > 0], lick_sim[clst_sim> 0])
    plot(clst_sim, lick_sim)
    
    clust_sim_all <- c(clust_sim_all, clst_sim)
    lick_sim_all <- c(lick_sim_all, lick_sim)
    pval_all <- c(pval_all, crt$p.value)
    cor_all <- c(cor_all, crt$estimate)
    N_trials_all <- c(N_trials_all, len(fa_trials))
    
    crt <- cor(clst_sim[lick_sim > 0], lick_sim[lick_sim> 0])
    print(crt)
  }
  
  cor_df <- data.frame(y=cor_all, x=pval_all)
  #gpoint_scatter <- 
    ggplot(cor_df, aes(x=x, y=y)) +
    geom_point() + 
    theme(legend.position="top") +
    big_text_base_plot_theme +
    ylab("Correlation (r)") +
    xlab("P-value") +
    geom_vline(xintercept=0.05, linetype="dashed", size=.25)
    #geom_hline(yintercept=-.15, linetype="dashed", size=.25)
  
  
  
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_point_scatter.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gpoint_scatter)
    
    dev.off()
    
  }
  
  
}
