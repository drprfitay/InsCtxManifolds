

get_signif <- function(p.values) { symnum(p.values,
                                          corr = FALSE, 
                                          na = FALSE, 
                                          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
                                          symbols = c("****"," ***", "**", "*", ".", " "))[[1]]}


output_path_f <- "~/../Desktop/Itay_group_meeting_dfs/"



get_thirsty_quenched_paths <- function() {
    mice <- list(`44`=
                 c(170518, 170519, 170523, 170524),
               `47`=
                 c(171121, 171213, 171214, 180102),
               `52`=
                 c(180404, 180405, 180329, 180522),
               `56`=
                 c(180628),
               `57`=
                 c(180621))
               
  thirsty_quenched_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    thirsty_quenched_paths <- c(thirsty_quenched_paths,
                                unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(thirsty_quenched_paths)
}

get_intrinsic_dimensionality <- function(paths, preset, original_mat=F, window_size=15, chunk=-1) {
  dimensions <- c()
  
  
  if (original_mat) {
    if (chunk == -1) { chunk <-  0 }
    
    for (p in paths) {
      
      
      mt <- get_reduced_mat_full_day(p, just_original_mat=T, window_size=window_size, chunk=chunk)
      dimensions <- c(dimensions,
                      twonn(mt, method = "linfit", c_trimmed = 0.05)$est[2])
    }   
    
    return(dimensions)
  }
  
  for (p in paths) {
    mt <- get_mat_with_preset(p, preset, chunk=chunk)
     dimensions <- c(dimensions,
                     twonn(mt, method = "linfit", c_trimmed = 0.05)$est[2])
  }
  
  return(dimensions)
}




get_colors_for_mat <- function(path, mat, chunk=-1, alpha=0.1, base_pal=viridis) {
  
  stim_master_mat <- 
    get_stim_mat(path,just_mat = T, window_size = 15, chunk=chunk)    
  behav <- get_behavior_variables(stim_master_mat, mat_size = nrow(mat), window_size = 15, seconds_after=1)
  cp_base <- base_pal(nrow(mat))
  reward_col <- adjustcolor(cp_base, alpha=alpha); reward_col[behav[,"Reward"] == 1] <- "red"
  tt3_col <- adjustcolor(cp_base, alpha=alpha); tt3_col[behav[,"TrialType 3"] == 1] <- "red"
  tt4_col <- adjustcolor(cp_base, alpha=alpha); tt4_col[behav[,"TrialType 4"] == 1] <- "red"
  tt5_col <- adjustcolor(cp_base, alpha=alpha); tt5_col[behav[,"TrialType 5"] == 1] <- "red"
  false_alarms_col <- adjustcolor(cp_base, alpha=alpha); false_alarms_col[behav[,"FalseAlarms"] == 1] <- "red"
  cumreward_col <- adjustcolor(rev(spec_cg(len(unique(behav[,"Cumulative reward"]))))[behav[,"Cumulative reward"]],
                               alpha)
  
  
  ent_trials = c("EntireTrials3", "EntireTrials4", "EntireTrials5",
                 "Trials3ITI", "Trials4ITI", "Trials5ITI")
  ent_trials_cols <- list()
  col <- rep("gray60", times=(nrow(mat)))
  
  for (et in ent_trials) {
    
    adj_col <- adjustcolor(col, alpha=alpha);
    
    for (i in 1:max(behav[,et])) {
      adj_col[which(behav[,et] == i)] <- adjustcolor(base_pal(max(behav[,et]))[i], alpha=0.5)
    }
    
    ent_trials_cols[[et]] <- adj_col
  }
  
  
  return(list(time=cp_base,
              reward=reward_col,
              tt3_col=tt3_col,
              tt4_col=tt4_col,
              tt5_col=tt5_col,
              false_alarms_col=false_alarms_col,
              cumreward_col=cumreward_col,
              EntireTrials3=ent_trials_cols$EntireTrials3,
              EntireTrials4=ent_trials_cols$EntireTrials4,
              EntireTrials5=ent_trials_cols$EntireTrials5,
              Trials3ITI=ent_trials_cols$Trials3ITI,
              Trials4ITI=ent_trials_cols$Trials4ITI,
              Trials5ITI=ent_trials_cols$Trials5ITI))
}


pair_plots_video <- function(path, prefix="") {
  
  dir.create("./tmp_input")
  dir.create("./tmp_output")
  
  mt <- get_mat_with_preset(path, "MilkaChoc", chunk=1)
  cmt <- get_colors_for_mat(path, mt, 1)
  
  
  dimensions <- combn(ncol(mt), 2)

  
  reward_vec <- as.character(as.numeric(cmt$reward == "red"))
  
  all_plots <- list()
  
  for (i in 1:ncol(dimensions)) {
    
    d1 <- dimensions[1,i]
    d2 <- dimensions[2,i]
    
    
    df <- data.frame(X=mt[,d1], Y=mt[,d2], REW=reward_vec)
    g <- ggplot(df, aes(x=X, y=Y)) + 
         geom_point(aes(color=REW)) +
         xlab(sprintf("Dim %d", d1)) +
         ylab(sprintf("Dim %d", d2)) + 
         theme_light() + 
         scale_color_manual(breaks=c("0", "1"),values=c(adjustcolor("gray50", alpha=0.1), 
                              adjustcolor("blue", alpha=0.5))) +
         theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               legend.position="NA")
    
    
    all_plots <- append(all_plots, list(g))
  }
  
  
  for (pt in 1:nrow(mt)) {
    all_plots_with_points <- list()
    for (i in 1:ncol(dimensions)) {
      
      d1 <- dimensions[1,i]
      d2 <- dimensions[2,i]
      
      pt_df <- data.frame(X=mt[pt,d1],
                          Y=mt[pt,d2])
      
      
      gf <- 
      all_plots[[i]] +
        geom_point(data=pt_df, col="red", size=3)
      
      all_plots_with_points <- append(all_plots_with_points, list(gf))
    } 
    
    
    
    all_plots_with_points$nrow = 5
    final_f <- do.call(arrangeGrob,all_plots_with_points)
    png(sprintf("./tmp_input/image_%d.png", pt), width = 500, height = 720, res = 58)
    plot(final_f)
    dev.off()
  }

  
  png_files <- sprintf("./tmp_input/image_%d.png", 1:nrow(mt))
  av::av_encode_video(png_files, sprintf('./tmp_output/output_%s.mp4', prefix), framerate = 20)
  

  
}


get_num_runs <- function(path, control=F) {
  
  if (control) {
    mt <- readMat(sprintf("%s%s", path, "dff.mat"))
    return(max(dim(mt[[1]])) / 57600)
  }
  runs <- list.files(path, recursive=F)
  runs <- runs[which(grepl("\\.R", runs))]
  return(len(runs))
}

cluster_map_centroid_correlation <- function(path_list) {
  
  output_path <- sprintf("%s//centroid_cluster_map//", output_path_f)
  dir.create(output_path)
  path_list <- get_thirsty_quenched_paths() 
  combination_matrix <- combn(len(path_list), 2)
  
  
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  annotation_df <- data.frame(Mice=mice_names)
  
  rownames(annotation_df) <- datasets_names
  
  indices  <- sample(1:len(mice_names), 3)
  
  while (!len(unique(mice_names[indices])) == 3) {
    indices  <- sample(1:len(mice_names), 3)
  }
  

  
  
  by_chunk_list <- c()
  
  for (chunk in c(1, -1)) {
    all_cluster_mats <- list()
    all_mats <- list()
  
  for (p in path_list) {
    mat <- get_mat_with_preset(p, "dalshtaim", chunk=chunk)
    cluster_mat <- get_clusters_mat(mat, min_frames_for_minc = 5, min_frames_for_outward = 5)
    
    all_cluster_mats <- append(all_cluster_mats,
                               list(cluster_mat))
    
    all_mats <- append(all_mats,
                       list(mat))
  }
  
  
  tables <- lapply(all_cluster_mats, function(clust_mat) {table(clust_mat)})
  tables_normed <- lapply(tables, function(tbl){tbl/sum(tbl)})
  max_cluster_num <- max(unlist(lapply(tables_normed, function(tab) {max(as.numeric(names(tab)))})))
  tables_normed_equal_length <- lapply(tables_normed,
                                       function(cluster_prob)
                                       { 
                                         cluster_names <- names(cluster_prob)
                                         missing <- as.character(0:max_cluster_num)[!as.character(0:max_cluster_num) %in% cluster_names]
                                         final_clusters <- c(cluster_prob, rep(0, times=(max_cluster_num + 1- len(cluster_prob))))
                                         names(final_clusters) <- c(cluster_names, missing)
                                         return(final_clusters[order(as.numeric(names(final_clusters)))])
                                         })
  
  clusters_table <- do.call(rbind, tables_normed_equal_length)
  binarized_table <- t(apply(clusters_table, 1, function(r) {as.numeric(r > 0)}))
  possible_configurations <- colSums(binarized_table) / nrow(binarized_table)
  names(possible_configurations) <- 0:(len(possible_configurations) - 1)
  
  possible_df <- data.frame(pct=possible_configurations * 100, 
                            x=factor(names(possible_configurations), 
                                     levels=as.character(0:(len(possible_configurations) - 1))))
  
  gposs_df <- ggplot(possible_df) + 
              geom_bar(aes(x=x, y=pct), stat="identity") + 
              xlab("Number of clusters (#)") +
              ylab("Fraction of datasets (%)") +
              ggtitle("Possible clusters configurations") +
              base_plot_theme + 
              theme(axis.text.x = element_text(angle = 90)) +
              scale_x_discrete(breaks=seq(0, max_cluster_num, 10))
    
    
  
  gposs_df_subset <- 
              ggplot(possible_df[2:10,]) + 
              geom_bar(aes(x=x, y=pct), stat="identity") + 
              xlab("Number of clusters (#)") +
              ylab("Fraction of datasets (%)") +
              ggtitle("Possible clusters configurations") +
              base_plot_theme
  
  prob_df <- data.frame(prob=colMeans(clusters_table), 
                        x=factor(names(possible_configurations), 
                                 levels=as.character(0:(len(possible_configurations) - 1))))
  
  gprob_df <- ggplot(prob_df[1:10,]) + 
                     geom_bar(aes(x=x, y=prob), stat="identity") + 
                     xlab("Number of clusters (#)") +
                     ylab("Probability") +
                     ggtitle("Probability for different cluster configurations") +
                     base_plot_theme
  
  gf <- grid.arrange(gprob_df, gposs_df, gposs_df_subset, nrow=1)
  
  pdf(sprintf("%s/clusters_configurations_%s.pdf", output_path, ifelse(chunk == 1, "First", "All")),
      height=3,
      width=3 * 3)
  plot(gf)
  dev.off()
  
  highest_3 <- sort(possible_df[possible_df$pct == 100,]$x, decreasing=T)[1:3]

  for (ind in indices) {
    tmp <- all_cluster_mats[ind][[1]]
    mat <- all_mats[ind][[1]]
    cent_min <- get_centroid_mat(mat, tmp, highest_3[1], min=T)
    cent_max <- get_centroid_mat(mat, tmp, highest_3[1], min=F)
    
    rownames(tmp) <- c()
    colnames(tmp) <- c()
    
    pdf(sprintf("%s/example_%s_%s_cluster_heatmap.pdf", output_path, datasets_names[ind], ifelse(chunk == 1, "First", "All")),
        height=3,
        width=3)
    pheatmap(tmp, col=spec_cg(500), border_color = NA, cluster_rows=F, cluster_cols=F,
             main=datasets_names[ind])
    dev.off()
    
    pairs_min <- pair_plots_gg_clusters(mat, labels =  cent_min[[3]], colors = cent_min[[2]])
    pairs_max <- pair_plots_gg_clusters(mat, labels =  cent_max[[3]], colors = cent_max[[2]])
    dev.off()

    
    png(sprintf("%s/example_%s_%s_min.png", output_path, datasets_names[ind], ifelse(chunk == 1, "First", "All")),
        units="in",
        res=1000,
        height=8,
        width=6)
    plot(pairs_min)
    dev.off()  
    
    
    png(sprintf("%s/example_%s_%s_max.png", output_path, datasets_names[ind], ifelse(chunk == 1, "First", "All")),
        units="in",
        res=1000,
        height=8,
        width=6)
    plot(pairs_max)
    dev.off() 
    
  }
  
  cluster_map_cors <- c()
  centroid_cors_1 <- c()
  centroid_cors_2 <- c()
  centroid_cors_3 <- c()
  
  centroid_cors_matrix_1 <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  centroid_cors_matrix_2 <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  centroid_cors_matrix_3 <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  cluster_map_cors_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  
  
  for (pair_idx in 1:ncol(combination_matrix)) {
    
    idx_1 <- combination_matrix[1,pair_idx]
    idx_2 <- combination_matrix[2,pair_idx]
    
    cluster_mat_1 <- all_cluster_mats[idx_1][[1]]
    cluster_mat_2 <- all_cluster_mats[idx_2][[1]]
    
    cluster_map_pair_cor <- cor(c(cluster_mat_1),  c(cluster_mat_2))
    cluster_map_cors <- c(cluster_map_cors, cluster_map_pair_cor)
    cluster_map_cors_matrix[idx_1, idx_2] <- cluster_map_pair_cor
    cluster_map_cors_matrix[idx_2, idx_1] <- cluster_map_pair_cor
    
    for (nc in highest_3) {
    
    mat_1 <- all_mats[combination_matrix[1,pair_idx]][[1]]
    mat_2 <- all_mats[combination_matrix[2,pair_idx]][[1]]
     
    centroid_mat_1 <- get_centroid_mat(mat_1, cluster_mat_1, nc)
    centroid_mat_2 <- get_centroid_mat(mat_2, cluster_mat_2, nc)
     
    print(nrow(centroid_mat_1[[1]]))
    print(nrow(centroid_mat_2[[1]]))
    
    dist_matrix_day1 <- as.matrix(dist(centroid_mat_1[[1]]))
    dist_matrix_day2 <- as.matrix(dist(centroid_mat_2[[1]]))
    
    
    h1 <- hclust(dist(dist_matrix_day1), method="single")
    h2 <- hclust(dist(dist_matrix_day2), method="single")

    centroid_pair_cor <- cor(c(dist_matrix_day1[h1$order, h1$order]),  c(dist_matrix_day2[h2$order, h2$order]))
    
    if (nc == highest_3[1]) {
      centroid_cors_1 <- c(centroid_cors_1, centroid_pair_cor)
      centroid_cors_matrix_1[idx_1, idx_2] <- centroid_pair_cor
      centroid_cors_matrix_1[idx_2, idx_1] <- centroid_pair_cor
    } else if (nc == highest_3[2]) {
      centroid_cors_2 <- c(centroid_cors_2, centroid_pair_cor)
      centroid_cors_matrix_2[idx_1, idx_2] <- centroid_pair_cor
      centroid_cors_matrix_2[idx_2, idx_1] <- centroid_pair_cor
    } else {
      centroid_cors_3 <- c(centroid_cors_3, centroid_pair_cor)
      centroid_cors_matrix_3[idx_1, idx_2] <- centroid_pair_cor
      centroid_cors_matrix_3[idx_2, idx_1] <- centroid_pair_cor
    }
    
    }

  }
  
  rownames(cluster_map_cors_matrix) <- datasets_names
  colnames(cluster_map_cors_matrix) <- datasets_names
  rownames(centroid_cors_matrix_3) <- datasets_names
  colnames(centroid_cors_matrix_3) <- datasets_names
  rownames(centroid_cors_matrix_2) <- datasets_names
  colnames(centroid_cors_matrix_2) <- datasets_names
  rownames(centroid_cors_matrix_1) <- datasets_names
  colnames(centroid_cors_matrix_1) <- datasets_names
  
  matrices_list <- list(`ClusterMap`=cluster_map_cors_matrix,
                        `CentroidCors1`=centroid_cors_matrix_1,
                        `CentroidCors2`=centroid_cors_matrix_2,
                        `CentroidCors3`=centroid_cors_matrix_3)
  
  for (matrix_name in names(matrices_list)) {
    pdf(sprintf("%s/%s_matrix_%s.pdf", output_path, matrix_name, ifelse(chunk == 1, "First", "All")),
        height=5,
        width=5.3)
    pheatmap(matrices_list[[matrix_name]], border_color = "NA", annotation_row = annotation_df, annotation_col = annotation_df)  
    dev.off()  
  }
  
  by_chunk_list <- append(by_chunk_list,
                          list(centroid_cors_1, centroid_cors_2, centroid_cors_3, cluster_map_cors))
  }
  
  
  centroid_cors_df <- data.frame(First_1=by_chunk_list[[1]],
                                 First_2=by_chunk_list[[2]],
                                 First_3=by_chunk_list[[3]],
                                 All_1=by_chunk_list[[5]],
                                 All_2=by_chunk_list[[6]],
                                 All_3=by_chunk_list[[7]])
  
  melted_centroid_df <- melt(centroid_cors_df)
  colnames(melted_centroid_df) <- c("X", "Y")
  
  centroid_box <- my_boxplot(melted_centroid_df, xl="Group", y="Corr (Pearson's r)", yrange=c(0,1))
  
  
  cluster_map_cors_df <- data.frame(First=by_chunk_list[[4]], All=by_chunk_list[[8]])
  melted_cluster_map_df <- melt(cluster_map_cors_df)
  colnames(melted_cluster_map_df) <- c("X", "Y")
  cluster_map_box <- my_boxplot(melted_cluster_map_df, xl="Group", y="Corr (Pearson's r)", yrange=c(0,1))
  
  
  pdf(sprintf("%s/centroid_boxplot.pdf", output_path),
      height=3, width=6)
  plot(centroid_box)
  dev.off()
  
  
  pdf(sprintf("%s/cluster_map_boxplot.pdf", output_path),
      height=3, width=3)
  plot(cluster_map_box)
  dev.off()
  
  
  names(by_chunk_list) <- c("First 1 cent", "First 2 cent", "First 3 cent", "First clustmap",
                            "All 1 cent",   "All 2 cent",   "All 3 cent", "All clustmap")
  
  final_df <- by_chunk_list
  save(file=sprintf("%s/centroid_clust_maps.Rda", output_path),
       final_df)
       
}

structure_distance_corr <- function(compare_shuffle = F, nclusters=100, chunk=1, nreps=5) {
  
  path_list <- get_thirsty_quenched_paths()
  combination_matrix <- combn(len(path_list), 2)
  
  
  results_all <- list()
  for (compare_shuffle in c(F ,T)) {
  
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
  
  for (pair_idx in 1:ncol(combination_matrix)) {
    day1 <- path_list[combination_matrix[1,pair_idx]]
    day2 <- path_list[combination_matrix[2,pair_idx]]
    
    s1 <- F
    s2 <- F
    
    presets <- c("dalshtaim", "dalshtaim")
    
    if (!compare_shuffle) {
      
      day1_mat <- get_mat_with_preset(day1, preset=presets[1], chunk=chunk)
      day2_mat <- get_mat_with_preset(day2, preset=presets[2], chunk=chunk)
      
    }
    
    struct_corr_pair <- list()
    
    for (metric_name in names(metric_functions)) {
      
      struct_corr_pair[[metric_name]] <- c()
      
    }
    

    for (i in 1:nreps) {
      
      
      
      if (compare_shuffle) {
        presets <- c("dalshtaim", "dalshtaim")
        presets[sample(1:2, 1)] <- "dalshtaimshuff"
        day1_mat <- get_mat_with_preset(day1, preset=presets[1], chunk=chunk)
        day2_mat <- get_mat_with_preset(day2, preset=presets[2], chunk=chunk)
        
      }
      
    
      km_day1 <- kmeans(apply(day1_mat, 2, scale), nclusters, iter.max=300)
      km_day2 <- kmeans(apply(day2_mat, 2, scale), nclusters, iter.max=300)
  
      
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
    }
  }
  
  }
}

intrinsic_dim_analysis <- function() {
  path_list <- get_thirsty_quenched_paths()
  
  idims_df <- c()
  
  for (p in path_list) {
    
    num_of_runs <- get_num_runs(p) 
    
    first_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=1, window_size = 15)
    last_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15)
    all_mat_original <-  get_reduced_mat_full_day(p, just_original_mat = T, window_size = 15)
    
    if (max(dim(last_run_mat_original)) < max(dim(first_run_mat_original))) {
      
      if (num_of_runs > 2) {
        print("Going to one run before final")
        num_of_runs = num_of_runs - 1
        last_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15)
      } else {
        print("Only two runs!")
      }
    }
    
    
    shuffled_first_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=1, window_size = 15, shuffled = T, time_shuffled = F)
    shuffled_last_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15, shuffled = T, time_shuffled = F)
    shuffled_all_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, window_size = 15, shuffled = T, time_shuffled = F)
    

    first_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, chunk=1, knn1=0.075, knn2 = 0)
    last_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, chunk=num_of_runs, knn1=0.075, knn2 = 0)
    all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, knn1=0.075, knn2 = 0)    
    first_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, chunk=1, knn1=0.075, knn2 = 0)
    last_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, chunk=num_of_runs, knn1=0.075, knn2 = 0)
    all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0)
    first_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, chunk=1, knn1=0.075, knn2 = 0)
    last_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, chunk=num_of_runs, knn1=0.075, knn2 = 0)
    all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0)
    
    

    shuffled_first_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, knn1=0.075, knn2 = 0, chunk=1, shuffled=T, time_shuffled=F)
    shuffled_last_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, knn1=0.075, knn2 = 0, chunk=num_of_runs, shuffled=T, time_shuffled=F)
    shuffled_all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=50, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)    
    shuffled_first_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0, chunk=1, shuffled=T, time_shuffled=F)
    shuffled_last_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0, chunk=num_of_runs, shuffled=T, time_shuffled=F)
    shuffled_all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)
    shuffled_first_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, chunk=1, shuffled=T, time_shuffled=F)
    shuffled_last_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, chunk=num_of_runs, shuffled=T, time_shuffled=F)
    shuffled_all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)
    
    
    
    
    
    
    matrices <- list(first_run_mat_original,
                     last_run_mat_original,
                     all_mat_original,
                     shuffled_first_mat_original,
                     shuffled_last_mat_original,
                     shuffled_all_mat_original,
                     first_run_mat_reduced,
                     last_run_mat_reduced,
                     all_mat_reduced,
                     shuffled_first_mat_reduced,
                     shuffled_last_mat_reduced,
                     shuffled_all_mat_reduced)
    
    
    idims <- unlist(lapply(matrices, function(mat) {return(twonn(mat, method = "linfit", c_trimmed = 0.01)$est[2])}))
    
    idims_df <- rbind(idims_df,
                      idims)
    
  }
}


intrinsic_dim_analysis_hungry_sated <- function() {
  path_list <- get_thirsty_quenched_paths()
  
  idims_df <- c()
  path_list <- c(insula_hunger_paths,
                 v1_paths,
                 por_paths)
  
  for (p in path_list) {
  
    control = T
    
    if (p %in% insula_hunger_paths) {
      control = F
    }
    
    num_of_runs <- get_num_runs(p, control = control) 
    
    first_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=1, window_size = 15, activity_threshold=0.25, control=control)
    last_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15, activity_threshold=0.25, control=control)
    all_mat_original <-  get_reduced_mat_full_day(p, just_original_mat = T, window_size = 15, activity_threshold=0.25, control=control)
    
    if (max(dim(last_run_mat_original)) < max(dim(first_run_mat_original))) {
      
      if (num_of_runs > 2) {
        print("Going to one run before final")
        num_of_runs = num_of_runs - 1
        last_run_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15)
      } else {
        print("Only two runs!")
      }
    }
    
    
    shuffled_first_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=1, window_size = 15, shuffled = T, time_shuffled = F, activity_threshold=0.25, control=control)
    shuffled_last_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, chunk=num_of_runs, window_size = 15, shuffled = T, time_shuffled = F, activity_threshold=0.25, control=control)
    shuffled_all_mat_original <- get_reduced_mat_full_day(p, just_original_mat = T, window_size = 15, shuffled = T, time_shuffled = F, activity_threshold=0.25, control=control)
    
    
    first_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, chunk=1, knn1=0.075, knn2 = 0, activity_threshold=0.25, control=control)
    last_run_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, chunk=num_of_runs, knn1=0.075, knn2 = 0, activity_threshold=0.25, control=control)
    all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, activity_threshold=0.25, control=control)    

    
    get_mat_with_preset(p, "dalshtaim", activity_threshold = 0.25, oldscope=control, chunk=1)
    get_mat_with_preset(p, "dalshtaim", activity_threshold = 0.25, oldscope=control, chunk=num_of_runs)
    get_mat_with_preset(p, "dalshtaim", activity_threshold = 0.25, oldscope=control)
    
    shuffled_first_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, chunk=1, shuffled=T, time_shuffled=F, activity_threshold=0.25, control=control)
    shuffled_last_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, chunk=num_of_runs, shuffled=T, time_shuffled=F, activity_threshold=0.25, control=control)
    shuffled_all_mat_reduced <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F, activity_threshold=0.25, control=control)    
    
    get_mat_with_preset(p, "dalshtaimshuff", activity_threshold = 0.25, oldscope=control, chunk=1)
    get_mat_with_preset(p, "dalshtaimshuff", activity_threshold = 0.25, oldscope=control, chunk=num_of_runs)
    get_mat_with_preset(p, "dalshtaimshuff", activity_threshold = 0.25, oldscope=control)

    
    
    
    
    
    
    matrices <- list(first_run_mat_original,
                     last_run_mat_original,
                     all_mat_original,
                     shuffled_first_mat_original,
                     shuffled_last_mat_original,
                     shuffled_all_mat_original,
                     first_run_mat_reduced,
                     last_run_mat_reduced,
                     all_mat_reduced,
                     shuffled_first_mat_reduced,
                     shuffled_last_mat_reduced,
                     shuffled_all_mat_reduced)
    
    
    idims <- unlist(lapply(matrices, function(mat) {return(twonn(mat, method = "linfit", c_trimmed = 0.01)$est[2])}))
    
    idims_df <- rbind(idims_df,
                      idims)
    
  }
}


structure_cluster_occupancy_analysis <- function(path_list, window_size=15, max_num_of_clusters=7) {
  output_path <- sprintf("%s//occupancy_analysis//", output_path_f)
  dir.create(output_path)
  path_list <- get_thirsty_quenched_paths()
  window_size=15
  max_num_of_clusters=7
  
  tt3_df <- c()
  tt4_df <- c()
  tt5_df <- c()
  reward_df <- c()
  
  
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  for (p in path_list) {
    


    # first_mat <- get_mat_with_preset(p, "dalshtaim", chunk=1)
    
    # 
    # color_palettes_first <- get_colors_for_mat(p, first_mat, chunk=1)
    
    
    # 
    # clust_mat_first <- get_clusters_mat(first_run_mat_reduced, min_frames_for_minc = 5, min_frames_for_outward = 5)
    # centroid_mat_first <- get_centroid_mat(first_run_mat_reduced, clust_mat_first, 8)
    
    #all_mat <- get_mat_with_preset(p, "dalshtaim")
    all_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0)    
    color_palettes <- get_colors_for_mat(p, all_mat, alpha=0.1)
    clust_mat_all <- get_clusters_mat(all_mat, min_frames_for_minc = 5, min_frames_for_outward = 5)
    centroid_mat_all <- get_centroid_mat(all_mat, clust_mat_all, max_num_of_clusters, min=T)
    
  
    timepoints <- 
      lapply(unique(centroid_mat[[3]][centroid_mat[[3]] != -1]), 
             function(clst) {as.numeric(which(centroid_mat_all[[3]] != -1 & centroid_mat_all[[3]] == clst) * window_size)})
    
    max_len <- max(unlist(lapply(timepoints, length)))
    
    for(idx in 1:length(timepoints)) {
      length(timepoints[[idx]]) <-  max_len
    }
    
    timepoints_df <- do.call(cbind, timepoints)
  
  
    neuronal_clust_labels <- centroid_mat_all[[3]]
    num_of_rewards <- c()
    num_of_tt3 <- c()
    num_of_tt4 <- c()
    num_of_tt5 <- c()
    all_cluster_labels <- sort(unique(neuronal_clust_labels[neuronal_clust_labels != - 1]))
    
    for (lbl in all_cluster_labels) {

      indices <- which(neuronal_clust_labels == lbl)
      
      if ("red" %in% names(table(color_palettes$reward[indices]))) {
        num_of_rewards_in_cluster <- table(color_palettes$reward[indices])["red"]
      } else {
        num_of_rewards_in_cluster <- 0
      }
      
      if ("red" %in% names(table(color_palettes$tt3_col[indices]))) {
        num_of_tt3_in_cluster <- table(color_palettes$tt3_col[indices])["red"]
      } else {
        num_of_tt3_in_cluster <- 0
      }

      if ("red" %in% names(table(color_palettes$tt4_col[indices]))) {
        num_of_tt4_in_cluster <- table(color_palettes$tt4_col[indices])["red"]
      } else {
        num_of_tt4_in_cluster <- 0
      }
      
      if ("red" %in% names(table(color_palettes$tt5_col[indices]))) {
        num_of_tt5_in_cluster <- table(color_palettes$tt5_col[indices])["red"]
      } else {
        num_of_tt5_in_cluster <- 0
      }      
      
      num_of_rewards <- c(num_of_rewards, num_of_rewards_in_cluster /len(indices))
      num_of_tt3 <- c(num_of_tt3, num_of_tt3_in_cluster /len(indices))
      num_of_tt4 <- c(num_of_tt4, num_of_tt4_in_cluster /len(indices))
      num_of_tt5 <- c(num_of_tt5, num_of_tt5_in_cluster /len(indices))
    }
    
    reward_df <- rbind(reward_df, sort(num_of_rewards, decreasing=T))
    tt3_df <- rbind(tt3_df, sort(num_of_tt3, decreasing=T))
    tt4_df <- rbind(tt4_df, sort(num_of_tt4, decreasing=T))
    tt5_df <- rbind(tt5_df, sort(num_of_tt5, decreasing=T))
    
  }
  
  
  occupancy_dfs <- list(reward_df, tt3_df, tt4_df, tt5_df)
  
  occupancy_dfs <- 
  lapply(occupancy_dfs, function(df) {
    df <- t(apply(df, 1, function(r) {r/sum(r)}))
    rownames(df) <- datasets_names
    colnames(df) <- sprintf("Activity cluster %d", 1:max_num_of_clusters)
    return(df)
  })
  
  names(occupancy_dfs) <- c("Reward",
                             "Trial type 3",
                             "Trial type 4",
                             "Trial type 5")
  
  for (df_name in names(occupancy_dfs)) {

    df <- occupancy_dfs[[df_name]]  
    pdf(sprintf("%s/%s_occupancy.pdf", output_path, df_name),
        height=4, width=4)
    pheatmap(df, cluster_rows=F, cluster_cols=F, col=spec_cg(500))
    dev.off()
    
  }
  
  mean_occupancy <- 
    lapply(occupancy_dfs,
           function(df) {colMeans(df) / sum(colMeans(df))})
  
  
  occupancy_plot_list <- list()
  names(mean_occupancy) <- c("Reward",
                             "Trial type 3",
                             "Trial type 4",
                             "Trial type 5")
  for (event_name in names(mean_occupancy)) {
    occupancy <- mean_occupancy[[event_name]]
    occ_df <- data.frame(Occ=occupancy, Clust=names(occupancy))
    
    yl_text = ifelse(event_name == "Reward", 
                     "Occupancy (%)",
                     "")
    g <- 
      ggplot(occ_df, aes(x=Clust, y=Occ)) + 
      geom_bar(stat="identity") + 
      base_plot_theme + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1.1)) + 
      xlab("") +
      ylab(yl_text) + 
      ylim(c(min(unlist(mean_occupancy)),
             max(unlist(mean_occupancy)))) +
      ggtitle(event_name)
    
    
    occupancy_plot_list <- append(occupancy_plot_list,
                                  list(g))
    
  }
  
  occupancy_plot_list$nrow <- 1
  gf <- do.call(grid.arrange, occupancy_plot_list)
  
  pdf(sprintf("%s/mean_occupancy.pdf", output_path),
      height=4, width=12)
  plot(gf)
  dev.off()
  dev.off()
  
  
  for (idx in 1:len(path_list)) {
    
    
      p <- path_list[idx]

      all_mat <- get_mat_with_preset(p, "dalshtaim", chunk=-1)
      color_palettes <- get_colors_for_mat(p, all_mat, alpha=0.05, chunk=-1)
      #color_palettes_2 <- get_colors_for_mat(p, all_mat, alpha=0.05, chunk=-1, base_pal=inferno)
      
      
      
      for (nm in names(color_palettes)[8:13]) {
      
      #g2 <- pair_plots_colored(all_mat, color_palettes_2$false_alarms_col)
      
      
      png(sprintf("%s\\%s\\%s\\%s.png", 
                  output_path_f,
                  "occupancy_analysis",
                  str_replace(datasets_names[idx]," ", "_"),
                  nm),
                  units="in",
                  res=1000,
                  height=8,
                  width=6)
      
      g <- pair_plots_colored(all_mat, color_palettes[[nm]])
      plot(g)
      dev.off()
      
   
  }
  }
}
# 
# 
# # Create some PNG images
# png("input%03d.png", width = 1280, height = 720, res = 58)
# for(i in 1:1000){
#   tmpcol <- c(rep("gray50", 3840), "red")
#   fmt <- rbind(mt, mt[i,])
#   plot(fmt[1:3840,4], fmt[1:3840,5], col=tmpcol[1:3840], pch=19)
#   points(fmt[3841,4], fmt[3841,5], col=tmpcol[3841], pch=19, cex=3)
# }
# dev.off()
# png_files <- sprintf("input%03d.png", 1:1000)
# av::av_encode_video(png_files, 'output.mp4', framerate = 20)
# utils::browseURL('output.mp4')




pair_plots_gg_clusters <- function(mat, labels, colors, main="") {
  dimensions <- combn(ncol(mat), 2)
  all_plots <- list()
  
  labels <- as.character(labels)
  unique_labels <- unique(labels)
  unique_colors <- unique(colors)
  unique_colors_border <- rep(adjustcolor("black", alpha=0.35), times=len(unique_labels))
  unique_colors_border[unique_labels == "-1"] <- adjustcolor("black", alpha=0.1)
  
  for (i in 1:ncol(dimensions)) {
    
    d1 <- dimensions[1,i]
    d2 <- dimensions[2,i]
    
    
    df <- data.frame(X=mat[,d1], Y=mat[,d2], col_lab=labels)
    g <- ggplot(df, aes(x=X, y=Y)) + 
      geom_point(aes(fill=col_lab, color=col_lab), shape=21, size=0.8) +
      xlab(sprintf("Dim %d", d1)) +
      ylab(sprintf("Dim %d", d2)) + 
      theme_light() + 
      scale_fill_manual(breaks=unique_labels, 
                         values=unique_colors) +
      scale_color_manual(breaks=unique_labels,
                         values=unique_colors_border) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="NA")
    
    
    all_plots <- append(all_plots, list(g))
  }
  
  all_plots$nrow <- 5
  
  a <- do.call(arrangeGrob, all_plots)
  
  return(a)
}


my_boxplot <- function(df, xl, yl, xaxis_mod=element_text(), xrange=NA, yrange=NA) {
  
  gbox <- 
    ggplot(df, aes(x=X, y=Y)) + 
    geom_jitter(data=df, aes(x=X, y=Y), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
    geom_boxplot(data=df, width=0.5, size=1, 
                 color=adjustcolor("gray65", alpha=1), 
                 fill=adjustcolor("gray80", alpha=0.8)) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "NA",
          panel.background = element_blank(),
          axis.text.x=xaxis_mod) + 
    xlab(xl) +
    ylab(yl) 
  
  if (all(!is.na(xrange))){
    gbox <- gbox + xlim(xrange)
  }
  
  if (all(!is.na(yrange))){
    gbox <- gbox + ylim(yrange)
  }
  
  return(gbox)
}


single_trial_heamtaps <- function() {
  output_path <- sprintf("%s//single_trial_analysis_shuffle//", output_path_f)
  dir.create(output_path)
  path_list <- get_thirsty_quenched_paths()
  window_size=15
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  for (idx in 1:len(path_list)) {
    
    p <- path_list[[idx]]
    
    mice_path <- sprintf("%s\\%s\\", 
                         output_path,
                         str_replace(datasets_names[idx]," ", "_"))
    dir.create(mice_path)

    #all_mat <- get_mat_with_preset(p, "dalshtaim")
    all_mat <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)
    stim_master_mat <- get_stim_mat(p,just_mat = T, window_size = 15, chunk=-1)    
    behav <- get_behavior_variables(stim_master_mat, mat_size = nrow(all_mat), window_size = 15, seconds_after=1)
    clust_mat_all <- get_clusters_mat_kmeans(all_mat, max_nc = 15)
    ph_cmat <- pheatmap(clust_mat_all$clust_mat, cluster_rows=F, cluster_cols=F, border_col=NA, cex=0.75)
    
    df_cols <- data.frame(mean=scale(clust_mat_all$clusters, center=F),
                          x=1:ncol(clust_mat_all$clust_mat),
                          sd=apply(clust_mat_all$clust_mat, 2, sd))
    
    df_rows <- data.frame(mean=scale(clust_mat_all$sphere_radii, center=F),
                          x=seq(500, nrow(all_mat), length.out=len(clust_mat_all$sphere_radii))/nrow(all_mat),
                          sd=apply(clust_mat_all$clust_mat, 1, sd))
    gcols <- 
      ggplot(df_cols, aes(x=x, y=mean)) + 
      geom_line(size=1) + 
      geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd), fill=adjustcolor("black", alpha=0.2), color=NA) + 
      theme_classic() + 
      xlab("Num of clusters (#)") + 
      ylab("Normalized MSE (AU)") +
      ggtitle("Cluster number selection") +
      geom_vline(xintercept=df_cols[clust_mat_all$nc,"x"], linetype="longdash", size=0.25, color=adjustcolor("black",alpha=0.5)) +
      geom_hline(yintercept=df_cols[clust_mat_all$nc,"mean"], linetype="longdash", size=0.25, color=adjustcolor("black",alpha=0.5))
    
    grows <- 
      ggplot(df_rows, aes(x=x, y=mean)) + 
      geom_line(size=1) + 
      geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd), fill=adjustcolor("black", alpha=0.2), color=NA) + 
      theme_classic() + 
      xlab("Sphere radii (%)") + 
      ylab("Normalized MSE (AU)") +
      ggtitle("Sphere size selection") +
      geom_vline(xintercept=df_rows[clust_mat_all$nr,"x"], linetype="longdash", size=0.25, color=adjustcolor("black",alpha=0.5)) +
      geom_hline(yintercept=df_rows[clust_mat_all$nr,"mean"], linetype="longdash", size=0.25, color=adjustcolor("black",alpha=0.5))

    
    gf <- arrangeGrob(ph_cmat[[4]], gcols, grows, nrow=1)
    png(sprintf("%s\\cluster_map.png",
                mice_path),
        units="in",
        res=500,
        height=4,
        width=12)
    plot(gf)
    dev.off()
    # 

    num_of_clusters <- clust_mat_all$nc
    cluster_labels <- clust_mat_all$labs
     
    p1 <- pair_plots_colored_clusters(all_mat, 
                                       cluster_labels = cluster_labels,
                                       to_label = 1)
    p2 <- pair_plots_colored_clusters(all_mat, 
                                       cluster_labels = cluster_labels,
                                       to_label = 0)

  
    chron_list <- list()
    occ_list <- list()
    prob_list <- list()
    chron_mat_list <- list()
    clust_order_prob_matrix <- list()
    
    
    for (ent in c("EntireTrials3", "EntireTrials4", "EntireTrials5")) {
      
      num_of_trials <- max(behav[,ent]) 
      occupancy_mat <- matrix(rep(0, times=num_of_trials * (num_of_clusters + 1)),
                              nrow=num_of_trials)
      chronological_mat <- matrix(rep(-1, times=num_of_trials * 20),
                                  nrow=num_of_trials)
      
      prob_matrix <- matrix(rep(0, times=(num_of_clusters + 1) ** 2),
                            nrow=(num_of_clusters + 1))
      
      colnames(prob_matrix) <- c(-1, 1:num_of_clusters)
      rownames(prob_matrix) <- c(-1, 1:num_of_clusters)
      
      colnames(occupancy_mat) <- c(-1, 1:num_of_clusters)
      rownames(occupancy_mat) <- 1:num_of_trials
      colnames(chronological_mat) <- 1:20
      rownames(chronological_mat) <- 1:num_of_trials
      
      trials_annot_df <- data.frame(trials=sprintf("Trial %d", 1:num_of_trials))
      clusters_annot_df <- data.frame(occupancy=sprintf("Cluster %d", c(-1,1:num_of_clusters)))
      second_annot_df <- data.frame(clusters=sprintf("Second %d", c(1:20)))
      
      trials_col <- bupu_cg(num_of_trials)
      names(trials_col) <-  sprintf("Trial %d", 1:num_of_trials)
      
      clusters_col <- ylgn_cg(num_of_clusters + 1)
      names(clusters_col) <- sprintf("Cluster %d", c(-1,1:num_of_clusters))
      
      chron_colors <- list(trials=trials_col)
      occ_colors <- list(trials=trials_col, 
                         occupancy=clusters_col)
      
      rownames(clusters_annot_df) <- c(-1, 1:num_of_clusters)
      rownames(trials_annot_df) <- 1:num_of_trials
      rownames(second_annot_df) <- 1:20
      

      
      colnames(chronological_mat) <- sprintf("t =  %.1f", 1:20 / 2)
      
      trials <- 1:max(behav[,ent]) 
      
      
      for (trial in trials) {
        trial_ind <- which(behav[,ent] == trial)
        trial_ind <- (trial_ind[1]):(trial_ind[1] + 19)
      
        clusters_in_trial <- cluster_labels[trial_ind]
        
        num_of_timepoints <- ifelse(len(trial_ind) >= 20,
                                    20,
                                    len(trial_ind))
        
        
        chronological_mat[trial,1:num_of_timepoints] <- clusters_in_trial[1:num_of_timepoints]
        
        tabled <- table(clusters_in_trial)
        
        for (tn in names(tabled)) {
          occupancy_mat[trial, tn] <- tabled[tn]
        }
      }
      
      prob_matrices <- list() 
      clust_annot <- list()
      tt <- substr(ent, unlist(gregexpr("[0-9]{1}", ent)), unlist(gregexpr("[0-9]{1}", ent)) + 1)
      title <- sprintf("%s - Trial Type %s (%d Clusters)", datasets_names[[idx]], tt, num_of_clusters)
      
      
      for (trial_idx in 1:len(trials)) {
        for (clust_idx in 1:(ncol(chronological_mat)  - 1)) {
          first_clust <- as.character(chronological_mat[trial_idx, clust_idx])
          second_clust <- as.character(chronological_mat[trial_idx, clust_idx + 1])
          #print(sprintf("From %s to %s", first_clust, second_clust))
          prob_matrix[first_clust, second_clust] <- prob_matrix[first_clust, second_clust] + 1
          
        }
        
        if (trial_idx %% 25 == 0 || trial_idx == len(trials)) {
          
          if (trial_idx == len(trials)) {
            name = "all"
          } else {
            name = as.character(trial_idx)
          }
          
          
          prb <- prob_matrix / sum(prob_matrix)
          ph_prob <- pheatmap(prb, cluster_rows=T, cluster_cols=T, col=rev(blues_cg(20)), breaks=seq(0,0.05, length.out=20),
                              main=sprintf("%s - Trial Type %s (%s)",datasets_names[[idx]], tt, name))
          

          prob_matrices[[name]] <- ph_prob[[4]]
          
          hc <- hclust(dist(prb), method="ward.D2")
          
          print(sprintf("For trial type %s, - %s", ent, name))
          print(rownames(prb)[hc$order])
          
          clust_annot[[name]]  <- rownames(prb)[hc$order]
        }
      }

      

      


          # ph_occ <- 
          # pheatmap(occupancy_mat[,include_main:(num_of_clusters + 1)], 
          #          cluster_rows=sorted, 
          #          cluster_cols=F, 
          #          annotation_row = trials_annot_df, 
          #          annotation_col = clusters_annot_df, 
          #          annotation_colors = occ_colors, 
          #          show_rownames = F, legend = T, 
          #          main = title,
          #          col=spec_cg(500))
          
          ph_chron <- 
          pheatmap(chronological_mat, 
                   cluster_rows=F, 
                   cluster_cols=F, 
                   annotation_row = trials_annot_df, 
                   annotation_colors = chron_colors, 
                   show_rownames = F, legend = T, 
                   main = title,
                   col=(spec_cg(num_of_clusters + 1)))
          # occ_list <- append(occ_list,
          #                    list(ph_occ[[4]]))
          
          chron_list <- append(chron_list,
                               list(ph_chron[[4]]))
          
          prob_list <- append(prob_list,
                              list(prob_matrices))
          
          clust_order_prob_matrix[[ent]] <- clust_annot
          
          chron_mat_list <- append(chron_mat_list,
                                   list(chronological_mat))
    }
    
    

    chron_list$nrow <- 1
    #chron_list$lay = layout_matrix
    # occ_p  <- do.call(grid.arrange, occ_list)
    chron_p <- do.call(arrangeGrob, chron_list)
    
    struct_p <- arrangeGrob(p1, p2, nrow=1)
    

    fp <- arrangeGrob(chron_p, struct_p, nrow=2)
    # png(sprintf("%s\\%s%sOCCUPANCY_ALL_%d.png",
    #             mice_path,
    #             ifelse(sorted, "sorted_", ""),
    #             ifelse(include_main == 1, "all_", "excluded_"),
    #             num_of_clusters),
    #     units="in",
    #     res=300,
    #     height=8,
    #     width=ifelse(sorted, 4 * 3, 5 * 3))
    #     
    # plot(occ_p)
    # dev.off()
    
    png(sprintf("%s\\trials_heatmap_clusters%d.png",
                mice_path,
                num_of_clusters),
        units="in",
        res=500,
        height=8 * 2 * 1.5,
        width=6 * 3 * 1.5)
    plot(fp)
    dev.off()
    
    for (p in names(prob_list[[1]])) {
      
      pl <- list()
      exit = F
      for (i in 1:len(prob_list)) {
        pl <- append(pl, list(prob_list[[i]][[p]]))
        
        if(is.null(prob_list[[i]][[p]])) {
          exit =T
        }
      }
      if (exit) {
        next
      }
      
      pl$nrow <- 1
      prob_p <- do.call(arrangeGrob, pl)
      png(sprintf("%s\\prob_matrices_%d_trials_%s.png",
                  mice_path,
                  num_of_clusters,
                  p),
          units="in",
          res=500,
          height=3,
          width=9)
      plot(prob_p)
      dev.off()
    }
    

    
    final_list <- list()
    final_list$cluster_mat <- clust_mat_all
    final_list$prob_annot <- clust_order_prob_matrix
    final_list$trial_mats <- chron_mat_list
    
    save(file=sprintf("%s\\structure_annotation_%d.Rda",
                      mice_path,
                      num_of_clusters),
         final_list)
  }
  
}


comp_v1_ins <- function() {
  
  
  v1_cluster_maps <- c()
  por_cluster_maps <- c()
  ins_thirst_cluster_maps <- c() 
  ins_hunger_cluster_maps <- c()
  
  
  corrs <- c()
  
  others <- c(v1_paths, por_paths)
  
  all_paths <- list(v1=v1_paths,
                 por=por_paths,
                 ins_thirst=insula_thirst_paths,
                 ins_hunger=insula_hunger_paths)
  
  cluster_maps <- list()
  for (group in names(all_paths)) {
    group_cluster_maps <- list()
    group_paths <- all_paths[[group]]
    
    for (p1 in group_paths) {
        mat1 <- get_mat_with_preset(p1, "dalshtaim", chunk=-1)
        cluster_mat1 <- get_clusters_mat(mat1, min_frames_for_minc = 5, min_frames_for_outward = 5)
        
        group_cluster_maps <- append(group_cluster_maps,
                                     list(cluster_mat1))
    }
    
    cluster_maps[[group]] <- group_cluster_maps
    
    if (group == "v1") {
      v1_cluster_maps <- group_cluster_maps
    } else if (group == "por") {
      por_cluster_maps <- group_cluster_maps
    } else if (group == "ins_thirst") {
      ins_thirst_cluster_maps <- group_cluster_maps
    } else if (group == "ins_hunger") {
      ins_hunger_cluster_maps <- group_cluster_maps
    }
  }
  
  comparisions <- list()
  comp_matrix <- combn(len(cluster_maps), 2)
  comp_matrix <- cbind(comp_matrix, matrix(rep(1:len(cluster_maps), each=2), nrow=2))
  
  comp_names <- c()
  
  for (idx in 1:ncol(comp_matrix)) {
   
    group_1_idx <- comp_matrix[1,idx]
    group_2_idx <- comp_matrix[2,idx]
    
    comparision <- sprintf("%s vs %s",
                           names(cluster_maps)[group_1_idx],
                           names(cluster_maps)[group_2_idx])
    
    comp_names <- c(comp_names, comparision)
    
    group_1 <- cluster_maps[[group_1_idx]]
    group_2 <- cluster_maps[[group_2_idx]]
    
    corrs <- c()
    
    for (cluster_map_group_1_idx in 1:len(group_1)) {
      for (cluster_map_group_2_idx in 1:len(group_2)) {
        
        if (group_1_idx == group_2_idx && cluster_map_group_1_idx == cluster_map_group_2_idx) {
          print(sprintf("Continuing (%d - %d), (%d - %d)",
                        group_1_idx,
                        group_2_idx,
                        cluster_map_group_1_idx,
                        cluster_map_group_2_idx))
          next
        }
            
        cluster_map_group_1 <-  group_1[[cluster_map_group_1_idx]]
        cluster_map_group_2 <- group_2[[cluster_map_group_2_idx]]
        
        corrs <- c(corrs,
                   euc_dist(c(cluster_map_group_1), c(cluster_map_group_2)))
      } 
    }
    
    print(sprintf("Comparing %s, mean - %f", comparision, mean(corrs)))
    
    comparisions <- append(comparisions,
                           list(corrs))
  }
  
  names(comparisions) <- comp_names
  
  final_df <- c()
  
  for (idx in 1:len(comparisions)) {
    
    corrs <- comparisions[[idx]]
    df <- data.frame(X=rep(comp_names[idx], times=len(corrs)),
                     Y=corrs)
    final_df <- rbind(final_df,
                      df)
  }
  
  my_boxplot(final_df, xl="", yl="Corr") + theme(axis.text.x = element_text(angle = 90))
}


across_mice_decoding <- function(window_size=15, trial_duration=20) {
  
  path_list <- insula_thirst_paths
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  for (trial_group in c("all")) {
  confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                              nrow=len(path_list))
  hit_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  miss_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  
  activity_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                             nrow=len(path_list))
  activity_hit_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  activity_miss_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
  
  for (idx_1 in 1:len(path_list)) {
    for (idx_2 in 1:len(path_list)) {
      if (idx_1 == idx_2) {
        next
      }
      
      p1 <- path_list[idx_1]
      p2 <- path_list[idx_2]
      
      orig_mat_1 <- get_reduced_mat_full_day(p1, window_size=15, just_original_mat = T)
      orig_mat_2 <- get_reduced_mat_full_day(p2, window_size=15, just_original_mat = T)
      
      if (nrow(orig_mat_2) > ncol(orig_mat_2)) {
        orig_mat_2 <- t(orig_mat_2)
      }
      
      if (nrow(orig_mat_1) > ncol(orig_mat_1)) {
        orig_mat_1 <- t(orig_mat_1)
      }
      mean_activity_1 <- colMeans(orig_mat_1)
      mean_activity_2 <- colMeans(orig_mat_2)
      
      stim_master_1 <- get_stim_mat(p1,just_mat = T, window_size = window_size, chunk=-1)
      stim_master_2 <- get_stim_mat(p2,just_mat = T, window_size = window_size, chunk=-1)
      
      annotated_struct_path_1 <-  sprintf("%s\\%s\\structure_annotation_%d.Rda", 
                                          output_path,
                                          str_replace(datasets_names[idx_1]," ", "_"),
                                          7)
      
      load(annotated_struct_path_1)
      annotated_mice_1 <- final_list
      
      annotated_struct_path_2 <-  sprintf("%s\\%s\\structure_annotation_%d.Rda", 
                                          output_path,
                                          str_replace(datasets_names[idx_2]," ", "_"),
                                          7)
      
      load(annotated_struct_path_2)
      annotated_mice_2 <- final_list
      
      labels_mice_1 <- annotated_mice_1$cluster_mat$labs
      labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
      mask_mice_1 <- annotated_mice_1$prob_annot$EntireTrials3[[trial_group]]
      mask_mice_2 <- annotated_mice_2$prob_annot$EntireTrials3[[trial_group]]
      
      reward_trials_mice_1 <- which(stim_master_1[,7] == 3)
      reward_trials_labels_mice_1 <- stim_master_1[which(stim_master_1[,7] == 3),8]
      reward_trials_mat_mice_1 <- c()
      activity_mat_trials_mice_1  <- c()
      
      for (r_trial in reward_trials_mice_1) {
        trial_binned_index <- as.numeric(get_binned_index(stim_master_1[r_trial,1], window_size))
        
        clustered_trial <- labels_mice_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
        
        if (sum(is.na(clustered_trial)) > 0 ) {
          clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
        }
        
        reward_trials_mat_mice_1 <- rbind(reward_trials_mat_mice_1,
                                          clustered_trial)
        
       
        activity_in_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
        activity_in_trial[is.na(activity_in_trial)] <- 0
        
        activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_in_trial)
        
        
      }
      
      annot_df_mice_1 <- data.frame(Result=reward_trials_labels_mice_1)
      rownames(reward_trials_mat_mice_1) <- 1:nrow(reward_trials_mat_mice_1)
      rownames(annot_df_mice_1) <- 1:nrow(reward_trials_mat_mice_1)
      
      annot_df_mice_1[annot_df_mice_1$Result == 0,"Result"] <- "Hit"
      annot_df_mice_1[annot_df_mice_1$Result == 1,"Result"] <- "Miss"
      
      
      reward_trials_mice_2 <- which(stim_master_2[,7] == 3)
      reward_trials_labels_mice_2 <- stim_master_2[which(stim_master_2[,7] == 3),8]
      reward_trials_mat_mice_2 <- c()
      activity_mat_trials_mice_2  <- c()
      
      for (r_trial in reward_trials_mice_2) {
        trial_binned_index <- as.numeric(get_binned_index(stim_master_2[r_trial,1], window_size))
        
        
        clustered_trial <- labels_mice_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
        
        if (sum(is.na(clustered_trial)) > 0 ) {
          clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
        }
        
        reward_trials_mat_mice_2 <- rbind(reward_trials_mat_mice_2,
                                          clustered_trial)
        
        
        activity_in_trial <- mean_activity_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
        activity_in_trial[is.na(activity_in_trial)] <- 0
        activity_mat_trials_mice_2 <- rbind(activity_mat_trials_mice_2, activity_in_trial)
      }
      
      annot_df_mice_2 <- data.frame(Result=reward_trials_labels_mice_2)
      rownames(reward_trials_mat_mice_2) <- 1:nrow(reward_trials_mat_mice_2)
      rownames(annot_df_mice_2) <- 1:nrow(reward_trials_mat_mice_2)
      
      annot_df_mice_2[annot_df_mice_2$Result == 0,"Result"] <- "Hit"
      annot_df_mice_2[annot_df_mice_2$Result == 1,"Result"] <- "Miss"
    
    
    
    translated_mat_mice_1_mice_2 <- c()
    
    
    for (trial_idx in 1:nrow(reward_trials_mat_mice_1)) {
      trial_before <- reward_trials_mat_mice_1[trial_idx,]
      trial_after <- trial_before
      
      for (cluster_idx in 1:len(mask_mice_1)) {
        trial_after[which(trial_after == mask_mice_1[cluster_idx])] <- mask_mice_2[cluster_idx] 
      }
      
      translated_mat_mice_1_mice_2 <- rbind(translated_mat_mice_1_mice_2,
                                            as.numeric(trial_after))
      
    }
    
    decoded_vec <- c()
    
    for (decoded_trial_idx in 1:nrow(reward_trials_mat_mice_2)) {
      decoded_trial <- reward_trials_mat_mice_2[decoded_trial_idx,]
      confidence <- 
        apply(translated_mat_mice_1_mice_2, 1, 
              function(trial) {euc_dist(trial, decoded_trial)})  
      
      
      is_hit <- mean(confidence[reward_trials_labels_mice_1 == 0]) 
      is_miss <- mean(confidence[reward_trials_labels_mice_1 == 1])
      decoded_vec <- c(decoded_vec,
                       ifelse(is_hit < is_miss, 0, 1))
    }
    
    
    decoding_statistics <-  sum(decoded_vec == reward_trials_labels_mice_2) / len(reward_trials_labels_mice_2)
    hit_decoding_statistics <- sum(decoded_vec[reward_trials_labels_mice_2 == 0] == 0) / sum(reward_trials_labels_mice_2 == 0)
    miss_decoding_statistics <- sum(decoded_vec[reward_trials_labels_mice_2 == 1] == 1) / sum(reward_trials_labels_mice_2 == 1)
    
    print(sprintf("Decoding statistics from %d to %d, %f (hits only: %f) (misses only :%f)", 
                  idx_1,
                  idx_2,
                  decoding_statistics,
                  hit_decoding_statistics,
                  miss_decoding_statistics))
    
    confusion_matrix[idx_1, idx_2] <- decoding_statistics
    miss_confusion_matrix[idx_1, idx_2] <- miss_decoding_statistics
    hit_confusion_matrix[idx_1, idx_2] <- hit_decoding_statistics
    
  
    activity_decoded_vec <- c()
    
    
    for (decoded_trial_idx in 1:nrow(activity_mat_trials_mice_2)) {
      decoded_trial <- activity_mat_trials_mice_2[decoded_trial_idx,]
      confidence <- 
        apply(activity_mat_trials_mice_1, 1, 
              function(trial) {cor(trial, decoded_trial)})  
      
      
      is_hit <- mean(confidence[reward_trials_labels_mice_1 == 0]) 
      is_miss <- mean(confidence[reward_trials_labels_mice_1 == 1])
      activity_decoded_vec <- c(activity_decoded_vec, ifelse(is_hit > is_miss, 0, 1))
    }
    
    activity_decoding_statistics <-  sum(activity_decoded_vec == reward_trials_labels_mice_2) / len(reward_trials_labels_mice_2)
    activity_hit_decoding_statistics <- sum(activity_decoded_vec[reward_trials_labels_mice_2 == 0] == 0) / sum(reward_trials_labels_mice_2 == 0)
    activity_miss_decoding_statistics <- sum(activity_decoded_vec[reward_trials_labels_mice_2 == 1] == 1) / sum(reward_trials_labels_mice_2 == 1)
    
    print(sprintf("Activity decoding statistics from %d to %d, %f (hits only: %f) (misses only :%f)", 
                  idx_1,
                  idx_2,
                  activity_decoding_statistics,
                  activity_hit_decoding_statistics,
                  activity_miss_decoding_statistics))
    
    activity_confusion_matrix[idx_1, idx_2] <- activity_decoding_statistics
    activity_miss_confusion_matrix[idx_1, idx_2] <- activity_miss_decoding_statistics
    activity_hit_confusion_matrix[idx_1, idx_2] <- activity_hit_decoding_statistics
    
    }
  }
  
  rownames(confusion_matrix) <- datasets_names
  colnames(confusion_matrix) <- datasets_names
  rownames(hit_confusion_matrix) <- datasets_names
  colnames(hit_confusion_matrix) <- datasets_names
  rownames(miss_confusion_matrix) <- datasets_names
  colnames(miss_confusion_matrix) <- datasets_names
  
  ph_confusion <- pheatmap(confusion_matrix, cluster_rows=F, cluster_cols=F, 
                           main="Hits & Misses", 
                           breaks=seq(0, 1, length.out=100),
                           col=turbo(100),
                           border_col=NA)
  ph_hits <- pheatmap(hit_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Hits",
                      breaks=seq(0, 1, length.out=100),
                      col=turbo(100),
                      border_col=NA)
  ph_miss <- pheatmap(miss_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Misses",
                      breaks=seq(0, 1, length.out=100),
                      col=turbo(100),
                      border_col=NA)
  
  
  boxdf <- data.frame()
  mt_list <- list(All=confusion_matrix,
                  Hits=hit_confusion_matrix,
                  Misses=miss_confusion_matrix)
  for (mt_name in names(mt_list)) {
    mt <- mt_list[[mt_name]]
    
    for (i in 1:nrow(mt)) {
      mt[i,i] <- NA
    }
    
    values <- c(mt)
    values <- values[!is.na(values)]
    
    boxdf <- rbind(boxdf,
                  data.frame(X=rep(mt_name, times=len(values)),
                             Y=values))
  }
  
  mboxplot <- my_boxplot(boxdf, xl="Group", yl="Classification (%)") + geom_hline(yintercept=0.5, linetype="longdash")
  
  gf <- grid.arrange(ph_confusion[[4]], ph_hits[[4]], ph_miss[[4]], mboxplot, top=sprintf("Trials used - %s", trial_group), nrow=1)
  
  output_path_2 <- sprintf("%s//across_mice_decoding//", output_path_f)
  dir.create(output_path_2)
  pdf(sprintf("%s//trials_%s.pdf", output_path_2, trial_group),
      height=3, width=12)
  
  plot(gf)
  dev.off()
  
  
  rownames(activity_confusion_matrix) <- datasets_names
  colnames(activity_confusion_matrix) <- datasets_names
  rownames(activity_hit_confusion_matrix) <- datasets_names
  colnames(activity_hit_confusion_matrix) <- datasets_names
  rownames(activity_miss_confusion_matrix) <- datasets_names
  colnames(activity_miss_confusion_matrix) <- datasets_names
  
  ph_confusion <- pheatmap(activity_confusion_matrix, cluster_rows=F, cluster_cols=F, 
                           main="Hits & Misses", 
                           breaks=seq(0, 1, length.out=100),
                           col=turbo(100),
                           border_col=NA)
  ph_hits <- pheatmap(activity_hit_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Hits",
                      breaks=seq(0, 1, length.out=100),
                      col=turbo(100),
                      border_col=NA)
  ph_miss <- pheatmap(activity_miss_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Misses",
                      breaks=seq(0, 1, length.out=100),
                      col=turbo(100),
                      border_col=NA)
  
  
  boxdf <- data.frame()
  mt_list <- list(All=activity_confusion_matrix,
                  Hits=activity_hit_confusion_matrix,
                  Misses=activity_miss_confusion_matrix)
  for (mt_name in names(mt_list)) {
    mt <- mt_list[[mt_name]]
    
    for (i in 1:nrow(mt)) {
      mt[i,i] <- NA
    }
    
    values <- c(mt)
    values <- values[!is.na(values)]
    
    boxdf <- rbind(boxdf,
                   data.frame(X=rep(mt_name, times=len(values)),
                              Y=values))
  }
  
  mboxplot <- my_boxplot(boxdf, xl="Group", yl="Classification (%)") + geom_hline(yintercept=0.5, linetype="longdash")
  
  gf <- grid.arrange(ph_confusion[[4]], ph_hits[[4]], ph_miss[[4]], mboxplot, top=sprintf("Trials used - %s", trial_group), nrow=1)
  
  output_path_2 <- sprintf("%s//across_mice_decoding//", output_path_f)
  dir.create(output_path_2)
  pdf(sprintf("%s//activity_decoding.pdf", output_path_2),
      height=3, width=12)
  
  plot(gf)
  dev.off()
  
  acm <- activity_confusion_matrix
  cm <- confusion_matrix
  for (i in 1:nrow(cm)) {acm[i,i] <- NA; cm[i,i] <- NA}
  cm <- c(cm)
  acm <- c(acm)
  cm <- cm[!is.na(cm)]
  acm <- acm[!is.na(acm)]
  decoding_df <- data.frame(values=c(cm,acm), groups=c(rep("By sequence", times=len(cm)), rep("By activity levels", times=len(acm))))
  g <- ggplot(data=decoding_df, aes(x=groups, y=values)) + geom_point()
  
  
  g2 <- g
  
  for (i in 1:len(cm)) {
    tmp_df <- data.frame(y=c(cm[i], acm[i]),
                         x=c(2, 1))
    
    g2 <- 
      g2 + geom_line(data=tmp_df, aes(x=x, y=y), color=adjustcolor("gray70", alpha=0.5))  
    
  }
  
  g2 <- g2 + theme_classic() + xlab("") + ylab("Classifcation (%)") + ylim(c(0,1))
  g3 <- g2 + geom_point(data=decoding_df, aes(x=groups, y=values)) + 
             geom_hline(yintercept=0.5, linetype="longdash") +
             geom_line(data=data.frame(x=c(1.25, 1.25, 1.75, 1.75),
                                       y=c(0.93, 0.95, 0.95, 0.93)),
                       aes(x=x,y=y)) +
              geom_text(x=1.5, y=0.98, label=sprintf("%s", 
                                                     get_signif(t.test(cm, acm, alternative="greater")$p.value)))
  
  }
  
}


across_mice_decoding_all_trials_build_prob_mat <- function() {
  output_path <- sprintf("%s//across_mice_decoding_all_trials//", output_path_f)
  dir.create(output_path)
  path_list <- get_thirsty_quenched_paths()
  window_size=15
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  for (idx in 6:len(path_list)) {
    
    p <- path_list[[idx]]
    
    mice_path <- sprintf("%s\\%s\\", 
                         output_path,
                         str_replace(datasets_names[idx]," ", "_"))
    dir.create(mice_path)
    
    all_mat <- get_mat_with_preset(p, "dalshtaim")
    stim_master_mat <- get_stim_mat(p,just_mat = T, window_size = 15, chunk=-1)    
    behav <- get_behavior_variables(stim_master_mat, mat_size = nrow(all_mat), window_size = 15, seconds_after=1)
    clust_mat_all <- get_clusters_mat_kmeans(all_mat)

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
      
    title <- sprintf("%s - %d Clusters", datasets_names[[idx]], num_of_clusters)
      
      
    for (trial_idx in 1:len(trials)) {
      for (clust_idx in 1:(ncol(chronological_mat)  - 1)) {
        first_clust <- as.character(chronological_mat[trial_idx, clust_idx])
        second_clust <- as.character(chronological_mat[trial_idx, clust_idx + 1])
        #print(sprintf("From %s to %s", first_clust, second_clust))
        prob_matrix[first_clust, second_clust] <- prob_matrix[first_clust, second_clust] + 1
          
      }
    }
      
    prb <- prob_matrix / sum(prob_matrix)
    ph_prob <- pheatmap(prb, cluster_rows=T, cluster_cols=T, col=rev(blues_cg(20)), breaks=seq(0,0.05, length.out=20),
                        main=sprintf("%s - All trials",datasets_names[[idx]]))
    
    
    prob_matrices[[name]] <- ph_prob[[4]]
    
    hc <- hclust(dist(prb), method="ward.D2")
    
    print(sprintf("For trial type %s, - %s", ent, name))
    print(rownames(prb)[hc$order])
    
    clust_annot[[name]]  <- rownames(prb)[hc$order]
      
    annot_df <- data.frame(trialtype=relevant_trials_mat[,"TrialType"], response=relevant_trials_mat[,"Response"] %% 2)
    rownames(annot_df) <- 1:num_of_trials
      

      
    ph_chron <- 
        pheatmap(chronological_mat, 
                 cluster_rows=F, 
                 cluster_cols=F, 
                 annotation_row = annot_df, 
                 show_rownames = F, legend = T, 
                 main = title,
                 col=(spec_cg(num_of_clusters + 1)))

    
    
  
    final_plot <- arrangeGrob(ph_chron[[4]], ph_prob[[4]], nrow=1)

    png(sprintf("%s\\heatmaps_%d.png",
                mice_path,
                num_of_clusters),
        units="in",
        res=500,
        height=4,
        width=8)
    plot(final_plot)
    dev.off()
    
    
    
    
    final_list <- list()
    final_list$cluster_mat <- clust_mat_all
    final_list$mask <- rownames(prb)[hc$order]
    final_list$trials_mat <- chronological_mat
    final_list$annot_df <- annot_df
    final_list$probability_mat <- prb
    
    save(file=sprintf("%s\\structure_annotation_all_trials_%d.Rda",
                      mice_path,
                      num_of_clusters),
         final_list)
  }
  
}


across_mice_decoding_all_trials_decode <- function()  {
  #output_path <- sprintf("%s//single_trial_analysis//", output_path_f)
  output_path <- sprintf("%s//across_mice_decoding_all_trials//", output_path_f)
  
  path_list <- insula_thirst_paths
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  for (trial_group in c("all")) {
    confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                               nrow=len(path_list))
    shuffle_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                                          nrow=len(path_list))
    
    cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                               nrow=len(path_list))
    shuffle_cosine_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                                       nrow=len(path_list))
    hit_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
    miss_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
    
    activity_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2),
                                        nrow=len(path_list))
    activity_hit_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
    activity_miss_confusion_matrix <- matrix(rep(1, times=len(path_list) ** 2), nrow=len(path_list))
    
    for (idx_1 in 1:len(path_list)) {
      for (idx_2 in 1:len(path_list)) {
        if (idx_1 == idx_2) {
          next
        }
        
        p1 <- path_list[idx_1]
        p2 <- path_list[idx_2]
        
        orig_mat_1 <- get_reduced_mat_full_day(p1, window_size=15, just_original_mat = T)
        orig_mat_2 <- get_reduced_mat_full_day(p2, window_size=15, just_original_mat = T)
        
        if (nrow(orig_mat_2) > ncol(orig_mat_2)) {
          orig_mat_2 <- t(orig_mat_2)
        }
        
        if (nrow(orig_mat_1) > ncol(orig_mat_1)) {
          orig_mat_1 <- t(orig_mat_1)
        }
        mean_activity_1 <- colMeans(orig_mat_1)
        mean_activity_2 <- colMeans(orig_mat_2)
        
        stim_master_1 <- get_stim_mat(p1,just_mat = T, window_size = window_size, chunk=-1)
        stim_master_2 <- get_stim_mat(p2,just_mat = T, window_size = window_size, chunk=-1)
        
        annotated_struct_path_1 <-  sprintf("%s\\%s\\structure_annotation_all_trials_%d.Rda", 
                                            output_path,
                                            str_replace(datasets_names[idx_1]," ", "_"),
                                            7)
        
        load(annotated_struct_path_1)
        annotated_mice_1 <- final_list
        
        annotated_struct_path_2 <-  sprintf("%s\\%s\\structure_annotation_all_trials_%d.Rda", 
                                            output_path,
                                            str_replace(datasets_names[idx_2]," ", "_"),
                                            7)
        
        load(annotated_struct_path_2)
        annotated_mice_2 <- final_list
        
        labels_mice_1 <- annotated_mice_1$cluster_mat$labs
        labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
      
        mask_mice_1 <- annotated_mice_1$prob_annot$EntireTrials3$all
        mask_mice_2 <- annotated_mice_2$prob_annot$EntireTrials3$all
        
        decoding_stats <- c()
        
        for (shuff_i in c(-1:200)) {
        
        shuffled_stim_mat <- stim_master_1
        shuffled_stim_mat[,"Response"]
        if (shuff_i != -1) {
          #print("NOT shuff")
          shuffled_stim_mat[,1] <- shuffle(shuffled_stim_mat[,1])
        }
        
        trials_mice_1 <- which(shuffled_stim_mat[,7] %in% c(3))
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
          
          
          activity_in_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          activity_in_trial[is.na(activity_in_trial)] <- 0
          
          activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_in_trial)
          
          
        }
        
        #annot_df_mice_1 <- data.frame(Result=paste(annotated_mice_1$annot_df[,1], annotated_mice_1$annot_df[,2]))
        rownames(trials_mat_mice_1) <- 1:nrow(trials_mat_mice_1)
        #rownames(annot_df_mice_1) <- 1:nrow(trials_mat_mice_1)

        
        
        trials_mice_2 <- which(stim_master_2[,7] %in% c(3))
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
          
          
          activity_in_trial <- mean_activity_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          activity_in_trial[is.na(activity_in_trial)] <- 0
          activity_mat_trials_mice_2 <- rbind(activity_mat_trials_mice_2, activity_in_trial)
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
        
        
        decoded_vec <- c()
        cosine_decoded_vec <- c()
        
        for (decoded_trial_idx in 1:nrow(trials_mat_mice_2)) {
          decoded_trial <- trials_mat_mice_2[decoded_trial_idx,]
          
          cosine_confidence <-   
              apply(translated_mat_mice_1_mice_2, 1, 
                    function(trial) {cos_dist_vec(trial, decoded_trial)})  
          confidence <- 
            apply(translated_mat_mice_1_mice_2, 1, 
                  function(trial) {euc_dist(trial, decoded_trial)})  
          
          
          conf_vec <- c()
          cosine_conf_vec <- c()
          trial_results <- sort(unique(trials_labels_mice_1))
          for (res in trial_results) {
            conf_vec <- c(conf_vec,
                          mean(confidence[trials_labels_mice_1 == res]))
            cosine_conf_vec <- c(cosine_conf_vec,
                                 mean(cosine_confidence[trials_labels_mice_1 == res]))
            
          }
          
          
          # is_hit <- mean(confidence[reward_trials_labels_mice_1 == 0]) 
          # is_miss <- mean(confidence[reward_trials_labels_mice_1 == 1])
          cosine_decoded_vec <- c(cosine_decoded_vec, trial_results[which.max(cosine_conf_vec)])
          decoded_vec <- c(decoded_vec, trial_results[which.min(conf_vec)])
        }
        
        
        decoding_statistics <-  sum(decoded_vec == trials_labels_mice_2) / len(trials_labels_mice_2)
        cosine_decoding_statistics <-  sum(cosine_decoded_vec == trials_labels_mice_2) / len(trials_labels_mice_2)
        # hit_decoding_statistics <- sum(decoded_vec[reward_trials_labels_mice_2 == 0] == 0) / sum(reward_trials_labels_mice_2 == 0)
        # miss_decoding_statistics <- sum(decoded_vec[reward_trials_labels_mice_2 == 1] == 1) / sum(reward_trials_labels_mice_2 == 1)
        
        print(sprintf("Decoding statistics from %d to %d, %f (cosine: %f)", 
                      idx_1,
                      idx_2,
                      decoding_statistics,
                      cosine_decoding_statistics))
        
        decoding_stats <- c(decoding_stats,
                           decoding_statistics)
        
        cosine_decoding_stats <- c(decoding_stats,
                                   cosine_decoding_statistics)
        }
        
        
        #print(table(trials_labels_mice_2[decoded_vec == trials_labels_mice_2]))
        
        decoding_statistics <- decoding_stats[1]
        cosine_decoding_statistics <- cosine_decoding_statistics[1]
        shuffled_decoding_statistics <- mean(decoding_stats[-1])
        shuffled_cosine <- mean(cosine_decoding_statistics[-1])
        
        print("######################## SHUFFLED DECODING of hits: ####")
        print(sprintf("%.3f vs %.3f", shuffled_decoding_statistics, decoding_statistics))
        print(sprintf("%.3f vs %.3f", shuffled_cosine, cosine_decoding_statistics))
        
        confusion_matrix[idx_1, idx_2] <- decoding_statistics
        shuffle_confusion_matrix[idx_1, idx_2] <- shuffled_decoding_statistics
        cosine_confusion_matrix[idx_1, idx_2] <- cosine_decoding_statistics
        shuffle_cosine_confusion_matrix[idx_1, idx_2] <- shuffled_cosine
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
    
    rownames(confusion_matrix) <- datasets_names
    colnames(confusion_matrix) <- datasets_names
    rownames(hit_confusion_matrix) <- datasets_names
    colnames(hit_confusion_matrix) <- datasets_names
    rownames(miss_confusion_matrix) <- datasets_names
    colnames(miss_confusion_matrix) <- datasets_names
    
    ph_confusion <- pheatmap(confusion_matrix, cluster_rows=F, cluster_cols=F, 
                             main="Hits & Misses", 
                             breaks=seq(0, 1, length.out=100),
                             col=turbo(100),
                             border_col=NA)
    ph_hits <- pheatmap(hit_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Hits",
                        breaks=seq(0, 1, length.out=100),
                        col=turbo(100),
                        border_col=NA)
    ph_miss <- pheatmap(miss_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Misses",
                        breaks=seq(0, 1, length.out=100),
                        col=turbo(100),
                        border_col=NA)
    
    
    boxdf <- data.frame()
    mt_list <- list(All=confusion_matrix,
                    Hits=hit_confusion_matrix,
                    Misses=miss_confusion_matrix)
    for (mt_name in names(mt_list)) {
      mt <- mt_list[[mt_name]]
      
      for (i in 1:nrow(mt)) {
        mt[i,i] <- NA
      }
      
      values <- c(mt)
      values <- values[!is.na(values)]
      
      boxdf <- rbind(boxdf,
                     data.frame(X=rep(mt_name, times=len(values)),
                                Y=values))
    }
    
    mboxplot <- my_boxplot(boxdf, xl="Group", yl="Classification (%)") + geom_hline(yintercept=0.5, linetype="longdash")
    
    gf <- grid.arrange(ph_confusion[[4]], ph_hits[[4]], ph_miss[[4]], mboxplot, top=sprintf("Trials used - %s", trial_group), nrow=1)
    
    output_path_2 <- sprintf("%s//across_mice_decoding//", output_path_f)
    dir.create(output_path_2)
    pdf(sprintf("%s//trials_%s.pdf", output_path_2, trial_group),
        height=3, width=12)
    
    plot(gf)
    dev.off()
    
    
    rownames(activity_confusion_matrix) <- datasets_names
    colnames(activity_confusion_matrix) <- datasets_names
    rownames(activity_hit_confusion_matrix) <- datasets_names
    colnames(activity_hit_confusion_matrix) <- datasets_names
    rownames(activity_miss_confusion_matrix) <- datasets_names
    colnames(activity_miss_confusion_matrix) <- datasets_names
    
    ph_confusion <- pheatmap(activity_confusion_matrix, cluster_rows=F, cluster_cols=F, 
                             main="Hits & Misses", 
                             breaks=seq(0, 1, length.out=100),
                             col=turbo(100),
                             border_col=NA)
    ph_hits <- pheatmap(activity_hit_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Hits",
                        breaks=seq(0, 1, length.out=100),
                        col=turbo(100),
                        border_col=NA)
    ph_miss <- pheatmap(activity_miss_confusion_matrix, cluster_rows=F, cluster_cols=F, main="Misses",
                        breaks=seq(0, 1, length.out=100),
                        col=turbo(100),
                        border_col=NA)
    
    
    boxdf <- data.frame()
    mt_list <- list(All=activity_confusion_matrix,
                    Hits=activity_hit_confusion_matrix,
                    Misses=activity_miss_confusion_matrix)
    for (mt_name in names(mt_list)) {
      mt <- mt_list[[mt_name]]
      
      for (i in 1:nrow(mt)) {
        mt[i,i] <- NA
      }
      
      values <- c(mt)
      values <- values[!is.na(values)]
      
      boxdf <- rbind(boxdf,
                     data.frame(X=rep(mt_name, times=len(values)),
                                Y=values))
    }
    
    mboxplot <- my_boxplot(boxdf, xl="Group", yl="Classification (%)") + geom_hline(yintercept=0.5, linetype="longdash")
    
    gf <- grid.arrange(ph_confusion[[4]], ph_hits[[4]], ph_miss[[4]], mboxplot, top=sprintf("Trials used - %s", trial_group), nrow=1)
    
    output_path_2 <- sprintf("%s//across_mice_decoding//", output_path_f)
    dir.create(output_path_2)
    pdf(sprintf("%s//activity_decoding.pdf", output_path_2),
        height=3, width=12)
    
    plot(gf)
    dev.off()
    
    acm <- activity_confusion_matrix
    cm <- confusion_matrix
    for (i in 1:nrow(cm)) {acm[i,i] <- NA; cm[i,i] <- NA}
    cm <- c(cm)
    acm <- c(acm)
    cm <- cm[!is.na(cm)]
    acm <- acm[!is.na(acm)]
    decoding_df <- data.frame(values=c(cm,acm), groups=c(rep("By sequence", times=len(cm)), rep("By activity levels", times=len(acm))))
    g <- ggplot(data=decoding_df, aes(x=groups, y=values)) + geom_point()
    
    
    g2 <- g
    
    for (i in 1:len(cm)) {
      tmp_df <- data.frame(y=c(cm[i], acm[i]),
                           x=c(2, 1))
      
      g2 <- 
        g2 + geom_line(data=tmp_df, aes(x=x, y=y), color=adjustcolor("gray70", alpha=0.5))  
      
    }
    
    g2 <- g2 + theme_classic() + xlab("") + ylab("Classifcation (%)") + ylim(c(0,1))
    g3 <- g2 + geom_point(data=decoding_df, aes(x=groups, y=values)) + 
      geom_hline(yintercept=0.5, linetype="longdash") +
      geom_line(data=data.frame(x=c(1.25, 1.25, 1.75, 1.75),
                                y=c(0.93, 0.95, 0.95, 0.93)),
                aes(x=x,y=y)) +
      geom_text(x=1.5, y=0.98, label=sprintf("%s", 
                                             get_signif(t.test(cm, acm, alternative="greater")$p.value)))
    
  }
  
}



get_all_cluster_numbers <- function()
{
  window_size=15
  mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = " ")
  
  res_df <- data.frame()
  
  params_mat <- as.data.frame(cross_df(list(time_shuffled = c(F),
                                            shuffled = c(T, F),
                                            ndim = c(20),
                                            max_nc = c(20))))
  
  
  for (param_idx in 1:nrow(params_mat)) {
    
    prms <- params_mat[param_idx,]
    
    if (prms[["shuffled"]] == F && prms[["time_shuffled"]] == "T") {
      next
    }
    
    for (idx in 1:len(path_list)) {
      all_mat <- get_reduced_mat_full_day(p, 
                                          "lem", 
                                          ndim=prms[["ndim"]], 
                                          window_size=15, 
                                          knn1=0.075, knn2 = 0, 
                                          shuffled=prms[["shuffled"]], 
                                          time_shuffled=prms[["time_shuffled"]])
      
      clust_mat_all <- get_clusters_mat_kmeans(all_mat, 
                                               max_nc = prms[["max_nc"]],
                                               new_method = T)
      
      
      
      result_of_params <- cbind(cbind(cbind(as.data.frame(prms), c(dataset=datasets_names[[idx]])), c(nc=clust_mat_all$nc)), c(nr=clust_mat_all$nr))
      colnames(result_of_params) <- c(colnames(prms), "dataset", "nc", "nr")
      
      res_df <- rbind(res_df,
                      result_of_params)
      
    }
  }
}

