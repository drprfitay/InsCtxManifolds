library(reshape2)
library(TDAstats)


figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"

figure_2_dimension_plots <- function() 
{
  
  write_path <- sprintf("%s\\figure_2\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\dimension_plots", write_path)
  dir.create(write_path)  
  
  
  sizes = list(big=c(width=2.5),
               medium=c(width=2),
               medium_2=c(width=1.5),
               medium_1=c(width=1.25),
               small=c(width=1))
  
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  dim_df <- data.frame()
  for (path_idx in 1:len(paths)) {
    
    p <- paths[[path_idx]]
    original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, just_original_mat = T)
    reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, time_shuffled=F)
    shuffled_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)
    #reg_id <- get_intrinsic_dim(reg_mat)
    
    reg_id <- twonn(reg_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
    shuffled_id <- twonn(shuffled_mat, c_trimmed = 0.05, method="linfit")
    
    dim_df <- rbind(dim_df,
                    data.frame(name=datasets_names[path_idx],
                               n_neurons=min(dim(original_mat)),
                               n_timepoints=max(dim(original_mat)),
                               reg_dim=reg_id$est[[2]],
                               shuffled_dim=shuffled_id$est[[2]]))
  }
  
  
  
  tmp_df <- dim_df[,c("reg_dim","shuffled_dim")]
  colnames(tmp_df) <- c("Regular", "Shuffle")
  
  melted_df <- melt(tmp_df)
  colnames(melted_df) <- c("x", "y")
  
  gbar <- 
  ggplot(melted_df, aes(x=x,y=y)) +
    geom_bar(stat="summary", width=.45,
             fill="gray70") +
    geom_errorbar(stat="summary", width = .3) +
    geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                size=1.5) + 
    theme_classic() +
    xlab("") +
    ylab("Estimated dimension") +
    theme(text=element_text(size=14, color="black"),
          legend.position="NA")  + 
    ylim(0,30) +
    geom_line(data=data.frame(x=c(1,1,2,2),
                              y=c(24.5,25,25,24.5)),
              aes(x=x,y=y)) +
    
    geom_point(data=data.frame(x=1.5,y=26.5),
               color="white",fill="white",stroke=0) + 
    geom_text(label=signif.num(wilcox.test(dim_df[,3], dim_df[,4], alternative="less", paired=T)$p.value),
              x=1.5,
              y=25.5) +
    scale_color_manual(values=c("#405EAB", "#CE3736")) +
    scale_y_continuous(expand=c(0,0)) 
    
    tmp_df_2 <- dim_df[,2:4]
    colnames(tmp_df_2) <- c("Neurons", "Regular", "Shuffle")
    
    melted_df_2 <- melt(tmp_df_2, id.vars = "Neurons")
    colnames(melted_df_2)[2:3] <- c("group", "Dim")
    
    gdim <- 
    ggplot(melted_df_2) +
      geom_point(aes(y=Dim, x=Neurons, color=group), alpha=.75,
                 size=3.5,
                 stroke=0,
                 fill=NA) + 
      theme_classic() +
      xlab("Number of neurons") +
      ylab("Estimated dimension") +
      theme(text=element_text(size=14, color="black"),
            legend.position="NA")  + 
      ylim(0,30)  +
      xlim(150,430) +
      scale_color_manual(values=c("#405EAB", "#CE3736"))
    
    
    gf <- arrangeGrob(gbar, gdim, nrow=1, widths = c(1.5,2))
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\%s",
                         write_path,
                         size_name))
      pdf(sprintf("%s\\%s\\dim_estimation_all.pdf",
                  write_path,
                  size_name),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]] * 1.75) 
      
      plot(gf)
      dev.off()
      
      
    }
}

figure_2_structure_similarity_plots <- function() 
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
  

  nclusters=100
  chunk=-1
  nreps=10
  #compare_shuffle = F
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
    plot_df <- data.frame(Similarity=c(rowMeans(results_all$Regular[[metric_name]]),
                                       rowMeans(results_all$Shuffle[[metric_name]])),
                          Group=c(regular_labels,
                                  shuffle_labels))
    
    ylim_max = max(plot_df$Similarity)
    ylim_min = min(plot_df$Similarity)
    
    wilx = wilcox.test(rowMeans(results_all$Regular[[metric_name]]),
                       rowMeans(results_all$Shuffle[[metric_name]]),
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
    
    plot_df_wthin_across <- 
      data.frame(Similarity=c(rowMeans(results_all$Regular[[metric_name]]),
                              rowMeans(results_all$Shuffle[[metric_name]])),
                 Group=factor(c(regular_within_across_labels, shuffle_within_across_labels),
                              levels=c(unique(regular_within_across_labels)[2],
                                       unique(shuffle_within_across_labels)[2],
                                       unique(regular_within_across_labels)[1],
                                       unique(shuffle_within_across_labels)[1])))
    
    ylim_max = max(plot_df$Similarity)
    
    wilx_within = wilcox.test(rowMeans(results_all$Regular[[metric_name]])[within_across_indices],
                              rowMeans(results_all$Shuffle[[metric_name]])[within_across_indices],
                              paired=T,
                              alternative=alternatives[[metric_name]])
    
    wilx_across = wilcox.test(rowMeans(results_all$Regular[[metric_name]])[!within_across_indices],
                              rowMeans(results_all$Shuffle[[metric_name]])[!within_across_indices],
                              paired=T,
                              alternative=alternatives[[metric_name]])
    
    wilx_across_within = wilcox.test(rowMeans(results_all$Regular[[metric_name]])[!within_across_indices],
                                     rowMeans(results_all$Regular[[metric_name]])[within_across_indices],
                                     paired=F)
    g_wa <- 
      ggplot(plot_df_wthin_across) +
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
            axis.ticks = element_line(color="black"),
            axis.text=element_text(size=14, color="black"),
            plot.title=element_text(size=10),
            legend.position="NA") +
      geom_line(data=data.frame(x=c(1,1,2,2),
                                y=c(ylim_max * 1.005,
                                    ylim_max * 1.007,
                                    ylim_max * 1.007,
                                    ylim_max * 1.005)),
                aes(x=x,y=y)) +
      geom_text(x=1.5, y=ylim_max * 1.008, label=signif.num(wilx_across$p.value)) +
      
      geom_line(data=data.frame(x=c(3,3,4,4),
                                y=c(ylim_max * 1.005,
                                    ylim_max * 1.007,
                                    ylim_max * 1.007,
                                    ylim_max * 1.005)),
                aes(x=x,y=y)) +
      geom_text(x=3.5, y=ylim_max * 1.008, label=signif.num(wilx_within$p.value)) +
      ggtitle(metric_name) + 
      
      geom_line(data=data.frame(x=c(1,1,3,3),
                                y=c(ylim_max * 1.010,
                                    ylim_max * 1.012,
                                    ylim_max * 1.012,
                                    ylim_max * 1.010)),
                aes(x=x,y=y)) +
      geom_text(x=2, y=ylim_max * 1.013, label=signif.num(wilx_across_within$p.value)) +
      ggtitle(metric_name) + 
      xlab("") + 
      ylab("Structure similarity")
    
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\across_within",
                         write_path))
      dir.create(sprintf("%s\\across_within\\%s",
                         write_path,
                         size_name))
      pdf(sprintf("%s\\across_within\\%s\\structure_similarity_%s.pdf",
                  write_path,
                  size_name,
                  metric_name),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]] * 1.65) 
      
      plot(g_wa)
      dev.off()
    }
    
    similarity_matrix <- matrix(rep(self_comp_values[[metric_name]], 
                                    times=(len(datasets_names) ** 2)),
                                nrow=len(datasets_names))
    
    for (comp_idx in 1:ncol(combination_matrix)) {
      idx1 <- combination_matrix[1,comp_idx]
      idx2 <- combination_matrix[2,comp_idx]
      similarity_matrix[idx1, idx2] <- rowMeans(results_all$Regular[[metric_name]])[comp_idx]
      similarity_matrix[idx2, idx1] <- rowMeans(results_all$Regular[[metric_name]])[comp_idx]
    }
    
    effective_range <- seq(ylim_min, ylim_max, length.out=50)
    
    colnames(similarity_matrix) <- label_datasets_names
    rownames(similarity_matrix) <- label_datasets_names
    
    mice_colors <- blues_cg(len(unique(mice_names)))
    names(mice_colors) <- unique(mice_names)
    
    annoCol<-list(Mice=mice_colors)
    
    ph_sim <- 
    pheatmap(similarity_matrix,
             breaks = effective_range,
             color=rdylbu_cg(50),
             annotation_row=annotation_df,
             annotation_col=annotation_df,
             annotation_colors =annoCol,
             border_col="NA")
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\similarity_matrices",
                         write_path))
      dir.create(sprintf("%s\\similarity_matrices\\%s",
                         write_path,
                         size_name))
      pdf(sprintf("%s\\similarity_matrices\\%s\\similarity_matrix_%s.pdf",
                  write_path,
                  size_name,
                  metric_name),
          height=sizes[[size_name]][["width"]],
          width=sizes[[size_name]][["width"]] * 1.3) 
      
      plot(ph_sim[[4]])
      dev.off()
    }
  }
}

figure_2_pairwise_topological_permutation_tests <- function()
{
  write_path <- sprintf("%s\\figure_2\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\topological_similarity_plots", write_path)
  dir.create(write_path)  
  
  
  shuffled_pariwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
  pairwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
  
  id_perm_test = F
  metric_functions <- list(KL=KL_dist,
                           JSD=JSD_dist,
                           wasserstein=wasserstein_dist)
  
  if (!id_perm_test) {
    pvalue_matrices_betti_0 <- list()
    pvalue_matrices_betti_1 <- list()
    pvalue_shuffle_matrics_betti_0 <- list()
    pvalue_shuffle_matrics_betti_1 <- list()
    
    for (metric_name in names(metric_functions)) {
      pvalue_matrices_betti_0[[metric_name]] <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
      pvalue_matrices_betti_1[[metric_name]] <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
      pvalue_shuffle_matrics_betti_0[[metric_name]] <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
      pvalue_shuffle_matrics_betti_1[[metric_name]] <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
    }
  }
  
  nreps = 200
  results_all <- data.frame()
  
  for (pair_idx in 1:ncol(combination_matrix)) {
    
    idx1 <- combination_matrix[1, pair_idx]
    idx2 <- combination_matrix[2, pair_idx]
    
    p1 <- paths[[idx1]]
    p2 <- paths[[idx2]]
    
    if (idx1 == idx2) {
      next
    }
    
    orig_reg_mat_1 <- get_mat_with_preset(p1,"dalshtaim")
    orig_reg_mat_2 <- get_mat_with_preset(p2, "dalshtaim")
    orig_shuffled_mat_1 <- get_mat_with_preset(p1,"dalshtaimshuff")
    orig_shuffled_mat_time <- get_mat_with_preset(p1,"dalshtaimtimeshuff")
    
    if (id_perm_test) {
      reg_mat_1 <- get_reduced_mat_full_day(p1, 
                                            "lem", 
                                            ndim=20, 
                                            window_size=15, 
                                            knn1=0.075, 
                                            knn2 = 0, 
                                            shuffled=F, 
                                            time_shuffled=F)
      
      reg_mat_2 <- get_reduced_mat_full_day(p2, 
                                            "lem", 
                                            ndim=20, 
                                            window_size=15, 
                                            knn1=0.075, 
                                            knn2 = 0, 
                                            shuffled=F, 
                                            time_shuffled=F)
      
      shuffled_mat_1 <- get_reduced_mat_full_day(p1, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, time_shuffled=F)
      
    } else {
      
      reg_mat_1 <- kmeans(orig_reg_mat_1, centers=80, iter.max=500)$centers
      reg_mat_2 <- kmeans(orig_reg_mat_2, centers=80, iter.max=500)$centers
      shuffled_mat_1 <- kmeans(orig_shuffled_mat_1, centers=80, iter.max=500)$centers
      time_shuffled_mat_1 <- kmeans(orig_shuffled_mat_time, centers=80, iter.max=500)$centers
      
    }
    
    mtall <- rbind(reg_mat_1,
                   reg_mat_2)
    
    mtshuffle <- rbind(shuffled_mat_1,
                       reg_mat_2)
    
    mttimeshuffle <- rbind(time_shuffled_mat_1,
                           reg_mat_2)
    
    
    if (id_perm_test) {
    id_1 <- twonn(reg_mat_1, c_trimmed = 0.04, method="linfit")
    id_2 <- twonn(reg_mat_2, c_trimmed = 0.04, method="linfit")
    id_shuffle_1 <- twonn(shuffled_mat_1, c_trimmed = 0.04, method="linfit")
    
    original_diff <- abs(id_1$est[[2]] - id_2$est[[2]])
    original_diff_shuffle <- abs(id_shuffle_1$est[[2]] - id_2$est[[2]])
    
    diff_dist <- c()
    diff_dist_shuffle <- c()
    
    } else {
      phomology_similarity <- list()
      phomology_similarity_shuffle <- list()
      phomology_similarity_time_shuffle <- list()
      
      permutation_phomology_similarity <- list()
      permutation_phomology_similarity_shuffle <- list()
      permutation_phomology_similarity_time_shuffle <- list()
      
      phom_reg_mat_1 <- calculate_homology(reg_mat_1)
      phom_reg_mat_2 <- calculate_homology(reg_mat_2)
      phom_shuffle <- calculate_homology(shuffled_mat_1)
      phom_time_shuffle <- calculate_homology(time_shuffled_mat_1)
      
      # Barcode 0 
      life_death_0_mat_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 0,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 0,2]
      life_death_0_mat_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 0,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 0,2]
      life_death_0_mat_1_shuffle <- phom_shuffle[phom_shuffle[,1] == 0,3] - phom_shuffle[phom_shuffle[,1] == 0,2]
      life_death_0_mat_1_time_shuffle <- phom_time_shuffle[phom_time_shuffle[,1] == 0,3] - phom_time_shuffle[phom_time_shuffle[,1] == 0,2]
      
      # Barcode 1
      life_death_1_mat_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 1,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 1, 2]
      life_death_1_mat_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 1, 3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 1, 2]
      life_death_1_mat_1_shuffle <- phom_shuffle[phom_shuffle[,1] == 1,3] - phom_shuffle[phom_shuffle[,1] == 1, 2]
      life_death_1_mat_1_time_shuffle <- phom_time_shuffle[phom_time_shuffle[,1] == 1,3] - phom_time_shuffle[phom_time_shuffle[,1] == 1,2]
      

      
      
      for (metric_name in names(metric_functions)) {
        
        phomology_similarity[[metric_name]] <- list()
        phomology_similarity_shuffle[[metric_name]] <- list()
        phomology_similarity_time_shuffle[[metric_name]] <- list()
        
        permutation_phomology_similarity[[metric_name]] <- list()
        permutation_phomology_similarity_shuffle[[metric_name]] <- list()
        permutation_phomology_similarity[[metric_name]][["betti_0"]] <- c()
        permutation_phomology_similarity[[metric_name]][["betti_1"]] <- c()
        permutation_phomology_similarity_shuffle[[metric_name]][["betti_0"]] <- c()
        permutation_phomology_similarity_shuffle[[metric_name]][["betti_1"]] <- c()
        
        
        phomology_similarity[[metric_name]][["betti_0"]] <- metric_functions[[metric_name]](life_death_0_mat_1,
                                                                                            life_death_0_mat_2)
        phomology_similarity[[metric_name]][["betti_1"]] <- metric_functions[[metric_name]](life_death_1_mat_1,
                                                                                            life_death_1_mat_2)
        
        
        phomology_similarity_shuffle[[metric_name]][["betti_0"]] <-
                    metric_functions[[metric_name]](life_death_0_mat_1_shuffle, life_death_0_mat_2)
        
        phomology_similarity_shuffle[[metric_name]][["betti_1"]] <-
                    metric_functions[[metric_name]](life_death_1_mat_1_shuffle, life_death_1_mat_2)
        
        print(sprintf("Dataset (%d) vs (%d) b0: %.3f, b1: %.3f [%s]",
                      idx1,
                      idx2,
                      phomology_similarity[[metric_name]][["betti_0"]],
                      phomology_similarity[[metric_name]][["betti_1"]],
                      metric_name))
        
        print(sprintf("Dataset (%d) - Shuffled vs (%d) b0: %.3f, b1: %.3f [%s]",
                      idx1,
                      idx2,
                      phomology_similarity_shuffle[[metric_name]][["betti_0"]],
                      phomology_similarity_shuffle[[metric_name]][["betti_1"]],
                      metric_name))
        
        print("------")
        print("------")
      }
    }
    
    for (i in 1:nreps) {
      
      logical_ind <- 1:nrow(mtall) %in% sample(1:nrow(mtall), nrow(reg_mat_1))
      ind1 <- which(logical_ind)
      ind2 <- which(!logical_ind)
      
      if (id_perm_test) {
        id_perm_1 <- twonn(mtall[ind1,], c_trimmed = 0.04, method="linfit")
        id_perm_2 <- twonn(mtall[ind2,], c_trimmed = 0.04, method="linfit")
        
        id_perm_shuffle_1 <- twonn(mtshuffle[ind1,], c_trimmed = 0.04, method="linfit")
        id_perm_shuffle_2 <- twonn(mtshuffle[ind2,], c_trimmed = 0.04, method="linfit")
        
        perm_diff <- abs(id_perm_1$est[[2]] - id_perm_2$est[[2]])
        perm_diff_shuffle <- abs(id_perm_shuffle_1$est[[2]] - id_perm_shuffle_2$est[[2]])
        
        print(perm_diff)
        print(perm_diff_shuffle)
        diff_dist <- c(diff_dist, perm_diff)
        diff_dist_shuffle <- c(diff_dist_shuffle, perm_diff_shuffle)
      } else {
        
        perm_phom_reg_mat_1 <- calculate_homology(mtall[ind1,])
        perm_phom_reg_mat_2 <- calculate_homology(mtall[ind2,])
        perm_phom_shuffle_1 <- calculate_homology(mtshuffle[ind1,])
        perm_phom_shuffle_2 <- calculate_homology(mtshuffle[ind2,])
        
        # Barcode 0 
        perm_life_death_0_mat_1 <- 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 0,3] - 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 0,2]
        
        perm_life_death_0_mat_2 <- 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 0,3] - 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 0,2]
        
        perm_life_death_0_shuffle_1 <- 
          perm_phom_shuffle_1[perm_phom_shuffle_1[,1] == 0,3] - 
          perm_phom_shuffle_1[perm_phom_shuffle_1[,1] == 0,2]
        
        perm_life_death_0_shuffle_2 <- 
          perm_phom_shuffle_2[perm_phom_shuffle_2[,1] == 0,3] - 
          perm_phom_shuffle_2[perm_phom_shuffle_2[,1] == 0,2]
        
        # Barcode 1
        perm_life_death_1_mat_1 <- 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 1,3] - 
          perm_phom_reg_mat_1[perm_phom_reg_mat_1[,1] == 1,2]
        
        perm_life_death_1_mat_2 <- 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 1,3] - 
          perm_phom_reg_mat_2[perm_phom_reg_mat_2[,1] == 1,2]
        
        perm_life_death_1_shuffle_1 <- 
          perm_phom_shuffle_1[perm_phom_shuffle_1[,1] == 1,3] - 
          perm_phom_shuffle_1[perm_phom_shuffle_1[,1] == 1,2]
        
        perm_life_death_1_shuffle_2 <- 
          perm_phom_shuffle_2[perm_phom_shuffle_2[,1] == 1,3] - 
          perm_phom_shuffle_2[perm_phom_shuffle_2[,1] == 1,2]
        
        for (metric_name in names(metric_functions)) {
          
          permutation_phomology_similarity[[metric_name]][["betti_0"]] <- 
            c(permutation_phomology_similarity[[metric_name]][["betti_0"]],
              metric_functions[[metric_name]](perm_life_death_0_mat_1, perm_life_death_0_mat_2))
          
          permutation_phomology_similarity[[metric_name]][["betti_1"]] <- 
            c(permutation_phomology_similarity[[metric_name]][["betti_1"]],
              metric_functions[[metric_name]](perm_life_death_1_mat_1,perm_life_death_1_mat_2))
          
          
          
          permutation_phomology_similarity_shuffle[[metric_name]][["betti_0"]] <- 
            c(permutation_phomology_similarity_shuffle[[metric_name]][["betti_0"]],
              metric_functions[[metric_name]](perm_life_death_0_shuffle_1, perm_life_death_0_shuffle_2))
          
          permutation_phomology_similarity_shuffle[[metric_name]][["betti_1"]] <- 
            c(permutation_phomology_similarity_shuffle[[metric_name]][["betti_1"]],
              metric_functions[[metric_name]](perm_life_death_1_shuffle_1, perm_life_death_1_shuffle_2))
        }
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
      
      
      results_all <- rbind(results_all,
                           data.frame(path1=idx1, 
                                      path2=idx2, 
                                      betti0=phomology_similarity_shuffle[[metric_name]]$betti_0,
                                      betti1=phomology_similarity_shuffle[[metric_name]]$betti_1,
                                      isShuffle="Shuffle",
                                      metric=metric_name))
      
      pval_func_betti_0 <- ecdf(permutation_phomology_similarity[[metric_name]]$betti_0)
      pval_func_betti_1 <- ecdf(permutation_phomology_similarity[[metric_name]]$betti_1)
      pval_func_shuffle_betti_0 <- ecdf(permutation_phomology_similarity_shuffle[[metric_name]]$betti_0)
      pval_func_shuffle_betti_1 <- ecdf(permutation_phomology_similarity_shuffle[[metric_name]]$betti_1)
      
      
      print(sprintf("Pvalues b0 (%.3f) b1(%.3f) b0[shuffle] (%.3f) b1[shuffle] (%.3f) %s",
                    1 - pval_func_betti_0(phomology_similarity[[metric_name]]$betti_0),
                    1 - pval_func_betti_1(phomology_similarity[[metric_name]]$betti_1),
                    1 - pval_func_shuffle_betti_0(phomology_similarity_shuffle[[metric_name]]$betti_0),
                    1 - pval_func_shuffle_betti_1(phomology_similarity_shuffle[[metric_name]]$betti_1),
                    metric_name))
      
      pvalue_matrices_betti_0[[metric_name]][idx1, idx2] <- 1 - pval_func_betti_0(phomology_similarity[[metric_name]]$betti_0)
      pvalue_matrices_betti_1[[metric_name]][idx1, idx2] <- 1 - pval_func_betti_1(phomology_similarity[[metric_name]]$betti_1)
      pvalue_shuffle_matrics_betti_0[[metric_name]][idx1, idx2] <- 1 - pval_func_shuffle_betti_0(phomology_similarity_shuffle[[metric_name]]$betti_0)
      pvalue_shuffle_matrics_betti_1[[metric_name]][idx1, idx2] <- 1 - pval_func_shuffle_betti_1(phomology_similarity_shuffle[[metric_name]]$betti_1)
    }
    
    pheatmap(pvalue_matrices_betti_0[["JSD"]], breaks=c(0,0.05,0.1,0.15,1), col=rdylbu_cg(4))

  }
  
  within_across_indices <- apply(combination_matrix, 2, function(c) {mice_names[c[1]] == mice_names[c[2]]})
  within_across_names <- rep("Across", times=len(within_across_indices))
  within_across_names[within_across_indices] <- "Within"
  
  regular_labels <- rep("Regular", times=ncol(combination_matrix))
  shuffle_labels <- rep("Shuffle", times=ncol(combination_matrix))
  
  regular_within_across_labels <- paste(regular_labels, within_across_names, sep=" - ")
  shuffle_within_across_labels <- paste(shuffle_labels, within_across_names, sep=" - ")
  barcodes = c(b0="betti0", b1="betti1")
  barcodes_labs = c(b0=TeX(sprintf("Topological similarity (%s)", "$\\beta_{0}$")),
                    b1=TeX(sprintf("Topological similarity (%s)", "$\\beta_{1}$")))
  
  for (barcode_name in names(barcodes)){
    barcode =  barcodes[[barcode_name]]
    
    for (metric_name in names(metric_functions)) {
      shuffle_df <- results_all[results_all$metric == metric_name &
                                  results_all$isShuffle == "Shuffle" ,]
      
      regular_df <- results_all[results_all$metric == metric_name &
                                  results_all$isShuffle == "Regular" ,]
      
      plot_df <- data.frame(phom=c(regular_df[[barcode]],shuffle_df[[barcode]]),
                            all_group=c(regular_labels, shuffle_labels))
      
      
      ylim_max = max(plot_df$phom)
      ylim_min = min(plot_df$phom)
      
      wilx = wilcox.test(regular_df[[barcode]],
                         shuffle_df[[barcode]],
                         paired=T,
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
                                  y=c(ylim_max * 1.005,
                                      ylim_max * 1.007,
                                      ylim_max * 1.007,
                                      ylim_max * 1.005)),
                  aes(x=x,y=y)) +
        geom_text(x=1.5, y=ylim_max * 1.008, label=signif.num(wilx$p.value)) +
        ggtitle(metric_name) + 
        xlab("") + 
        ylab(barcodes_labs[[barcode_name]])
      
      
      for (size_name in names(sizes)) {
        
        dir.create(sprintf("%s\\all",
                           write_path))
        dir.create(sprintf("%s\\all\\%s",
                           write_path,
                           size_name))
        pdf(sprintf("%s\\all\\%s\\topological_similarity_%s_%s.pdf",
                    write_path,
                    size_name,
                    barcode_name,
                    metric_name),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]]) 
        
        plot(g)
        dev.off()
      }
      
      plot_df_wthin_across <- 
        data.frame(Similarity=c(regular_df[[barcode]],shuffle_df[[barcode]]),
                   Group=factor(c(regular_within_across_labels, shuffle_within_across_labels),
                                levels=c(unique(regular_within_across_labels)[2],
                                         unique(shuffle_within_across_labels)[2],
                                         unique(regular_within_across_labels)[1],
                                         unique(shuffle_within_across_labels)[1])))
      
      #ylim_max = max(plot_df$Similarity)
      
      wilx_within = wilcox.test(regular_df[[barcode]][within_across_indices],
                                shuffle_df[[barcode]][within_across_indices],
                                paired=T,
                                alternative="less")
      
      wilx_across = wilcox.test(regular_df[[barcode]][!within_across_indices],
                                shuffle_df[[barcode]][!within_across_indices],
                                paired=T,
                                alternative="less")
      
      wilx_across_within = wilcox.test(regular_df[[barcode]][!within_across_indices],
                                regular_df[[barcode]][within_across_indices],
                                paired=F)
      g_wa <- 
        ggplot(plot_df_wthin_across) +
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
              axis.ticks = element_line(color="black"),
              axis.text=element_text(size=14, color="black"),
              plot.title=element_text(size=10),
              legend.position="NA") +
        geom_line(data=data.frame(x=c(1,1,2,2),
                                  y=c(ylim_max * 1.005,
                                      ylim_max * 1.007,
                                      ylim_max * 1.007,
                                      ylim_max * 1.005)),
                  aes(x=x,y=y)) +
        geom_text(x=1.5, y=ylim_max * 1.008, label=signif.num(wilx_across$p.value)) +
        
        geom_line(data=data.frame(x=c(3,3,4,4),
                                  y=c(ylim_max * 1.005,
                                      ylim_max * 1.007,
                                      ylim_max * 1.007,
                                      ylim_max * 1.005)),
                  aes(x=x,y=y)) +
        geom_text(x=3.5, y=ylim_max * 1.008, label=signif.num(wilx_within$p.value)) +
        ggtitle(metric_name) + 
        
        geom_line(data=data.frame(x=c(1,1,3,3),
                                  y=c(ylim_max * 1.010,
                                      ylim_max * 1.012,
                                      ylim_max * 1.012,
                                      ylim_max * 1.010)),
                  aes(x=x,y=y)) +
        geom_text(x=2, y=ylim_max * 1.013, label=signif.num(wilx_across_within$p.value)) +
        ggtitle(metric_name) + 
        xlab("") + 
        ylab(barcodes_labs[[barcode_name]])
      
      
      for (size_name in names(sizes)) {
        
        dir.create(sprintf("%s\\across_within",
                           write_path))
        dir.create(sprintf("%s\\across_within\\%s",
                           write_path,
                           size_name))
        pdf(sprintf("%s\\across_within\\%s\\topological_similarity_%s_%s.pdf",
                    write_path,
                    size_name,
                    barcode_name,
                    metric_name),
            height=sizes[[size_name]][["width"]],
            width=sizes[[size_name]][["width"]] * 1.65) 
        
        plot(g_wa)
        dev.off()
      }
      
    }
  }
}

figure_2_pairwise_topological_permutation <- function()
{
  
  
  b0_shuffled_pariwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
  b0_pairwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
  b1_shuffled_pariwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))
  b1_pairwise <- matrix(rep(0, times=len(paths) ** 2), nrow=len(paths))

  nreps = 200
  results_all <- data.frame()
  
  for (pair_idx in 1:ncol(combination_matrix)) {
    
    path_idx_1 <- combination_matrix[1, pair_idx]
    path_idx_2 <- combination_matrix[2, pair_idx]
    p1 <- path_list[path_idx_1]
    p2 <- path_list[path_idx_2]
      
      if (path_idx_1 == path_idx_2) {
        next
      }
      
      orig_reg_mat_1 <- get_mat_with_preset(p1,"dalshtaim")
      orig_reg_mat_2 <- get_mat_with_preset(p2, "dalshtaim")
      orig_shuffled_mat_1 <- get_mat_with_preset(p1,"dalshtaimshuff")
      reg_mat_1 <- kmeans(orig_reg_mat_1, centers=80, iter.max=500)$centers
      reg_mat_2 <- kmeans(orig_reg_mat_2, centers=80, iter.max=500)$centers
      shuffled_mat_1 <- kmeans(orig_shuffled_mat_1, centers=80, iter.max=500)$centers
      
      prmt <- permutation_test(reg_mat_1, reg_mat_2, iterations=nreps)
      prmt_shuffle <- permutation_test(shuffled_mat_1, reg_mat_2, iterations=nreps)
      
      b0_pairwise[path_idx_1,path_idx_2] <- prmt[[1]]$pvalue
      b1_pairwise[path_idx_1,path_idx_2] <- prmt[[2]]$pvalue
      b0_shuffled_pariwise[path_idx_1,path_idx_2] <- prmt_shuffle[[1]]$pvalue
      b1_shuffled_pariwise[path_idx_1,path_idx_2] <- prmt_shuffle[[2]]$pvalue
      
      b0_pairwise[path_idx_2,path_idx_1] <- prmt[[1]]$pvalue
      b1_pairwise[path_idx_2,path_idx_1] <- prmt[[2]]$pvalue
      b0_shuffled_pariwise[path_idx_2,path_idx_1] <- prmt_shuffle[[1]]$pvalue
      b1_shuffled_pariwise[path_idx_2,path_idx_1] <- prmt_shuffle[[2]]$pvalue

      results_all <- rbind(results_all,
                           data.frame(path1=path_idx_1, 
                                      path2=path_idx_2, 
                                      betti0=prmt[[1]]$wasserstein,
                                      betti1=prmt[[2]]$wasserstein,
                                      isShuffle="Regular"))
      
      results_all <- rbind(results_all,
                           data.frame(path1=path_idx_1, 
                                      path2=path_idx_2, 
                                      betti0=prmt_shuffle[[1]]$wasserstein,
                                      betti1=prmt_shuffle[[2]]$wasserstein,
                                      isShuffle="Shuffle"))
      
      print(sprintf("Pvalues b0 (%.3f) b1(%.3f) b0[shuffle] (%.3f) b1[shuffle] (%.3f)",
                    prmt[[1]]$pvalue,
                    prmt[[2]]$pvalue,
                    prmt_shuffle[[1]]$pvalue,
                    prmt_shuffle[[2]]$pvalue))
      
      
  }
}

figure_2_example_dimensionality_dataset <- function() {
  
  write_path <- sprintf("%s\\figure_2\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\structure_examples", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=2,
                     height=2),
               medium=c(width=1.5,
                        height=1.5),
               medium_2=c(width=1.25,
                          height=1.25),
               medium_1=c(width=1,
                          height=1),
               small=c(width=.75,
                       height=.75))
  paths <- get_thirsty_quenched_paths()
  
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  for (i in 1:len(paths)) {
  reduced_mat <- get_mat_with_preset(paths[[i]], "dalshtaim")  
  
  p_lem2_3 <- pair_plots_colored(reduced_mat,
                                 rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                 plot_nrow = 3)
  
  p_lem2_5 <- pair_plots_colored(reduced_mat,
                                 rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                 plot_nrow = 5)
  
  
  p_lem2_list <- pair_plots_colored(reduced_mat, 
                                    rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                    return_grid = F,
                                    pt_size=.55)  
  
  p_lem2_list <- lapply(1:len(p_lem2_list), 
                        function(i) {return(p_lem2_list[[i]] + 
                                              theme(line=element_blank(),
                                                    #rect=element_rect(color="white"),
                                                    plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                    panel.border=element_rect(color="white"),
                                                    plot.background = element_rect("white"),
                                                    axis.title=element_text(color="white")))})
  
  
  for (size_name in names(sizes)) {
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    png(sprintf("%s\\%s\\lem2_%s_3_rows.png",
                write_path,
                size_name,
                datasets_names[i]),
        height=sizes[[size_name]][["height"]] * 3,
        width=sizes[[size_name]][["width"]] * 5,
        unit="in",
        res=1500)
    
    plot(p_lem2_3)
    dev.off()
    
    png(sprintf("%s\\%s\\lem2_%s_5_rows.png",
                write_path,
                size_name,
                datasets_names[i]),
        height=sizes[[size_name]][["height"]] * 5,
        width=sizes[[size_name]][["width"]] * 3,
        unit="in",
        res=1500)
    
    plot(p_lem2_5)
    dev.off()
    
    p_lem2_5_no_annot <- p_lem2_list
    p_lem2_3_no_annot <- p_lem2_list
    
    p_lem2_5_no_annot$nrow <- 5
    pf_5 <- do.call(arrangeGrob, p_lem2_5_no_annot)
    
    p_lem2_3_no_annot$nrow <- 3
    pf_3 <- do.call(arrangeGrob, p_lem2_3_no_annot)
    
    png(sprintf("%s\\%s\\lem2_%s_3_rows_no_annot.png",
                write_path,
                size_name,
                datasets_names[i]),
        height=sizes[[size_name]][["height"]] * 3,
        width=sizes[[size_name]][["width"]] * 5,
        unit="in",
        res=1500)
    
    plot(pf_3)
    dev.off()
    
    
    png(sprintf("%s\\%s\\lem2_%s_5_rows_no_annot.png",
                write_path,
                size_name,
                datasets_names[i]),
        height=sizes[[size_name]][["height"]] * 5,
        width=sizes[[size_name]][["width"]] * 3,
        unit="in",
        res=1500)
    
    plot(pf_5)
    dev.off()
    
  }

  } 
}



figure_2_pairwise_topological_permutation_generic <- function(metadata_1, metadata_2, remove_diagonal = F, niter = 10, single_example=F)
{

  metric_functions <- list(KL=KL_dist,
                           JSD=JSD_dist,
                           wasserstein=wasserstein_dist)
  
  pvalue_matrices_betti_0 <- list()
  pvalue_matrices_betti_1 <- list()
  iter_pvalue_matrices_betti_0 <- list()
  iter_pvalue_matrices_betti_1 <- list()
  pvalue_shuffle_matrics_betti_0 <- list()
  pvalue_shuffle_matrics_betti_1 <- list()
    
  for (metric_name in names(metric_functions)) {
      pvalue_matrices_betti_0[[metric_name]] <- matrix(rep(0, times=len(metadata_1) * len(metadata_2)), nrow=len(metadata_1))
      pvalue_matrices_betti_1[[metric_name]] <- matrix(rep(0, times=len(metadata_1) * len(metadata_2)), nrow=len(metadata_1))
      iter_pvalue_matrices_betti_0[[metric_name]] <- c()
      iter_pvalue_matrices_betti_1[[metric_name]] <- c()
  }
  
  nreps = 200
  results_all <- data.frame()
  
  
  
  for (i1 in 1:len(metadata_1)) {
    for (i2 in 1:len(metadata_2)) {
      
      if (remove_diagonal && i1 == i2) {
        next
      }
      
      iter_results_all <- data.frame()
      iter_metric_pvalue_betti0 <- list()
      iter_metric_pvalue_betti1 <- list()
      
      for (metric_name in names(metric_functions)) {
       iter_metric_pvalue_betti0[[metric_name]] <- c()
       iter_metric_pvalue_betti0[[metric_name]] <- c()
      }
      
      for (ni_idx in 1:niter) { 

        
        orig_reg_mat_1 <- metadata_1[[i1]]$red_mat
        orig_reg_mat_2 <- metadata_2[[i2]]$red_mat
        reg_mat_1 <- kmeans(orig_reg_mat_1, centers=80, iter.max=500)$centers
        reg_mat_2 <- kmeans(orig_reg_mat_2, centers=80, iter.max=500)$centers
        
        mtall <- rbind(reg_mat_1, reg_mat_2)
        
        phomology_similarity <- list()
        permutation_phomology_similarity <- list()
        
        phom_reg_mat_1 <- calculate_homology(reg_mat_1)
        phom_reg_mat_2 <- calculate_homology(reg_mat_2)
        
        if (single_example) {
          
          betti0_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 0,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 0,2]
          betti0_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 0,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 0,2]
          
          
          plot_df_b0 <- data.frame(life=c(betti0_1, betti0_2),
                                   grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(betti0_1), len(betti0_2)))))
          betti1_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 1,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 1,2]
          betti1_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 1,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 1,2]
          
          
          plot_df_b1 <- data.frame(life=c(betti1_1, betti1_2),
                                   grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(betti1_1), len(betti1_2)))))        
          
          
          geaxmple_b0 <- 
            ggplot(plot_df_b0, aes(x = life, fill=grp)) + 
            geom_density(lwd = 1.2,
                         linetype = 2,
                         colour = "black", alpha=.5) +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA")  +
            xlab(TeX(sprintf("Feature lifespan (%s)", "$\\beta_{0}$"))) +
            ylab("Density") +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA") +
            scale_y_continuous(expand=c(0,0))+
            scale_x_continuous(expand=c(0,0)) +
            scale_fill_manual(breaks=c("Dataset 1", "Dataset 2"),
                              values=c("#8430ab", "#ffa200"))
          
          geaxmple_b1 <- 
            ggplot(plot_df_b1, aes(x = life, fill=grp)) + 
            geom_density(lwd = 1.2,
                         linetype = 2,
                         colour = "black", alpha=.5) +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA")  +
            xlab(TeX(sprintf("Feature lifespan (%s)", "$\\beta_{1}$"))) +
            ylab("Density") +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA") +
            scale_y_continuous(expand=c(0,0))+
            scale_x_continuous(expand=c(0,0)) +
            scale_fill_manual(breaks=c("Dataset 1", "Dataset 2"),
                              values=c("#76ba00", "#ff0059"))
          
          
          phom_reg_mat_1 <- as.data.frame(phom_reg_mat_1)
          phom_reg_mat_2 <- as.data.frame(phom_reg_mat_2)
          phom_reg_mat_1$nrow <- 1:nrow(phom_reg_mat_1)
          phom_reg_mat_2$nrow <- 1:nrow(phom_reg_mat_2)
          colnames(phom_reg_mat_1) <- c("Dim", "Birth", "Death", "Row")
          colnames(phom_reg_mat_2) <- c("Dim", "Birth", "Death", "Row")
          gexample_dataset_1 <- 
            ggplot(as.data.frame(phom_reg_mat_1)) +
            geom_line(data=data.frame(x=c(0, max(phom_reg_mat_1[,c("Birth", "Death")])),
                                      y=c(0, max(phom_reg_mat_1[,c("Birth", "Death")]))),
                      aes(x=x,y=y)) +
            geom_point(aes(x=Birth, y=Death, color=factor(Dim)), alpha=.6) +
            xlab("Feature birth") +
            ylab("Featue death") +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA") +
            scale_color_manual(breaks=c(0, 1),
                               values=c("#8430ab", "#76ba00"))
          
          gexample_dataset_2 <- 
            ggplot(as.data.frame(phom_reg_mat_2)) +
            geom_line(data=data.frame(x=c(0, max(phom_reg_mat_1[,c("Birth", "Death")])),
                                      y=c(0, max(phom_reg_mat_1[,c("Birth", "Death")]))),
                      aes(x=x,y=y)) +
            geom_point(aes(x=Birth, y=Death, color=factor(Dim)), alpha=.6) +
            xlab("Feature birth") +
            ylab("Featue death") +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA") +
            scale_color_manual(breaks=c(0, 1),
                               values=c("#ffa200", "#ff0059"))
          
          
          
          gbarcode_dataset_1 <- 
            ggplot() +
            xlab("Search radius") +
            ylab("") +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA",
                  axis.line.y = element_blank(),
                  axis.text.y  = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank()) + 
            scale_x_continuous(expand=c(0,0))
          
          for (rowidx in 1:nrow(phom_reg_mat_1)) {
            row <- phom_reg_mat_1[rowidx,]
            
            col = ifelse(row[["Dim"]] == 0, "#8430ab", "#76ba00")
            gbarcode_dataset_1 <- 
              gbarcode_dataset_1 + 
              geom_line(data=data.frame(x=c(row[["Birth"]], row[["Death"]]),
                                        y=c(row[["Row"]], row[["Row"]])),
                        aes(x=x, y=y),
                        color=col)
            
            
          }
          
          gbarcode_dataset_2 <- 
            ggplot() +
            xlab("Search radius") +
            ylab("") +
            theme_classic() +
            theme(text=element_text(size=14, color="black"),
                  legend.position="NA",
                  axis.line.y = element_blank(),
                  axis.text.y  = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank()) + 
            scale_x_continuous(expand=c(0,0))
          
          for (rowidx in 1:nrow(phom_reg_mat_2)) {
            row <- phom_reg_mat_2[rowidx,]
            
            col = ifelse(row[["Dim"]] == 0,"#ffa200", "#ff0059")
            gbarcode_dataset_2 <- 
              gbarcode_dataset_2 + 
              geom_line(data=data.frame(x=c(row[["Birth"]], row[["Death"]]),
                                        y=c(row[["Row"]], row[["Row"]])),
                        aes(x=x, y=y),
                        color=col)
          }
          
          
          for (size_name in names(sizes)) {
            write_path <- sprintf("%s\\%s", figures_base_path, "topo_illustrations")
            dir.create(write_path)
            write_path <- sprintf("%s\\%s", write_path, size_name)
            dir.create(write_path)
            
            pdf(sprintf("%s\\dataset_1_barcode.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(gbarcode_dataset_1)
            dev.off()
            
            pdf(sprintf("%s\\dataset_2_barcode.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(gbarcode_dataset_2)
            dev.off()
            
            pdf(sprintf("%s\\dataset_1_persist.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(gexample_dataset_1)
            dev.off()
            
            pdf(sprintf("%s\\dataset_2_persist.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(gexample_dataset_2)
            dev.off()
            
            pdf(sprintf("%s\\distribution_b0.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(geaxmple_b0)
            dev.off()
            
            pdf(sprintf("%s\\distribution_b1.pdf",
                        write_path),
                height=sizes[[size_name]],
                width=sizes[[size_name]])
            
            plot(geaxmple_b1)
            dev.off()
            
          }
          
          single_example = F
          
        }
        
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
                        i1,
                        i2,
                        phomology_similarity[[metric_name]][["betti_0"]],
                        phomology_similarity[[metric_name]][["betti_1"]],
                        metric_name))
          
          
          print("------")
          print("------")
        }
        
        
        for (i in 1:nreps) {
          
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
            
            
            if (len(perm_life_death_1_mat_1) == 0 || len(perm_life_death_1_mat_2) == 0) {
              next
            }
            
            permutation_phomology_similarity[[metric_name]][["betti_1"]] <- 
              c(permutation_phomology_similarity[[metric_name]][["betti_1"]],
                metric_functions[[metric_name]](perm_life_death_1_mat_1,perm_life_death_1_mat_2))
            
          }
        }
        
        
        for (metric_name in names(metric_functions)) {
          iter_results_all <- rbind(iter_results_all,
                               data.frame(path1=i1, 
                                          path2=i2, 
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
          
         
          iter_metric_pvalue_betti0[[metric_name]] <- c(iter_metric_pvalue_betti0[[metric_name]], 1 - pval_func_betti_0(phomology_similarity[[metric_name]]$betti_0))
          iter_metric_pvalue_betti1[[metric_name]] <- c(iter_metric_pvalue_betti1[[metric_name]], 1 - pval_func_betti_1(phomology_similarity[[metric_name]]$betti_1))
           
        }
      }
      
      for (metric_name in names(metric_functions)) {
        pvalue_matrices_betti_0[[metric_name]][i1, i2] <- median(iter_metric_pvalue_betti0[[metric_name]])
        pvalue_matrices_betti_1[[metric_name]][i1, i2] <- median(iter_metric_pvalue_betti1[[metric_name]])
        
        iter_pvalue_matrices_betti_0[[metric_name]] <- rbind(iter_pvalue_matrices_betti_0[[metric_name]],
                                                             iter_metric_pvalue_betti0[[metric_name]])
        iter_pvalue_matrices_betti_1[[metric_name]] <- rbind(iter_pvalue_matrices_betti_1[[metric_name]],
                                                             iter_metric_pvalue_betti1[[metric_name]])
        
      }
      
      results_all <- rbind(results_all,
                           ddply(iter_results_all, .(metric), function(met_df) {c(colMeans(met_df[,1:4]), met_df[1,5])}))
      
      pheatmap(pvalue_matrices_betti_0[["wasserstein"]], breaks=seq(0,1, length.out=100), col=rdylbu_cg(100))
      
    }
    }
 
  
  final_res = list(topological_similarity=results_all,
                  pvalue_b0=pvalue_matrices_betti_0,
                  pvalue_b1=pvalue_matrices_betti_1,
                  all_iters_pval_b0=iter_pvalue_matrices_betti_0,
                  all_iters_pval_b1=iter_pvalue_matrices_betti_1)
  
  save(file=sprintf("%s\\figure_2\\data\\%s",
                    figures_base_path,
                    output_name),
       final_res)
  
  return(final_res)
}

v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")

por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
               "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
               "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\")

control_paths <-  c(v1_paths, por_paths)


### generate_topological_similarities







metadata_14 <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
m14_m14 <- figure_2_pairwise_topological_permutation_generic(metadata_14, metadata_14)

metadata_14_shuffle <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths(), preset="dalshtaimshuff")
m14_m14shuffle <- figure_2_pairwise_topological_permutation_generic(metadata_14, metadata_14_shuffle)

metadata_SFO <- across_mice_decoding_build_metadata(path_list=get_SFO_paths())
m14_mSFO <- figure_2_pairwise_topological_permutation_generic(metadata_14, metadata_SFO)

metadata_hunger <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
m14_mHunger <- figure_2_pairwise_topological_permutation_generic(metadata_14, metadata_hunger)

metadata_v1_por <- across_mice_decoding_build_metadata(v1_paths, control=T)
mHunger_mv1POR <- figure_2_pairwise_topological_permutation_generic(metadata_hunger, metadata_v1_por)
mThirst_mv1POR <- figure_2_pairwise_topological_permutation_generic(metadata_14, metadata_v1_por)


mHunger_mAGRP <- figure_2_pairwise_topological_permutation_generic(metadata_hunger,
                                                                    metadata_agrp)

mHunger_mHunger <- figure_2_pairwise_topological_permutation_generic(metadata_hunger, metadata_hunger)


mAGRP_mAGRP <- figure_2_pairwise_topological_permutation_generic(metadata_agrp, metadata_agrp)

# all_results <- list(mHunger_mv1POR=mHunger_mv1POR,
#                     mThirst_mv1POR=mThirst_mv1POR,
#                     m14_mHunger=m14_mHunger,
#                     m14_mSFO=m14_mSFO,
#                     m14_m14shuffle=m14_m14shuffle,
#                     m14_m14=m14_m14)
all_results <- list(mHunger_mHunger=mHunger_mHunger,
                    mHunger_mAGRP=mHunger_mAGRP,
                    mAGRP_mAGRP=mAGRP_mAGRP)
                    
save(file=sprintf("%s/figure_2_comp.Rda", figures_base_path), all_results)



mHunger_mv1POR = all_results$mHunger_mv1POR
mThirst_mv1POR = all_results$mThirst_mv1POR
m14_mHunger = all_results$m14_mHunger
m14_mSFO = all_results$m14_mSFO
m14_m14shuffle = all_results$m14_m14shuffle
m14_m14 = all_results$m14_m14



write_path <- sprintf("%s\\figure_2\\", figures_base_path)
dir.create(write_path)
write_path <- sprintf("%s\\structure_similarity_statistics", write_path)
dir.create(write_path)  


sizes = list(big=c(width=3),
             medium=c(width=2.75),
             medium_2=c(width=2.5),
             medium_1=c(width=2),
             small=c(width=1.75))


#unique(m14_m14$topological_similarity$metric)
for (metric in c("KL","JSD","wasserstein")) {
  for (barc in c("betti0", "betti1")) {
    
    # m14_m14_vec <- ddply(m14_m14$topological_similarity[m14_m14$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    # m14_mSFO_vec <- ddply(m14_mSFO$topological_similarity[m14_mSFO$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    # m14_mHunger_vec <- ddply(m14_mHunger$topological_similarity[m14_mHunger$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    # m14_m14shuffle_vec <- ddply(m14_m14shuffle$topological_similarity[m14_m14shuffle$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    # mThirst_mv1POR_vec <- ddply(mThirst_mv1POR$topological_similarity[mThirst_mv1POR$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    # mHunger_mv1POR_vec <- ddply(mHunger_mv1POR$topological_similarity[mHunger_mv1POR$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    mHunger_mHunger_vec <- ddply(mHunger_mHunger$topological_similarity[mHunger_mHunger$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    mHunger_mAGRP_vec <- ddply(mHunger_mAGRP$topological_similarity[mHunger_mAGRP$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    mAGRP_mAGRP_vec <- ddply(mAGRP_mAGRP$topological_similarity[mAGRP_mAGRP$topological_similarity$metric == metric,], .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    
    # 
    # 
    # vals <- list(`1`=m14_m14_vec,
    #              `2`=m14_mSFO_vec,
    #              `3`=m14_mHunger_vec,
    #              `4`=m14_m14shuffle_vec,
    #              `5`=mThirst_mv1POR_vec,
    #              `6`=mHunger_mv1POR_vec)
    
    vals <- list(`7`=mHunger_mHunger_vec,
                 `8`=mHunger_mAGRP_vec,
                 `9`=mAGRP_mAGRP_vec)
    
    pval_name <- ifelse(barc == "betti0", "pvalue_b0", "pvalue_b1")
    
    pvals <- list(#`1`=c(m14_m14[[pval_name]][[metric]]),
                  #`2`=c(m14_mSFO[[pval_name]][[metric]]),
                  #`3`=c(m14_mHunger[[pval_name]][[metric]]),
                  #`4`=c(m14_m14shuffle[[pval_name]][[metric]]),
                  #`5`=c(mThirst_mv1POR[[pval_name]][[metric]]),
                  #`6`=c(mHunger_mv1POR[[pval_name]][[metric]]),
                  `7`=c(mHunger_mHunger[[pval_name]][[metric]]),
                  `8`=c(mHunger_mAGRP[[pval_name]][[metric]]),
                  `9`=c(mAGRP_mAGRP[[pval_name]][[metric]]))
    # 
    comparisions <- 
              list(#`wt_sfo`=list(comp_ind=list(c(1,2)), is_paired=T),
                   #`wt_hunger`=list(comp_ind=list(c(1,3)), is_paired=T),
                   #`wt_shuffle`=list(comp_ind=list(c(1,4)), is_paired=T),
                   #`wt_porv1`=list(comp_ind=list(c(1,5)), is_paired=T),
                   #`hunger_porv1`=list(comp_ind=list(c(3,6)), is_paired=F),
                   #`all`=list(comp_ind=list(c(1,2), c(1,3), c(1,4), c(1,5)), is_paired=T),
                  "hunger_agrp"=list(comp_ind=list(c(7,8)), is_paired=T),
                  "agrp_agrp"=list(comp_ind=list(c(9,8)), is_paired=T))
    
    labels <-  
      list(`1`="m14_m14",
           `2`="m14_mSFO",
           `3`="m14_mHunger",
           `4`="m14_m14shuffle",
           `5`="mThirst_mv1POR",
           `6`="mHunger_mv1POR_vec",
           `7`="mHunger_mHunger",
           `8`="mHunger_mAGRP",
           `9`="mAGRP_mAGRP")         
             
             
    
    comp_barcodes_labs = c(betti0=TeX(sprintf("Topological similarity (%s)", "$\\beta_{0}$")),
                           betti1=TeX(sprintf("Topological similarity (%s)", "$\\beta_{1}$")))
    
    pval_barcodes_labs = c(betti0=TeX(sprintf("%s (%s)", "$\\P_{value}$",  "$\\beta_{0}$")),
                           betti1=TeX(sprintf("%s (%s)", "$\\P_{value}$", "$\\beta_{1}$")))
    
    
    for (comp_conf_name in names(comparisions)) {
      
      comp_conf = comparisions[[comp_conf_name]]
      comps_df <- data.frame()
      pvalues_df <- data.frame()
      
      for (comp in unique(unlist(comp_conf$comp_ind))){
        comps_df <- rbind(comps_df,
                          data.frame(similarity=vals[[as.character(comp)]],
                                     x=rep(labels[[as.character(comp)]], times=len(vals[[as.character(comp)]]))))
        
        pvalues_df <- rbind(pvalues_df,
                            data.frame(pval=pvals[[as.character(comp)]],
                                       x=rep(labels[[as.character(comp)]], times=len(pvals[[as.character(comp)]]))))
      }
      
      
      gpval_bar <- 
      ggplot(pvalues_df, aes(x=x,y=pval)) + 
        #geom_point() + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 
        geom_bar(stat="summary", alpha=.8, width=.45, fill="gray70") +
        geom_errorbar(stat="summary", width=.25, size=.75) + 
        scale_y_continuous(expand=c(0,0)) + 

        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(pval_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      gpval_box <- 
        ggplot(pvalues_df, aes(x=x,y=pval)) + 
        #geom_point() + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 

        scale_y_continuous(expand=c(0,0)) + 
        geom_boxplot(aes(),
                     width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9),
                     outlier.shape = NA) + 
        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(pval_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      
      gcomp_bar <- 
        ggplot(comps_df, aes(x=x,y=similarity)) + 
        #geom_point() + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 
        geom_bar(stat="summary", alpha=.8, width=.45, fill="gray70") +
        geom_errorbar(stat="summary", width=.25, size=.75) + 
        scale_y_continuous(expand=c(0,0)) + 
        
        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(comp_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      gcomp_box <- 
        ggplot(comps_df, aes(x=x,y=similarity)) + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 
        geom_boxplot(aes(),
                     width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9),
                     outlier.shape = NA) + 
        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(comp_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      # 
      # for (comp in comp_conf$comp_ind) {
      #   
      #   vals_1 <- vals[[(comp[[1]])]]
      #   vals_2 <- vals[[(comp[[2]])]]
      #   
      #   pvals_1 <-   pvals[[(comp[[1]])]]
      #   pvals_2 <-   pvals[[(comp[[2]])]]
      #   
      #   wlk <- wilcox.test(vals_1, vals_2, paired=comp_conf$is_paired, alternative="less")
      #   
      #   max_val <- max(vals_1, vals_2, na.rm=T)
      #   
      #   if (len(comp_conf$comp_ind) == 1) { 
      #     line_df <- data.frame(x=c(1,1,2,2),
      #                         y=c(max_val, max_val * 1.05, max_val * 1.05, max_val))
      #     signif_df <- data.frame(x=1.5, y=max_val, label=signif.num(wlk$p.value))
      #   } else {
      #     line_df <- data.frame(x=c(comp[[1]],comp[[1]],comp[[2]],comp[[2]]),
      #                           y=c(max_val, max_val * 1.05, max_val * 1.05, max_val))
      #     signif_df <- data.frame(x=(comp[[1]] + comp[[2]]) * .5, y=max_val, label=signif.num(wlk$p.value * len(comp_conf$comp_ind)))
      #   }
      #   
      #   gcomp_bar <- gcomp_bar + geom_line(data=line_df, aes(x=x,y=y)) +
      #                            geom_text(data=signif_df, aes(x=x,y=y,label=label))
      #   
      #   gcomp_box <- gcomp_box + geom_line(data=line_df, aes(x=x,y=y)) +
      #     geom_text(data=signif_df, aes(x=x,y=y,label=label))
      #   
      #   
      #   
      #   wlk <- wilcox.test(pvals_1, pvals_2, paired=F, alternative="greater")
      #   
      #   max_val <- max(pvals_1, pvals_2, na.rm=T)
      #   
      #   if (len(comp_conf$comp_ind) == 1) { 
      #     line_df <- data.frame(x=c(1,1,2,2),
      #                           y=c(max_val, max_val * 1.05, max_val * 1.05, max_val))
      #     signif_df <- data.frame(x=1.5, y=max_val, label=signif.num(wlk$p.value))
      #   } else {
      #     line_df <- data.frame(x=c(comp[[1]],comp[[1]],comp[[2]],comp[[2]]),
      #                           y=c(max_val, max_val * 1.05, max_val * 1.05, max_val))
      #     signif_df <- data.frame(x=(comp[[1]] + comp[[2]]) * .5, y=max_val, label=signif.num(wlk$p.value * len(comp_conf$comp_ind)))
      #   }
      #   
      #   gpval_bar <- gpval_bar + geom_line(data=line_df, aes(x=x,y=y)) +
      #     geom_text(data=signif_df, aes(x=x,y=y,label=label))
      #   
      #   gpval_box <- gpval_box + geom_line(data=line_df, aes(x=x,y=y)) +
      #     geom_text(data=signif_df, aes(x=x,y=y,label=label))
      #                            
      # 
      # }
      
      for (size_name in names(sizes)) {
        tmp_write_path <- sprintf("%s\\%s", write_path, size_name)
        dir.create(tmp_write_path)
        tmp_write_path <- sprintf("%s\\comparisions_of_pvalues_topological", tmp_write_path)
        dir.create(tmp_write_path)
        
        pdf(sprintf("%s\\boxplot_comp_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_conf_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gcomp_box)
        dev.off()
        
        pdf(sprintf("%s\\barplot_comp_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_conf_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gcomp_bar)
        dev.off()
        
        pdf(sprintf("%s\\pvalue_boxplot_comp_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_conf_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gpval_box)
        dev.off()
        
        pdf(sprintf("%s\\pvalue_barplot_comp_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_conf_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gpval_bar)
        dev.off()
      }
    }
  }
}


datasets <- get_datasets_names(get_thirsty_quenched_paths())

mice_groups <- unlist(lapply(str_split(datasets, " "), function(str) {str[[1]]}))


for (metric in unique(m14_m14$topological_similarity$metric)) {
  for (barc in c("betti0", "betti1")) {
    

    
    
    
    vals <- list(`1`=m14_m14,
                 `4`=m14_m14shuffle)
    
    reg_working_df <- m14_m14$topological_similarity[m14_m14$topological_similarity$metric == metric,]
    shuffle_working_df <- m14_m14shuffle$topological_similarity[m14_m14shuffle$topological_similarity$metric == metric,]
    
    
    within <- 
      apply(working_df[,c("path1", "path2")], 
            1,
            function(df) {
              idx1 <- df[1]
              idx2 <- df[2]
              mice_groups[idx1] == mice_groups[idx2]

            })
    
    across <- !within
    
    reg_across_working_df <- reg_working_df[across,]
    shuffle_across_working_df <- shuffle_working_df[across,]
    reg_within_working_df <- reg_working_df[within,]
    shuffle_within_working_df <- shuffle_working_df[within,]
    
    
    assert(all(reg_across_working_df[,1] - shuffle_across_working_df[,1] == 0))
    assert(all(reg_across_working_df[,2] - shuffle_across_working_df[,2] == 0))
    assert(all(reg_within_working_df[,1] - shuffle_within_working_df[,1] == 0))
    assert(all(reg_within_working_df[,2] - shuffle_within_working_df[,2] == 0))
    
    
    reg_across_vec <- ddply(reg_across_working_df, .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    shuffle_across_vec <- ddply(shuffle_across_working_df, .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    reg_within_vec <- ddply(reg_within_working_df, .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    shuffle_within_vec <- ddply(shuffle_within_working_df, .(path1), function(df) {colMeans(df[,c("betti0", "betti1")])})[,barc]
    
    
    boxplot(reg_across_vec,
    shuffle_across_vec,
    reg_within_vec,
    shuffle_within_vec)
    
    

      
      comps_df <- data.frame(similarity=c(reg_across_vec,
                                          shuffle_across_vec,
                                          reg_within_vec,
                                          shuffle_within_vec),
                             x=c(rep("Across", times=len(reg_across_vec)),
                                 rep("Across - Shuffle", times=len(shuffle_across_vec)),
                                 rep("Within", times=len(reg_within_vec)),
                                 rep("Within - Shuffle", times=len(shuffle_within_vec))))

      
      
      
      gcomp_bar <- 
        ggplot(comps_df, aes(x=x,y=similarity)) + 
        #geom_point() + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 
        geom_bar(stat="summary", alpha=.8, width=.45, fill="gray70") +
        geom_errorbar(stat="summary", width=.25, size=.75) + 
        scale_y_continuous(expand=c(0,0)) + 
        
        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(comp_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      gcomp_box <- 
        ggplot(comps_df, aes(x=x,y=similarity)) + 
        geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
                    size=1.5, alpha=.75) + 
        geom_boxplot(aes(),
                     width=0.5, 
                     size=1, 
                     color=adjustcolor("gray65", alpha=1), 
                     fill=adjustcolor("gray80", alpha=0.9),
                     outlier.shape = NA) + 
        theme_classic() +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")  +
        xlab("") +
        ylab(comp_barcodes_labs[[barc]]) +
        theme(text=element_text(size=14, color="black"),
              legend.position="NA")
      
      
      
      across_across <- wilcox.test(reg_across_vec, shuffle_across_vec, paired=T, alternative = "less", correct = F)
      within_within <- wilcox.test(reg_within_vec, shuffle_within_vec, paired=T, alternative = "less", correct = F)
      within_across <- wilcox.test(reg_within_vec, reg_across_vec, paired=F, alternative = "less", correct = F)
      
      max_val <- max(comps_df$similarity)
      
      signif_df <- data.frame(x=c(1.5, 3.5, 2), 
                              y=c(max_val, max_val, 0.35 * max_val), 
                              label=c(signif.num(across_across$p.value),
                                      signif.num(within_within$p.value),
                                      signif.num(within_across$p.value)))
      
      gcomp_box <- 
        gcomp_box + geom_text(data=signif_df, aes(x=x,y=y,label=label))
      
      gcomp_bar <- 
        gcomp_bar + geom_text(data=signif_df, aes(x=x,y=y,label=label))
      
      for (size_name in names(sizes)) {
        tmp_write_path <- sprintf("%s\\%s", write_path, size_name)
        dir.create(tmp_write_path)
        tmp_write_path <- sprintf("%s\\across_within_comparisions", tmp_write_path)
        dir.create(tmp_write_path)
        
        pdf(sprintf("%s\\boxplot_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gcomp_box)
        dev.off()
        
        pdf(sprintf("%s\\barplot_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]])
        
        plot(gcomp_bar)
        dev.off()
        
      }
    }
}



get_mice_names <- function(paths,control=F) {unlist(lapply(str_split(get_datasets_names(paths, control=control), " "), function(s) {s[[1]]}))}


v1por_names <- paste(c(rep("V1", 4), rep("POR", 3)), v1por_names, sep="_")
thirst_names <- get_mice_names(get_thirsty_quenched_paths())
hunger_names <- get_mice_names(get_hungry_sated_paths())
agrp_names <- get_mice_names(get_agrp_paths())
sfo_names <- get_mice_names(get_SFO_paths())
ylgn_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "YlGn")))

all_results_annot <- list(#mHunger_mv1POR=list(n1=hunger_names, n2=v1por_names),
                          #mThirst_mv1POR=list(n1=thirst_names, n2=v1por_names),
                          #m14_mHunger=list(n1=thirst_names, n2=hunger_names),
                          #m14_mSFO=list(n1=thirst_names, n2=sfo_names),
                          #m14_m14shuffle=list(n1=thirst_names, n2=thirst_names),
                          #m14_m14=list(n1=thirst_names, n2=thirst_names),
                          mHunger_mHunger=list(n1=hunger_names, hunger_names),
                          mHunger_mAGRP=list(n1=hunger_names, agrp_names),
                          mAGRP_mAGRP=list(n1=agrp_names, agrp_names))


for (metric in unique(mHunger_mHunger$topological_similarity$metric)) {
  for (barc in c("betti0", "betti1")) {
    for (comp_name in names(all_results)) {
      paths1 <- unique(all_results[[comp_name]]$topological_similarity$path1)
      paths2 <- unique(all_results[[comp_name]]$topological_similarity$path2)
      
      sim_mat <- matrix(rep(0, times=len(paths1) * len(paths2)), nrow=len(paths1))
      
      for (idx_i in 1:len(paths1)) {
        for(idx_j in 1:len(paths2)) {
          sim_mat[idx_i, idx_j] <- all_results[[comp_name]]$topological_similarity[all_results[[comp_name]]$topological_similarity$path1 == idx_i &
                                                                                     all_results[[comp_name]]$topological_similarity$path2 == idx_j &
                                                                                     all_results[[comp_name]]$topological_similarity$metric == metric, barc]
        }
      }
    
    pval_mat <- all_results[[comp_name]][[pval_name]][[metric]]
      
    row_annotation_df <- data.frame(Mice1=all_results_annot[[comp_name]]$n1)
    col_annotation_df <- data.frame(Mice2=all_results_annot[[comp_name]]$n2)
    rownames(row_annotation_df) <- 1:len(all_results_annot[[comp_name]]$n1)
    rownames(col_annotation_df) <- 1:len(all_results_annot[[comp_name]]$n2)
    
    rownames(sim_mat) <- 1:len(all_results_annot[[comp_name]]$n1)
    colnames(sim_mat) <- 1:len(all_results_annot[[comp_name]]$n2)
    rownames(pval_mat) <- 1:len(all_results_annot[[comp_name]]$n1)
    colnames(pval_mat) <- 1:len(all_results_annot[[comp_name]]$n2)
    
    
    col_mice_1 <- ylgn_cg(len(unique(all_results_annot[[comp_name]]$n1)))
    names(col_mice_1) <- unique(all_results_annot[[comp_name]]$n1)
    col_mice_2 <- ylgn_cg(len(unique(all_results_annot[[comp_name]]$n2)))
    names(col_mice_2) <- unique(all_results_annot[[comp_name]]$n2)
    
    annotation_colors <- list(Mice1=col_mice_1, Mice2=col_mice_2)

    
    pval_name <- ifelse(barc == "betti0", "pvalue_b0", "pvalue_b1")
    
    
    ph_similarity <- 
      pheatmap(sim_mat, cluster_rows=T, cluster_cols=T, border_col=NA,
               annotation_row = row_annotation_df,
               annotation_col = col_annotation_df,
               show_rownames = F,
               show_colnames = F,
               treeheight_row = 20,
               treeheight_col = 20,
               annotation_colors = annotation_colors)
    
    ph_pval <- 
      pheatmap(pval_mat, cluster_rows=T, cluster_cols=T, border_col=NA,
               annotation_row = row_annotation_df,
               annotation_col = col_annotation_df,
               show_rownames = F,
               show_colnames = F,
               annotation_colors = annotation_colors,
               treeheight_row = 20,
               treeheight_col = 20,
               breaks=c(0,0.05,0.1,0.15,0.25,0.5,0.75,1), col=rev(rdylbu_cg(7)))
    
    
    ph_similarity_unc <- 
      pheatmap(sim_mat, cluster_rows=F, cluster_cols=F, border_col=NA,
               annotation_row = row_annotation_df,
               annotation_col = col_annotation_df,
               show_rownames = F,
               show_colnames = F,
               treeheight_row = 20,
               treeheight_col = 20,
               annotation_colors = annotation_colors)
    
    ph_pval_unc <- 
      pheatmap(pval_mat, cluster_rows=F, cluster_cols=F, border_col=NA,
               annotation_row = row_annotation_df,
               annotation_col = col_annotation_df,
               show_rownames = F,
               show_colnames = F,
               annotation_colors = annotation_colors,
               treeheight_row = 20,
               treeheight_col = 20,
               breaks=c(0,0.05,0.1,0.15,0.25,0.5,0.75,1), col=rev(rdylbu_cg(7)))
    

      for (size_name in names(sizes)) {
        tmp_write_path <- sprintf("%s\\%s", write_path, size_name)
        dir.create(tmp_write_path)
        tmp_write_path <- sprintf("%s\\heatmaps", tmp_write_path)
        dir.create(tmp_write_path)
  
        pdf(sprintf("%s\\similarity_heatmap_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]] * 1.25)
  
        plot(ph_similarity[[4]])
        dev.off()
  
        pdf(sprintf("%s\\pval_heatmap_%s_metric_%s_barc_%s.pdf",
                    tmp_write_path,
                    comp_name,
                    metric,
                    barc),
            height=sizes[[size_name]],
            width=sizes[[size_name]] * 1.25)
  
        plot(ph_pval[[4]])
        dev.off()

    }
    }
  }
}

output_name = "hunger_v1_topological_comparision.Rda"
load(sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), verbose=T)
hv1 <- final_res

metric="wasserstein"
dfa <- rbind(#th$topological_similarity[th$topological_similarity$metric == metric,],
             tv1$topological_similarity[tv1$topological_similarity$metric == metric,],
             tt$topological_similarity[tt$topological_similarity$metric == metric,])
             #hv1$topological_similarity[hv1$topological_similarity$metric == metric,])

cmps <- c(#nrow(th$topological_similarity[th$topological_similarity$metric == metric,]),
          nrow(tv1$topological_similarity[tv1$topological_similarity$metric == metric,]),
          nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]))
          #nrow(hv1$topological_similarity[hv1$topological_similarity$metric == metric,]))

dfa$group = rep(c("T-H", "T-V1", "T-T", "H-V1"), cmps)
dfa <- as.data.frame(dfa)
dfa$betti0 <- as.numeric(dfa$betti0)
dfa$betti1 <- as.numeric(dfa$betti1)
dff <- dfa[!(dfa$path1 == dfa$path2 & dfa$group == "T-T"),]

ggplot(dff, aes(x=group, y=as.numeric(betti0)), group=group) +
  #geom_jitter(alpha=.1) + 
  geom_violin(aes(fill=group), color=NA) + 
  stat_summary(size=1) + 
  theme_light() + 
  base_plot_theme


fraction_hist_df <- 
ddply(dff, .(group), 
      function(group_df) {
        hb0 <-  hist(group_df$betti0, breaks=seq(min(dff$betti0),
                                                 max(dff$betti0),
                                                 length.out=31), plot=F)
        
        hb1 <-  hist(group_df$betti1, breaks=seq(min(dff$betti1),
                                                 max(dff$betti1),
                                                 length.out=31), plot=F)
        
        
        hist_df <- 
        data.frame(frac_betti0=hb0$counts / sum(hb0$counts), 
                   breaks_betti0=hb0$breaks[-1], 
                   frac_betti1=hb1$counts / sum(hb1$counts), 
                   breaks_betti1=hb1$breaks[-1], 
                   group=rep(group_df$group[1],len(hb0$counts)))
        print(sum(hist_df$frac_betti0))
        
        return(hist_df)
      })

# 
 ggplot(fraction_hist_df, aes(y=frac_betti0, x=breaks_betti0, group=group)) + 
   geom_bar(aes(fill=group), stat="identity", size=2, alpha=.5,  position=position_nudge()) + 
   #geom_line(aes(color=group), size=.2, alpha=1, position=position_dodge(.0001)) + 
   theme_light() +
   scale_y_continuous(expand=c(0,0)) +
   base_plot_theme

 
 
 

 
 output_name = "hunger_v1_topological_comparision.Rda"
 load(sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), verbose=T)
 hv1 <- final_res
 
 metric="wasserstein"
 dfa <- rbind(th$topological_similarity[th$topological_similarity$metric == metric,],
              tv1$topological_similarity[tv1$topological_similarity$metric == metric,],
              tt$topological_similarity[tt$topological_similarity$metric == metric,],
              hv1$topological_similarity[hv1$topological_similarity$metric == metric,])
 
 cmps <- c(nrow(th$topological_similarity[th$topological_similarity$metric == metric,]),
           nrow(tv1$topological_similarity[tv1$topological_similarity$metric == metric,]),
           nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]),
           nrow(hv1$topological_similarity[hv1$topological_similarity$metric == metric,]))
 
 dfa$group = rep(c("T-H", "T-V1", "T-T", "H-V1"), cmps)
 dfa <- as.data.frame(dfa)
 dfa$betti0 <- as.numeric(dfa$betti0)
 dfa$betti1 <- as.numeric(dfa$betti1)
 dff <- dfa[!(dfa$path1 == dfa$path2 & dfa$group == "T-T"),]