
figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\figures_structure_dynamics\\", base_output_path) 

dir.create(base_output_path)
dir.create(output_path)

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_HS <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
metadata_control <- across_mice_decoding_build_metadata(get_control_paths(), control=T)
metadata_SFO <- across_mice_decoding_build_metadata(get_SFO_paths())
metadata_AGRP <- across_mice_decoding_build_metadata(get_agrp_paths())
metadata_hypertonicsaline <- across_mice_decoding_build_metadata(get_hypertonic_saline_paths())

plot_TT_trial_similarities <- function()
{
  base_figure_path <- sprintf("%s\\TT_trial_similarities", output_path)
  dir.create(base_figure_path)
  
  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
  metadata <- metadata_TT

  trial_annot <- c("Hit"="3_0",
                   "Miss"="3_1",
                   "CR (A)"="4_0",
                   "FA (A)"="4_1",
                   "CR (N)"="5_0",
                   "FA (N)"="5_1")
  
  
  final_all_df <- data.frame()
  
  all_pooled_similarity_df <- list()
  
  for (i in 1:1) {
  
  #mpc <- -sample(1:7, 1)
  pooled_similarity_df <- list()
  
  for (idx in c(1:len(paths))) {
    m_obj <- metadata[[idx]]
    dataset_trial_annot <- paste(m_obj$annot_df[,1], m_obj$annot_df[,2], sep="_")
    
    #trans_to_use <- c(0,sample(101:107, 7))
    trans_to_use <- c(mpc, sample(c(1:7), 7))
    
    for (annot_name in names(trial_annot)) {
      tt_annot_sig <- trial_annot[[annot_name]]
      tt_indices <- which(dataset_trial_annot == tt_annot_sig)
      tt_mat <- m_obj$trials_mat[tt_indices,]
      #sw <- swap_labels(tt_mat, trans = trans_to_use)
      
      sim <- cosine(t(tt_mat))
      #sim <- jaccard_mat_sim(tt_mat)
      
      pooled_similarity_df[[annot_name]] <- c(pooled_similarity_df[[annot_name]], sim[lower.tri(sim)])
      all_pooled_similarity_df[[annot_name]] <- c(all_pooled_similarity_df[[annot_name]], sim[lower.tri(sim)])
      
    }
    
    print(idx)
  }
  
  N_breaks <- 20
  
  hist_df <- 
  lapply(names(pooled_similarity_df),
         function(lst_name) {
           lst <- pooled_similarity_df[[lst_name]]
           
           h <- hist(lst, plot=F, seq(min(unlist(pooled_similarity_df), na.rm=T),
                                      max(unlist(pooled_similarity_df), na.rm=T),
                                      length.out=(N_breaks + 1)))
           
           df <- 
           data.frame(breaks=h$breaks,
                      freq=c(0, h$counts/sum(h$counts)),
                      group=rep(lst_name, times=(N_breaks + 1)))
           
           return(df)
         })
  
  all_df <- do.call(rbind,hist_df)
  
  final_all_df <- rbind(final_all_df,
                        all_df)
  }
  
    #g <- 
  ggplot(final_all_df, aes(x=breaks, y=freq)) + 
    #geom_bar(aes(fill=group), color=adjustcolor("gray30", alpha=.85), size=.01, alpha=.2, width=.01,  position=position_nudge(), stat="summary") + 
    geom_line(aes(color=group),  size=.85, alpha=.85, stat="summary") + 
    #geom_ribbon(aes(fill=group),  color=NA, alpha=.25, stat="summary") +
    big_text_base_plot_theme + 
    scale_y_continuous(expand=c(0,0)) + 
    xlab("Cosine similarity")  +
    ylab("Fraction") + 
    scale_color_brewer(palette="Spectral") +
    scale_fill_brewer(palette="Spectral")
  
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_trial_similarity_hists.pdf",base_figure_path, size_name), height=size, width=size)
    plot(g)
    dev.off()
    
    pdf(sprintf("%s\\pooled_%s_trial_similarity_hists_legend.pdf",base_figure_path, size_name), height=size, width=size)
    plot(g +  theme(legend.position = "top"))
    dev.off()
  }
  
  statistics_df <- 
        lapply(names(all_pooled_similarity_df)[-1],
                function(lst_name) {
                  print(lst_name)
                  
                  ks <- wilcox.test(pooled_similarity_df$Hit, pooled_similarity_df[[lst_name]], alternative="greater")
                  
                  res_df <- 
                  data.frame(Comp=sprintf("%s - %s", "Hit", lst_name),
                             Statistics=ks$statistic,
                             Pvalue=ks$p.value,
                             Corrected=ks$p.value * (len(pooled_similarity_df) - 1),
                             Method=ks$method,
                             Alternative=ks$alternative,
                             DN=ks$data.name,
                             Mean_Hit=mean(pooled_similarity_df$Hit),
                             Mean_compared=mean(pooled_similarity_df[[lst_name]]),
                             Median_Hit=median(pooled_similarity_df$Hit),
                             Median_compared=median(pooled_similarity_df[[lst_name]]),
                             SD_Hit=sd(pooled_similarity_df$Hit),
                             SD_compared=sd(pooled_similarity_df[[lst_name]]),
                             SEM_Hit=sem(pooled_similarity_df$Hit),
                             SEM_compared=sem(pooled_similarity_df[[lst_name]]),
                             N_Hit=len(pooled_similarity_df$Hit),
                             N_compared=len(pooled_similarity_df[[lst_name]]))
                  
                  res_df$Signif <- signif.num(res_df$Pvalue)
                  res_df$Corrected_signif <- signif.num(res_df$Corrected)
                  
                  return(res_df)
                  
                })
  
  statistics_df <- do.call(rbind, statistics_df)
  write.csv(file=sprintf("%s//wilcox_pooled_statistics.csv", statistics_figure_path), statistics_df)
}

swap_labels <- function(mt, trans=NA) {
  
  original <- c(-1,1:7)
  
  if (all(is.na(trans))) {
    trans <- sample(original, len(original))
  } else {
    print(trans)
  }
  
  mt <- 
  apply(mt, 1, 
        function(r) {
          
          old_r <- r
          new_r <- old_r
          
          for (ix in 1:len(original)) {
            new_r[old_r == original[ix]] <- trans[ix] 
          }
    
      return(new_r)
  })
  
  
  return(t(mt))
}

plot_control_trial_dynamics <- function()
{
  paths <- get_control_paths() 
  datasets_names <- get_datasets_names(paths, sep="_", control = T)
  metadata <- metadata_control
  
  base_figure_path <- sprintf("%s\\control_trial_dynamics", output_path)
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

plot_HS_trial_dynamics <- function()
{
  paths <- get_hungry_sated_paths() 
  datasets_names <- get_datasets_names(paths, sep="_")
  metadata <- clean_metadata_HS(metadata_HS)
  
  base_figure_path <- sprintf("%s\\HS_trial_dynamics", output_path)
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

plot_TT_trial_dynamics <- function()
{
  paths <- get_thirsty_quenched_paths() 
  datasets_names <- get_datasets_names(paths, sep="_")
  metadata <- metadata_TT
  
  base_figure_path <- sprintf("%s\\TT_trial_dynamics", output_path)
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

plot_TT_decoding_translated <- function()
{
  base_figure_path <- sprintf("%s\\TT_trial_dynamics_translated", output_path)
  dir.create(base_figure_path)
  datasets_names <- get_datasets_names(paths, sep="_")
  metadata <- metadata_TT
  
  names(metadata) <- datasets_names
  
  for (m1_name in names(metadata)) {
    for (m2_name in names(metadata)) {
      
      if (m1_name == m2_name) {
        next
      }
      
      m1 <- metadata[[m1_name]]
      m2 <- metadata[[m2_name]]
      
      cluster_color_label <- spec_cg(len(unique(m1$cluster_mat$labs)))
      names(cluster_color_label) <- c(-1, 1:(len(unique(m1$cluster_mat$labs)) - 1))
      m1_m2 <- matrix(rep(0, times=len(m1$trials_mat)), nrow=nrow(m1$trials_mat))
      
      for (trial_idx in 1:nrow(m1$trials_mat)) {
        trial <- m1$trials_mat[trial_idx,]
        
        for (clust_idx in 1:len(m1$mask)){
          m1_m2[trial_idx, trial == m1$mask[clust_idx]] = as.numeric(m2$mask[clust_idx])
        }
      }
      
      ord_m1 <- order(paste(m1$annot_df$trialtype, m1$annot_df$response))
      ord_m2 <- order(paste(m2$annot_df$trialtype, m2$annot_df$response))
      
      ph1 <- pheatmap(m1$trials_mat[ord_m1,][m1$annot_df[ord_m1,1] == 3,], 
                      annotation_row = m1$annot_df[ord_m1,][m1$annot_df[ord_m1,1] == 3,c(2,1)],
                      cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F,
                      border_col=NA,
                      legend=F,
                      annotation_names_row = F,
                      annotation_legend = F,
                      col=cluster_color_label)
      
      rownames(m1_m2) <- rownames(m1$annot_df)
      ph2 <- pheatmap(m1_m2[ord_m1,][m1$annot_df[ord_m1,1] == 3,], 
                      annotation_row = m1$annot_df[ord_m1,][m1$annot_df[ord_m1,1] == 3,c(2,1)],
                      cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F,
                      border_col=NA,
                      legend=F,
                      annotation_names_row = F,
                      annotation_legend = F,
                      col=cluster_color_label)
      
      ph3 <- pheatmap(m2$trials_mat[ord_m2,][m2$annot_df[ord_m2,1] == 3,], 
                      annotation_row = m2$annot_df[ord_m2,][m2$annot_df[ord_m2,1] == 3, c(2,1)],
                      cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F,
                      border_col=NA,
                      legend = F,
                      annotation_names_row = F,
                      annotation_legend = F,
                      col=cluster_color_label)
      
      
      gall <- do.call(plot_grid, list(ph1[[4]], ph2[[4]], ph3[[4]], nrow=1))
      
      
      png(sprintf("%s//%s_to_%s.png",
                  base_figure_path,
                  m1_name,
                  m2_name),
          unit="in",
          res=500,
          height=2.5,
          width=7.5)
      plot(gall)
      dev.off()
      
      pdf(sprintf("%s//%s_to_%s.pdf",
                  base_figure_path,
                  m1_name,
                  m2_name),
          height=2.5,
          width=7.5)
      plot(gall)
      dev.off()
      
    }
  }
}

plot_TT_decoding_results <- function()
{
  base_figure_path <- sprintf("%s\\TT_decoding_results", output_path)
  dir.create(base_figure_path)
  
  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path) 
  
  average_across_datasets = F
  
  
  load(file=sprintf("%s\\data\\figure_3\\confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\pvalues_decoding.Rda", base_output_path), verbose=T)
  
  confusion_matrix_list <- 
  lapply(confusion_matrix_list,
         function(mt) {
           for (row_idx in 1:nrow(mt)) {
             mt[row_idx,row_idx] <- NA
           }
           return(mt)
         })
  
  shuffle_matrix_list <- 
  lapply(shuffle_matrix_list,
         function(mt) {
           for (row_idx in 1:nrow(mt)) {
             mt[row_idx,row_idx] <- NA
           }
           return(mt)
         })
  
  mice_names_for_path <- unlist(lapply(str_split(get_datasets_names(get_thirsty_quenched_paths()), " "), function(lst) {lst[[1]]}))
  names(mice_names_for_path) <- 1:14
  
  
  statistics_df <- data.frame()
  metrics_to_use <- c("cosine")  
    
  for (met in metrics_to_use) { 
    matrices_to_use <- names(confusion_matrix_list)[grepl(metrics_to_use[1], names(confusion_matrix_list))]
    matrices_to_use <- unlist(str_split(matrices_to_use, "_"))[1:(len(matrices_to_use) * 2) %% 2 == 1]
    
    
    for (mat_name in matrices_to_use) {
      
      shuff_mat <- shuffle_matrix_list[[sprintf("%s_%s", mat_name, met)]]
      data_mat <- confusion_matrix_list[[sprintf("%s_%s", mat_name, met)]]
      
      melted_shuff_mat <- melt(shuff_mat)
      melted_data_mat <- melt(data_mat)
      melted_shuff_mat$Mice1 <- mice_names_for_path[melted_shuff_mat$Var1]
      melted_shuff_mat$Mice2 <- mice_names_for_path[melted_shuff_mat$Var2]
      melted_data_mat$Mice1 <- mice_names_for_path[melted_data_mat$Var1]
      melted_data_mat$Mice2 <- mice_names_for_path[melted_data_mat$Var2]
      
      melted_data_mat$original_group <- "Data"
      melted_shuff_mat$original_group <- "Shuffle"
      
      mat_final <- rbind(melted_data_mat, melted_shuff_mat)
      mat_final$group <- paste(mat_final$original_group, ifelse(as.numeric(mat_final$Mice1 == mat_final$Mice2) > 0, "Within", "Across"), sep="_")
      mat_final <- mat_final[mat_final$Var1 != mat_final$Var2,]
      
      colnames(mat_final) <- c("Path1", "Path2", "Accuracy", "Mice1", "Mice2", "Original_group", "Group")
      
      if (average_across_datasets) {
       df_final <-
         ddply(mat_final,
      
               .(Group),
               function(group_df) {
                 res <-
                   ddply(group_df,
                         .(Path1),
                         function(mice_path_df) {
                           return(mean(mice_path_df[,"Accuracy"], na.rm=T))
                         })
      
               })

      colnames(df_final) <- c("Group", "Idx", "Accuracy")
      
      } else {
        df_final <- mat_final
      }
      
      plot_groups <- c("Data_Within", "Shuffle_Within", "Data_Across", "Shuffle_Across")
      
      g <- 
      ggplot(df_final, 
             aes(x=factor(Group, levels=plot_groups), 
                 y=Accuracy), 
             group=Group)
      
      if (average_across_datasets) {
        g <- g + geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)
      } else {
        g <- g + geom_jitter(position=position_jitterdodge(.25), aes(fill=Group))
      }
        g <- g +
        scale_shape_manual(values=18) +
        geom_boxplot(fill=NA, width=.5, outlier.shape=NA)  + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Decoding accuracy (%)") + 
        ylim(0,1) +
        xlab("")  +
        ggtitle(mat_name) +
        theme(plot.title = element_text(size=9))
      
      for (size_name in names(a4_sizes)) {
        
        if (average_across_datasets) {
          dir.create(sprintf("%s\\AVERAGED_OVER_MICE", base_figure_path))  
        } else {
          dir.create(sprintf("%s\\ALL_PAIRWISE", base_figure_path))  
        }
        
        size = a4_sizes[[size_name]]
        pdf(sprintf("%s\\%s%s_%s_%s.pdf",base_figure_path, ifelse(average_across_datasets, "AVERAGED_OVER_MICE\\", "ALL_PAIRWISE\\"),size_name, met, mat_name), height=size, width=size)
        plot(g)
        dev.off()
  
      }
      
      
      ##### Statistics
      Data_Within <- df_final[df_final$Group == "Data_Within",]
      Shuffle_Within <- df_final[df_final$Group == "Shuffle_Within",]
      Data_Across <- df_final[df_final$Group == "Data_Across",]
      Shuffle_Across <- df_final[df_final$Group == "Shuffle_Across",]
      
      
      if (average_across_datasets) { 
        data_within_order <-  order(Data_Within[,"Idx"])
        data_across_order <-  order(Data_Across[,"Idx"])
        shuffle_within_order <-  order(Shuffle_Within[,"Idx"])
        shuffle_across_order <-  order(Shuffle_Across[,"Idx"])
      } else {
        data_within_order <-  order(paste(Data_Within[,"Path1"], Data_Within[,"Path2"]))
        data_across_order <-  order(paste(Data_Across[,"Path1"], Data_Across[,"Path2"]))
        shuffle_within_order <-  order(paste(Shuffle_Within[,"Path1"], Shuffle_Within[,"Path2"]))
        shuffle_across_order <-  order(paste(Shuffle_Across[,"Path1"], Shuffle_Across[,"Path2"]))
      }
      
      comp_configurations <- list(list(name = "within_wilc", 
                                       group_A = Data_Within[data_within_order,"Accuracy"],
                                       group_B = Shuffle_Within[shuffle_within_order ,"Accuracy"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "across_wilc", 
                                       group_A = Data_Across[data_across_order,"Accuracy"],
                                       group_B = Shuffle_Across[shuffle_across_order,"Accuracy"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "across_within", 
                                       group_A = Data_Within[,"Accuracy"],
                                       group_B = Data_Across[,"Accuracy"],
                                       paired= F,
                                       alternative="two.sided")) 

      for (stat_test_configuration in comp_configurations) {
       
        wilc <- wilcox.test(stat_test_configuration$group_A,
                            stat_test_configuration$group_B,
                            paired=stat_test_configuration$paired,
                            correct = F,
                            alternative=stat_test_configuration$alternative) 
        
        
        stat_df <- data.frame(conf_name=stat_test_configuration$name,
                              statistic=wilc$statistic,
                              pval=wilc$p.value,
                              method=wilc$method,
                              alternative=wilc$alternative,
                              mean_group_A=mean(stat_test_configuration$group_A),
                              mean_group_B=mean(stat_test_configuration$group_B),
                              sd_group_A=sd(stat_test_configuration$group_A),
                              sd_group_B=sd(stat_test_configuration$group_B),
                              sem_group_A=sem(stat_test_configuration$group_A),
                              sem_group_B=sem(stat_test_configuration$group_B),
                              N_group_A=len(stat_test_configuration$group_A),
                              N_group_B=len(stat_test_configuration$group_B),
                              corrected_pvalue=wilc$p.value * len(comp_configurations),
                              signif=signif.num(wilc$p.value),
                              signif_corrected=signif.num(wilc$p.value * len(comp_configurations)),
                              metric_name=met,
                              matrix_name=mat_name)
        
        statistics_df <- rbind(statistics_df,
                               stat_df)
      }
    }
  }
  
  
  write.csv(file=sprintf("%s//statistics//%sstatistics.csv", base_figure_path, ifelse(average_across_datasets, "AVERAGED_OVER_MICE_", "")),
            statistics_df)
}

plot_TH_HC_decoding_results <- function()
{
  base_figure_path <- sprintf("%s\\TH_HC_decoding_results", output_path)
  dir.create(base_figure_path)
  
  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path) 

  
  
  load(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_pvalues_decoding.Rda", base_output_path), verbose=T)
  
  th_confusion_matrix_list <- confusion_matrix_list
  th_shuffle_matrix_list <- shuffle_matrix_list
  th_pvalues_all <- pvalues_all
  
  load(file=sprintf("%s\\data\\figure_3\\control_hunger_new_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\control_hunger_new_shuffle_confusion_matrices.Rda", base_output_path), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\control_hunger_new_pvalues_decoding.Rda", base_output_path), verbose=T)
  
  ch_confusion_matrix_list <- confusion_matrix_list
  ch_shuffle_matrix_list <- shuffle_matrix_list
  ch_pvalues_all <- pvalues_all
  
  confusion_matrix_list <- 
    lapply(confusion_matrix_list,
           function(mt) {
             for (row_idx in 1:nrow(mt)) {
               mt[row_idx,row_idx] <- NA
             }
             return(mt)
           })
  
  shuffle_matrix_list <- 
    lapply(shuffle_matrix_list,
           function(mt) {
             for (row_idx in 1:nrow(mt)) {
               mt[row_idx,row_idx] <- NA
             }
             return(mt)
           })
  
  mice_names_for_path <- unlist(lapply(str_split(get_datasets_names(get_thirsty_quenched_paths()), " "), function(lst) {lst[[1]]}))
  names(mice_names_for_path) <- 1:14
  
  
  all_matrices <- list(ch_confusion_matrix_list$Hit_cosine,
                       ch_shuffle_matrix_list$Hit_cosine,
                       th_confusion_matrix_list$Hit_cosine,
                       th_shuffle_matrix_list$Hit_cosine)
  
  final_mat <- 
    do.call(rbind,
            lapply(all_matrices, function(mt) {
            melted_mat <- melt(mt)
            ddply(melted_mat, .(Var1), function(vrdf) {mean(vrdf[,"value"])})}))
  
  final_mat <- as.data.frame(final_mat)
  final_mat$group <- rep( c("V-I", "V-I (shuff)", "T-H", "T-H (shuff)"), times=c(7,7,14,14))
  
  colnames(final_mat) <- c("#", "Accuracy", "Group")
  df_final <- final_mat
  
  group_levels <- c("T-H", "T-H (shuff)", "V-I", "V-I (shuff)")
  g <- 
    ggplot(df_final, 
           aes(x=factor(Group, levels=group_levels), 
               y=Accuracy), 
           group=Group) + 
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75) + 
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape=NA)  + 
    theme_light() + 
    big_text_base_plot_theme + 
    ylab("Decoding accuracy (%)") + 
    ylim(0,1) +
    xlab("")  +
    ggtitle("Hits") +
    theme(plot.title = element_text(size=9))
  
  for (size_name in names(a4_sizes)) {
    
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_th_hc_hit.pdf",base_figure_path, size_name), height=size, width=1.5 * size)
    plot(g)
    dev.off()
    
  }
  
    
    dt <- dunn.test::dunn.test(df_final$Accuracy, df_final$Group, method="bonferroni")
    stat_df <- as.data.frame(dt)
    stat_df$signif.adj <- signif.num(stat_df$P.adjusted)
    stat_df$signif <- signif.num(stat_df$P)
  
    values_df <- ddply(df_final, .(Group), function(group_df) {acc_vec <- group_df[,"Accuracy"];
                                                               c(mean(acc_vec), sd(acc_vec), sem(acc_vec), len(acc_vec))})
    colnames(values_df) <- c("mean", "sd", "sem", "N")
    
  write.csv(file=sprintf("%s//statistics//statistics.csv", base_figure_path),
            stat_df)
  write.csv(file=sprintf("%s//statistics//values.csv", base_figure_path),
            values_df)
}

plot_TT_decoding_pvalue_distributions <- function()
{
  base_figure_path <- sprintf("%s\\TT_decoding_pvalue_distributions", output_path)
  dir.create(base_figure_path)

  load(file=sprintf("%s\\data\\figure_3\\pvalues_decoding.Rda", base_output_path), verbose=T)
  
  mt <- matrix(1:14 ** 2, nrow=14)
  rownames(mt) <- 1:14
  colnames(mt) <- 1:14
  melted_mt <- melt(t(mt))
  mice_names_for_path <- unlist(lapply(str_split(get_datasets_names(get_thirsty_quenched_paths()), " "), function(lst) {lst[[1]]}))
  names(mice_names_for_path) <- 1:14
  
  
  melted_mt$Mice1 <- mice_names_for_path[melted_mt$Var1]
  melted_mt$Mice2 <- mice_names_for_path[melted_mt$Var2]
  
  melted_mt <- melted_mt[melted_mt$Var1 != melted_mt$Var2,]
  across_days_indices <- melted_mt$Mice1 == melted_mt$Mice2
  across_mice_indices <- melted_mt$Mice1 != melted_mt$Mice2
  
  metrics_to_use <- c("cosine")  
  
  
  values_df_all <- data.frame()
  
  for (met in metrics_to_use) { 
    
    across_days <- pvalues_all[[met]][across_days_indices,]
    across_mice <- pvalues_all[[met]][across_mice_indices,]
    
    
    # VALUES & STATISTICS
    values_df  <- data.frame()
    values_df <- rbind(values_df,
                       c(colMeans(1-across_days), 
                         apply(1-across_days, 2, sd), 
                         apply(1-across_days, 2, sem)))
    
    values_df <- rbind(values_df,
                       c(colMeans(1-across_mice), 
                         apply(1-across_mice, 2, sd), 
                         apply(1-across_mice, 2, sem)))
    
    values_df <- cbind(values_df,
                    `N_Datasets`=c(len(unique(melted_mt[across_days_indices,"Var1"])),
                       len(unique(melted_mt[across_mice_indices,"Var1"]))))
    values_df <- cbind(values_df,
                       `N_Pairwise`=c(nrow(across_days), nrow(across_mice)))
    
    values_df$metric = met
    
    colnames(values_df)[1:(ncol(across_days) * 3)] <- paste(rep(c("Mean", "SD", "Sem"), each=ncol(across_days)), colnames(across_days))
    
    values_df_all <- rbind(values_df_all,
                           values_df)
    
    melted_across_days <- melt(across_days)
    melted_across_mice <- melt(across_mice)
    
    colnames(melted_across_days) <- c("#", "Trial", "Pvalue")
    colnames(melted_across_mice) <- c("#", "Trial", "Pvalue")
    
    gviolinacrossmice <- 
    ggplot(melted_across_mice, aes(x=Trial, y=Pvalue, group=Trial)) +
      geom_violin(aes(fill=Trial),color=NA, fill="gray60") + 
      geom_jitter(aes(fill=Trial), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
      #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
      stat_summary(size=1) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Pvalue") + 
      xlab("") 
    
    gviolinacrossdays <- 
      ggplot(melted_across_days, aes(x=Trial, y=Pvalue, group=Trial)) +
      geom_violin(aes(fill=Trial),color=NA, fill="gray60") + 
      geom_jitter(aes(fill=Trial), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
      #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
      stat_summary(size=1) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Pvalue") + 
      xlab("") 
    
      
      
      
      across_mice_pvalue_hist_df <- 
        ddply(melted_across_mice, .(Trial), 
              function(group_df) {
                
                
                h <-  hist(group_df$Pvalue, breaks=seq(0,1,length.out=20), plot=F)

                hist_df <- 
                  data.frame(frac=c(0, h$counts / sum(h$counts)), 
                             breaks=h$breaks)
                
                
                return(hist_df)
              })
      
      
      glineacrossmice <- 
      ggplot(across_mice_pvalue_hist_df, aes(x=breaks, y=frac)) + 
        geom_line(aes(color=Trial),  size=.85, alpha=.85) + 
        big_text_base_plot_theme + 
        scale_y_continuous(expand=c(0,0)) + 
        xlab("Pvalue")  +
        ylab("Fraction") + 
        scale_color_brewer(palette="Spectral") +
        scale_fill_brewer(palette="Spectral")
      
      
      across_days_pvalue_hist_df <- 
        ddply(melted_across_days, .(Trial), 
              function(group_df) {
                
                
                h <-  hist(group_df$Pvalue, breaks=seq(0,1,length.out=20), plot=F)
                
                hist_df <- 
                  data.frame(frac=c(0, h$counts / sum(h$counts)), 
                             breaks=h$breaks)
                
                
                return(hist_df)
              })
      
      
      glineacrossdays <- 
        ggplot(across_days_pvalue_hist_df, aes(x=breaks, y=frac)) + 
        geom_line(aes(color=Trial),  size=.85, alpha=.85) + 
        big_text_base_plot_theme + 
        scale_y_continuous(expand=c(0,0)) + 
        xlab("Pvalue")  +
        ylab("Fraction") + 
        scale_color_brewer(palette="Spectral") +
        scale_fill_brewer(palette="Spectral")
      
      
      
      for (size_name in names(a4_sizes)) {
        size = a4_sizes[[size_name]]
        pdf(sprintf("%s\\%s_across_mice_violin.pdf",base_figure_path, size_name), height=size, width=size); plot(gviolinacrossmice); dev.off()
        pdf(sprintf("%s\\%s_across_days_violin.pdf",base_figure_path, size_name), height=size, width=size); plot(gviolinacrossdays); dev.off()
        
        pdf(sprintf("%s\\%s_across_mice_line.pdf",base_figure_path, size_name), height=size, width=size); plot(glineacrossmice); dev.off()
        pdf(sprintf("%s\\%s_across_days_line.pdf",base_figure_path, size_name), height=size, width=size); plot(glineacrossdays); dev.off()
        
        pdf(sprintf("%s\\%s_across_mice_line_legend.pdf",base_figure_path, size_name), height=size, width=size); plot(glineacrossmice  + theme(legend.position = "top")); dev.off()
        pdf(sprintf("%s\\%s_across_days_line_legend.pdf",base_figure_path, size_name), height=size, width=size); plot(glineacrossdays  + theme(legend.position = "top")); dev.off()
      }
      
      
      write.csv(file=sprintf("%s//values.csv", base_figure_path), values_df)
  }
}

across_mice_decoding_all_trials_decode <- function()
{
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
  
  Metrics_objective <- list(cosine=which.max,
                            euc=which.min,
                            lv=which.min,
                            jaccard=which.max,
                            cosineqg=which.max,
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
  true_vs_decoded_all <- list()
  for (metric_name in names(Metrics)) {
    pvalues_all[metric_name] <- c()
    hits_decoding[[metric_name]] <- list()
    true_vs_decoded_all[[metric_name]] <- list()
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
          
          if (shuff_i == -1) {
            true_vs_decoded_all[[metric_name]] <- 
              append(true_vs_decoded_all[[metric_name]],
                     list(list(train=idx_1, 
                               test=idx_2,
                               decoded=decoded_vec[[metric_name]],
                               gt=trials_labels_mice_2)))
            
            rs3 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 3]
            rs3f <- sum(rs3==0)/len(rs3)
            
            rs5 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 5]
            rs5f <- sum(rs5==0)/len(rs5)
            print(sprintf("%.3f - %.3f (FAs)", rs3f, rs5f))
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
  
  save(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_confusion_matrices.Rda", base_output_path), confusion_matrix_list)
  save(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_shuffle_confusion_matrices.Rda", base_output_path), shuffle_matrix_list)
  save(file=sprintf("%s\\data\\figure_3\\thirst_hunger_new_pvalues_decoding.Rda", base_output_path), pvalues_all)
  save(file=sprintf("%s\\data\\figure_3\\thirst_thirst_new_decoded_all.Rda", base_output_path), true_vs_decoded_all)
}

new_version_across_mice_decoding <- function()
{

    
    control_datasets=F
    paths1 <- get_thirsty_quenched_paths()
    paths2 <- get_hypertonic_saline_paths()
    metadata_1 <- across_mice_decoding_build_metadata(paths1, control=control_datasets)
    metadata_2 <- across_mice_decoding_build_metadata(paths2)
    prefix = "thirst_hypertonic"
    is_mode = F
    remove_main_diagonal=F
    
    
    
    if (all(paths2 == get_hungry_sated_paths())) {
      print("Cleaning hs metadata due to licking trials that arent real")
      metadata_2 <- clean_metadata_HS(metadata_2)
      clusters_found <- unlist(lapply(1:len(paths2), function(pidx) {len(unique(metadata_2[[pidx]]$cluster_mat$labs))}))
      ds_to_remove <- -which(clusters_found != 8)
      metadata_2 <-  metadata_2[ds_to_remove] # Remove the only dataset that wasn't automatically detected with 7 clusters
      paths2 <- paths2[ds_to_remove]
    }
    
    if (all(paths1 == get_hungry_sated_paths())) {
      print("Cleaning hs metadata due to licking trials that arent real")
      metadata_1 <- clean_metadata_HS(metadata_2)
      clusters_found <- unlist(lapply(1:len(paths1), function(pidx) {len(unique(metadata_1[[pidx]]$cluster_mat$labs))}))
      ds_to_remove <- -which(clusters_found != 8)
      metadata_1 <-  metadata_1[ds_to_remove] # Remove the only dataset that wasn't automatically detected with 7 clusters
      paths1 <- paths1[ds_to_remove]
    }
    
    cos_dist_vec <- function(a,b) {cosine(a,b)[1,1]}
    shuffle <- function(vec) {sample(vec, len(vec))}
    
    confusion_matrix_list <- list()
    shuffle_matrix_list <- list()
    f1_confusion_matrix_list <- list()
    f1_shuffle_matrix_list <- list()
    weird_confusion_matrix_list <- list()
    weird_shuffle_matrix_list <- list()
    weird_all_confusion_matrix_list <- list()
    weird_all_shuffle_matrix_list <- list()
    
    activity_confusion_matrix_list <- list()
    activity_shuffle_matrix_list <- list()
    f1_activity_confusion_matrix_list <- list()
    f1_activity_shuffle_matrix_list <- list()
    weird_activity_confusion_matrix_list <- list()
    weird_activity_shuffle_matrix_list <- list()
    weird_all_activity_confusion_matrix_list <- list()
    weird_all_activity_shuffle_matrix_list <- list()    
    
    pvalues_all <- list()
    f1_pvalues_all <- list()
    weird_pvalues_all <- list()
    weird_all_pvalues_all <- list()
    activity_pvalues_all <- list()
    f1_activity_pvalues_all <- list()
    weird_activity_pvalues_all <- list()
    weird_all_activity_pvalues_all <- list()
    
    
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
    
    Metrics_objective <- list(cosine=which.max,
                              euc=which.min,
                              lv=which.min,
                              jaccard=which.max,
                              cosineqg=which.max,
                              hamming=which.min,
                              manhattan=which.min)
    
    for (resp_name in Responses) {
      for (metric_name in names(Metrics)) {
        confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        activity_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        activity_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]]<- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        f1_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        f1_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        weird_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        weird_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        weird_all_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        weird_all_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1))
        
        
        f1_activity_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
        f1_activity_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
        weird_activity_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
        weird_activity_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
        weird_all_activity_confusion_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
        weird_all_activity_shuffle_matrix_list[[sprintf("%s_%s", resp_name, metric_name)]] <- 
          matrix(rep(1, times=len(metadata_1) * len(metadata_2)),nrow=len(metadata_1)) 
        
      }
    }
    
    hits_decoding <- list()
    true_vs_decoded_all <- list()
    
  
    activity_true_vs_decoded_all <- list()
    for (metric_name in names(Metrics)) {
      pvalues_all[metric_name] <- c()
      true_vs_decoded_all[[metric_name]] <- list()
      
      activity_pvalues_all[metric_name] <- c()
      activity_true_vs_decoded_all[[metric_name]] <- list()
      
      f1_pvalues_all <- list()
      weird_pvalues_all <- list()
      weird_all_pvalues_all <- list()

      f1_activity_pvalues_all <- list()
      weird_activity_pvalues_all <- list()
      weird_all_activity_pvalues_all <- list()
    }

    
    window_size= 15
    trial_duration = 16

    for (idx_1 in 1:len(metadata_1)) {
      for (idx_2 in 1:len(metadata_2)) {
        
        if (remove_main_diagonal) {
          if (idx_1 == idx_2) {
            next
          }
        }
        
        
        annotated_mice_1 <- metadata_1[[idx_1]]
        annotated_mice_2 <- metadata_2[[idx_2]]
        
        
        if (is_mode) {
          mean_activity_1 = get_mode(annotated_mice_1, 
                                     path=paths[[idx_1]],
                                     axis=T,
                                     projection_only=T)
          
          mean_activity_2 = get_mode(annotated_mice_2, 
                                     path=paths[[idx_2]],
                                     axis=T,
                                     projection_only=T)
          
          mean_activity_1 <- as.vector(mean_activity_1)
          mean_activity_2 <- as.vector(mean_activity_2)
          
        } else {
          orig_mt_1 <- (get_reduced_mat_full_day(paths1[[idx_1]], window_size=15, just_original_mat = T, normed = F, control = control_datasets))
          orig_mt_2 <- (get_reduced_mat_full_day(paths2[[idx_2]], window_size=15, just_original_mat = T, normed = F, control = control_datasets))
          
          if(ncol(orig_mt_1) > nrow(orig_mt_1)) {orig_mt_1 <- t(orig_mt_1)}
          if(ncol(orig_mt_2) > nrow(orig_mt_2)) {orig_mt_2 <- t(orig_mt_2)}
          
          mean_activity_1 <- rowMeans(orig_mt_1)
          mean_activity_2 <- rowMeans(orig_mt_2)
        }
        

        
        stim_master_1 <- annotated_mice_1$stim_master_mat
        stim_master_2 <- annotated_mice_2$stim_master_mat 
        
        labels_mice_1 <- annotated_mice_1$cluster_mat$labs
        labels_mice_2 <- annotated_mice_2$cluster_mat$labs 
        
        mask_mice_1 <- annotated_mice_1$mask
        mask_mice_2 <- annotated_mice_2$mask
        
        decoding_stats <- list()
        activity_decoding_stats <- list()
        f1_decoding_stats <- list()
        f1_activity_decoding_stats <- list()
        weird_decoding_stats <- list()
        weird_activity_decoding_stats <- list()
        weird_all_decoding_stats <- list()
        weird_all_activity_decoding_stats <- list()
        
        for (metric_name in names(Metrics)) {
          decoding_stats[[metric_name]] <- c()
          activity_decoding_stats[[metric_name]] <- c()
          f1_decoding_stats[[metric_name]] <- c()
          f1_activity_decoding_stats[[metric_name]] <- c()
          weird_decoding_stats[[metric_name]] <- c()
          weird_activity_decoding_stats[[metric_name]] <- c()
          weird_all_decoding_stats[[metric_name]] <- c()
          weird_all_activity_decoding_stats[[metric_name]] <- c()
        }
        
        
        trials_mice_1 <- which(stim_master_1[,"Response"] %in% c(0,1,2,3,4,5))
        trials_labels_mice_1 <- stim_master_1[trials_mice_1,"Response"]
        trials_mat_mice_1 <- c()
        activity_mat_trials_mice_1  <- c()
        
        for (r_trial in trials_mice_1) {
          trial_binned_index <- as.numeric(get_binned_index(stim_master_1[r_trial,"Frames"], window_size))
          
          clustered_trial <- labels_mice_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          activity_trial <- mean_activity_1[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
            activity_trial[is.na(activity_trial)] <- activity_trial[which(!is.na(activity_trial))[len(which(!is.na(activity_trial)))]]
          }
          
          trials_mat_mice_1 <- rbind(trials_mat_mice_1, clustered_trial)
          activity_mat_trials_mice_1 <- rbind(activity_mat_trials_mice_1, activity_trial)
          
          
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
          activity_trial <- mean_activity_2[trial_binned_index:(trial_binned_index + trial_duration - 1)]
          
          if (sum(is.na(clustered_trial)) > 0 ) {
            clustered_trial[is.na(clustered_trial)] <- clustered_trial[which(!is.na(clustered_trial))[len(which(!is.na(clustered_trial)))]]
            activity_trial[is.na(activity_trial)] <- activity_trial[which(!is.na(activity_trial))[len(which(!is.na(activity_trial)))]]
          }
          
          trials_mat_mice_2 <- rbind(trials_mat_mice_2,clustered_trial)
          
          activity_mat_trials_mice_2 <- rbind(activity_mat_trials_mice_2, activity_trial)
          
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
        
        
        similarity_mt <- cosine(t(rbind(translated_mat_mice_1_mice_2, trials_mat_mice_2)))
        similarity_mt <- similarity_mt[1:nrow(translated_mat_mice_1_mice_2),
                                       (nrow(translated_mat_mice_1_mice_2) + 1):(nrow(translated_mat_mice_1_mice_2) + nrow(trials_mat_mice_2))]
        
        
        activity_similarity_mt <- cosine(t(rbind(activity_mat_trials_mice_1, activity_mat_trials_mice_2)))
        
        activity_similarity_mt <- activity_similarity_mt[1:nrow(activity_mat_trials_mice_1),
                                            (nrow(activity_mat_trials_mice_1) + 1):(nrow(activity_mat_trials_mice_1) + nrow(activity_mat_trials_mice_2))]
      
        
        
        for (shuff_i in c(-1:100)) {
          
          
          if (shuff_i != -1) {
            trials_labels_mice_1 <- shuffle(stim_master_1[trials_mice_1,"Response"])
            trials_labels_mice_2 <- shuffle(stim_master_2[trials_mice_2,"Response"])
          } else {
            trials_labels_mice_1 <- stim_master_1[trials_mice_1,"Response"]
            trials_labels_mice_2 <- stim_master_2[trials_mice_2,"Response"]
          }
          
          
          activity_decoded_vec <- list()
          decoded_vec <- list()
          for (metric_name in names(Metrics)) {
            decoded_vec[[metric_name]] <- c()
            activity_decoded_vec[[metric_name]] <- c()
          }
          
          for (decoded_trial_idx in 1:nrow(trials_mat_mice_2)) {
            
            
            confidence <- similarity_mt[,decoded_trial_idx]
            conf_vec <- c()
            
            activity_confidence <- activity_similarity_mt[,decoded_trial_idx]
            activity_conf_vec <- c()
            
            conf_thresh <- mean(confidence) + -0.2 * sd(confidence)
            act_conf_thresh <- mean(activity_confidence) + -0.2 * sd(activity_confidence)
            #conf_thresh <- 0
            #act_conf_thresh <- 0.2

            #which.max(table(trials_labels_mice_1[order(confidence, decreasing=T)[1:50]]))
            trial_results <- sort(unique(trials_labels_mice_1))

            for (res in trial_results) {
              conf_vec <- rbind(conf_vec,mean(confidence[trials_labels_mice_1 == res & confidence > conf_thresh]))
              activity_conf_vec <- rbind(activity_conf_vec, mean(activity_confidence[trials_labels_mice_1 == res & activity_confidence > act_conf_thresh]))

            }

            colnames(conf_vec) <- names(Metrics)
            colnames(activity_conf_vec) <- names(Metrics)
            for (metric_name in names(Metrics)) {
              decoded_vec[[metric_name]] <- c(decoded_vec[[metric_name]],
                                              trial_results[Metrics_objective[[metric_name]](conf_vec[,metric_name])])

              activity_decoded_vec[[metric_name]] <- c(activity_decoded_vec[[metric_name]],
                                              trial_results[Metrics_objective[[metric_name]](activity_conf_vec[,metric_name])])


              # decoded_vec[[metric_name]] <- c(decoded_vec[[metric_name]],
              #                                 as.numeric(names(which.max(table(trials_labels_mice_1[order(confidence, decreasing=T)[1:10]])))))
              # 
              # activity_decoded_vec[[metric_name]] <- c(activity_decoded_vec[[metric_name]],
              #                                          as.numeric(names(which.max(table(trials_labels_mice_1[order(activity_confidence, decreasing=T)[1:10]])))))
              # 
              # itu <- 1:10
              # 
              # dec <-
              #   unlist(lapply(unique(trials_labels_mice_1[order(confidence, decreasing=T)[itu]]),
              #                 function(tt) {
              #                   ind <- trials_labels_mice_1[order(confidence, decreasing=T)[itu]] == tt;
              #                   return(sum(confidence[order(confidence, decreasing=T)[itu]][ind]))
              #                 }))
              # names(dec) <- unique(trials_labels_mice_1[order(confidence, decreasing=T)[itu]])
              # 
              # act_dec <-
              # unlist(lapply(unique(trials_labels_mice_1[order(activity_confidence, decreasing=T)[itu]]),
              #               function(tt) {
              #                 ind <- trials_labels_mice_1[order(activity_confidence, decreasing=T)[itu]] == tt;
              #                 return(sum(activity_confidence[order(activity_confidence, decreasing=T)[itu]][ind]))
              #               }))
              # 
              # names(act_dec) <- unique(trials_labels_mice_1[order(activity_confidence, decreasing=T)[itu]])
              # 
              # 
              # decoded_vec[[metric_name]] <- c(decoded_vec[[metric_name]],
              #                                 as.numeric(names(which.max(dec))))
              # 
              # activity_decoded_vec[[metric_name]] <- c(activity_decoded_vec[[metric_name]],
              #                                 as.numeric(names(which.max(act_dec))))
            }
          }    
          
          
          # 
          # for (tt1 in c(0:5)) {
          #   for (tt2 in c(0:5)) {
          # 
          #     work_mt <- activity_similarity_mt[trials_labels_mice_1 == tt1,trials_labels_mice_2 == tt2]
          #     pht <-
          #     pheatmap(work_mt, border_col=NA,
          #              col=rdylbu_cg(200), breaks=seq(-1,1,length.out=200),
          #              cluster_cols=F,
          #              show_rownames=F,
          #              show_colnames=F,
          #              main = sprintf("Rows (train - %d), Cols (test - %d)", tt1, tt2))
          # 
          # 
          #     pdf(sprintf("Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\trash\\tmp\\activity\\pheatmap_%d_%d.pdf", tt1, tt2),
          #         height=4,
          #         width=4)
          #     plot(pht[[4]])
          #     dev.off()
          #   }
          # }



          
          decoding_statistics <- list()
          activity_decoding_statistics <- list()
          for (metric_name in names(Metrics)) { 
            decoding_statistics[[metric_name]] <- c()
            activity_decoding_statistics[[metric_name]] <- c()
          }
          
          
          for (metric_name in names(Metrics)) { 
            for (response_type in c(0,1,2,3,4,5)) {
              response_indices <- trials_labels_mice_2 == response_type
              
              decoding_accuracy = sum(decoded_vec[[metric_name]][response_indices] == trials_labels_mice_2[response_indices]) / 
                len(trials_labels_mice_2[response_indices])
              
              activity_decoding_accuracy =  sum(activity_decoded_vec[[metric_name]][response_indices] == trials_labels_mice_2[response_indices]) / 
                                len(trials_labels_mice_2[response_indices])
              
              decoding_statistics[[metric_name]] <- 
                c(decoding_statistics[[metric_name]], decoding_accuracy)
              
              activity_decoding_statistics[[metric_name]] <- 
                c(activity_decoding_statistics[[metric_name]], activity_decoding_accuracy)
            }
            
            if (shuff_i == -1) {
              true_vs_decoded_all[[metric_name]] <- 
                append(true_vs_decoded_all[[metric_name]],
                       list(list(train=idx_1, 
                                 test=idx_2,
                                 decoded=decoded_vec[[metric_name]],
              
                                                    gt=trials_labels_mice_2)))
              
              
              
              rs2 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 2]
              rs2f <- sum(rs2==0)/len(rs2)
              
              rs3 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 3]
              rs3f <- sum(rs3==0)/len(rs3)
              
              rs4 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 4]
              rs4f <- sum(rs4==0)/len(rs4)
            
              rs5 <- decoded_vec[[metric_name]][trials_labels_mice_2 == 5]
              rs5f <- sum(rs5==0)/len(rs5)
              
              print(sprintf("%.3f - %.3f (FAs) %.3f - %.3f (CRs)", rs3f, rs5f, rs2f, rs4f))
              
              
              activity_true_vs_decoded_all[[metric_name]] <- 
                append(activity_true_vs_decoded_all[[metric_name]],
                       list(list(train=idx_1, 
                                 test=idx_2,
                                 decoded=activity_decoded_vec[[metric_name]],
                                 gt=trials_labels_mice_2)))
              
              
              rs2 <- activity_decoded_vec[[metric_name]][trials_labels_mice_2 == 2]
              rs2f <- sum(rs2==0)/len(rs2)
              
              rs3 <- activity_decoded_vec[[metric_name]][trials_labels_mice_2 == 3]
              rs3f <- sum(rs3==0)/len(rs3)
              
              rs4 <- activity_decoded_vec[[metric_name]][trials_labels_mice_2 == 4]
              rs4f <- sum(rs4==0)/len(rs4)
              
              rs5 <- activity_decoded_vec[[metric_name]][trials_labels_mice_2 == 5]
              rs5f <- sum(rs5==0)/len(rs5)
              
              print(sprintf("%.3f - %.3f (FAs) %.3f - %.3f (CRs)", rs3f, rs5f, rs2f, rs4f))
            }
            
            
            f1_clust <- get_f1_scores(trials_labels_mice_2, decoded_vec[[metric_name]])
            f1_activity <- get_f1_scores(trials_labels_mice_2, activity_decoded_vec[[metric_name]])
            
            f1_decoding_stats[[metric_name]] <- rbind(f1_decoding_stats[[metric_name]],
                                                      f1_clust$f1)
            f1_activity_decoding_stats[[metric_name]] <- rbind(f1_activity_decoding_stats[[metric_name]],
                                                               f1_activity$f1)
            weird_decoding_stats[[metric_name]] <- rbind(weird_decoding_stats[[metric_name]],
                                                         f1_clust$weird)
            weird_activity_decoding_stats[[metric_name]] <- rbind(weird_activity_decoding_stats[[metric_name]],
                                                                  f1_activity$weird)
            weird_all_decoding_stats[[metric_name]] <- rbind(weird_all_decoding_stats[[metric_name]],
                                                             f1_clust$weird_all)
            weird_all_activity_decoding_stats[[metric_name]] <- rbind(weird_all_activity_decoding_stats[[metric_name]],
                                                                      f1_activity$weird_all)
            
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
            
            activity_decoding_stats[[metric_name]] <- rbind(activity_decoding_stats[[metric_name]], activity_decoding_statistics[[metric_name]])
          }
        }
        
        for (metric_name in names(Metrics)) { 
          
          
          colnames(decoding_stats[[metric_name]]) <- Responses
          decoding_statistics <- decoding_stats[[metric_name]][1,]
          shuffled_statistics <- colMeans(decoding_stats[[metric_name]][-1,])
          
          colnames(activity_decoding_stats[[metric_name]]) <- Responses
          activity_decoding_statistics <- activity_decoding_stats[[metric_name]][1,]
          activity_shuffled_statistics <- colMeans(activity_decoding_stats[[metric_name]][-1,])
          
          colnames(f1_decoding_stats[[metric_name]]) <- Responses
          f1_decoding_statistics <- f1_decoding_stats[[metric_name]][1,]
          f1_shuffled_statistics <- colMeans(f1_decoding_stats[[metric_name]][-1,])
          
          colnames(f1_activity_decoding_stats[[metric_name]]) <- Responses
          f1_activity_decoding_statistics <- f1_activity_decoding_stats[[metric_name]][1,]
          f1_activity_shuffled_statistics <- colMeans(f1_activity_decoding_stats[[metric_name]][-1,])
          
          colnames(weird_decoding_stats[[metric_name]]) <- Responses
          weird_decoding_statistics <- weird_decoding_stats[[metric_name]][1,]
          weird_shuffled_statistics <- colMeans(weird_decoding_stats[[metric_name]][-1,])
          
          colnames(weird_activity_decoding_stats[[metric_name]]) <- Responses
          weird_activity_decoding_statistics <- weird_activity_decoding_stats[[metric_name]][1,]
          weird_activity_shuffled_statistics <- colMeans(weird_activity_decoding_stats[[metric_name]][-1,])
          
          colnames(weird_all_decoding_stats[[metric_name]]) <- Responses
          weird_all_decoding_statistics <- weird_all_decoding_stats[[metric_name]][1,]
          weird_all_shuffled_statistics <- colMeans(weird_all_decoding_stats[[metric_name]][-1,])
          
          colnames(weird_all_activity_decoding_stats[[metric_name]]) <- Responses
          weird_all_activity_decoding_statistics <- weird_all_activity_decoding_stats[[metric_name]][1,]
          weird_all_activity_shuffled_statistics <- colMeans(weird_all_activity_decoding_stats[[metric_name]][-1,])
          
          ecdfs <- apply(decoding_stats[[metric_name]][-1,], 2, ecdf)
          pvalues <-  unlist(lapply(1:len(decoding_stats[[metric_name]][1,]), 
                                    function(acc_idx) {
                                      ecdfs[[acc_idx]](decoding_stats[[metric_name]][1,acc_idx])
                                    }))
          
          activity_ecdfs <- apply(activity_decoding_stats[[metric_name]][-1,], 2, ecdf)
          activity_pvalues <-  unlist(lapply(1:len(activity_decoding_stats[[metric_name]][1,]), 
                                    function(acc_idx) {
                                      activity_ecdfs[[acc_idx]](activity_decoding_stats[[metric_name]][1,acc_idx])
                                    }))
          
          
          f1_ecdfs <- apply(f1_decoding_stats[[metric_name]][-1,], 2, ecdf)
          f1_pvalues <- unlist(lapply(1:len(f1_decoding_stats[[metric_name]][1,]),
                                      function(acc_idx) {
                                        f1_ecdfs[[acc_idx]](f1_decoding_stats[[metric_name]][1,acc_idx])
                                      }))
          weird_ecdfs <- apply(weird_decoding_stats[[metric_name]][-1,], 2, ecdf)
          weird_pvalues <- unlist(lapply(1:len(weird_decoding_stats[[metric_name]][1,]),
                                         function(acc_idx) {
                                           weird_ecdfs[[acc_idx]](weird_decoding_stats[[metric_name]][1,acc_idx])
                                         }))
          
          weird_all_ecdfs <- apply(weird_all_decoding_stats[[metric_name]][-1,], 2, ecdf)
          weird_all_pvalues <- unlist(lapply(1:len(weird_all_decoding_stats[[metric_name]][1,]),
                                             function(acc_idx) {
                                               weird_all_ecdfs[[acc_idx]](weird_all_decoding_stats[[metric_name]][1,acc_idx])
                                             }))
          f1_activity_ecdfs <- apply(f1_activity_decoding_stats[[metric_name]][-1,], 2, ecdf)
          f1_activity_pvalues <- unlist(lapply(1:len(f1_activity_decoding_stats[[metric_name]][1,]),
                                               function(acc_idx) {
                                                 f1_activity_ecdfs[[acc_idx]](f1_activity_decoding_stats[[metric_name]][1,acc_idx])
                                               }))
          weird_activity_ecdfs <- apply(weird_activity_decoding_stats[[metric_name]][-1,], 2, ecdf)
          weird_activity_pvalues <- unlist(lapply(1:len(weird_activity_decoding_stats[[metric_name]][1,]),
                                                  function(acc_idx) {
                                                    weird_activity_ecdfs[[acc_idx]](weird_activity_decoding_stats[[metric_name]][1,acc_idx])
                                                  }))
          weird_all_activity_ecdfs <- apply(weird_all_activity_decoding_stats[[metric_name]][-1,], 2, ecdf)
          weird_all_activity_pvalues <- unlist(lapply(1:len(weird_all_activity_decoding_stats[[metric_name]][1,]),
                                                      function(acc_idx) {
                                                        weird_all_activity_ecdfs[[acc_idx]](weird_all_activity_decoding_stats[[metric_name]][1,acc_idx])
                                                      }))
          
          names(pvalues) <- Responses
          names(activity_pvalues) <- Responses
          names(f1_pvalues) <- Responses
          names(weird_pvalues) <- Responses
          names(weird_all_pvalues) <- Responses
          names(f1_activity_pvalues) <- Responses
          names(weird_activity_pvalues) <- Responses
          names(weird_all_activity_pvalues) <- Responses
          
          
          
          f1_pvalues_all[[metric_name]] <- rbind(f1_pvalues_all[[metric_name]],f1_pvalues)
          weird_pvalues_all[[metric_name]] <- rbind(weird_pvalues_all[[metric_name]],weird_pvalues)
          weird_all_pvalues_all[[metric_name]] <- rbind(weird_all_pvalues_all[[metric_name]],weird_all_pvalues)
          f1_activity_pvalues_all[[metric_name]] <- rbind(f1_activity_pvalues_all[[metric_name]],f1_activity_pvalues)
          weird_activity_pvalues_all[[metric_name]] <- rbind(weird_activity_pvalues_all[[metric_name]],weird_activity_pvalues)
          weird_all_activity_pvalues_all[[metric_name]] <- rbind(weird_all_activity_pvalues_all[[metric_name]],weird_all_activity_pvalues)
          
          pvalues_all[[metric_name]] <- rbind(pvalues_all[[metric_name]],pvalues)
          activity_pvalues_all[[metric_name]] <- rbind(activity_pvalues_all[[metric_name]],pvalues)
          
          print(sprintf("######################## SHUFFLED DECODING: (Shuffle vs decoder) - %d vs  %d ####", idx_1, idx_2))
          
          for (resp_name in Responses)  {
            
            mat_name <- sprintf("%s_%s", resp_name, metric_name)
            confusion_matrix_list[[mat_name]][idx_1, idx_2] <- decoding_statistics[resp_name]
            shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- shuffled_statistics[resp_name]
            
            activity_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- activity_decoding_statistics[resp_name]
            activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- activity_shuffled_statistics[resp_name]
            
            
            
            f1_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- f1_decoding_statistics[resp_name]
            f1_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- f1_shuffled_statistics[resp_name]
            weird_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- weird_decoding_statistics[resp_name]
            weird_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- weird_shuffled_statistics[resp_name]
            weird_all_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- weird_all_decoding_statistics[resp_name]
            weird_all_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- weird_all_shuffled_statistics[resp_name]
            
            
            f1_activity_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- f1_activity_decoding_statistics[resp_name]
            f1_activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- f1_activity_shuffled_statistics[resp_name]
            weird_activity_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- weird_activity_decoding_statistics[resp_name]
            weird_activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- weird_activity_shuffled_statistics[resp_name]
            weird_all_activity_confusion_matrix_list[[mat_name]][idx_1, idx_2] <- weird_all_activity_decoding_statistics[resp_name]
            weird_all_activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2] <- weird_all_activity_shuffled_statistics[resp_name]
            
            
            #print(sprintf("%.3f vs %.3f", shuffled_decoding_statistics, decoding_statistics))
            print(sprintf("(%d). %s: %.3f [A: %.3f] vs %.3f [A: %.3f] [Pvalue: %.3f, A_Pvalue: %.3f] -> %s", 
                          nrow(pvalues_all[[metric_name]]),
                          resp_name, 
                          shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          activity_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          pvalues[[resp_name]],
                          activity_pvalues[[resp_name]],
                          metric_name))
            print(sprintf("--> F1(C[%.3f  - Shf: %.3f] vs A([%.3f - Shf: %.3f]) F1 Pvalue (C[%.3f] vs A[%.3f])",
                          f1_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          f1_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          f1_activity_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          f1_activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          f1_pvalues[[resp_name]],
                          f1_activity_pvalues[[resp_name]]))


            print(sprintf("--> W(C[%.3f  - Shf: %.3f] vs A[%.3f  - Shf: %.3f]) W Pvalue (C[%.3f] vs A[%.3f])",
                          weird_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          weird_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          weird_activity_confusion_matrix_list[[mat_name]][idx_1, idx_2],
                          weird_activity_shuffle_matrix_list[[mat_name]][idx_1, idx_2],
                          weird_pvalues[[resp_name]],
                          weird_activity_pvalues[[resp_name]]))
            
            
            
            
            
          }
          
          print("#######################")
          
        }
      }
    }
    
    rownames(confusion_matrix) <- datasets_names
    # colnames(confusion_matrix) <- datasets_names
    
    dir.create(sprintf("%s\\data\\figure_3\\new_decoder", base_output_path))
    dir.create(sprintf("%s\\data\\figure_3\\new_decoder\\%s\\", base_output_path, prefix))
    
    if (is_mode) {
      dir.create(sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\", base_output_path, prefix))
      
      mode_confusion_matrix_list <- activity_confusion_matrix_list
      mode_shuffle_matrix_list <- activity_shuffle_matrix_list
      f1_mode_confusion_matrix_list <- f1_activity_confusion_matrix_list
      f1_mode_shuffle_matrix_list <- f1_activity_shuffle_matrix_list
      weird_mode_confusion_matrix_list <- weird_activity_confusion_matrix_list
      weird_mode_shuffle_matrix_list <- weird_activity_shuffle_matrix_list
      weird_all_mode_confusion_matrix_list <- weird_all_activity_confusion_matrix_list
      weird_all_mode_shuffle_matrix_list <- weird_all_activity_shuffle_matrix_list
      mode_pvalues_all <- activity_pvalues_all
      f1_mode_pvalues_all <- f1_activity_pvalues_all
      weird_mode_pvalues_all <- weird_activity_pvalues_all
      weird_all_mode_pvalues_all <- weird_all_activity_pvalues_all
      mode_true_vs_decoded_all <- activity_true_vs_decoded_all
      
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\confusion_matrix_list.Rda", base_output_path, prefix), confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\shuffle_matrix_list.Rda", base_output_path, prefix), shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\pvalues_all.Rda", base_output_path, prefix), pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_confusion_matrix_list.Rda", base_output_path, prefix), f1_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_shuffle_matrix_list.Rda", base_output_path, prefix), f1_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_pvalues_all.Rda", base_output_path, prefix), f1_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_confusion_matrix_list.Rda", base_output_path, prefix), weird_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_shuffle_matrix_list.Rda", base_output_path, prefix), weird_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_pvalues_all.Rda", base_output_path, prefix), weird_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_confusion_matrix_list.Rda", base_output_path, prefix), weird_all_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_shuffle_matrix_list.Rda", base_output_path, prefix), weird_all_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_pvalues_all.Rda", base_output_path, prefix), weird_all_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\mode_confusion_matrix_list.Rda", base_output_path, prefix), mode_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\mode_shuffle_matrix_list.Rda", base_output_path, prefix), mode_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\mode_pvalues_all.Rda", base_output_path, prefix), mode_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_mode_confusion_matrix_list.Rda", base_output_path, prefix), f1_mode_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_mode_shuffle_matrix_list.Rda", base_output_path, prefix), f1_mode_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\f1_mode_pvalues_all.Rda", base_output_path, prefix), f1_mode_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_mode_confusion_matrix_list.Rda", base_output_path, prefix), weird_mode_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_mode_shuffle_matrix_list.Rda", base_output_path, prefix), weird_mode_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_mode_pvalues_all.Rda", base_output_path, prefix), weird_mode_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_mode_confusion_matrix_list.Rda", base_output_path, prefix), weird_all_mode_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_mode_shuffle_matrix_list.Rda", base_output_path, prefix), weird_all_mode_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\weird_all_mode_pvalues_all.Rda", base_output_path, prefix), weird_all_mode_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\true_vs_decoded_all.Rda", base_output_path, prefix), true_vs_decoded_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\mode\\mode_true_vs_decoded_all.Rda", base_output_path, prefix), mode_true_vs_decoded_all)      
    } else {
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\confusion_matrix_list.Rda", base_output_path, prefix), confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\shuffle_matrix_list.Rda", base_output_path, prefix), shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\pvalues_all.Rda", base_output_path, prefix), pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_confusion_matrix_list.Rda", base_output_path, prefix), f1_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_shuffle_matrix_list.Rda", base_output_path, prefix), f1_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_pvalues_all.Rda", base_output_path, prefix), f1_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_confusion_matrix_list.Rda", base_output_path, prefix), weird_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_shuffle_matrix_list.Rda", base_output_path, prefix), weird_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_pvalues_all.Rda", base_output_path, prefix), weird_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_confusion_matrix_list.Rda", base_output_path, prefix), weird_all_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_shuffle_matrix_list.Rda", base_output_path, prefix), weird_all_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_pvalues_all.Rda", base_output_path, prefix), weird_all_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\activity_confusion_matrix_list.Rda", base_output_path, prefix), activity_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\activity_shuffle_matrix_list.Rda", base_output_path, prefix), activity_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\activity_pvalues_all.Rda", base_output_path, prefix), activity_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_activity_confusion_matrix_list.Rda", base_output_path, prefix), f1_activity_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_activity_shuffle_matrix_list.Rda", base_output_path, prefix), f1_activity_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_activity_pvalues_all.Rda", base_output_path, prefix), f1_activity_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_activity_confusion_matrix_list.Rda", base_output_path, prefix), weird_activity_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_activity_shuffle_matrix_list.Rda", base_output_path, prefix), weird_activity_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_activity_pvalues_all.Rda", base_output_path, prefix), weird_activity_pvalues_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_activity_confusion_matrix_list.Rda", base_output_path, prefix), weird_all_activity_confusion_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_activity_shuffle_matrix_list.Rda", base_output_path, prefix), weird_all_activity_shuffle_matrix_list)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_activity_pvalues_all.Rda", base_output_path, prefix), weird_all_activity_pvalues_all)    
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\true_vs_decoded_all.Rda", base_output_path, prefix), true_vs_decoded_all)
      save(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\activity_true_vs_decoded_all.Rda", base_output_path, prefix), activity_true_vs_decoded_all)
    }
}