



sizes = list(big=c(width=3.5),
             medium=c(width=3),
             medium_2=c(width=2.5),
             medium_1=c(width=2),
             small=c(width=1.75))
write_path <- sprintf("%s\\figure_SFO", figures_base_path)
dir.create(write_path)



similarity_df <- data.frame()

metadata_list <- list(SFO=across_mice_decoding_build_metadata(get_SFO_paths(), preset = "dalshtaim"),
                      SFO_scam=across_mice_decoding_build_metadata(get_SFO_paths(), preset = "dalshtaimscam"),
                      AGRP=across_mice_decoding_build_metadata(get_agrp_paths(), preset = "dalshtaim"),
                      HS=across_mice_decoding_build_metadata(get_hungry_sated_paths(), preset = "dalshtaim"),
                      TQ=across_mice_decoding_build_metadata(get_thirsty_quenched_paths(), preset = "dalshtaim"))

SFO_chunks = list(`1`=4:5,
                  `2`=5:6,
                  `3`=6:8,
                  `4`=7:9)


AGRP_chunks = list(`1`=5:9,
                  `2`=6:10,
                  `3`=5:9,
                  `4`=3:4)

pre_SFO_chunks = list(`1`=c(1,3),
                  `2`=c(1,4),
                  `3`=1:3,
                  `4`=1:3)

pre_AGRP_chunks = list(`1`=1:4,
                   `2`=1:5,
                   `3`=1:4,
                   `4`=1:2)

similarity_dist <-list(SFO=list(),
                       SFO_scam=list(),
                       AGRP=list(),
                       HS=list(),
                       TQ=list())


first_similarity_dist <-list(SFO=list(),
                             SFO_scam=list(),
                             AGRP=list(),
                             HS=list(),
                             TQ=list())

second_similarity_dist <-list(SFO=list(),
                              SFO_scam=list(),
                              AGRP=list(),
                              HS=list(),
                              TQ=list())


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
  
  
# 
#   for (size_name in names(sizes)) {
#     tmp_write_path <- sprintf("%s\\example_cluster_map", write_path)
#     dir.create(tmp_write_path)
#     tmp_write_path <- sprintf("%s\\%s", tmp_write_path, size_name)
#     dir.create(tmp_write_path)
#     
#     if (size_name == "medium") {
#     ph <- pheatmap(metadata[[i]]$cluster_mat$clust_mat, cluster_rows=F, cluster_cols=F, border_color=NA, show_colnames = F, show_rownames = F)
#     } else {
#       ph <- pheatmap(metadata[[i]]$cluster_mat$clust_mat, cluster_rows=F, cluster_cols=F, border_color=NA)
#     }
#     pdf(sprintf("%s\\%s_%d.pdf",
#                 tmp_write_path,
#                 metadata_name,
#                 i),
#         height=sizes[[size_name]],
#         width=sizes[[size_name]])
#     
#     plot(ph[[4]])
#     dev.off()
#     
#     
#   }
  

  
  if (metadata_name == "SFO" || metadata_name == "SFO_scam") {
    first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],15] %in% SFO_chunks[[as.character(i)]]]
    second_half_trials <- reward_trials[stim_master_mat[relevant_trials[reward_trials],15] %in% SFO_chunks[[as.character(i)]]]
    print("Using this SFO")
    print(SFO_chunks[[as.character(i)]])
  } else if (metadata_name == "AGRP") {
    
    
    if (i==1) {
      first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],9] %in% AGRP_chunks[[as.character(i)]]]#[1:5]  
    } else {
      first_half_trials <- reward_trials[!stim_master_mat[relevant_trials[reward_trials],9] %in% AGRP_chunks[[as.character(i)]]]
    }
    
    second_half_trials <- reward_trials[stim_master_mat[relevant_trials[reward_trials],9] %in% AGRP_chunks[[as.character(i)]]]
    print("Using this AGRP")
    print(AGRP_chunks[[as.character(i)]])
  } else {

    first_half_trials <- reward_trials[1:round(len(reward_trials) * .5)]
    second_half_trials <- reward_trials[(round(len(reward_trials) * .5) + 1):len(reward_trials)]
  }
  
  if (metadata_name == "SFO") {

    if (i == 4) {
      relevant_hits =  hits_decoding$cosine[1:len(hits_decoding$cosine) %% 4 == 0]
    } else {
      relevant_hits =  hits_decoding$cosine[1:len(hits_decoding$cosine) %% 4 == i]
    }
    
    print(sprintf("%d - %d - %d",
                  len(reward_trials),
                  len(relevant_hits[[1]]),
                  len(reward_trials) - len(relevant_hits[[1]])))
    
    
    tmpo <- 
    lapply(relevant_hits,
           
           function(hits) {
             
             first_half_hits <- hits[1:len(first_half_trials)]
             second_half_hits <- hits[(len(first_half_trials) + 1):len(relevant_hits)]
             
             return(c(sum(first_half_hits) / len(first_half_hits),
                      sum(second_half_hits) / len(second_half_hits)))
             
           })
    
    hits_dec_all <- rbind(hits_dec_all,
                          cbind(cbind(do.call(rbind, tmpo), 1:len(tmpo)), i))
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
                                                                                                                  len(first_half_trials):nrow(similarity_matrix)]))
  
  print(sprintf("all %f %f %d", median(similarity_matrix, na.rm=T), mean_sim, i))
  
  half_similarity <- c(half_similarity,
                       mean_sim)
    
  }
  
  similarity_df <- rbind(similarity_df,
                         data.frame(similarity=half_similarity,
                                    x=rep(metadata_name, times=len(half_similarity))))
}


 wilk <- wilcox.test(similarity_df[similarity_df$x == "AGRP",]$similarity,
                    similarity_df[similarity_df$x == "SFO",]$similarity,
                    alternative="greater",
                    correct=F)
gcomp_box <- 
  ggplot(similarity_df, aes(x=x,y=similarity)) + 

  geom_boxplot(aes(),
               width=0.5, 
               size=1, 
               color=adjustcolor("gray65", alpha=1), 
               fill=adjustcolor("gray80", alpha=0.9),
               outlier.shape = NA) + 
  geom_jitter(position=position_jitterdodge(.25), aes(fill=x),
              size=1.5, alpha=.75) + 
  theme_classic() +
  theme(text=element_text(size=14, color="black"),
        legend.position="NA")  +
  xlab("") +
  ylab("Hit trial cosine similarity") +
  theme(text=element_text(size=14, color="black"),
        legend.position="NA")# + 
  geom_text(data=data.frame(x=1.5, y=max(similarity_df$similarity), label=signif.num(wilk$p.value)),
            aes(x=x,y=y,label=label),
            size=5)
  
  
  
  fdf <- 
  lapply(names(similarity_dist),
         function(nm) {
           
           pooled_trial_sim_dist <- 
            unlist(lapply(similarity_dist[[nm]], c))
           
            h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=50), plot=F) 
            
            return(data.frame(Frac=h$counts / sum(h$counts), Breaks=h$breaks[-1], Group=rep(nm, times=len(h$counts))))
         })
  
  hist_df <- do.call(rbind, fdf)
  
  first_half_fdf <- lapply(names(first_similarity_dist),
                           function(nm) {
                             
                             pooled_trial_sim_dist <- 
                               unlist(lapply(first_similarity_dist[[nm]], c))
                             
                             h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=20), plot=F) 
                             
                             return(data.frame(Frac=h$counts / sum(h$counts), Breaks=h$breaks[-1], Group=rep(nm, times=len(h$counts))))
                           })
  
  first_hist_df <- do.call(rbind, first_half_fdf)
  
  
  second_half_fdf <- lapply(names(second_similarity_dist),
                           function(nm) {
                             
                             pooled_trial_sim_dist <- 
                               unlist(lapply(second_similarity_dist[[nm]], c))
                             
                             h <- hist(pooled_trial_sim_dist, breaks=seq(-1,1, length.out=20), plot=F) 
                             
                             return(data.frame(Frac=h$counts  / sum(h$counts), Breaks=h$breaks[-1], Group=rep(nm, times=len(h$counts))))
                           })
  
  second_hist_df <- do.call(rbind, second_half_fdf)
    
# 
# # sizes = list(big=c(width=3),
#              medium=c(width=2.75),
#              medium_2=c(width=2.5),
#              medium_1=c(width=2),
#              small=c(width=1.75))

write_path <- sprintf("%s\\figure_SFO", figures_base_path)
dir.create(write_path)

for (size_name in names(sizes)) {
  tmp_write_path <- sprintf("%s\\%s", write_path, size_name)
  dir.create(tmp_write_path)

  
  pdf(sprintf("%s\\trial_similarity_SFO_wt.pdf",
              tmp_write_path),
      height=sizes[[size_name]],
      width=sizes[[size_name]])
  
  plot(gcomp_box)
  dev.off()
  
  
}





