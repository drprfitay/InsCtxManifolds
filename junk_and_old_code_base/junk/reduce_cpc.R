

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_HS <- across_mice_decoding_build_metadata(get_hungry_sated_paths())
metadata_control <- across_mice_decoding_build_metadata(get_control_paths(), control=T)
metadata_SFO <- across_mice_decoding_build_metadata(get_SFO_paths())
metadata_AGRP <- across_mice_decoding_build_metadata(get_agrp_paths())
metadata_hypertonicsaline <- across_mice_decoding_build_metadata(get_hypertonic_saline_paths())[-4]

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



reduce_cpc <- function(paths) {
  
  metadata <- across_mice_decoding_build_metadata(paths)
  for (idx in 1:len(paths)) {
    p <- paths[idx]
    m_obj <- metadata[[idx]]
    original <- get_reduced_mat_full_day(p, window_size = 15, just_original_mat = T, activity_threshold = .5)
    cpc_ind <- m_obj$cluster_mat$labs == -1
    cpc_original <- original[,cpc_ind]
    
    rmat <- reduce_dim(t(cpc_original), "lem", 20, knn1=.075,knn2=0)
    reg_id <- twonn(rmat, c_trimmed = 0.04, method="linfit")
    rmat_f <- reduce_dim(rmat, "lem", 6, knn1=.025,knn2=0)
  }
}


cluster_cumulative <- function(mt) {
  
  cl1_all <- c()
  cl2_all <- c()
  reward_cl1_all <- c()
  reward_cl2_all <- c()
  
  cmsm_all <- c()
  
  for (idx in 1:len(metadata_TT)) {
    m_obj  <- metadata_TT[[idx]]
    
    
    used_trials <- m_obj$stim_master_mat[,"TrialType"] %in% c(1,3,4,5)
    cr_used_trials <- m_obj$stim_master_mat[,"Response"] %in% c(2,4)
    stim_master_mat <- m_obj$stim_master_mat
    cumsum_all <- rep(0, times=nrow(stim_master_mat))
    drink_ind <- which(!is.nan(stim_master_mat[,6]))
    cumsum_all[drink_ind] <- 1
    
    cumsum_f <- cumsum(cumsum_all)
    
    used_cumsum <- cumsum_f[used_trials] / max(cumsum_f[used_trials])
    cr_used_cumsum <- cumsum_f[cr_used_trials] / max(cumsum_f[cr_used_trials])
  
    cmsm_all <- rbind(cmsm_all,
                      used_cumsum[floor(seq(1, len(used_cumsum), length.out=200))])

    
    
    cr_zero_q <- 1
    cr_twenty_q <- which(cr_used_cumsum >= .2)[1]
    cr_forty_q <- which(cr_used_cumsum >=  .4)[1]
    cr_sixty_q <- which(cr_used_cumsum >= .6)[1]
    cr_eighty_q<- which(cr_used_cumsum >= .8)[1]
    cr_oneh_q <- sum(cr_used_trials)
    
    zero_q <- 1
    twenty_q <- which(used_cumsum >= .2)[1]
    forty_q <- which(used_cumsum >=  .4)[1]
    sixty_q <- which(used_cumsum >= .6)[1]
    eighty_q<- which(used_cumsum >= .8)[1]
    oneh_q <- sum(used_trials)
    
    
    cr_ind <- (m_obj$annot_df[,1] == 4 | m_obj$annot_df[,1] == 5) & (m_obj$annot_df[,2] == 0)
    # thirst_f <- m_obj$trials_mat[cr_ind,][cr_zero_q:cr_twenty_q,15:20]
    # second_f <- m_obj$trials_mat[cr_ind,][(cr_twenty_q + 1):cr_forty_q,15:20]
    # third_f <-  m_obj$trials_mat[cr_ind,][(cr_forty_q + 1):cr_sixty_q,15:20]
    # fourth_f <- m_obj$trials_mat[cr_ind,][(cr_sixty_q + 1):cr_eighty_q,15:20]
    # fifth_f <- m_obj$trials_mat[cr_ind,][(cr_eighty_q + 1):cr_oneh_q,15:20]


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
  
  
  gclust <- 
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
  
}


cluster_cumulative_agrp <- function(mt) {
  
  cl1_all <- c()
  cl2_all <- c()
  cl3_all <- c()
  reward_cl1_all <- c()
  reward_cl2_all <- c()
  
  cmsm_all <- c()
  
  manipulation_chunks <- AGRP_chunks
  metadata_manipulation <- metadata_AGRP
  for (idx in 1:len(metadata_manipulation)) {
    m_obj  <- metadata_manipulation[[idx]]
    
    
    used_trials <- m_obj$stim_master_mat[,"TrialType"] %in% c(1,3,4,5)
    cr_used_trials <- m_obj$stim_master_mat[,"Response"] %in% c(2,4)
    stim_master_mat <- m_obj$stim_master_mat
    cumsum_all <- rep(0, times=nrow(stim_master_mat))
    drink_ind <- which(!is.nan(stim_master_mat[,6]))
    cumsum_all[drink_ind] <- 1
    
    pre_manipulation_trials_all <- !stim_master_mat[,ncol(stim_master_mat)] %in% manipulation_chunks[[idx]]
    post_manipulation_trials_all <- stim_master_mat[,ncol(stim_master_mat)] %in% manipulation_chunks[[idx]]
    
    
    pre_manipulation_trials <- !stim_master_mat[used_trials,ncol(stim_master_mat)] %in% manipulation_chunks[[idx]]
    post_manipulation_trials <- stim_master_mat[used_trials,ncol(stim_master_mat)] %in% manipulation_chunks[[idx]]
    

    pre_trials_mat <- m_obj$trials_mat[pre_manipulation_trials,]
    post_trials_mat <- m_obj$trials_mat[post_manipulation_trials,]
    

    pre_stim_mat <- stim_master_mat[used_trials,][pre_manipulation_trials,]
    post_stim_mat <- stim_master_mat[used_trials,][post_manipulation_trials,]

    
    cumsum_f <- cumsum(cumsum_all)
    
    pre_used_cumsum <- cumsum_f[used_trials]
    pre_used_cumsum <- pre_used_cumsum[pre_manipulation_trials]
    pre_used_cumsum <- (pre_used_cumsum - min(pre_used_cumsum)) / (max(pre_used_cumsum) -min(pre_used_cumsum))
    post_used_cumsum <- cumsum_f[used_trials]
    post_used_cumsum <- post_used_cumsum[post_manipulation_trials]
    post_used_cumsum <- (post_used_cumsum - min(post_used_cumsum)) / (max(post_used_cumsum) - min(post_used_cumsum))
    
    
    
    cmsm_all <- rbind(cmsm_all,
                      used_cumsum[floor(seq(1, len(used_cumsum), length.out=200))])
    
    pre_zero_q <- 1
    pre_twenty_q <- which(pre_used_cumsum >= .25)[1]
    pre_forty_q <- which(pre_used_cumsum >=  .5)[1]
    pre_sixty_q <- which(pre_used_cumsum >= .75)[1]
    pre_eighty_q<- sum(pre_manipulation_trials)#pre_eighty_q<- which(pre_used_cumsum >= .8)[1]
    #pre_oneh_q <- sum(pre_manipulation_trials)
    
    post_zero_q <- 1
    post_twenty_q <- which(post_used_cumsum >= .25)[1]
    post_forty_q <- which(post_used_cumsum >=  .5)[1]
    post_sixty_q <- which(post_used_cumsum >= .75)[1]
    post_eighty_q <- sum(post_manipulation_trials)#which(post_used_cumsum >= .8)[1]
    #post_oneh_q <- sum(post_manipulation_trials)
    
    
    pre_thirst_f <- pre_trials_mat[pre_zero_q:pre_twenty_q,15:20]                     [pre_stim_mat[pre_zero_q:pre_twenty_q,"Response"] %in% c(2,4),]
    pre_second_f <- pre_trials_mat[(pre_twenty_q + 1):pre_forty_q,15:20]                     [pre_stim_mat[(pre_twenty_q + 1):pre_forty_q,"Response"] %in% c(2,4),]
    pre_third_f <-  pre_trials_mat[(pre_forty_q + 1):pre_sixty_q,15:20]                     [pre_stim_mat[(pre_forty_q + 1):pre_sixty_q,"Response"] %in% c(2,4),]
    pre_fourth_f <- pre_trials_mat[(pre_sixty_q + 1):pre_eighty_q,15:20]                     [pre_stim_mat[(pre_sixty_q + 1):pre_eighty_q,"Response"] %in% c(2,4),]
    #pre_fifth_f <- pre_trials_mat[(pre_eighty_q + 1):pre_oneh_q,15:20]                     [pre_stim_mat[(pre_eighty_q + 1):pre_oneh_q,"Response"] %in% c(2,4),]
    post_thirst_f <- post_trials_mat[post_zero_q:post_twenty_q,15:20]                     [post_stim_mat[post_zero_q:post_twenty_q,"Response"] %in% c(2,4),]
    post_second_f <- post_trials_mat[(post_twenty_q + 1):post_forty_q,15:20]                     [post_stim_mat[(post_twenty_q + 1):post_forty_q,"Response"] %in% c(2,4),]
    post_third_f <-  post_trials_mat[(post_forty_q + 1):post_sixty_q,15:20]                     [post_stim_mat[(post_forty_q + 1):post_sixty_q,"Response"] %in% c(2,4),]
    post_fourth_f <- post_trials_mat[(post_sixty_q + 1):post_eighty_q,15:20]                     [post_stim_mat[(post_sixty_q + 1):post_eighty_q,"Response"] %in% c(2,4),]
    #post_fifth_f <- post_trials_mat[(post_eighty_q + 1):post_oneh_q,15:20]                     [post_stim_mat[(post_eighty_q + 1):post_oneh_q,"Response"] %in% c(2,4),]
    
    
    pre_tabled_thirst_f <- rep(0, times=8); names(pre_tabled_thirst_f) <- c(-1,1:7); pre_tabled_thirst_f[names(table(pre_thirst_f))] <- table(pre_thirst_f)
    pre_tabled_second_f <- rep(0, times=8); names(pre_tabled_second_f) <- c(-1,1:7); pre_tabled_second_f[names(table(pre_second_f))] <- table(pre_second_f)
    pre_tabled_third_f <-  rep(0, times=8); names(pre_tabled_third_f) <- c(-1,1:7);  pre_tabled_third_f[names(table(pre_third_f))] <- table(pre_third_f)
    pre_tabled_fourth_f <- rep(0, times=8); names(pre_tabled_fourth_f) <- c(-1,1:7); pre_tabled_fourth_f[names(table(pre_fourth_f))] <- table(pre_fourth_f)
    #pre_tabled_fifth_f <-  rep(0, times=8); names(pre_tabled_fifth_f) <- c(-1,1:7); pre_tabled_fifth_f[names(table(pre_fifth_f))] <- table(pre_fifth_f)
    post_tabled_thirst_f <- rep(0, times=8); names(post_tabled_thirst_f) <- c(-1,1:7); post_tabled_thirst_f[names(table(post_thirst_f))] <- table(post_thirst_f)
    post_tabled_second_f <- rep(0, times=8); names(post_tabled_second_f) <- c(-1,1:7); post_tabled_second_f[names(table(post_second_f))] <- table(post_second_f)
    post_tabled_third_f <-  rep(0, times=8); names(post_tabled_third_f) <- c(-1,1:7);  post_tabled_third_f[names(table(post_third_f))] <- table(post_third_f)
    post_tabled_fourth_f <- rep(0, times=8); names(post_tabled_fourth_f) <- c(-1,1:7); post_tabled_fourth_f[names(table(post_fourth_f))] <- table(post_fourth_f)
    #post_tabled_fifth_f <-  rep(0, times=8); names(post_tabled_fifth_f) <- c(-1,1:7); post_tabled_fifth_f[names(table(post_fifth_f))] <- table(post_fifth_f)
    
    
    
    tmp <- rbind(pre_tabled_thirst_f,
                 pre_tabled_second_f,
                 pre_tabled_third_f,
                 pre_tabled_fourth_f,
                 #pre_tabled_fifth_f,
                 post_tabled_thirst_f,
                 post_tabled_second_f,
                 post_tabled_third_f,
                 post_tabled_fourth_f)
                 #post_tabled_fifth_f)
    
    
    pre_reward_thirst_f <- pre_trials_mat[pre_zero_q:pre_twenty_q,8:13]  [pre_stim_mat[pre_zero_q:pre_twenty_q,"Response"] == 0, ]
    pre_reward_second_f <- pre_trials_mat[(pre_twenty_q + 1):pre_forty_q,8:13]                [pre_stim_mat[(pre_twenty_q + 1):pre_forty_q,"Response"] == 0, ]
    pre_reward_third_f <-  pre_trials_mat[(pre_forty_q + 1):pre_sixty_q,8:13]                [pre_stim_mat[(pre_forty_q + 1):pre_sixty_q,"Response"] == 0, ]
    pre_reward_fourth_f <- pre_trials_mat[(pre_sixty_q + 1):pre_eighty_q,8:13]                [pre_stim_mat[(pre_sixty_q + 1):pre_eighty_q,"Response"] == 0, ]
    #pre_reward_fifth_f <- pre_trials_mat[(pre_eighty_q + 1):pre_oneh_q,8:13]                [pre_stim_mat[(pre_eighty_q + 1):pre_oneh_q,"Response"] == 0, ]
    
    post_reward_thirst_f <- post_trials_mat[post_zero_q:post_twenty_q,8:13]                [post_stim_mat[post_zero_q:post_twenty_q,"Response"] == 0, ]
    post_reward_second_f <- post_trials_mat[(post_twenty_q + 1):post_forty_q,8:13]                [post_stim_mat[(post_twenty_q + 1):post_forty_q,"Response"] == 0, ]
    post_reward_third_f <-  post_trials_mat[(post_forty_q + 1):post_sixty_q,8:13]                [post_stim_mat[(post_forty_q + 1):post_sixty_q,"Response"] == 0, ]
    post_reward_fourth_f <- post_trials_mat[(post_sixty_q + 1):post_eighty_q,8:13]                [post_stim_mat[(post_sixty_q + 1):post_eighty_q,"Response"] == 0, ]
    #post_reward_fifth_f <- post_trials_mat[(post_eighty_q + 1):post_oneh_q,8:13]                [post_stim_mat[(post_eighty_q + 1):post_oneh_q,"Response"] == 0, ]
    
    
    pre_reward_tabled_thirst_f <- rep(0, times=8); names(pre_reward_tabled_thirst_f) <- c(-1,1:7); pre_reward_tabled_thirst_f[names(table(pre_reward_thirst_f))] <- table(pre_reward_thirst_f)
    pre_reward_tabled_second_f <- rep(0, times=8); names(pre_reward_tabled_second_f) <- c(-1,1:7); pre_reward_tabled_second_f[names(table(pre_reward_second_f))] <- table(pre_reward_second_f)
    pre_reward_tabled_third_f <-  rep(0, times=8); names(pre_reward_tabled_third_f) <- c(-1,1:7);  pre_reward_tabled_third_f[names(table(pre_reward_third_f))] <- table(pre_reward_third_f)
    pre_reward_tabled_fourth_f <- rep(0, times=8); names(pre_reward_tabled_fourth_f) <- c(-1,1:7); pre_reward_tabled_fourth_f[names(table(pre_reward_fourth_f))] <- table(pre_reward_fourth_f)
    #pre_reward_tabled_fifth_f <-  rep(0, times=8); names(pre_reward_tabled_fifth_f) <- c(-1,1:7); pre_reward_tabled_fifth_f[names(table(pre_reward_fifth_f))] <- table(pre_reward_fifth_f)
    post_reward_tabled_thirst_f <- rep(0, times=8); names(post_reward_tabled_thirst_f) <- c(-1,1:7); post_reward_tabled_thirst_f[names(table(post_reward_thirst_f))] <- table(post_reward_thirst_f)
    post_reward_tabled_second_f <- rep(0, times=8); names(post_reward_tabled_second_f) <- c(-1,1:7); post_reward_tabled_second_f[names(table(post_reward_second_f))] <- table(post_reward_second_f)
    post_reward_tabled_third_f <-  rep(0, times=8); names(post_reward_tabled_third_f) <- c(-1,1:7);  post_reward_tabled_third_f[names(table(post_reward_third_f))] <- table(post_reward_third_f)
    post_reward_tabled_fourth_f <- rep(0, times=8); names(post_reward_tabled_fourth_f) <- c(-1,1:7); post_reward_tabled_fourth_f[names(table(post_reward_fourth_f))] <- table(post_reward_fourth_f)
    #post_reward_tabled_fifth_f <-  rep(0, times=8); names(post_reward_tabled_fifth_f) <- c(-1,1:7); post_reward_tabled_fifth_f[names(table(post_reward_fifth_f))] <- table(post_reward_fifth_f)
    

    
    reward_tmp <- 
    rbind(pre_reward_tabled_thirst_f,
          pre_reward_tabled_second_f,
          pre_reward_tabled_third_f,
          pre_reward_tabled_fourth_f,
          #pre_reward_tabled_fifth_f,
          post_reward_tabled_thirst_f,
          post_reward_tabled_second_f,
          post_reward_tabled_third_f,
          post_reward_tabled_fourth_f)
          #post_reward_tabled_fifth_f)

    tabled_reward <- table(pre_trials_mat[m_obj$annot_df[pre_manipulation_trials,1] == 3 & m_obj$annot_df[pre_manipulation_trials,2] == 0,8:13])
    reward_cl1 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[1]
    reward_cl2 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[2]

  
    
    
    
    cl1_vec <- c()
    cl2_vec <- c()
    cl3_vec <- c()
    
    
    contingency_list <- list(pre_tabled_thirst_f,
                             pre_tabled_second_f,
                             pre_tabled_third_f,
                             pre_tabled_fourth_f,
                             #pre_tabled_fifth_f,
                             post_tabled_thirst_f,
                             post_tabled_second_f,
                             post_tabled_third_f,
                             post_tabled_fourth_f)
                             #post_tabled_fifth_f)
    
     cl1 <- names(sort(post_tabled_fourth_f, decreasing=T)[names(sort(post_tabled_fourth_f, decreasing=T)) != -1])[1]
     cl2 <- names(sort(post_tabled_fourth_f, decreasing=T)[names(sort(post_tabled_fourth_f, decreasing=T)) != -1])[2]
     cl3 <- names(sort(post_tabled_fourth_f, decreasing=T)[names(sort(post_tabled_fourth_f, decreasing=T)) != -1])[3]
    
    for (tabled in contingency_list) {
      
      cl1_vec <- c(cl1_vec, tabled[cl1])
      cl2_vec <- c(cl2_vec, tabled[cl2])
      cl3_vec <- c(cl3_vec, tabled[cl3])
      
    }
    
    reward_cl1_vec <- c()
    reward_cl2_vec <- c()
    
    reward_contingency_list <- 
      list(pre_reward_tabled_thirst_f,
            pre_reward_tabled_second_f,
            pre_reward_tabled_third_f,
            pre_reward_tabled_fourth_f,
            #pre_reward_tabled_fifth_f,
            post_reward_tabled_thirst_f,
            post_reward_tabled_second_f,
            post_reward_tabled_third_f,
            post_reward_tabled_fourth_f)
            #post_reward_tabled_fifth_f)

    for (tabled in reward_contingency_list) {

      reward_cl1_vec <- c(reward_cl1_vec, tabled[reward_cl1])
      reward_cl2_vec <- c(reward_cl2_vec, tabled[reward_cl2])

    }


     
    cl1_all <- rbind(cl1_all, cl1_vec)
    cl2_all <- rbind(cl2_all, cl2_vec)
    cl3_all <- rbind(cl3_all, cl3_vec)
    
    reward_cl1_all <- rbind(reward_cl1_all, reward_cl1_vec)
    reward_cl2_all <- rbind(reward_cl2_all, reward_cl2_vec)
  }
  
  
  # cl1_all <- t(apply(cl1_all, 1, function(r) {c(r[1:5]/sum(r[1:5]),r[6:10]/sum(r[6:10]))}))
  # cl2_all <- t(apply(cl2_all, 1, function(r) {c(r[1:5]/sum(r[1:5]),r[6:10]/sum(r[6:10]))}))
  # cl3_all <- t(apply(cl2_all, 1, function(r) {c(r[1:5]/sum(r[1:5]),r[6:10]/sum(r[6:10]))}))
  # 
  # reward_cl1_all <- t(apply(reward_cl1_all, 1, function(r) {c(r[1:5]/sum(r[1:5]),r[6:10]/sum(r[6:10]))}))
  # reward_cl2_all <- t(apply(reward_cl2_all, 1, function(r) {c(r[1:5]/sum(r[1:5]),r[6:10]/sum(r[6:10]))}))
  
  cl1_all <- t(apply(cl1_all, 1, function(r) {c(r[1:4]/sum(r[1:4]),r[5:8]/sum(r[5:8]))}))
  cl2_all <- t(apply(cl2_all, 1, function(r) {c(r[1:4]/sum(r[1:4]),r[5:8]/sum(r[5:8]))}))
  cl3_all <- t(apply(cl2_all, 1, function(r) {c(r[1:4]/sum(r[1:4]),r[5:8]/sum(r[5:8]))}))
  
  reward_cl1_all <- t(apply(reward_cl1_all, 1, function(r) {c(r[1:4]/sum(r[1:4]),r[5:8]/sum(r[5:8]))}))
  reward_cl2_all <- t(apply(reward_cl2_all, 1, function(r) {c(r[1:4]/sum(r[1:4]),r[5:8]/sum(r[5:8]))}))
  
  colnames(cl1_all) <- 1:ncol(cl1_all)
  colnames(cl2_all) <- 1:ncol(cl2_all)
  colnames(reward_cl1_all) <- 1:ncol(reward_cl1_all)
  colnames(reward_cl2_all) <- 1:ncol(reward_cl2_all)
  
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
                    mdf_reward_cl1)
                    #mdf_reward_cl2)
  
  
  #gclust <- 
    ggplot(df_final) + 
    geom_vline(xintercept=4.5, linetype="dashed", size=.5, col="gray50") +
     # geom_vline(xintercept=.4, linetype="dashed", size=.5, col="gray50") +
     # geom_vline(xintercept=.6, linetype="dashed", size=.5, col="gray50") +
     # geom_vline(xintercept=.8, linetype="dashed", size=.5, col="gray50") +
     # geom_vline(xintercept=1, linetype="dashed", size=.5, col="gray50") +
    geom_point(aes(x=Thirst, y=Fraction, color=Group), stat="summary") +
    geom_line(aes(x=Thirst, y=Fraction, color=Group), stat="summary") +
    geom_ribbon(aes(x=Thirst, y=Fraction, fill=Group), stat="summary", alpha=.5, color=NA) + 
    big_text_base_plot_theme_wl +
    theme(legend.position="top") + 
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
  
  plot_grid(gbehav, gclust)
  
}
