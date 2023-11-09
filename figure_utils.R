library("philentropy")
library("RColorBrewer")
library("tidyverse")
library("intRinsic")
library("dbscan")
library("cowplot")
library("pheatmap")
library("reticulate")
library("testit")
library("dplyr")
library("plyr")
library("OneR")
library("stringr")
library("R.matlab")
library("ggplot2")
library('latex2exp')
library("gridExtra") 
library("dimRed")
library("scatterplot3d")
library("Rtsne")
library("viridis")
library("umap")
library("lsa")
library("reshape2")
library("zoo")

get_stim_mat <- function(path, window_size=15, just_mat=F, verbose=F, chunk=-1)
{
  behav_files <- list.files(sprintf("%s\\behavior\\",path))
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  
  runs_ind <- 1:len(actual_runs)
  
  if (all(chunk != -1)) {
    runs_ind <- runs_ind[chunk]
  }
  runs <- len(runs_ind)
  
  # cumulative_frames_per_mat <- cumsum(unlist(frames_per_mat))
  # cumulative_frames_per_mat <- cumulative_frames_per_mat - cumulative_frames_per_mat[1]
  cumulative_frames_per_mat <-   cumsum(c(0,unlist(frames_per_mat[runs_ind])))[1:runs]
  names(cumulative_frames_per_mat) <- runs_ind
  
  licking_files <- behav_files[grepl("f2p", behav_files)]
  behav_files <- behav_files[!grepl("f2p", behav_files)]
  behavior_mat_list <- lapply(behav_files,
                              function(tv_mat) {
                                return(readMat(sprintf("%s\\behavior\\%s",
                                                       path,
                                                       tv_mat))$Stim.Master)
                              })
  
  # Good luck reading this lol
  behavior_indices <- 
    unlist(lapply(unlist(lapply(str_split(behav_files, ".TrialVar"), 
                                function(sp) {
                                  sp[[1]][1]
                                })),
                  function(p) {
                    tmp <- strsplit(p, "")[[1]]; 
                    return(as.numeric(paste(tmp[(length(tmp) - 1):(length(tmp))], collapse="")))
                  }))
  
  fm <- behavior_mat_list
  ret_list <- list()
  
  
  for (i in runs_ind) {
    if (!i %in% behavior_indices) {
      ret_list <- append(ret_list, list())
      
      if (verbose) {  
        print(sprintf("Skipping behavior mat %d (spont)",i)) 
      }
      next
    }
    
    j <- which(i == behavior_indices)
    
    if (verbose) {
      print(sprintf("Adding %d frames to behavior mat %d", cumulative_frames_per_mat[as.character(i)], i))
    }
    
    fm[[j]] <- cbind(fm[[j]], c(`RunNum`=i))

    fm[[j]][,1] <- fm[[j]][,1] + cumulative_frames_per_mat[as.character(i)]
    fm[[j]][,6] <- fm[[j]][,6] +  cumulative_frames_per_mat[as.character(i)]
    
    ret_list <- append(ret_list, list(fm[[j]]))
    
  }
  
  
  # if (just_mat) {
    
    # if (chunk != -1) {
    #   stim_master_mat <- do.call(rbind,ret_list[chunk])  
    # } else {
    #   
    # }
    stim_master_mat <- do.call(rbind,ret_list)
    colnames(stim_master_mat) <- c(1:ncol(stim_master_mat))
    colnames(stim_master_mat)[c(1,3,6,7,8)] <- c("Frames",
                                                 "Grating",
                                                 "Reward",
                                                 "TrialType",
                                                 "Response")
    return(stim_master_mat)
  # }
  
  # response_ind <- unlist(lapply(ret_list, function(sdf) {
  #   if (!all(is.na(sdf))) { return(sdf[,6])}
  # }))
  # 
  # response_ind <- response_ind[!is.nan(response_ind)]
  # binned_responses <- get_binned_index(response_ind, window_size)
  # 
  # bound <- do.call(rbind, ret_list)
  # binned_by_trial_type <- 
  #   lapply(sort(unique(bound[,7])), 
  #          function(t) {get_binned_index(bound[which(bound[,7] == t),1],window_size)})
  # 
  # names(binned_by_trial_type) <- sort(unique(bound[,7]))
  # 
  # return(list(binned_response=binned_responses,
  #             binned_by_trial_type=binned_by_trial_type))
  
}

get_stim_mat_control <- function(path, just_mat=F, window_size=15)
{
  behavior_mat <- readMat(sprintf("%s\\behavior.mat", path))
  
  bound <- c()
  for(d in 1:dim(behavior_mat$mystruct)[1]) {
    
    bound <- cbind(bound,
                   unlist(behavior_mat$mystruct[[d]]))
  }
  colnames(bound) <- c("TrialType", "Frames_Broke", "Frames", "Reward_Broke", "Reward", "Response", "Grating")
  
  if(just_mat) {
    return(bound)
  }
  
  response_ind <- bound[,"Reward"]
  response_ind <- response_ind[!is.nan(response_ind)]
  binned_responses <- get_binned_index(response_ind, window_size)
  
  binned_by_trial_type <- 
    lapply(sort(unique(bound[,"TrialType"])), 
           function(t) {get_binned_index(bound[which(bound[,"TrialType"] == t),"Frames"],window_size)})
  
  names(binned_by_trial_type) <- sort(unique(bound[,"TrialType"]))
  
  return(list(binned_response=binned_responses,
              binned_by_trial_type=binned_by_trial_type))
}

 
get_pupil_files <- function(path, window_size=30, load_all=T, save=F)
{
  
  if (load_all) {
    load(file=sprintf("%s\\pupil_files_all\\all_pupil.Rda",path))
    return(pupil_res)
  }
  
  pupil_files <- list.files(sprintf("%s\\pupil_files\\",path))
  pupil_vec <- lapply(pupil_files,
                              function(pf) {
                                pm <- readMat(sprintf("%s\\pupil_files\\%s",
                                                       path,
                                                       pf))$pupil.mat
                                raw_vec <- pm[,3]
                                
                                # if (len(raw_vec) > 58000 ) {
                                #   raw_vec <- raw_vec[1:len(raw_vec) %% 2]
                                # } else if (len(raw_vec) < 30000)  {
                                #   raw_vec <- rep(raw_vec, each=2)
                                # }
                                
                   
                                
                                return(raw_vec)
                              })
    
  
  raw_pupil_vec <- unlist(pupil_vec)
  binned_vec <- time_bin_average_vec(raw_pupil_vec, window_size)
  
  if (save) {
    pupil_res = list(raw=raw_pupil_vec,
               smoothed=binned_vec)
    
    dir.create(sprintf("%s\\pupil_files_all\\", path))
    save(file=sprintf("%s\\pupil_files_all\\all_pupil.Rda",path),
         pupil_res)
  }

  
  # 
  
  # lapply(pup)
  # ### HORRIBLE PATCH
  # if (path == get_thirsty_quenched_paths()[7]) {
  #   raw_pupil_mat_list <- lapply(pupil_mat_list,
  #                                function(mt) {return(mt[,3])})
  #   
  #   raw_pupil_mat_list[[1]] <- rep(raw_pupil_mat_list[[1]], each=2)
  #   raw_pupil_vec <- unlist(raw_pupil_mat_list)
  #   binned_vec <- time_bin_average_vec(raw_pupil_vec, window_size)
  # }
  # #####
  # 
  # 
  # pupil_mat <- do.call(rbind, pupil_mat_list)
  # 
  # raw_pupil_vec <- pupil_mat[,3]
  # pupil_vec_d <- rep(raw_pupil_vec, each=2)
  # binned_vec <- time_bin_average_vec(pupil_vec_d, window_size)
  
  
  return (list(raw=raw_pupil_vec, smoothed=binned_vec))
  
}

get_behavior_variables <- function(path, window_size=15, begin=0, end=2, pupil_q=0.1, control=F)
{
  
  if (control) {
    stim_master_mat <- get_stim_mat_control(path)
  } else {
    stim_master_mat <- get_stim_mat(path)
  }
  
  frame_rate = 30
  
  mt <- get_mat_with_preset(path, "dalshtaim")
  mat_size <- max(dim(mt))
  pupil <- get_pupil_files(path, window_size = window_size)
  
  trials_ind <- sapply(stim_master_mat[,"Frames"], function(i) {get_binned_index(i, window_size)})
  non_na_ind <- !is.na(trials_ind)
  stim_master_mat <- stim_master_mat[non_na_ind,]
  
  base <- rep(0, times=mat_size)
  
  #responses <- unique(stim_master_mat[,"Response"])
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm",
                "8"="PavHit", "9"="PavMiss")
  
  
  all_variables <- list()
  
  for(resp in as.numeric(names(Responses))) {
    trials <- stim_master_mat[stim_master_mat[,"Response"] == resp,"Frames"]
    
    trial_indices <- lapply(trials,
                            function(trial_start) 
                              {(trial_start + begin * frame_rate):(trial_start + end * frame_rate - 1)})
    
    
    trial_indices_all <- lapply(trial_indices,
                                function(trial_ind) {
                                  return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, window_size)})))
                                })
    
    trial_indices_all <- unique(unlist(trial_indices_all))
    
    binnarized_trial_ind <- base
    binnarized_trial_ind[trial_indices_all] <- 1
    all_variables[[Responses[as.character(resp)]]] <- binnarized_trial_ind
  }
  
  pupil_lower_q <- quantile(pupil$smoothed, pupil_q)
  pupil_upper_q <- quantile(pupil$smoothed, 1 - pupil_q)
  
  pupil_vec <- base
  pupil_vec[pupil$smoothed >= pupil_upper_q] <- 2
  pupil_vec[pupil$smoothed <= pupil_lower_q] <- 1
  

  rewards_vec <- base
  
  rewardsAll_trials <- stim_master_mat[!is.nan(stim_master_mat[,"Reward"]),"Frames"]
  
  rewardsAll_indices <- lapply(rewardsAll_trials,
                          function(trial_start) 
                          {(trial_start + begin * frame_rate):(trial_start + end * frame_rate - 1)})
  
  rewardsAll_indices_all <- lapply(rewardsAll_indices,
                              function(trial_ind) {
                                return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, window_size)})))
                              })
  
  rewardsAll_indices_all <- unique(unlist(rewardsAll_indices_all))
  
  rewards_vec[rewardsAll_indices_all] <- 1
  
  faAll_vec <- base
  
  faAll_trials <- stim_master_mat[stim_master_mat[,"Response"] %in% c("3", "5"),"Frames"]
  
  faAll_indices <- lapply(faAll_trials,
                          function(trial_start) 
                          {(trial_start + begin * frame_rate):(trial_start + end * frame_rate - 1)})
  
  faAll_indices_all <- lapply(faAll_indices,
                              function(trial_ind) {
                                return(unlist(lapply(trial_ind, function(i) {get_binned_index(i, window_size)})))
                              })
  
  faAll_indices_all <- unique(unlist(faAll_indices_all))
  faAll_vec[faAll_indices_all] <- 1
  
  all_variables[["Pupil"]] <- pupil_vec
  all_variables[["AllRewards"]] <- rewards_vec
  all_variables[["AllFalseAlarms"]] <- faAll_vec
  
  return(list(variables=all_variables,
              stim_mater_mat=stim_master_mat))
}

pair_plots_colored_clusters <- function(mat, 
                                        cluster_labels, 
                                        col_palette=spec_cg, 
                                        stroke_alpha=0.5, 
                                        to_label=1, 
                                        sample_ind=-1,
                                        specific_dims=NA,
                                        remove_outliers=F,
                                        plot_labels=T,
                                        return_grid=T)
{
  
  
  
  clust_labels_unique <- sort(unique(cluster_labels))
  
  if (remove_outliers) {
    mat <- mat[cluster_labels != -1,]
    cluster_labels <- cluster_labels[cluster_labels != -1]
  }
  
  centroid_mat <- c()
  
  for (lbl in clust_labels_unique) {
    
    centroid_mat <- rbind(centroid_mat,
                          colMeans(mat[which(cluster_labels == lbl),]))
  }
  
  rownames(centroid_mat) <- clust_labels_unique
  
  
  clust_col <- col_palette(nrow(centroid_mat))
  names(clust_col) <- clust_labels_unique
  colors <- clust_col[as.character(cluster_labels)]
  colors[names(colors) == "-1"] <- adjustcolor(colors[names(colors) == "-1"], 0.1)
  colors[names(colors) != "-1"] <- adjustcolor(colors[names(colors) != "-1"], 0.5)
  
  if (sum(!is.na(specific_dims)) != 0 ) {
    dimensions <- combn(specific_dims, 2)
  } else {
    dimensions <- combn(ncol(mat), 2)
  }
  all_plots <- list()
  
  
  if (sample_ind != -1) {
    ind <- sort(sample(1:nrow(mat), sample_ind))
  } else {
    ind <- 1:nrow(mat)
  }
  
  for (i in 1:ncol(dimensions)) {
    
    d1 <- dimensions[1,i]
    d2 <- dimensions[2,i]
    
    
    labels=as.character(1:nrow(mat[ind,]))
    rewards_ind = which(colors[ind] == "red")
    stroke <- rep("0", len(ind))
    stroke[names(colors) != "-1"] <- "1"
    
    
    df <- data.frame(X=mat[ind,d1], Y=mat[ind, d2], col_lab=labels, pt_fill=colors[ind], stroke_color=stroke)
    

    
    gt <- ggplot(df, aes(x=X, y=Y)) + 
      geom_point(aes(fill=pt_fill, color=stroke_color), shape=21, size=1.5) +
      xlab(sprintf("Dim %d", d1)) +
      ylab(sprintf("Dim %d", d2)) + 
      theme_light() + 
      scale_color_manual(breaks=c("0", "1"), values=c(adjustcolor("red", 0.001), adjustcolor("black",  stroke_alpha))) +
      scale_fill_identity() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="NA")  
    
    if (plot_labels) {
      ctm <- as.data.frame(centroid_mat[,c(d1, d2)])
      ctm$id <- rownames(centroid_mat)
      colnames(ctm) <- c("x", "y", "id")
      
      if (i %% 2 == to_label) {
        g <- gt + geom_label(data=ctm, 
                             aes(x=x, y=y, label=id),
                             size=2.5,
                             fill=adjustcolor("white", alpha=0.75))
      } else {
        g <- gt
      }
    } else {
      g <- gt
    }
    
    
    all_plots <- append(all_plots, list(g))
  }
  
  if (return_grid) {
  
  if (sum(!is.na(specific_dims)) != 0 ) {
    all_plots$nrow <- 1
  } else {
    all_plots$nrow <- 5
  }
  
  
  a <- do.call(arrangeGrob, all_plots)
  
  return(a) }
  else {
    return(all_plots)
  }
}

pair_plots_colored <- function(mat, 
                               colors, 
                               sample_ind=-1, 
                               stroke_all=F, 
                               stroke_alpha=0.5, 
                               plot_nrow=5, 
                               return_grid=T, 
                               pt_size=1.5,
                               other_col_val="#ffc400",
                               just_annotation=F)
{
  dimensions <- combn(ncol(mat), 2)
  all_plots <- list()
  
  
  if (sample_ind != -1) {
    ind <- sort(sample(1:nrow(mat), sample_ind))
  } else {
    ind <- 1:nrow(mat)
  }
  
  for (i in 1:ncol(dimensions)) {
    
    d1 <- dimensions[1,i]
    d2 <- dimensions[2,i]
    
    
    labels=as.character(1:nrow(mat[ind,]))
    rewards_ind = which(colors[ind] == "red")
    other_inds = which(colors[ind] == other_col_val)
    stroke <- rep("0", len(ind))
    stroke[rewards_ind] <- "1"
    stroke[other_inds] <- "2"
    
    if (stroke_all) {
      stroke <- rep("1", len(ind))
    }
    
    df <- data.frame(X=mat[ind,d1], Y=mat[ind, d2], col_lab=labels, pt_fill=colors[ind], stroke_color=stroke)
    
    
    #reward_df = df[rewards_ind,]
    
    g <- ggplot(df, aes(x=X, y=Y))
    
    if (!just_annotation) {
       g <- g + geom_point(aes(fill=pt_fill, color=stroke_color), shape=21, size=pt_size)
    }
      
    g <- g + 
      xlab(sprintf("Dim %d", d1)) +
      ylab(sprintf("Dim %d", d2)) + 
      theme_light() + 
      scale_color_manual(breaks=c("0", "1", "2"),
                         values=c(adjustcolor("red", 0.001), 
                                  adjustcolor("black",  stroke_alpha),
                                  adjustcolor("black",  stroke_alpha))) +
      scale_fill_identity() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="NA") 
    
    #geom_point(data=reward_df, fill=adjustcolor("red", alpha=0.8), color=adjustcolor("black", alpha=0.35), shape=21, size=1.5)
    
    
    
    all_plots <- append(all_plots, list(g))
  }
  
  if (return_grid) {
    
      all_plots$nrow <- plot_nrow
    
    
    a <- do.call(arrangeGrob, all_plots)
    
    return(a) }
  else {
    return(all_plots)
  }
}

get_intrinsic_dim <- function(mat, top=0.96)  {
  dist_mat <- as.matrix(dist(mat))
  mus <- apply(dist_mat, 
               1, 
               function(r) {
                 mus_temp <- r[order(r)[2:3]]
                 return(mus_temp[2] / mus_temp[1]) })
  
  mus <- sort(mus)
  
  linreg_x <- log(mus[1:(floor(len(mus)) * top)])
  linreg_y <- -log(1- (0:(len(linreg_x) -1) / len(linreg_x)))
  linreg <- lm(linreg_y ~ linreg_x - 1)
  
  slope = linreg$coefficients[[1]]
  
  return(slope)
}

get_datasets_names <- function(path_list, sep = " ", control=F){
  
  if (control) {
    mice_name_indices <- unlist(gregexpr("fov[0-9]{1}", path_list))
  } else {
    mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
  }
  
  mice_names <- unlist(lapply(1:len(mice_name_indices), 
                              function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
  days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
  days <- unlist(lapply(1:len(days_indices),
                        function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
  datasets_names <- paste(mice_names, days, sep = sep)
  
  return(datasets_names)
}

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
         symbols = c("****", "***", "**", "*", "NS"))
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard_mat_sim <- function(mt) {
  
  mtrows <- nrow(mt)
  sim_mt <- matrix(rep(0,  times=mtrows ** 2), nrow=mtrows)
  cmbns <- combn(mtrows, 2)
  cmbns <- cbind(cmbns)#, rbind(1:mtrows, 1:mtrows))
  
  for (cmb_idx in 1:ncol(cmbns)) {
      col <- cmbns[,cmb_idx]
      jcd <- stringdist(paste(mt[col[1],], collapse=""), paste(mt[col[2],], collapse=""), "jw");
      #jcd <- hamming_dist(mt[col[1],], mt[col[2],]);
      sim_mt[col[2], col[1]] <- jcd
  }
  
  return(sim_mt)

}

wasserstein_dist <- function(a, b, p=1, wa=NULL, wb=NULL) {
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
    return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  stopifnot(is.null(wa) || length(wa) == m)
  stopifnot(is.null(wb) || length(wb) == n)
  if (is.null(wa)) {
    wa <- rep(1,m)
  } else { # remove points with zero weight
    wha <- wa > 0
    wa <- wa[wha]
    a <- a[wha]
    m <- length(a)
  }
  if (is.null(wb)) {
    wb <- rep(1,n)
  } else { # remove points with zero weight
    whb <- wb > 0
    wb <- wb[whb]
    b <- b[whb]
    n <- length(b)
  }
  
  orda <- order(a)
  ordb <- order(b)
  a <- a[orda]
  b <- b[ordb]
  wa <- wa[orda]
  wb <- wb[ordb]
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))  
  cub <- c(cumsum(ub))  
  arep <- hist(cub, breaks = c(-Inf, cua, Inf), plot = FALSE)$counts + 1
  brep <- hist(cua, breaks = c(-Inf, cub, Inf), plot = FALSE)$counts + 1
  # we sum over rectangles with cuts on the vertical axis each time one of the two ecdfs makes a jump
  # arep and brep tell us how many times each of the a and b data have to be repeated in order to get the points on the horizontal axis
  # note that sum(arep)+sum(brep) = m+n-1 (we do not count the height-zero final rectangle where both ecdfs jump to 1)
  
  aa <- rep(a, times=arep)
  bb <- rep(b, times=brep)
  
  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
  #  print(rbind(uu1-uu0, pmax(aa,bb)-pmin(aa,bb)))
  return(areap)
}

JSD_dist <- function(a, b) {
  nbins=20
  a_b_breaks <- seq(min(c(a,b)), max(c(a,b)), length.out=nbins+1)
  a_counts <-  hist(a, breaks=a_b_breaks, plot=F)
  b_counts <-  hist(b, breaks=a_b_breaks, plot=F)
  
  prob_counts_a <- a_counts$counts / sum(a_counts$counts)
  prob_counts_b <- b_counts$counts / sum(b_counts$counts)
  
  prob_mat <- rbind(prob_counts_a,
                    prob_counts_b)
  
  d <- JSD(prob_mat)
  names(d) <- c()
  return(d)  
}

KL_dist <- function(a, b) {
  nbins=20
  a_b_breaks <- seq(min(c(a,b)), max(c(a,b)), length.out=nbins+1)
  a_counts <-  hist(a, breaks=a_b_breaks, plot=F)
  b_counts <-  hist(b, breaks=a_b_breaks, plot=F)
  
  prob_counts_a <- a_counts$counts / sum(a_counts$counts)
  prob_counts_b <- b_counts$counts / sum(b_counts$counts)
  
  prob_mat <- rbind(prob_counts_a,
                    prob_counts_b)
  
  d <- KL(prob_mat)
  
  names(d) <- c()
  return(d)  
}

dprime_dist <- function(a , b) {
  return(abs(mean(a) - mean(b)) / (0.5 * sd(a) + 0.5 * sd(b)))
}

mse <- function(a, b) {
  return(mean((a - b) ** 2))
}

sem <- function(r) {
  ln <- sum(!is.na(r))
  if (ln <= 1) {
    return(0)
  }
  return(sqrt(var(as.numeric(r), na.rm=T))/sqrt(ln))
}

spearman_cor <- function(a,b) {cor(a,b, method="spearman")}

get_path <- function(mice_id, day) {
  base_path = "Y:\\livneh\\itayta\\data"
  
  return(sprintf("%s\\IC%d\\day_%d", base_path, mice_id, day))
  
  
}

get_control_paths <- function() {
  
  v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")
  
  por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\")
  
  return(c(v1_paths, por_paths))
}
get_SFO_paths <- function() {
  mice <- list(
               `47`=
                 c(171122),
               `49`=
                 c(180216),
               `56`=
                 c(180629),
               `57`=
                 c(180622))
  
  SFO_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    SFO_paths <- c(SFO_paths,
                                unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(SFO_paths)
}
get_satiation_paths <- function() {
  mice <- list(
    `13`=
      c(150407),
    `17`=
      c(150619),
    `19`=
      c(150916),
    `26`=
      c(151119),
    `47`=
      c(171215))
  
  hungry_sated_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    hungry_sated_paths <- c(hungry_sated_paths,
                            unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(hungry_sated_paths)
}
get_agrp_paths <- function() {
  mice <- list(
               `13`=
                 c(150508),
               `16`=
                 c(150624),
               `17`=
                 c(150620),
               `26`=
                 c(151217))
  
  hungry_sated_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    hungry_sated_paths <- c(hungry_sated_paths,
                            unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(hungry_sated_paths)
}
get_hungry_sated_paths <- function() {
  mice <- list(`19`=
                 c(150911),
               `17`=
                 c(150615, 150619),
               `16`=
                 c(150623),
               `13`=
                 c(150406, 150407, 150507),
               `26`=
                 c(151216),
               `32`=
                 c(161214),
               `42`=
                 c(161117))
  
  hungry_sated_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    hungry_sated_paths <- c(hungry_sated_paths,
                                unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(hungry_sated_paths)
}
get_hypertonic_saline_paths <- function() {
  mice <- list(`71`=
                 c(190906),
               `77`=
                 c(190908),
               `67`=
                 c(190909))
  
  thirsty_quenched_paths <- c()
  
  for (mice_id in names(mice)) {
    days <- mice[[mice_id]]
    
    thirsty_quenched_paths <- c(thirsty_quenched_paths,
                                unlist(lapply(days, function(day) {get_path(as.numeric(mice_id), day)})))
  }
  
  return(thirsty_quenched_paths)
}
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

ph <- function(mt, ...) {pheatmap(mt,cluster_rows=F, cluster_cols=F, ...)}

across_mice_decoding_build_metadata <- function(path_list=get_thirsty_quenched_paths(), chunk_list=list(), control=F, preset="dalshtaim", reward_trials_only=F)
{
  
  #path_list <- get_thirsty_quenched_paths()
  results <- list()
  
  for (idx in 1:len(path_list)) {
    
    if (len(chunk_list) > 0) {
      chunk = chunk_list[[idx]]
    } else {
      chunk = -1
    }
    
    p <- path_list[[idx]]
    
    if (preset == "Isomap") {
      all_mat <- get_reduced_mat_full_day(p, ndim = 6, type = "isomap", knn1=100, knn2=0, window_size=15)
    } else {
      all_mat <- get_mat_with_preset(p, preset=preset, chunk=chunk)
    }
    if (control) {
      stim_master_mat <- get_stim_mat_control(p,just_mat = T, window_size = 15)
    } else {
      stim_master_mat <- get_stim_mat(p,just_mat = T, window_size = 15, chunk=chunk)
    }
    clust_mat_all <- get_clusters_mat_kmeans(p, chunk=chunk, red_preset = preset)
    
    num_of_clusters <- clust_mat_all$nc
    cluster_labels <- clust_mat_all$labs
    
    relevant_trials_mat <- stim_master_mat[stim_master_mat[,"TrialType"] %in% c(1,3,4,5),]
    
    num_of_trials <- nrow(relevant_trials_mat)
    
    chronological_mat <- matrix(rep(-1, times=num_of_trials * 20), nrow=num_of_trials)
    extended_chronological_mat <- matrix(rep(-1, times=num_of_trials * 40), nrow=num_of_trials)
    prob_matrix <- matrix(rep(0, times=(num_of_clusters + 1) ** 2), nrow=(num_of_clusters + 1))
    
    colnames(prob_matrix) <- c(-1, 1:num_of_clusters)
    rownames(prob_matrix) <- c(-1, 1:num_of_clusters)
    
  
    rownames(chronological_mat) <- 1:num_of_trials
    colnames(chronological_mat) <- sprintf("t =  %.1f", 1:20 / 2)
    
    
    rownames(extended_chronological_mat) <- 1:num_of_trials
    colnames(extended_chronological_mat) <- sprintf("t =  %.1f", 1:40 / 2)
    
    trials <- 1:num_of_trials
    
    
    for (trial in trials) {
      chronological_matrices_conf = c(reg=19, ext=39)
      
      for (conf_name in names(chronological_matrices_conf)) {
        num_of_frames <- chronological_matrices_conf[[conf_name]]  
        
        trial_ind <- get_binned_index(relevant_trials_mat[trial,"Frames"], 15)
        trial_ind <- (trial_ind):(trial_ind + num_of_frames)
        
        clusters_in_trial <- cluster_labels[trial_ind]
        
        if (sum(is.na(clusters_in_trial)) > 0) {
          clusters_in_trial[is.na(clusters_in_trial)] <- 
            clusters_in_trial[!is.na(clusters_in_trial)][len(clusters_in_trial[!is.na(clusters_in_trial)])]
        }
        
        
        if (conf_name == "reg") {
          chronological_mat[trial,] <- clusters_in_trial
        } else {
          extended_chronological_mat[trial,] <- clusters_in_trial
        }
      }
    }
    
    #title <- sprintf("%s - %d Clusters", datasets_names[[idx]], num_of_clusters)
    
    
    non_re_trial = T
    for (trial_idx in 1:len(trials)) {
      
      if (reward_trials_only && relevant_trials_mat[trial_idx,"TrialType"] != 3) {
        
        if (non_re_trial){
          print("Non reward")
          non_re_trial = F
        }
        
        next
      }
      
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
    final_list$extended_trials_mat <- extended_chronological_mat
    final_list$annot_df <- annot_df
    final_list$probability_mat <- prb
    final_list$stim_master_mat <- stim_master_mat
    final_list$red_mat <- all_mat
    
    results[[idx]] <- final_list
    
    
  }
  
  return(results)
}

getmode <- function(v)
{
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get_licking_mat <- function(metadata_obj,  reward_trials, path)
{
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  run_indices <- unique(reward_trials[,ncol(reward_trials)])
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {print(run_ind);get_licking_vec(path, run_ind)})
  
  names(lick_vectors_all) <- run_indices
  
  
  lick_mat <- c()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    trial_run <- relevant_trial[len(relevant_trial)]
    relevant_lick_vec <- lick_vectors_all[[as.character(trial_run)]]
    frame_in_run <- relevant_trial["Frames"] - frames_to_add[trial_run]
    
    
    lv <- relevant_lick_vec[frame_in_run:(frame_in_run + 300)]
    lick_mat <- rbind(lick_mat, lv)
    
    
  }
  
  return(lick_mat)
}

get_licking_vec <-  function(path, lick_file_idx, verbose=F)
{
  
  search_path <- sprintf("%s//licking",path)
  licking_files = list.files(search_path)[(grep("IC.*.licking", list.files(search_path)))] 
  
  if(verbose){
    print(licking_files)
    print(lick_file_idx)
  }
  
  regexp_str = "00[%d]{1}"
  
  # Unsafe
  if (lick_file_idx > 9) {
    regexp_str = "01[%d]{1}"
    lick_file_idx = lick_file_idx - 10
  }
  
  requested_file <- licking_files[which(regexpr(sprintf(regexp_str, lick_file_idx), licking_files) != -1)]
  
  if (verbose) {
    print(requested_file)
  }
  
  lick_vec <- readMat(sprintf("%s//%s", search_path, requested_file))
  
  return(as.vector(lick_vec$lick.vec))
}

clean_metadata_HS <- function(metadata_HS)
{
  
  for (i in 1:8) {
    reward_mat <- metadata_HS[[i]]$trials_mat[metadata_HS[[i]]$annot_df[,1] == 3 & metadata_HS[[i]]$annot_df[,2] == 0,]
    phr <- pheatmap(reward_mat, cluster_rows=F, cluster_cols=F, legend = F, border_col=NA, show_colnames = F, show_rownames = F)
    
    path <- get_hungry_sated_paths()[i]
    metadata_obj <- metadata_HS[[i]]
    tt_mat_ind <- which(metadata_obj$stim_master_mat[,"TrialType"] %in% c(1,3,4,5))
    reward_indices <- which(metadata_obj$stim_master_mat[,"TrialType"] %in% c(3,4,5) & metadata_obj$stim_master_mat[,"Response"] == 0)
    reward_trials <- metadata_obj$stim_master_mat[reward_indices,]
    licking_mat <- get_licking_mat(metadata_obj,reward_trials, path)
    
    faulty_reward_trials <- reward_indices[rowSums(licking_mat[,61:120]) == 0]
    to_keep <- !(1:nrow(metadata_obj$stim_master_mat) %in% faulty_reward_trials)
    which_to_keep <- which(to_keep)
    
    tt_mat_filtered_ind <- tt_mat_ind %in% which_to_keep
    
    metadata_HS[[i]]$stim_master_mat <- metadata_HS[[i]]$stim_master_mat[to_keep,]
    metadata_HS[[i]]$trials_mat <- metadata_HS[[i]]$trials_mat[which(tt_mat_filtered_ind), ]
    metadata_HS[[i]]$annot_df <- metadata_HS[[i]]$annot_df[which(tt_mat_filtered_ind), ]
    
    
    reward_mat <- metadata_HS[[i]]$trials_mat[metadata_HS[[i]]$annot_df[,1] == 3 & metadata_HS[[i]]$annot_df[,2] == 0,]
    phr <- pheatmap(reward_mat, cluster_rows=F, cluster_cols=F, legend = F, border_col=NA, show_colnames = F, show_rownames = F)
  }
  
  return(metadata_HS)
}

get_lick_bouts <- function(lick_vec, window_size=15, interval=30, unbinned_bouts=F)
{
  
  licks_indices <- which(lick_vec > 0)
  licks_interval <- diff(licks_indices)
  lick_bouts <- licks_indices[which(licks_interval > interval)  + 1]
  
  if (!unbinned_bouts) {
    binned_bouts = unlist(lapply(lick_bouts, function(i) {get_binned_index(i, window_size)}))
  } else {
    binned_bouts <- unlist(lick_bouts)
  }
  
  return(binned_bouts)  
  
}

get_task_vs_freely_lick_distribution <- function(metadata, path, satiation_run_idx, extra_s=3, bouts_interval=30)
{
  
  
  actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
  actual_runs <- sort(actual_runs)
  
  frames_per_mat <- lapply(actual_runs,
                           function(r_path)
                           {
                             load(sprintf("%s\\%s", path, r_path))
                             return(ncol(fmat))
                           })
  
  orig_mt <- get_reduced_mat_full_day(path, window_size=15, just_original_mat = T, activity_threshold = .5)
  
  
  reward_trials <- metadata$stim_master_mat[metadata$stim_master_mat[,"TrialType"] %in% c(1,3,4,5) & metadata$stim_master_mat[,"Response"] == 0,]
  run_indices <- unique(reward_trials[,ncol(reward_trials)])
  
  frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
  
  binned_frames_per_mat <- unlist(frames_per_mat) / 15
  binned_frames_to_add <- cumsum(binned_frames_per_mat) - binned_frames_per_mat[1]
  
  lick_vectors_all <- 
    lapply(run_indices, function(run_ind) {get_licking_vec(path, run_ind)})
  
  
  names(lick_vectors_all) <- run_indices
  
  
  activity_task_mat <- c()
  task_lick_bouts_mat <- c()
  task_licking_mat <- c()
  neuronal_cluster_task_licks_dist <- list()
  lick_mat_anwywas <- c()
  
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    
    relevant_trial <- reward_trials[reward_trial_idx, ]
    run <- relevant_trial[len(relevant_trial)]
    trial_frame_offset = relevant_trial["Frames"] - frames_to_add[run]
    licks_in_trial_unbinned <- lick_vectors_all[[as.character(run)]][trial_frame_offset:(trial_frame_offset+(20 * 15 - 1))]
    
    lick_mat_anwywas <- rbind(lick_mat_anwywas, licks_in_trial_unbinned)
    
    if( sum(licks_in_trial_unbinned > 0) <= 0 ) {
      print("No licks in trial continue")
      next
    }
    
    if(all((!(which(licks_in_trial_unbinned > 0)) > 60) | (which(licks_in_trial_unbinned > 0) > 120))) {
      if (verbose) {
        print(sprintf("%d. Only anticipatory licks?", reward_trial_idx))
      }
      next
    }
    
    
    lick_in_trial_binned = which(time_bin_average_vec(licks_in_trial_unbinned, 15) > 0)
    
    lick_in_trial_binned <- unique(c(lick_in_trial_binned, max(lick_in_trial_binned): (max(lick_in_trial_binned) + extra_s)))
    
    lick_in_trial_binned <- lick_in_trial_binned[lick_in_trial_binned <= 20]
    
    cluster_labels_of_trial <- metadata$trials_mat[as.character(which(metadata$annot_df[,1] == 3 & metadata$annot_df[,2] == 0)[reward_trial_idx]),]
    
    neuronal_cluster_task_licks_dist <- append(neuronal_cluster_task_licks_dist, list(cluster_labels_of_trial[lick_in_trial_binned]))
    
    
    
    unbinned_first_lick_in_respone <- which(licks_in_trial_unbinned > 0)[which(which(licks_in_trial_unbinned > 0) > 60)[1]]
    binned_first_lick_in_response <- get_binned_index(unbinned_first_lick_in_respone, 15)
    
    if (binned_first_lick_in_response + 3 > 20) {
      binned_first_lick_in_response = 20 - 3
    }
    
    task_lick_bouts_mat <- rbind(task_lick_bouts_mat,
                                 cluster_labels_of_trial[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 8)])
    
    task_licking_mat <- rbind(task_licking_mat,
                              as.numeric(licks_in_trial_unbinned[(unbinned_first_lick_in_respone - 60):(unbinned_first_lick_in_respone + 120)]>0))
    
    global_frame <- get_binned_index(relevant_trial["Frames"], 15)
    
    global_activity <- (rowMeans(orig_mt[((global_frame):(global_frame+19))[(binned_first_lick_in_response-3):(binned_first_lick_in_response + 8)],]))
    activity_task_mat <- rbind(activity_task_mat,
                               global_activity)
  }
  
  
  
  
  freely_licking_vec <- get_licking_vec(path, satiation_run_idx)
  freely_licking_indices <- which(freely_licking_vec > 0) + frames_to_add[satiation_run_idx]
  binned_freely_licking <- unique(unlist(lapply(freely_licking_indices, function(i) {get_binned_index(i, 15)})))
  #binned_freely_licking <- unique(unlist(lapply(binned_freely_licking, function(i) {i:(i+extra_s)})))
  
  neuronal_cluster_task_licks_dist <- table(unlist(neuronal_cluster_task_licks_dist))
  neuronal_cluster_freely_licks_dist <- table(metadata$cluster_mat$labs[binned_freely_licking])
  
  
  freely_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  task_licks_cluster <- rep(0, times=len(unique(metadata$cluster_mat$labs)))
  
  names(freely_licks_cluster) <- c(-1, 1:(len(freely_licks_cluster) - 1))
  names(task_licks_cluster) <- c(-1, 1:(len(task_licks_cluster) - 1))
  
  task_licks_cluster[names(neuronal_cluster_task_licks_dist)] <- neuronal_cluster_task_licks_dist
  freely_licks_cluster[names(neuronal_cluster_freely_licks_dist)] <- neuronal_cluster_freely_licks_dist
  
  freely_binned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval)
  freely_unbinned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval, unbinned_bouts = T)
  freely_bouts_indices <-  freely_binned_bouts + binned_frames_to_add[satiation_run_idx]
  
  freely_lick_bouts_mat <- c()
  freely_licking_mat <- c()
  freely_indices_mat <- c()
  freely_activity_mat <- c()
  
  for (idx in 1:len(freely_binned_bouts)) {
    bout <- freely_binned_bouts[idx]
    bout_trial_ind <- (bout - 3):(bout + 8)
    freely_lick_bouts_mat <- rbind(freely_lick_bouts_mat,
                                   metadata$cluster_mat$labs[bout_trial_ind])
    
    freely_licking_mat <- rbind(freely_licking_mat, 
                                as.numeric(freely_licking_vec[(freely_unbinned_bouts[idx] - 60):(freely_unbinned_bouts[idx] + 120)] > 0))
    
    freely_indices_mat <- rbind(freely_indices_mat,
                                (freely_unbinned_bouts[idx] - 60):(freely_unbinned_bouts[idx] + 120))    
    
    
    
    global_free_activity <- (rowMeans(orig_mt[bout_trial_ind + binned_frames_to_add[satiation_run_idx],]))
    freely_activity_mat <- rbind(freely_activity_mat,
                                 global_free_activity)
  }
  
  
  is_done = FALSE
  tmp_freely_indices_mat = freely_indices_mat
  rownames(tmp_freely_indices_mat) <- 1:nrow(freely_indices_mat)
  ri = 1
  
  while (!is_done) {
    
    if (ri == nrow(tmp_freely_indices_mat)) {
      is_done = TRUE
      print("Done!")
      next
    }
    
    overlap = sum(tmp_freely_indices_mat[ri,] %in% tmp_freely_indices_mat[ri+1,])
    if (overlap > 0) {
      print(sprintf("Removing row %d due to %d overlap to row %d", ri+1, overlap, ri))
      tmp_freely_indices_mat <- tmp_freely_indices_mat[-(ri+1),]
    } else {
      ri <- ri +1
    }
  }
  
  #keep = rowSums(apply(freely_indices_mat, 1, function(r) {apply(freely_indices_mat, 1, function(r2) {sum(max(r) %in% r2)})})) <= 1
  keep = as.numeric(rownames(tmp_freely_indices_mat))
  
  freely_licking_mat <- freely_licking_mat[keep,]
  freely_lick_bouts_mat <- freely_lick_bouts_mat[keep,]
  freely_activity_mat <- freely_activity_mat[keep,]
  
  return(list(freely=freely_licks_cluster,
              task=task_licks_cluster,
              freely_licking_mat=freely_licking_mat,
              task_licking_mat=task_licking_mat,
              freely_cluster_mat=freely_lick_bouts_mat,
              task_cluster_mat=task_lick_bouts_mat,
              unaligned_licks=lick_mat_anwywas,
              activity_freely=freely_activity_mat,
              activity_task=activity_task_mat))
}

combined_max <- function(vec, combs = c(1:2))
{
  vec <- c(vec)
  ret <- max(unlist(lapply(combs,
                           function(nc) {
                             return(max(colSums(combn(vec, nc))))
                           })))
  
  return(ret)
}

cross_tabulate <- function(vec1, vec2)
{
  
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

calculate_MI_from_cont <- function(tabulated)
{
  
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

get_mode <- function(m_obj, path, axis=T, projection_only=F)
{
  
  org_mt <- get_reduced_mat_full_day(path, just_original_mat = T, window_size = 15)
  if(ncol(org_mt) > nrow(org_mt)) {org_mt <- t(org_mt)}
  
  red_mat <- m_obj$red_mat
  sm <- m_obj$stim_master_mat
  cr_sm <- sm[sm[,"Response"] == 2 | sm[,"Response"] == 4,]
  first <- cr_sm[cr_sm[,ncol(cr_sm)] == 1,]
  last <- cr_sm[cr_sm[,ncol(cr_sm)] == max(cr_sm[,ncol(cr_sm)]),]
  
  last_trials <- which(sm[,"Frames"] %in% last[,"Frames"]) + 1
  last_trials <- last_trials[last_trials <= nrow(sm)]
  
  last_final<- sm[last_trials,]
  first_final <- sm[which(sm[,"Frames"] %in% first[,"Frames"]) + 1,]
  
  last_indices <- lapply(last_final[,"Frames"], function(i) {(get_binned_index(i,15) - 6):(get_binned_index(i,15) - 1)})
  first_indices <- lapply(first_final[,"Frames"], function(i) {(get_binned_index(i,15) - 6):(get_binned_index(i,15) - 1)})
  
  last_indices <- unlist(last_indices)
  first_indices <- unlist(first_indices)
  
  
  
  #N_cells_to_remove  <- floor(len(first_avg_vec) * .15)
  
  
  first_mat <- org_mt[first_indices,]
  last_mat <- org_mt[last_indices,]
  
  first_avg_vec <- colMeans(org_mt[first_indices,])
  last_avg_vec <- colMeans(org_mt[last_indices,])
  mode_vec <- (last_avg_vec - first_avg_vec)
  
  #reward_trials <- sm[sm[,"Response"] == 0,]
  trials <- sm[sm[,"TrialType"] %in% c(3,4,5),]
  projected <- unlist(org_mt[,] %*% mode_vec) 
  
  if (axis) {
    projected <- (projected - c(last_avg_vec %*% mode_vec)) / (c(first_avg_vec %*% mode_vec) - c(last_avg_vec %*% mode_vec))
  } else {
    projected <- rowMeans(org_mt)
  }
  
  if (projection_only) {
    return(projected)
  }
  
  proj_reward_ind <- lapply(trials[,"Frames"], function(i) {projected[(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)]})
  #reward_ind <- lapply(reward_trials[trials_to_use,"Frames"], function(i) {(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)})
  
  
  all_mt <- 
    do.call(rbind,
            proj_reward_ind)
  
  return(list(proj=all_mt,
              sm=trials))
}

get_f1_scores <- function(gt, decoded_vec)
{
  
  
  mt <- matrix(rep(0,36), nrow=6)
  colnames(mt) <- 0:5
  rownames(mt) <- 0:5
  for (true in sort(unique(gt))) {
    true_labels = gt == true
    for (decoded in sort(unique(decoded_vec)))
    {
      dec = sum(decoded_vec[true_labels] == decoded)
      mt[as.character(decoded), as.character(true)] <- dec
    }
  }
  
  fp_fn_mt <- mt
  for(i in 1:6) {fp_fn_mt[i,i] <- 0}
  
  all_fp <- rowSums(fp_fn_mt)
  all_fn <- colSums(fp_fn_mt)
  all_tp <- unlist(lapply(1:6, function(i) {mt[i,i]}))
  all_precision <- all_tp / (all_tp + all_fp  + 10^-30)
  all_recall <- all_tp / (all_tp + all_fn + 10 ^-30)
  all_f1 = 2*(all_precision * all_recall) / (all_precision  + all_recall + 10^-30)
  weird_stat <- ((all_tp / colSums(mt)) + (1- (sum(mt[c("2", "4"),1]) / sum(sum(mt[c("2", "4"),]))))) / 2
  weird_stat_all <- ((all_tp / colSums(mt)) + (1- (sum(mt[c(as.character(1:5)),1]) / sum(sum(mt[as.character(1:5),]))))) / 2
  
  return(list(f1=all_f1, weird=weird_stat, weird_all=weird_stat_all))
  # mltd <- melt(confusion_matrix_list$Hit_cosine)[,1:2]
  # act <- cbind(mltd, f1_all[,1], weird_all[,1])
  # clst <- cbind(mltd, f1_all_clust[,1], weird_all_clust[,1])
  # act <- act[act$Var1 != act$Var2,]
  # clst <- clst[clst$Var1 != clst$Var2,]
  # colnames(act) <- c("Var1", "Var2", "F1", "W")
  # colnames(clst) <- c("Var1", "Var2", "F1", "W")
  # actm <- ddply(act, .(Var1), function(mdf) {colMeans(mdf[,c("F1", "W")])})
  # clstm <- ddply(clst, .(Var1), function(mdf) {colMeans(mdf[,c("F1", "W")])})
}


spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
bupu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "BuPu")))
ylgn_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "YlGn")))
prgn_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "PRGn")))
puor_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "PuOr")))
blues_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "Blues")))
rdylbu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "RdYlBu")))
ylorbr_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "YlOrRd")))
rdbu_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "RdBu")))
rdylbu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "RdYlBu")))
blues_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "Blues")))

base_plot_theme <- theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black"),
                         legend.position="NA",
                         panel.border = element_blank(),
                         panel.background = element_blank(),
                         axis.ticks.length.y=unit(-.3, "lines"),
                         axis.ticks.length.x=unit(.3, "lines"),
                         axis.ticks=element_line(colour="black"),
                         text=element_text(colour="black"))


big_text_base_plot_theme <- theme(panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  legend.position="NA",
                                  panel.border = element_blank(),
                                  panel.background = element_blank(),
                                  axis.ticks.length.y=unit(-.3, "lines"),
                                  axis.ticks.length.x=unit(.3, "lines"),
                                  axis.ticks=element_line(colour="black"),
                                  text=element_text(size=13, colour="black"),
                                  axis.text = element_text(size=13, colour="black"))

big_text_base_plot_theme_wl <- theme(panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     axis.line = element_line(colour = "black"),
                                     panel.border = element_blank(),
                                     panel.background = element_blank(),
                                     axis.ticks.length.y=unit(-.3, "lines"),
                                     axis.ticks.length.x=unit(.3, "lines"),
                                     axis.ticks=element_line(colour="black"),
                                     text=element_text(size=13, colour="black"),
                                     axis.text = element_text(size=13, colour="black"))

a4_sizes <- c(half=8.3/2,
              third=8.3/3,
              three_tenths=3 * 8.3/10,
              fourth=8.3/4,
              fifth=8.3/5,
              sixth=8.3/6)

small_a4_sizes <- c(fourth=8.3/4,
                    fifth=8.3/5,
                    sixth=8.3/6,
                    seventh=8.3/8,
                    eigth=8.3/8)