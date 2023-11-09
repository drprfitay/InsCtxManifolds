
figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\supplemental_figrures\\", base_output_path) 
dir.create(output_path)

metadata_TT <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
metadata_TT_isomap <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths(), preset="Isomap")

pair_plots_video <- function(m_obj, prefix="")
{
  
  base_figure_path <- sprintf("%s\\Reward_movie_dynamics", output_path)
  dir.create(base_figure_path)
  
  input_path <- sprintf("%s\\input\\", base_figure_path)
  output_path <- sprintf("%s\\output\\", base_figure_path)
  dir.create(input_path)
  dir.create(output_path)
  
  metadata <- across_mice_decoding_build_metadata()
  mt <- m_obj$red_mat
  #cmt <- get_colors_for_mat(path, mt, 1)
  
  
  dimensions <- combn(ncol(mt), 2)
  clpl <- spec_cg(max(m_obj$cluster_mat$labs) + 1)
  names(clpl) <- c(-1, 1:max(m_obj$cluster_mat$labs))
  f_col=clpl[as.character(m_obj$cluster_mat$labs)]
  
  #reward_vec <- as.character(as.numeric(cmt$reward == "red"))
  
  all_plots_cluster <- list()
  all_plots_gray <- list()
  
  for (i in 1:ncol(dimensions)) {
    
    
    d1 <- dimensions[1,i]
    d2 <- dimensions[2,i]
    
    
    df <- data.frame(X=mt[,d1], Y=mt[,d2], clust=f_col)
    g <- ggplot(df, aes(x=X, y=Y)) + 
      geom_point(aes(color=clust), size=.25, alpha=.15) +
      xlab(sprintf("Dim %d", d1)) +
      ylab(sprintf("Dim %d", d2)) + 
      theme_light() + 
      scale_color_identity() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="NA")
    
    
    
    g2 <- ggplot(df, aes(x=X, y=Y)) + 
      geom_point(color=adjustcolor("gray50", alpha=0.1), size=.25) +
      xlab(sprintf("Dim %d", d1)) +
      ylab(sprintf("Dim %d", d2)) + 
      theme_light() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="NA") 
    
    
    all_plots_cluster <- append(all_plots_cluster, list(g))
    all_plots_gray <- append(all_plots_gray, list(g2))
  }
  
  
  reward_trials <- m_obj$stim_master_mat[m_obj$stim_master_mat[,"Response"] == 0,]
  annoCol = list(tri=c("0"="gray30", "1"="red"))
  reward_trial_mat_clusters <- m_obj$trials_mat[m_obj$annot_df[,1] == 3 & m_obj$annot_df[,2] == 0,]
  rownames(reward_trial_mat_clusters) <- 1:nrow(reward_trial_mat_clusters)
  
  ongoing_counter = 1
  for (reward_trial_idx in 1:nrow(reward_trials)) {
    print(reward_trial_idx)
    first_reward_ind <- get_binned_index(reward_trials[reward_trial_idx,"Frames"], 15)
    trial_indices <- first_reward_ind:(first_reward_ind + 19)
    
    
    
    
    trial_annot_df <- data.frame(tri=rep(0, times=nrow(reward_trial_mat_clusters)))
    rownames(trial_annot_df) <- 1:nrow(reward_trial_mat_clusters)
    trial_annot_df[reward_trial_idx,] <- 1
    
    
    ph_trial <- ph(reward_trial_mat_clusters, 
                   col=clpl, 
                   border_col=NA, 
                   legend=F, 
                   show_rownames=F, 
                   show_colnames=F, 
                   annotation_row=trial_annot_df, 
                   annotation_names_row = F, 
                   annotation_legend=F, 
                   annotation_color=annoCol)
    
    for (ti in trial_indices) {
      
      all_plots_cluster_with_points <- list()
      all_plots_gray_with_points <- list()
      
      for (i in 1:ncol(dimensions)) {
      
      d1 <- dimensions[1,i]
      d2 <- dimensions[2,i]
      
      pt_df <- data.frame(X=mt[ti,d1],
                          Y=mt[ti,d2])
      
      
      gf <- 
        all_plots_cluster[[i]] +
        geom_point(data=pt_df, col="red", size=3)
      
      gf2 <- 
        all_plots_gray[[i]] +
        geom_point(data=pt_df, col="red", size=3)
      
      all_plots_cluster_with_points <- append(all_plots_cluster_with_points, list(gf))
      all_plots_gray_with_points <- append(all_plots_gray_with_points, list(gf2))
      } 
    
    
    
    
    all_plots_cluster_with_points$nrow = 5
    all_plots_gray_with_points$nrow = 5
    final_f_cluster <- do.call(plot_grid,all_plots_cluster_with_points)
    final_f_gray <- do.call(plot_grid,all_plots_gray_with_points)
    
    final_all <- plot_grid(ph_trial[[4]], final_f_cluster, final_f_gray, rel_widths = c(1,1.5,1.5), nrow=1)
    png(sprintf("%s\\image_%d.png", input_path, ongoing_counter), width = 1250, height = 720, res = 58)
    ongoing_counter <- ongoing_counter + 1
    plot(final_all)
    dev.off()
    
    }
  }
  
  
  png_files <- sprintf("%s\\image_%d.png", input_path, 1:((nrow(reward_trials) * 20)))
  av::av_encode_video(png_files, sprintf('%s\\output_5hz.mp4', output_path), framerate = 5)
  av::av_encode_video(png_files[101:181], sprintf('%s\\output_2hz.mp4', output_path), framerate = 2)
  av::av_encode_video(png_files[501:521], sprintf('%s\\output_1hz.mp4', output_path), framerate = 1)
  av::av_encode_video(png_files, sprintf('%s\\output_10hz.mp4', output_path), framerate = 10)
  
}

modes_analysis <- function()
{
  
  
  meta <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
  all_all_mt <- c()
  all_all_mt_b <- c()
  all_proj <- list()
  all_proj_ind <- list()
  for (idx in 1:14) {
    
    org_mt <- get_reduced_mat_full_day(get_thirsty_quenched_paths()[idx], just_original_mat = T, window_size = 15)
    ti <- 1:ncol(org_mt)
    if(ncol(org_mt) > nrow(org_mt)) {org_mt <- t(org_mt)}
    
    m_obj <- meta[[idx]]
    
    sm <- m_obj$stim_master_mat
    cr_sm <- sm[sm[,"Response"] == 2 | sm[,"Response"] == 4,]
    first <- cr_sm[cr_sm[,ncol(cr_sm)] == 1,]
    last <- cr_sm[cr_sm[,ncol(cr_sm)] == max(cr_sm[,ncol(cr_sm)]),]
    
    last_trials <- which(sm[,"Frames"] %in% last[,"Frames"]) + 1
    last_trials <- last_trials[last_trials <= nrow(sm)]
    
    last_final<- sm[last_trials,]
    first_final <- sm[which(sm[,"Frames"] %in% first[,"Frames"]) + 1,]
    
    #last_indices <- lapply(last_final[,"Frames"], function(i) { (i-91):(i - 1) }) # 
    #first_indices <- lapply(first_final[,"Frames"], function(i) { (i-91):(i - 1) }) # 
    
    last_indices <- lapply(last_final[,"Frames"], function(i) {(get_binned_index(i,15) - 6):(get_binned_index(i,15) - 1)})
    first_indices <- lapply(first_final[,"Frames"], function(i) {(get_binned_index(i,15) - 6):(get_binned_index(i,15) - 1)})
    
    last_indices <- unlist(last_indices)
    first_indices <- unlist(first_indices)
    #org_mt <- t(org_mt)
    
    
    #N_cells_to_remove  <- floor(len(first_avg_vec) * .15)
    
    
    first_mat <- org_mt[first_indices,]
    last_mat <- org_mt[last_indices,]
    
    first_avg_vec <- colMeans(org_mt[first_indices,ti])
    last_avg_vec <- colMeans(org_mt[last_indices,ti])
    
    #first_to_keep <- order(first_avg_vec)[N_cells_to_remove:(len(first_avg_vec) - N_cells_to_remove)]
    #last_to_keep <- order(last_avg_vec)[N_cells_to_remove:(len(last_avg_vec) - N_cells_to_remove)]
    
    #tmp <- rep(0, times=ncol(org_mt))
    #tmp[first_to_keep] <- first_avg_vec[first_to_keep]
    #first_avg_vec <- tmp
    #first_avg_vec[first_to_keep] <- 0
    #tmp <- rep(0, times=ncol(org_mt))
    #tmp[last_to_keep] <- last_avg_vec[last_to_keep]
    #last_avg_vec <- tmp
    #last_avg_vec[last_to_keep] <- 0
    
    
    mode_vec <- (last_avg_vec - first_avg_vec)
    mode_vec_b <- (first_avg_vec - last_avg_vec)
    
    
    reward_trials <- sm[sm[,"Response"] == 0,]
    
    
    #projected <- unlist(org_mt %*% mode_vec) / c(mode_vec %*% mode_vec)
    
    projected <- unlist(org_mt[,ti] %*% mode_vec) 
    
    projected <- (projected - c(last_avg_vec %*% mode_vec)) / (c(first_avg_vec %*% mode_vec) - c(last_avg_vec %*% mode_vec))
    
    projected_b <- unlist(org_mt[,ti] %*% mode_vec_b) 
    
    projected_b <- (projected_b - c(last_avg_vec %*% mode_vec_b)) / (c(first_avg_vec %*% mode_vec_b) - c(last_avg_vec %*% mode_vec_b))
    

    #proj_reward_ind <- lapply(reward_trials[,"Frames"], function(i) {projected[(i-60):(i+240)]}) #
    #proj_reward_ind_b <- lapply(reward_trials[,"Frames"], function(i) {projected_b[(i-60):(i+240)]}) #
    proj_reward_ind <- lapply(reward_trials[,"Frames"], function(i) {projected[(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)]})
    proj_reward_ind_b <- lapply(reward_trials[,"Frames"], function(i) {projected_b[(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)]})
    
    all_mt <- 
    do.call(rbind,
            proj_reward_ind)
    
    all_mt_b <- 
      do.call(rbind,
              proj_reward_ind_b)
    plot(colMeans(all_mt))
    
    all_proj <- append(all_proj, list(projected))
    all_proj_ind <- append(all_proj_ind, list((lapply(reward_trials[1:45,"Frames"], function(i) {(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)}))))
    
    plot(unlist(proj_reward_ind), type="l")
    plot(colMeans(all_mt), type="l")
    all_all_mt <- rbind(all_all_mt, all_mt)
    all_all_mt_b <- rbind(all_all_mt_b, all_mt_b)
  }
  
  
  
}

dimensions_estimation <- function()
{
  
  base_figure_path <- sprintf("%s\\Full_dimensionality_estimation", output_path)
  dir.create(base_figure_path)
  
  output_name = "TT_calculated_dimensionality_estimation_all.Rda"
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  
  if (create) {
  paths <- get_thirsty_quenched_paths()
  
  dim_df <- data.frame()
  for (path_idx in 1:len(paths)) {
    
    p <- paths[[path_idx]]
    original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, just_original_mat = T, control = F)
    shuffled_original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, just_original_mat = T, control = F, shuffled = T, time_shuffled = F)
    reg_mat_20 <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=F, time_shuffled=F)
    shuffled_mat_20 <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=T, time_shuffled=F)
    reg_mat_30 <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=F, time_shuffled=F)
    shuffled_mat_30 <- get_reduced_mat_full_day(p, "lem", ndim=30, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=T, time_shuffled=F)
    reg_mat_40 <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=F, time_shuffled=F)
    shuffled_mat_40 <- get_reduced_mat_full_day(p, "lem", ndim=40, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=T, time_shuffled=F)
    reg_mat_100 <- get_reduced_mat_full_day(p, "lem", ndim=100, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=F, time_shuffled=F)
    shuffled_mat_100 <- get_reduced_mat_full_day(p, "lem", ndim=100, window_size=15, knn1=0.075, knn2 = 0, control=F, shuffled=T, time_shuffled=F)
    #reg_id <- get_intrinsic_dim(reg_mat)
    
    shuff_org_id <- twonn(shuffled_original_mat, c_trimmed = 0.04, method="linfit")
    org_id <- twonn(original_mat, c_trimmed = 0.04, method="linfit")
    reg_id_20 <- twonn(reg_mat_20, c_trimmed = 0.04, method="linfit")
    shuffled_id_20 <- twonn(shuffled_mat_20, c_trimmed = 0.04, method="linfit")
    reg_id_30 <- twonn(reg_mat_30, c_trimmed = 0.04, method="linfit")
    shuffled_id_30 <- twonn(shuffled_mat_30, c_trimmed = 0.04, method="linfit")
    reg_id_40 <- twonn(reg_mat_40, c_trimmed = 0.04, method="linfit")
    shuffled_id_40 <- twonn(shuffled_mat_40, c_trimmed = 0.04, method="linfit")
    reg_id_100 <- twonn(reg_mat_100, c_trimmed = 0.04, method="linfit")
    shuffled_id_100 <- twonn(shuffled_mat_100, c_trimmed = 0.04, method="linfit")
    
    dim_df <- rbind(dim_df,
                    data.frame(name=datasets_names[path_idx],
                               n_neurons=min(dim(original_mat)),
                               n_timepoints=max(dim(original_mat)),
                               org_id=unlist(org_id$est[[2]]),
                               reg_id_100=unlist(reg_id_100$est[[2]]),
                               reg_id_40=unlist(reg_id_40$est[[2]]),
                               reg_id_30=unlist(reg_id_30$est[[2]]),
                               reg_id_20=unlist(reg_id_20$est[[2]]),
                               shuff_org_id=unlist(shuff_org_id$est[[2]]),
                               shuffled_id_100=unlist(shuffled_id_100$est[[2]]),
                               shuffled_id_40=unlist(shuffled_id_40$est[[2]]),
                               shuffled_id_30=unlist(shuffled_id_30$est[[2]]),
                               shuffled_id_20=unlist(shuffled_id_20$est[[2]])))
    
    
    plot(c(unlist(dim_df[path_idx,4:8])), xlab="", ylab="Dim", ylim=c(min(dim_df[path_idx,4:13]), max(dim_df[path_idx,4:13])), type="l")
    lines(c(unlist(dim_df[path_idx,9:13])), xlab="", ylab="Dim", ylim=c(min(dim_df[path_idx,4:13]), max(dim_df[path_idx,4:13])), type="l")
  }
  save(file=sprintf("%s\\data\\supplemental\\%s", base_output_path, output_name), dim_df)
  }
  else {
    load(sprintf("%s\\data\\supplemental\\%s", base_output_path, output_name), verbose=T)
  }
  
  reg_df <- dim_df[,4:8]
  shuff_df <- dim_df[,9:13]
  colnames(reg_df) <- c("Full", "100", "40", "30", "20")
  colnames(shuff_df) <- c("Full", "100", "40", "30", "20")
  
  melted_reg <- melt(reg_df)
  melted_shuff <- melt(shuff_df)
  
  melted_reg$Group <- "Data"
  melted_shuff$Group <- "Shuffle"
  
  final_df <- as.data.frame(rbind(melted_reg, melted_shuff))
  colnames(final_df)[1:2] <- c("Emb", "Intr")
  
  gdim <- 
  ggplot(final_df) + 
    geom_line(aes(x=Emb, y=Intr, group=Group, color=Group), stat="summary") + 
    geom_ribbon(aes(x=Emb, y=Intr, group=Group, fill=Group), fun.data="mean_sdl", stat="summary", alpha=.1, color=NA) +
    big_text_base_plot_theme_wl +
    geom_hline(yintercept=6, linetype="dashed") +
    theme(legend.position="top") +
    ylab("Estimated dimension") +
    xlab("Embedding dimension")
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_dimensionality_all.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gdim)
    
    dev.off()
    
  }
}

dimensions_estimation_simulation <- function()
{
  
  base_figure_path <- sprintf("%s\\Dimensionality_estimation_simulation", output_path)
  dir.create(base_figure_path)

  
  tmp2d <- cbind(runif(1000),runif(1000));
  twonn(tmp2d, c_trimmed = 0.01, method="linfit")
  
  tmp1d <- cbind(rep(0.5, 250) + rnorm(250, sd=.01), runif(250, min=0, max=100))
  twonn(tmp1d, c_trimmed = 0.005, method="linfit")
  
  
  simulations <- list(list(name="1D",
                           plot_slope=1,
                           data=tmp1d,
                           Nsamp=200,
                           q=0.995),
                      list(name="2D",
                           plot_slope=2,
                           data=tmp2d,
                           Nsamp=300,
                           q=.99))
  
  for(sim_conf in simulations){
    
    data_to_use = sim_conf$data
    name = sim_conf$name
    plot_slope = sim_conf$plot_slop
    Nsamp = sim_conf$Nsamp
    qt <- sim_conf$q
    
    
    dist_mat <- as.matrix(dist(data_to_use))
    mus <- apply(dist_mat, 
                 1, 
                 function(r) {
                   mus_temp <- r[order(r)[2:3]]
                   return(mus_temp[2] / mus_temp[1]) })
    
    mus <- sort(mus)
    
    linreg_x <- log(mus[1:(floor(len(mus)) * qt)])
    linreg_y <- -log(1- (0:(len(linreg_x) -1) / len(linreg_x)))
    linreg <- lm(linreg_y ~ linreg_x - 1)
    
    slope = linreg$coefficients[[1]]
    
    ind <- sort(sample(1:len(linreg_x), Nsamp))
    
    plot_df <- data.frame(x=linreg_x[ind],
                          y=linreg_y[ind],
                          i=rep(idx, Nsamp))
    
    
    x_to_use <- seq(min(plot_df$x[!is.infinite(plot_df$x)]),
                    max(plot_df$x[!is.infinite(plot_df$x)]),
                    length.out=10)
    
    y_to_use <- seq(min(plot_df$y[!is.infinite(plot_df$y)]),
                    max(plot_df$y[!is.infinite(plot_df$y)]),
                    length.out=10)
    
    max_x <- max(x_to_use)
    min_x <- min(x_to_use)
    max_y <- max(y_to_use)
    min_y <- min(y_to_use)
    
    
    
    
    b = max_y - plot_slope * min_x
    old_min_x <- min_x
    get_x_from_y <- function(y) {(y - b) / plot_slope}
    min_x <- get_x_from_y(min_y)
    x_to_use <- seq(min_x,
                    max_x,
                    length.out=20)
    
    gdim_all <- ggplot() +
      theme_classic() +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0))
    
    
    
    for(intercept_x in x_to_use) {
      b = min_y -(plot_slope * intercept_x )
      get_x_from_y <- function(y) {(y - b) / plot_slope}
      get_y_from_x <- function(x) {(plot_slope * x) + b}
      
      if (intercept_x < old_min_x) {
        xt <- c(old_min_x, get_x_from_y(max_y))
        yt <- c(get_y_from_x(old_min_x), max_y)
        
        if(xt[2] > max_x) {
          xt[2] <- max_x
          yt[2] <- get_y_from_x(max_x)
        }
      } else {
        
        xt <- c(intercept_x, max_x)
        yt <- c(min_y, get_y_from_x(max_x))
        
        if (yt[2] > max_y) {
          yt[2] <- max_y
          xt[2]  <- get_x_from_y(max_y)
        }
        print(xt)
        print(yt)
      }
      
      
      slope_lines_df <- data.frame(x=xt,
                                   y=yt)
      
      
      gdim_all <- gdim_all +
        geom_line(data=slope_lines_df, aes(x=x,y=y), linetype="dashed", col="gray80")
      
    }
    
    colnames(data_to_use) <- c("y", "x")
    gscatter <- 
      ggplot(as.data.frame(data_to_use), aes(x=x,y=y)) + 
      geom_point(color="#931A1D", alpha=.5, stroke=0, size=2) +
      big_text_base_plot_theme_wl +
      xlab("Dim 1") + 
      ylab("Dim 2")
      if (name == "1D") {
        gscatter <- 
          gscatter + ylim(c(0,1))  
      }
      
    
    gdim_f <- 
      gdim_all +
      #geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
      geom_point(data=plot_df, aes(x=x, y=y), alpha=.5, stroke=0, size=2) +
      xlab("log(distance)") + 
      ylab(sprintf("-log(1-quantile) - %.3f", slope)) + 
      big_text_base_plot_theme +
      theme(plot.title=element_text(size=10))
    
    
    gf <- plot_grid(gscatter, gdim_f, nrow=1)
    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_dim_estimation_simulation_example.pdf",
                  base_figure_path,
                  size_name,
                  name),
          height=size,
          width=2 * size)
      
      plot(gf)
      
      dev.off()
      
    }
  }
}

axes_vs_structure <- function()
{
  
  base_figure_path <- sprintf("%s\\axes_vs_structure", output_path)
  dir.create(base_figure_path)
  
  datasets_names <- get_datasets_names(get_thirsty_quenched_paths(), sep="_")
  
  confs <- 
  list(list(mouse_idx=1,
            name=datasets_names[[1]],
            plane_idx=3,
            dim1=1,
            dim2=4,
            trials_to_use=5:14),
       list(mouse_idx=5,
            name=datasets_names[[5]],
            plane_idx=5,
            dim1=1,
            dim2=6,
            trials_to_use=58:67),
       list(mouse_idx=10,
            name=datasets_names[[10]],
            plane_idx=4,
            dim1=1,
            dim2=5,
            trials_to_use=1:9))
  
  meta <- across_mice_decoding_build_metadata(get_thirsty_quenched_paths())
  for (conf in confs) {
    
    idx <- conf$mouse_idx
    name <- conf$name
    dim1 <- conf$dim1
    dim2 <- conf$dim2
    trials_to_use <- conf$trials_to_use
    plane_idx <- conf$plane_idx
    
    org_mt <- get_reduced_mat_full_day(get_thirsty_quenched_paths()[idx], just_original_mat = T, window_size = 15)
    
    if(ncol(org_mt) > nrow(org_mt)) {org_mt <- t(org_mt)}
    
    m_obj <- meta[[idx]]
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
  
    reward_trials <- sm[sm[,"Response"] == 0,]
    projected <- unlist(org_mt[,] %*% mode_vec) 
    

    projected <- (projected - c(last_avg_vec %*% mode_vec)) / (c(first_avg_vec %*% mode_vec) - c(last_avg_vec %*% mode_vec))
    proj_reward_ind <- lapply(reward_trials[,"Frames"], function(i) {projected[(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)]})
    reward_ind <- lapply(reward_trials[trials_to_use,"Frames"], function(i) {(get_binned_index(i, 15)):(get_binned_index(i, 15) + 19)})
    
    
    all_mt <- 
      do.call(rbind,
              proj_reward_ind)
    
    p_lem2_list <- pair_plots_colored(red_mat, 
                                      rep(adjustcolor("gray50", alpha=.3), times=nrow(red_mat)),
                                      return_grid = F,
                                      pt_size=.55, 
                                      just_annotation = F)
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    
    selected_plane <- p_lem2_list[[plane_idx]]
    
    axis_plot <- 
      ggplot() +
      big_text_base_plot_theme + 
      xlab("Time (s)") +
      ylab("Thirsty-Quenched axis")
    
    for (r_idx in 1:len(reward_ind)) {
      
      global_ind <- reward_ind[[r_idx]]
      trial_df <- data.frame(x=red_mat[global_ind,dim1],
                             y=red_mat[global_ind,dim2],
                             id=viridis(20),
                             val=all_mt[r_idx,])
      
      
      for (tp in 1:(nrow(trial_df) - 1)) {
        tp_df <- trial_df[tp:(tp+1),]
        selected_plane <- 
        selected_plane + 
          geom_line(data=tp_df, aes(x=x,y=y), color=adjustcolor("gray20", alpha=.5))
         
      }
      selected_plane <- 
      selected_plane + 
        geom_point(data=trial_df, aes(x=x,y=y, color=id)) +
        scale_color_identity()
      
      
      axis_plot <- 
        axis_plot +
        geom_point(data=trial_df,
                   aes(x=1:20 * 0.5, y=val, color=id)) +
        geom_line(data=trial_df,
                  aes(x=1:20 * 0.5, y=val),
                  color=adjustcolor("gray20", alpha=.25)) + 
        
        scale_color_identity()
    }
    
    all_df <- data.frame(val=colMeans(all_mt))
    
    axis_plot <- 
      axis_plot +
      geom_line(data=all_df,
                aes(x=1:20 * 0.5, y=val),
                color=adjustcolor("gray20", alpha=.75),
                size=1.5) + 
      
      scale_color_identity()
    global_ind <- reward_ind[[r_idx]]
    
    thirsty_df <- data.frame(x=red_mat[first_indices,dim1],
                             y=red_mat[first_indices,dim2])
    avg_thirsty_vec <- data.frame(x=colMeans(thirsty_df)[1], y=colMeans(thirsty_df)[2])
    
    quenched_df <- data.frame(x=red_mat[last_indices,dim1],
                              y=red_mat[last_indices,dim2])
    avg_quenched_vec <- data.frame(x=colMeans(quenched_df)[1], y=colMeans(quenched_df)[2])
    
    
    thirsty_plane <- 
      p_lem2_list[[plane_idx]] + 
      geom_point(data=thirsty_df, aes(x=x,y=y), color="blue", alpha=.15) +
      geom_point(data=avg_thirsty_vec, aes(x=x,y=y),color="black", alpha=.5, size=5)
      
    
    quenched_plane <-   
      p_lem2_list[[plane_idx]] + 
        geom_point(data=quenched_df, aes(x=x,y=y), color="red", alpha=.15) +
        geom_point(data=avg_quenched_vec, aes(x=x,y=y),color="black", alpha=.5, size=5)
    
    
    #gf <- plot_grid(thirsty_plane, quenched_plane, selected_plane, axis_plot, nrow=1)
    
    
    cluster_color_label <- spec_cg(len(unique(m_obj$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(m_obj$cluster_mat$labs)) - 1))
    col_clusters <- cluster_color_label[as.character(m_obj$cluster_mat$labs)]

    
    
    p_lem2_list <- pair_plots_colored(red_mat, 
                                      col_clusters,
                                      return_grid = F,
                                      pt_size=.55, 
                                      stroke_alpha = 0.0001,
                                      just_annotation = F)  
    
    
    
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    sat_clusters <- identify_satiety_clusters(m_obj)
    centroid <- colMeans(red_mat)
    satiety_1_ind <- which(m_obj$cluster_mat$labs == sat_clusters$satieity_1)
    satiety_2_ind <- which(m_obj$cluster_mat$labs == sat_clusters$satiety_2)
    reward_1_ind <- which(m_obj$cluster_mat$labs == sat_clusters$reward_1)
    reward_2_ind <- which(m_obj$cluster_mat$labs == sat_clusters$reward_2)
    distances_from_centroid <- apply(red_mat, 1, function(rl) {euc_dist(rl, centroid)})
    max_tp_satiety_1 <- satiety_1_ind[which.max(distances_from_centroid[satiety_1_ind])]
    max_tp_satiety_2 <- satiety_2_ind[which.max(distances_from_centroid[satiety_2_ind])]
    max_tp_reward_1 <- reward_1_ind[which.max(distances_from_centroid[reward_1_ind])]
    max_tp_reward_2 <- reward_2_ind[which.max(distances_from_centroid[reward_2_ind])]
    
    
    tps_of_interest <- c(max_tp_satiety_1,
                         max_tp_satiety_2,
                         max_tp_reward_1,
                         max_tp_reward_2)
    
    clustered_plane <- p_lem2_list[[plane_idx]]
    
    annot_df <- data.frame(x=red_mat[tps_of_interest, dim1],
                           y=red_mat[tps_of_interest, dim2],
                           label=c("Satiety 1", "Satiety2", "Reward 1", "Reward 2"))
    # clustered_plane <- 
    #   clustered_plane + 
    #   geom_label(data=annot_df, aes(x=x, y=y, label=label))
    
    gf <- plot_grid(clustered_plane, thirsty_plane, quenched_plane, selected_plane, axis_plot, nrow=1)
    
    
    for (size_name in names(a4_sizes)) {
      
      dir.create(sprintf("%s\\%s\\", base_figure_path, name))
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s\\%s_all.pdf",
                  base_figure_path,
                  name,
                  size_name),
          height=size,
          width=size * 4)
      
      plot(gf)
      
      dev.off()
      
      
      pdf(sprintf("%s\\%s\\%s_clustered_plane.pdf",base_figure_path,name,size_name),height=size,width=size)
      plot(clustered_plane)
      dev.off()
      
      png(sprintf("%s\\%s\\%s_clustered_plane.png",base_figure_path,name,size_name),height=size,width=size, unit="in", res=300)
      plot(clustered_plane)
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s_thirsty_plane.pdf",base_figure_path,name,size_name),height=size,width=size)
      plot(thirsty_plane)
      dev.off()      
      
      png(sprintf("%s\\%s\\%s_thirsty_plane.png",base_figure_path,name,size_name),height=size,width=size, unit="in", res=300)
      plot(thirsty_plane)
      dev.off()      
      
      pdf(sprintf("%s\\%s\\%s_quenched_plane.pdf",base_figure_path,name,size_name),height=size,width=size)
      plot(quenched_plane)
      dev.off()
      
      png(sprintf("%s\\%s\\%s_quenched_plane.png",base_figure_path,name,size_name),height=size,width=size, unit="in", res=300)
      plot(quenched_plane)
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s_selected_plane.pdf",base_figure_path,name,size_name),height=size,width=size)
      plot(selected_plane)
      dev.off()      
      
      png(sprintf("%s\\%s\\%s_selected_plane.png",base_figure_path,name,size_name),height=size,width=size, unit="in", res=300)
      plot(selected_plane)
      dev.off()      
      
      pdf(sprintf("%s\\%s\\%s_axis_plot.pdf",base_figure_path,name,size_name),height=size,width=size)
      plot(axis_plot)
      dev.off()            
      
    }

  }

}

identify_satiety_clusters <- function(m_obj)
{
  used_trials <- m_obj$stim_master_mat[,"TrialType"] %in% c(1,3,4,5)
  stim_master_mat <- m_obj$stim_master_mat
  cumsum_all <- rep(0, times=nrow(stim_master_mat))
  drink_ind <- which(!is.nan(stim_master_mat[,6]))
  cumsum_all[drink_ind] <- 1
  
  cumsum_f <- cumsum(cumsum_all)
  
  used_cumsum <- cumsum_f[used_trials] / max(cumsum_f[used_trials])

  eighty_q<- which(used_cumsum >= .8)[1]
  oneh_q <- sum(used_trials)
  
  fifth_f <- m_obj$trials_mat[(eighty_q + 1):oneh_q,15:20]
  tabled_fifth_f <-  rep(0, times=8)
  names(tabled_fifth_f) <- c(-1,1:7)
  tabled_fifth_f[names(table(fifth_f))] <- table(fifth_f)
  tabled_reward <- table(m_obj$trials_mat[m_obj$annot_df[,1] == 3 & m_obj$annot_df[,2] == 0,8:13])
  reward_cl1 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[1]
  reward_cl2 <- names(sort(tabled_reward, decreasing=T)[names(sort(tabled_reward, decreasing=T)) != -1])[2]
  cl1 <- names(sort(tabled_fifth_f, decreasing=T)[names(sort(tabled_fifth_f, decreasing=T)) != -1])[1]
  cl2 <- names(sort(tabled_fifth_f, decreasing=T)[names(sort(tabled_fifth_f, decreasing=T)) != -1])[2]
  
  
  return (list(satieity_1=cl1,
               satiety_2=cl2,
               reward_1=reward_cl1,
               reward_2=reward_cl2))
  
}

licking_rate_lick_trials <- function()
{
  
  base_figure_path <- sprintf("%s\\licking_rate_lick_trials", output_path)
  dir.create(base_figure_path)
  
  metadata_satiation <- across_mice_decoding_build_metadata(get_satiation_paths())
  satiation_run_indices_all = c(5, 4, 3, 3, 3)
  satiation_run_indices_interval <- c(30, 16, 30, 16,30)
  satiation_paths <- get_satiation_paths()
  datasets_names <- get_datasets_names(satiation_paths, sep="_")
  for (idx in 1:len(metadata_satiation)) {
    path <- satiation_paths[idx]
    satiation_run_idx <- satiation_run_indices_all[idx]
    bouts_interval <- satiation_run_indices_interval[idx]
    
    actual_runs=list.files(path)[(grep("IC.*.R", list.files(path)))]
    actual_runs <- sort(actual_runs)
    
    frames_per_mat <- lapply(actual_runs,
                             function(r_path)
                             {
                               load(sprintf("%s\\%s", path, r_path))
                               return(ncol(fmat))
                             })
    
    frames_to_add <- cumsum(unlist(frames_per_mat)) - frames_per_mat[[1]]
    
    task_licking_vec <-  get_licking_vec(path, 1)
    #task_unbinned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval, unbinned_bouts = T)
    task_unbinned_bouts <- metadata_satiation[[idx]]$stim_master_mat[metadata_satiation[[idx]]$stim_master_mat[,"Response"] == 0,"Reward"]
    freely_licking_vec <- get_licking_vec(path, satiation_run_idx)
    freely_unbinned_bouts <- get_lick_bouts(freely_licking_vec, interval=bouts_interval, unbinned_bouts = T)
    first_idx <- freely_unbinned_bouts[10] - 60
    last_idx <- freely_unbinned_bouts[20] + 200
    
    task_first_idx <- task_unbinned_bouts[1] - 60
    task_last_idx <- task_unbinned_bouts[5] + 200
    
    
    lick_df <-
      data.frame(x=(first_idx - first_idx + 1):(last_idx - first_idx + 1),
                 y=rollsum(as.numeric(freely_licking_vec >0), k=60)[first_idx:last_idx])
    
    task_lick_df <-
      data.frame(x=(task_first_idx - task_first_idx + 1):(task_last_idx - task_first_idx + 1),
                 y=rollsum(as.numeric(task_licking_vec >0), k=60)[task_first_idx:task_last_idx])
    
    glick <- 
      ggplot(lick_df) + geom_line(aes(x=x/30,y=y/2))
    
    gtasklick <- 
      ggplot(task_lick_df) + geom_line(aes(x=x/30,y=y/2))
    
    bouts <- (freely_unbinned_bouts[10:20] - first_idx +1)/30
    task_bouts <- (task_unbinned_bouts[1:5] - task_first_idx +1)/30
    
    for (bout_idx in 1:(len(bouts))) {
    
      xintercept <-  bouts[bout_idx]
      task_xintercept <- task_bouts[bout_idx]
      poly_df <- data.frame(x=c(xintercept - 2, xintercept - 2, xintercept + 4, xintercept +4),
                            y=c(0, max(lick_df$y), max(lick_df$y), 0))      
      
      task_poly_df <- data.frame(x=c(task_xintercept - 2, task_xintercept - 2, task_xintercept + 4, task_xintercept +4),
                                    y=c(0, max(task_lick_df$y), max(task_lick_df$y), 0))      
      glick <- 
      glick +
        geom_vline(xintercept =xintercept, linetype="dashed") +
        geom_polygon(data=poly_df, aes(x=x,y=y),
                     color=NA,
                     fill="#FFBE2C",
                     alpha=.2)
      
      gtasklick <- 
        gtasklick +
        geom_vline(xintercept =task_xintercept, linetype="dashed") +
        geom_polygon(data=task_poly_df, aes(x=x,y=y),
                     color=NA,
                     fill="#FFBE2C",
                     alpha=.2)
    }

    glick <- 
      glick + 
      xlab("Time (seconds)") + 
      ylab("Lick rate (Hz)") +
      big_text_base_plot_theme +
      scale_y_continuous(expand=c(0,0))
    
    gtasklick <- 
      gtasklick + 
      xlab("Time (seconds)") + 
      ylab("Lick rate (Hz)") +
      big_text_base_plot_theme +
      scale_y_continuous(expand=c(0,0))

    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_%s_lick_trials.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[idx]),
          height=size,
          width=3 * size)
      
      plot(glick)
      
      dev.off()
      pdf(sprintf("%s\\%s_%s_hit_lick_trials.pdf",
                  base_figure_path,
                  size_name,
                  datasets_names[idx]),
          height=size,
          width=3 * size)
      
      plot(gtasklick)
      
      dev.off()
      
    }
  
  }
}

freely_behaving_activity_levels <- function()
{
  
  base_figure_path <- sprintf("%s\\freely_behaving_activity_levels", output_path)
  dir.create(base_figure_path)
  
  metadata_satiation <- across_mice_decoding_build_metadata(get_satiation_paths())
  satiation_run_indices_all = c(5, 4, 3, 3, 3)
  satiation_run_indices_interval <- c(30, 16, 30, 16,30)
  satiation_paths <- get_satiation_paths()
  datasets_names <- get_datasets_names(satiation_paths, sep="_")
  
  for (i in 1:len(metadata_satiation)) {
    res = get_task_vs_freely_lick_distribution(metadata_satiation[[i]], 
                                               satiation_paths[i], 
                                               satiation_run_indices_all[i], 
                                               bouts_interval = satiation_run_indices_interval[i])
    
    
    activity_mt_all <- rbind(res$activity_task, res$activity_freely)
    annot_df <- data.frame(rep(c("Task", "Freely"), times=c(nrow(res$activity_task), nrow(res$activity_freely))))
    rownames(activity_mt_all) <- 1:nrow(activity_mt_all)
    rownames(annot_df) <- 1:nrow(activity_mt_all)
    phactall <- 
            ph(activity_mt_all, 
               col=rdylbu_cg(20), 
               breaks=seq(-.5,.85,length.out=20),
               show_rownames = F,
               annotation_row =annot_df,
               annotation_names_row = F,
               annotation_names_col = F,
               annotation_legend = F,
               border_col=NA)
    
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_%s_activity_levels_mt.pdf",
                base_figure_path,
                size_name,
                datasets_names[i]),
        height=size,
        width=size)
    
    plot(phactall[[4]])
    
    dev.off()
    
  }
  }
}

topo_illustrations <- function()
{
  library(TDAstats)
  base_figure_path <- sprintf("%s\\topo_illustrations", output_path)
  dir.create(base_figure_path)
    
  paths <- get_thirsty_quenched_paths()
  orig_reg_mat_1 <- get_mat_with_preset(paths[1], preset="dalshtaim")
  orig_reg_mat_2 <- get_mat_with_preset(paths[5], preset="dalshtaim")
  reg_mat_1 <- kmeans(orig_reg_mat_1, centers=80, iter.max=500)$centers
  reg_mat_2 <- kmeans(orig_reg_mat_2, centers=80, iter.max=500)$centers
  
  mtall <- rbind(reg_mat_1, reg_mat_2)
  
  
    
  logical_ind <- 1:nrow(mtall) %in% sample(1:nrow(mtall), nrow(reg_mat_1))
  ind1 <- which(logical_ind)
  ind2 <- which(!logical_ind)
  
  perm_phom_reg_mat_1 <- calculate_homology(mtall[ind1,])
  perm_phom_reg_mat_2 <- calculate_homology(mtall[ind2,])
  
  p_lem2_list <- pair_plots_colored(reg_mat_1, 
                                    rep(adjustcolor("gray20", alpha=1), times=nrow(reg_mat_1)),
                                    return_grid = F,
                                    pt_size=1, 
                                    just_annotation = F)
  
  p_lem2_list <- lapply(1:len(p_lem2_list), 
                        function(i) {return(p_lem2_list[[i]] + 
                                              theme(line=element_blank(),
                                                    #rect=element_rect(color="white"),
                                                    plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                    panel.border=element_rect(color="white"),
                                                    plot.background = element_rect("white")))})
  
  p_lem2_list$nrow <- 3
  p_mat_1 <- do.call(plot_grid, p_lem2_list)
  
  p_lem2_list <- pair_plots_colored(reg_mat_2, 
                                    rep(adjustcolor("gray20", alpha=1), times=nrow(reg_mat_2)),
                                    return_grid = F,
                                    pt_size=1, 
                                    just_annotation = F)
  
  p_lem2_list <- lapply(1:len(p_lem2_list), 
                        function(i) {return(p_lem2_list[[i]] + 
                                              theme(line=element_blank(),
                                                    #rect=element_rect(color="white"),
                                                    plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                    panel.border=element_rect(color="white"),
                                                    plot.background = element_rect("white")))})
  
  p_lem2_list$nrow <- 3
  p_mat_2 <- do.call(plot_grid, p_lem2_list)
  
  
  
  p_lem2_list <- pair_plots_colored(mtall[ind1,], 
                                    rep(adjustcolor("gray20", alpha=1), times=nrow(reg_mat_2)),
                                    return_grid = F,
                                    pt_size=1, 
                                    just_annotation = F)
  
  p_lem2_list <- lapply(1:len(p_lem2_list), 
                        function(i) {return(p_lem2_list[[i]] + 
                                              theme(line=element_blank(),
                                                    #rect=element_rect(color="white"),
                                                    plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                    panel.border=element_rect(color="white"),
                                                    plot.background = element_rect("white")))})
  
  p_lem2_list$nrow <- 3
  p_mat_perm_1 <- do.call(plot_grid, p_lem2_list)
  
  
  p_lem2_list <- pair_plots_colored(mtall[ind2,], 
                                    rep(adjustcolor("gray20", alpha=1), times=nrow(reg_mat_2)),
                                    return_grid = F,
                                    pt_size=1, 
                                    just_annotation = F)
  
  p_lem2_list <- lapply(1:len(p_lem2_list), 
                        function(i) {return(p_lem2_list[[i]] + 
                                              theme(line=element_blank(),
                                                    #rect=element_rect(color="white"),
                                                    plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                    panel.border=element_rect(color="white"),
                                                    plot.background = element_rect("white")))})
  
  p_lem2_list$nrow <- 3
  p_mat_perm_2 <- do.call(plot_grid, p_lem2_list)
  
  
  phom_reg_mat_1 <- calculate_homology(reg_mat_1)
  phom_reg_mat_2 <- calculate_homology(reg_mat_2)
  
  phom_perm_mat_1 <- calculate_homology(mtall[ind1,])
  phom_perm_mat_2 <- calculate_homology(mtall[ind2,])
    
  betti0_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 0,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 0,2]
  betti0_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 0,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 0,2]
    
    
  plot_df_b0 <- data.frame(life=c(betti0_1, betti0_2),
                             grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(betti0_1), len(betti0_2)))))
  betti1_1 <- phom_reg_mat_1[phom_reg_mat_1[,1] == 1,3] - phom_reg_mat_1[phom_reg_mat_1[,1] == 1,2]
  betti1_2 <- phom_reg_mat_2[phom_reg_mat_2[,1] == 1,3] - phom_reg_mat_2[phom_reg_mat_2[,1] == 1,2]
    
    
  plot_df_b1 <- data.frame(life=c(betti1_1, betti1_2),
                             grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(perm_betti1_1), len(perm_betti1_2)))))        
    
    
  
  perm_betti0_1 <- phom_perm_mat_1[phom_perm_mat_1[,1] == 0,3] - phom_perm_mat_1[phom_perm_mat_1[,1] == 0,2]
  perm_betti0_2 <- phom_perm_mat_2[phom_perm_mat_2[,1] == 0,3] - phom_perm_mat_2[phom_perm_mat_2[,1] == 0,2]
  
  
  perm_plot_df_b0 <- data.frame(life=c(perm_betti0_1, perm_betti0_2),
                           grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(betti0_1), len(betti0_2)))))
  perm_betti1_1 <- phom_perm_mat_1[phom_perm_mat_1[,1] == 1,3] - phom_perm_mat_1[phom_perm_mat_1[,1] == 1,2]
  perm_betti1_2 <- phom_perm_mat_2[phom_perm_mat_2[,1] == 1,3] - phom_perm_mat_2[phom_perm_mat_2[,1] == 1,2]
  
  
  perm_plot_df_b1 <- data.frame(life=c(perm_betti1_1, perm_betti1_2),
                           grp=c(rep(c("Dataset 1", "Dataset 2"), times=c(len(perm_betti1_1), len(perm_betti1_2)))))        
  
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
    
    geaxmple_perm_b0 <- 
      ggplot(perm_plot_df_b0, aes(x = life, fill=grp)) + 
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
    
    geaxmple_perm_b1 <- 
      ggplot(perm_plot_df_b1, aes(x = life, fill=grp)) + 
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
    
    
    gpermb0 <- ggplot(data.frame(x=rnorm(mean=-9, sd=.5, 800))) + geom_histogram(aes(x=x), color="black", binwidth=.1) + big_text_base_plot_theme_wl + scale_y_continuous(expand=c(0,0))
    gpermb1 <- ggplot(data.frame(x=rnorm(mean=-9, sd=.5, 800))) + geom_histogram(aes(x=x), color="black", binwidth=.1) + big_text_base_plot_theme_wl + scale_y_continuous(expand=c(0,0))
    for (size_name in names(a4_sizes)) {
      write_path <- sprintf("%s\\%s", base_figure_path, size_name)
      dir.create(write_path)
      
      pdf(sprintf("%s\\dataset_1_barcode.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gbarcode_dataset_1)
      dev.off()
      
      pdf(sprintf("%s\\dataset_2_barcode.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gbarcode_dataset_2)
      dev.off()
      
      pdf(sprintf("%s\\dataset_1_persist.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gexample_dataset_1)
      dev.off()
      
      pdf(sprintf("%s\\dataset_2_persist.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gexample_dataset_2)
      dev.off()
      
      pdf(sprintf("%s\\distribution_b0.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(geaxmple_b0)
      dev.off()
      
      pdf(sprintf("%s\\distribution_b1.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(geaxmple_b1)
      dev.off()
      
      pdf(sprintf("%s\\perm_distribution_b0.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(geaxmple_perm_b0)
      dev.off()
      
      pdf(sprintf("%s\\perm_distribution_b1.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(geaxmple_perm_b1)
      dev.off()
      
      pdf(sprintf("%s\\perm_dist_all_b0.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gpermb0)
      dev.off()
      
      pdf(sprintf("%s\\perm_dist_all_b1.pdf",
                  write_path),
          height=a4_sizes[[size_name]],
          width=a4_sizes[[size_name]])
      
      plot(gpermb1)
      dev.off()
      
      
      
      pdf(sprintf("%s\\dataset_1_mat.pdf",
                  write_path),
          height=a4_sizes[[size_name]] * 1.5,
          width=a4_sizes[[size_name]] * 2.5)
      
      plot(p_mat_1)
      dev.off()
      
      pdf(sprintf("%s\\dataset_2_mat.pdf",
                  write_path),
          height=a4_sizes[[size_name]] * 1.5,
          width=a4_sizes[[size_name]] * 2.5)
      
      plot(p_mat_2)
      dev.off()
      
      
      pdf(sprintf("%s\\perm_dataset_1_mat.pdf",
                  write_path),
          height=a4_sizes[[size_name]] * 1.5,
          width=a4_sizes[[size_name]] * 2.5)
      
      plot(p_mat_perm_1)
      dev.off()
      
      pdf(sprintf("%s\\perm_dataset_2_mat.pdf",
                  write_path),
          height=a4_sizes[[size_name]] * 1.5,
          width=a4_sizes[[size_name]] * 2.5)
      
      plot(p_mat_perm_2)
      dev.off()
      
    }
}

structure_example_ultra_small <- function()
{
  
  base_figure_path <- sprintf("%s\\structure_example_ultra_small", output_path)
  dir.create(base_figure_path)
  metadata <- metadata_TT
  datasets_names <- get_datasets_names(get_thirsty_quenched_paths(), sep="_")
  
  for (idx in 1:len(metadata)) {
    m_obj <- metadata[[idx]]
    reduced_mat <- m_obj$red_mat
    
    cluster_color_label <- spec_cg(len(unique(m_obj$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(m_obj$cluster_mat$labs)) - 1))
    clustered_col <- cluster_color_label[as.character(m_obj$cluster_mat$labs)]
    
    p_lem2_list <- pair_plots_colored(reduced_mat, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                      return_grid = F,
                                      pt_size=.1, 
                                      just_annotation = F)  
    
    p_lem2_color_list <- pair_plots_colored(reduced_mat, 
                                            clustered_col,
                                            return_grid = F,
                                            pt_size=.1, 
                                            just_annotation = F)  
    
    p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                            rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                            return_grid = F,
                                            pt_size=.1, 
                                            just_annotation = T)  
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    p_lem2_color_list <- lapply(1:len(p_lem2_color_list), 
                          function(i) {return(p_lem2_color_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    p_lem2_no_annot <- p_lem2_list
    p_lem2_no_annot$nrow <- 3
    pf_structure <- do.call(plot_grid, p_lem2_no_annot)
    
    p_lem2_colored_struct <- p_lem2_color_list
    p_lem2_colored_struct$nrow <- 3
    pf_colored_structure <- do.call(plot_grid, p_lem2_colored_struct)
    
    p_lem2_annot_only <- p_lem2_list_annot
    p_lem2_annot_only$nrow <- 3
    pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
    
    
  dir.create(sprintf("%s\\%s\\",
                     base_figure_path,
                     datasets_names[idx]))
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    
    png(sprintf("%s\\%s\\%s_lem_struct.png",
                base_figure_path,
                datasets_names[idx],
                size_name),
        height=size * 1,
        width=size * 1.667,
        unit="in",
        res=1500)
    
    plot(pf_structure)
    
    dev.off()
    
    png(sprintf("%s\\%s\\%s_colored_lem_struct.png",
                base_figure_path,
                datasets_names[idx],
                size_name),
        height=size * 1,
        width=size * 1.667,
        unit="in",
        res=1500)
    
    plot(pf_colored_structure)
    
    dev.off()
    
    pdf(sprintf("%s\\%s\\%s_annot.pdf",
                base_figure_path,
                datasets_names[idx],
                size_name),
        height=size * 1,
        width=size * 1.667)
    
    plot(pf_annot_only)
    
    dev.off()
    
  }
  }
}

structure_example_ultra_small_isomap <- function()
{
  
  base_figure_path <- sprintf("%s\\structure_example_ultra_small_isomap", output_path)
  dir.create(base_figure_path)
  metadata <- metadata_TT_isomap
  datasets_names <- get_datasets_names(get_thirsty_quenched_paths(), sep="_")
  
  for (idx in 1:len(metadata)) {
    m_obj <- metadata[[idx]]
    reduced_mat <- m_obj$red_mat
    
    cluster_color_label <- spec_cg(len(unique(m_obj$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(m_obj$cluster_mat$labs)) - 1))
    clustered_col <- cluster_color_label[as.character(m_obj$cluster_mat$labs)]
    
    p_lem2_list <- pair_plots_colored(reduced_mat, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                      return_grid = F,
                                      pt_size=.1, 
                                      just_annotation = F)  
    
    p_lem2_color_list <- pair_plots_colored(reduced_mat, 
                                            clustered_col,
                                            return_grid = F,
                                            pt_size=.1, 
                                            just_annotation = F)  
    
    p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                            rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                            return_grid = F,
                                            pt_size=.1, 
                                            just_annotation = T)  
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    p_lem2_color_list <- lapply(1:len(p_lem2_color_list), 
                                function(i) {return(p_lem2_color_list[[i]] + 
                                                      theme(line=element_blank(),
                                                            #rect=element_rect(color="white"),
                                                            plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                            panel.border=element_rect(color="white"),
                                                            plot.background = element_rect("white"),
                                                            axis.title=element_text(color="white")))})
    
    
    p_lem2_no_annot <- p_lem2_list
    p_lem2_no_annot$nrow <- 3
    pf_structure <- do.call(plot_grid, p_lem2_no_annot)
    
    p_lem2_colored_struct <- p_lem2_color_list
    p_lem2_colored_struct$nrow <- 3
    pf_colored_structure <- do.call(plot_grid, p_lem2_colored_struct)
    
    p_lem2_annot_only <- p_lem2_list_annot
    p_lem2_annot_only$nrow <- 3
    pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
    
    
    dir.create(sprintf("%s\\%s\\",
                       base_figure_path,
                       datasets_names[idx]))
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      
      png(sprintf("%s\\%s\\%s_isomap_struct.png",
                  base_figure_path,
                  datasets_names[idx],
                  size_name),
          height=size * 1,
          width=size * 1.667,
          unit="in",
          res=1500)
      
      plot(pf_structure)
      
      dev.off()
      
      png(sprintf("%s\\%s\\%s_colored_isomap_struct.png",
                  base_figure_path,
                  datasets_names[idx],
                  size_name),
          height=size * 1,
          width=size * 1.667,
          unit="in",
          res=1500)
      
      plot(pf_colored_structure)
      
      dev.off()
      
      pdf(sprintf("%s\\%s\\%s_annot.pdf",
                  base_figure_path,
                  datasets_names[idx],
                  size_name),
          height=size * 1,
          width=size * 1.667)
      
      plot(pf_annot_only)
      
      dev.off()
      
    }
  }
}

TT_trial_dynamics_isomap <- function()
{
  paths <- get_thirsty_quenched_paths() 
  datasets_names <- get_datasets_names(paths, sep="_")
  metadata <- metadata_TT_isomap
  
  base_figure_path <- sprintf("%s\\TT_trial_dynamics_isomap", output_path)
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
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata[[idx]]$cluster_mat$labs)) - 1))
    
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

cluster_map_illustrations <- function()
{
  
  base_figure_path <- sprintf("%s\\cluster_map_illustrations", output_path)
  dir.create(base_figure_path)
  
  
  metadata <- metadata_TT
  m_obj <- metadata[[1]]
  
  
  mat <- m_obj$red_mat
  cm <- colMeans(mat)
  out_n_arr <- floor(seq(500, nrow(mat), length.out=40))
  
  distance_mat <- apply(mat, 1, function(r) {euc_dist(r, cm)})
  ordered_distance <- order(distance_mat, decreasing = T)
  
  
  
  
  for (np in out_n_arr) {
    outmost_n <- ordered_distance[1:np]
    working_mat <- mat[outmost_n,]
    
    
    for (num_of_clusters in 1:20) {
      
      print(sprintf("Running for %d: %d", np, num_of_clusters))
      reps_mse <- c()
      labels <- c()
      
      for (i in 1:20) {
        km <- kmeans(working_mat, centers=num_of_clusters, iter.max=100)
        reps_mse <- c(reps_mse,sum(km$withinss))
        labels <- rbind(labels, km$cluster)
      }
    
      
      hc <- hclust(dist(t(labels)), method = "ward.D2")
      flabs <- cutree(hc, k=num_of_clusters)
      final_labels <- rep(-1, times=nrow(mat))
      final_labels[outmost_n] <- flabs
      
      cluster_color_label <- spec_cg(num_of_clusters + 1)
      names(cluster_color_label) <- c(-1, 1:(num_of_clusters))
      clustered_col <- cluster_color_label[as.character(final_labels)]
      
      
      p_lem2_color_list <- pair_plots_colored(mat, 
                                              clustered_col,
                                              return_grid = F,
                                              pt_size=.1, 
                                              just_annotation = F)  
      
      p_lem2_color_list <- lapply(1:len(p_lem2_color_list), 
                                  function(i) {return(p_lem2_color_list[[i]] + 
                                                        theme(line=element_blank(),
                                                              #rect=element_rect(color="white"),
                                                              plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                              panel.border=element_rect(color="white"),
                                                              plot.background = element_rect("white"),
                                                              axis.title=element_text(color="white")))})
      
      
    
      for (size_name in names(a4_sizes)) {
        size = a4_sizes[[size_name]]
        png(sprintf("%s\\%s_tp%.3f_nc%d_plane_1_example.png",
                    base_figure_path,
                    size_name,
                    np/nrow(mat),
                    num_of_clusters),
            height=size,
            width=size,
            unit="in",
            res=1500)
        
        plot(p_lem2_color_list[[8]])
        
        dev.off()
        
        png(sprintf("%s\\%s_tp%.3f_nc%d_plane_2_example.png",
                    base_figure_path,
                    size_name,
                    np/nrow(mat),
                    num_of_clusters),
            height=size,
            width=size,
            unit="in",
            res=1500)
        
        plot(p_lem2_color_list[[13]])
        
        dev.off()
        
      }
      
    }
    
    
    
  }
  
  
  cmt_ph <- ph(m_obj$cluster_mat$clust_mat, border_col=NA, show_rownames=F, show_colnames=F)
  
  radii_df <- data.frame(y=rowMeans(m_obj$cluster_mat$clust_mat),
                         x=floor(seq(500, nrow(mat), length.out=40)) / nrow(mat))
  
  cluster_df <- data.frame(y=colMeans(m_obj$cluster_mat$clust_mat),
                           x=1:20)
                           
                
  
  
  gop_clust <- 
  ggplot(cluster_df, aes(x=x, y=y)) + 
    geom_point() + 
    geom_vline(xintercept=m_obj$cluster_mat$nc,
               linetype="dashed") + 
    ylab("MSE") + 
    xlab("Number of clusters") + 
    big_text_base_plot_theme
  
  gop_radii <- 
    ggplot(radii_df, aes(x=1-x, y=y)) + 
    geom_point() + 
    geom_vline(xintercept=1-radii_df$x[m_obj$cluster_mat$nr],
               linetype="dashed") + 
    ylab("MSE") + 
    xlab("Central point cloud radius (%)") + 
    big_text_base_plot_theme
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_cluster_mat.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(cmt_ph[[4]])
    
    dev.off()
    
    pdf(sprintf("%s\\%s_optimal_radii.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gop_radii)
    
    dev.off()
    
    pdf(sprintf("%s\\%s_optimal_clust.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gop_clust)
    
    dev.off()    
    
  }
}

mutual_information_illustration <- function()
{
  
  base_figure_path <- sprintf("%s\\mutual_information_illustration", output_path)
  dir.create(base_figure_path)
  
  metadata <- metadata_TT
  max_window_size <- 3
  offset <- 6
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  calculate_entropy <- function(pmf) {- sum(pmf * log2(pmf + 10 ** -30))}
  frame_rate=30
  binning_window_size=15
  thirsty_q <- .15
  cont_breaks = 7
  n_shuffles = 150
  n_subsamples = 100
  entropy_window_size = max_window_size - offset
  Responses = c("0"="Hit", "1"="Miss", "2"="CR (N)", "3"="FA (N)", "4"="CR (A)", "5"="FA (A)",  "15"="Pupil", "16"="Consumed water")
      
  for (path_idx in c(1,5,12)) {
    
    dataset_nm <- datasets_names[[path_idx]]
    reduced_mat <- metadata[[path_idx]]$red_mat
    stim_master_mat <- metadata[[path_idx]]$stim_master_mat
  
    base <- rep(0, times=nrow(reduced_mat))
  
    
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
    
    ongoing_trials_frames = metadata[[path_idx]]$stim_master_mat[ongoing_trials + 1,1]
    ongoing_ind <- unique(unlist(lapply(ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - max_window_size * 2):(get_binned_index(i, 15)-1)})))
    
    ongoing_ind_logic <- 1:nrow(reduced_mat) %in% ongoing_ind
    ongoing_ind_no_main <- ongoing_ind_logic & no_main_indices
    
    
    MI_all <- c()
    MI_main_all <- c()
    
    iteration_MI_all <- c()
    iteration_MI_main_all <- c()
    
    work_stim_mat <- stim_master_mat
    shuffle_ongoing_trials_frames = work_stim_mat[ongoing_trials + 1,1]
    shuffle_ongoing_ind <- unique(unlist(lapply(shuffle_ongoing_trials_frames, function(i) {(get_binned_index(i, 15) - max_window_size * 2):(get_binned_index(i, 15)-1)})))
    
    
    shuffle_ongoing_ind_logic <- 1:nrow(reduced_mat) %in% shuffle_ongoing_ind
    shuffle_ongoing_ind_no_main <- shuffle_ongoing_ind_logic & no_main_indices
    ind_to_use <- 1:nrow(reduced_mat)
    ind_to_use_no_main <- ind_to_use[no_main_indices]
    ind_to_use_ongoing <- ongoing_ind
    ind_to_use_ongoing_no_main <- ongoing_ind_no_main
    
    
    indices_all <- list()
    
    for(resp in as.numeric(names(Responses))) {
      resp_name <- Responses[as.character(resp)]
      
      indices_all[[resp_name]] <- base
      
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

        
        tabulated_main <- cross_tabulate(cluster_labels,
                                         variable_indices)
        
        tabulated_no_main <- cross_tabulate(cluster_labels_no_main,
                                            variable_indices_no_main)
        

        iteration_MI_all <- c(iteration_MI_all, calculate_MI_from_cont(tabulated_no_main))
        iteration_MI_main_all <- c(iteration_MI_main_all, calculate_MI_from_cont(tabulated_main))
        
        phmain <- 
        ph(tabulated_main / sum(tabulated_main), 
           display_numbers=T, breaks=seq(0,0.2,length.out=100),
           col=rev(colorspace::heat_hcl(100)),
           show_rownames=F,
           show_colnames=F,
           main=sprintf("%.4f", log2(calculate_MI_from_cont(tabulated_main))))
        
        
        bar_col_df <- data.frame(y=rowSums(tabulated_main) / sum(tabulated_main),
                                 x=sprintf(c("%s", "not %s"), resp_name))
        
        bar_row_df <- data.frame(y=as.numeric(colSums(tabulated_main) / sum(tabulated_main)),
                                 x=as.character(c(-1,1:7)))
        
        gbar_row <- 
          ggplot(bar_row_df, aes(x=x, y=y)) + 
          geom_bar(stat="identity") + 
          scale_y_continuous(expand=c(0,0)) + 
          big_text_base_plot_theme
          
        gbar_col <- 
          ggplot(bar_col_df, aes(x=x, y=y)) + 
            geom_bar(stat="identity") + 
            scale_y_continuous(expand=c(0,0)) + 
            big_text_base_plot_theme
        
        for (size_name in names(a4_sizes)) {
          size = a4_sizes[[size_name]]
          pdf(sprintf("%s\\%s_%s_%s_prob_table.pdf",
                      base_figure_path,
                      size_name,
                      datasets_names[path_idx],
                      resp_name),
              height=size,
              width=size)
          
          plot(phmain[[4]])
          
          dev.off()
          
          pdf(sprintf("%s\\%s_%s_%s_marginal_col.pdf",
                      base_figure_path,
                      size_name,
                      datasets_names[path_idx],
                      resp_name),
              height=size,
              width=size)
          
          plot(gbar_col)
          
          dev.off()
          pdf(sprintf("%s\\%s_%s_%s_marginal_row.pdf",
                      base_figure_path,
                      size_name,
                      datasets_names[path_idx],
                      resp_name),
              height=size,
              width=size)
          
          plot(gbar_row)
          
          dev.off()
          
        }
        
      }
    }
    
  }
}

across_mice_decoding_all_trials_decode_axis_or_activity <- function()
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
  axis = F
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
  
  paths <- get_thirsty_quenched_paths()
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
      
      
      mode_1 <- get_mode(annotated_mice_1, paths[[idx_1]], axis = axis)
      mode_2 <- get_mode(annotated_mice_2, paths[[idx_2]], axis = axis)
      work_mode_1 <- mode_1$proj[which(apply(mode_1$proj, 1, function(r) {sum(is.na(r))})==0),]
      work_mode_2 <- mode_2$proj[which(apply(mode_2$proj, 1, function(r) {sum(is.na(r))})==0),]

      
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
          shuffled_stim_mat[,"Response"] <- shuffle(shuffled_stim_mat[,"Response"])
        }
        
        trials_mice_1 <- which(shuffled_stim_mat[,"Response"] %in% c(0,1,2,3,4,5))[which(apply(mode_1$proj, 1, function(r) {sum(is.na(r))})==0)]
        trials_labels_mice_1 <- shuffled_stim_mat[trials_mice_1,"Response"][which(apply(mode_1$proj, 1, function(r) {sum(is.na(r))})==0)]

        
        
        trials_mice_2 <- which(stim_master_2[,"Response"] %in% c(0,1,2,3,4,5))[which(apply(mode_2$proj, 1, function(r) {sum(is.na(r))})==0)]
        trials_labels_mice_2 <- stim_master_2[trials_mice_2,"Response"][which(apply(mode_2$proj, 1, function(r) {sum(is.na(r))})==0)]

  
        
        
        decoded_vec <- list()
        for (metric_name in names(Metrics)) {
          decoded_vec[[metric_name]] <- c()
        }
        
        for (decoded_trial_idx in 1:nrow(work_mode_2)) {
          decoded_trial <- work_mode_2[decoded_trial_idx,]
          
          confidence <-   
            apply(work_mode_1, 1, 
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
          
          
          true_vs_decoded_all[[metric_name]] <- 
              append(true_vs_decoded_all[[metric_name]],
                     list(list(train=idx_1, 
                          test=idx_2,
                          decoded=decoded_vec[[metric_name]],
                          gt=trials_labels_mice_2)))
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
  
  dir.create(sprintf("%s\\data\\supplemental\\", base_output_path))
  
  save(file=sprintf("%s\\data\\supplemental\\axes_thirst_thirst_new_confusion_matrices.Rda", base_output_path), confusion_matrix_list)
  save(file=sprintf("%s\\data\\supplemental\\axes_thirst_thirst_new_shuffle_confusion_matrices.Rda", base_output_path), shuffle_matrix_list)
  save(file=sprintf("%s\\data\\supplemental\\axes_thirst_thirst_new_pvalues_decoding.Rda", base_output_path), pvalues_all)
  save(file=sprintf("%s\\data\\supplemental\\activity_axes_thirst_thirst_new_decoded.Rda", base_output_path), mv)
  #lapply(mv, function(item) {res <- c(); for (s in sort(unique(item$gt))) {ind <- item$gt == s; res <- c(res, sum(item$decoded[ind] == item$gt[ind]) / sum(ind))}; return(res)})
}

plot_decoding_results <- function(prefix = "thirst_thirst", is_mode=F)
{
  if (is_mode) {
    prefix = sprintf("%s\\mode", prefix)
  }
  
  base_figure_path <- sprintf("%s\\general_decoding_results", output_path)
  dir.create(base_figure_path)
  base_figure_path <- sprintf("%s\\general_decoding_results\\%s", output_path, prefix)
  dir.create(base_figure_path)
  
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\confusion_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\shuffle_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\pvalues_all.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_confusion_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_shuffle_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_pvalues_all.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_confusion_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_shuffle_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_pvalues_all.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_confusion_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_shuffle_matrix_list.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_pvalues_all.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\%s_confusion_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\%s_shuffle_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\%s_pvalues_all.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_%s_confusion_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_%s_shuffle_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\f1_%s_pvalues_all.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_%s_confusion_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_%s_shuffle_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_%s_pvalues_all.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_%s_confusion_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_%s_shuffle_matrix_list.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\weird_all_%s_pvalues_all.Rda", base_output_path, prefix, ifelse(is_mode, "mode", "activity")), verbose=T)
  
  
  if (!is_mode) {
  decoder_scores_list <- 
    list(`F1`=
           list(`clusters`=f1_confusion_matrix_list,
                `clusters_shuffle`=f1_shuffle_matrix_list,
                `activity`=f1_activity_confusion_matrix_list,
                `activity_shuffle`=f1_activity_shuffle_matrix_list,
                `clusters_pvalue`=f1_pvalues_all,
                `activity_pvalue`=f1_activity_pvalues_all),
           
         `Weird`=
           list(`clusters`=weird_confusion_matrix_list,
                `clusters_shuffle`=weird_shuffle_matrix_list,
                `activity`=weird_activity_confusion_matrix_list,
                `activity_shuffle`=weird_activity_shuffle_matrix_list,
                `clusters_pvalue`=weird_pvalues_all,
                `activity_pvalue`=weird_activity_pvalues_all),
           
         `Weird_all`=
           list(`clusters`=weird_all_confusion_matrix_list,
                `clusters_shuffle`=weird_all_shuffle_matrix_list,
                `activity`=weird_all_activity_confusion_matrix_list,
                `activity_shuffle`=weird_all_activity_shuffle_matrix_list,
                `clusters_pvalue`=weird_all_pvalues_all,
                `activity_pvalue`=weird_all_activity_pvalues_all),
         
         `Accuracy`=
           list(`clusters`=confusion_matrix_list,
                `clusters_shuffle`=shuffle_matrix_list,
                `activity`=activity_confusion_matrix_list,
                `activity_shuffle`=activity_shuffle_matrix_list,
                `clusters_pvalue`=pvalues_all,
                `activity_pvalue`=activity_pvalues_all)
           )
  } else {
    decoder_scores_list <- 
      list(`F1`=
             list(`clusters`=f1_confusion_matrix_list,
                  `clusters_shuffle`=f1_shuffle_matrix_list,
                  `activity`=f1_mode_confusion_matrix_list,
                  `activity_shuffle`=f1_mode_shuffle_matrix_list,
                  `clusters_pvalue`=f1_pvalues_all,
                  `activity_pvalue`=f1_mode_pvalues_all),
           
           `Weird`=
             list(`clusters`=weird_confusion_matrix_list,
                  `clusters_shuffle`=weird_shuffle_matrix_list,
                  `activity`=weird_mode_confusion_matrix_list,
                  `activity_shuffle`=weird_mode_shuffle_matrix_list,
                  `clusters_pvalue`=weird_pvalues_all,
                  `activity_pvalue`=weird_mode_pvalues_all),
           
           `Weird_all`=
             list(`clusters`=weird_all_confusion_matrix_list,
                  `clusters_shuffle`=weird_all_shuffle_matrix_list,
                  `activity`=weird_all_mode_confusion_matrix_list,
                  `activity_shuffle`=weird_all_mode_shuffle_matrix_list,
                  `clusters_pvalue`=weird_all_pvalues_all,
                  `activity_pvalue`=weird_all_mode_pvalues_all),
           
           `Accuracy`=
             list(`clusters`=confusion_matrix_list,
                  `clusters_shuffle`=shuffle_matrix_list,
                  `activity`=mode_confusion_matrix_list,
                  `activity_shuffle`=mode_shuffle_matrix_list,
                  `clusters_pvalue`=pvalues_all,
                  `activity_pvalue`=mode_pvalues_all)
      )
  }
  
 
  processed <- 
  lapply(decoder_scores_list,
         
         function(decoder_scores) {
           
           all_decoders_results <- 
           lapply(c("clusters","clusters_shuffle","activity", "activity_shuffle"),
                  function(matrix_name) {
                      scores = decoder_scores[[matrix_name]]
                      
                      
                      all_trialtypes_scored <- 
                      lapply(names(scores),
                             function(trial_name) {
                               
                               trial_matrix <- scores[[trial_name]]
                               melted <- melt(trial_matrix)
                               melted <- melted[melted$Var1 != melted$Var2,]
                               averaged <- ddply(melted, .(Var1), function(mdf) {mean(mdf[,"value"])})
                               colnames(averaged) <- c("DatasetIdx", "Score")
                               averaged$TrialType = trial_name
                               
                               return(averaged)
                             })
                      
                      all_trialtypes_scored <- do.call(rbind, all_trialtypes_scored)
                      all_trialtypes_scored$Class <- matrix_name
                    
                      return(all_trialtypes_scored)  
                    }
                  )
           
           all_decoders_results <- do.call(rbind, all_decoders_results)
           
           activity_melted <- melt(decoder_scores$activity_pvalue)
           clusters_melted <- melt(decoder_scores$clusters_pvalue)
           colnames(activity_melted) <- c("#", "TrialType", "Pvalue", "##")
           colnames(clusters_melted) <- c("#", "TrialType", "Pvalue", "##")
           
           activity_melted$Class = "activity"
           clusters_melted$Class = "clusters"
           pvalues_df_all <- rbind(activity_melted, clusters_melted)
           
           final_res = list(decoder=all_decoders_results,
                            pvalues=pvalues_df_all)
           return(final_res)
         })
  
  
  
  metrics_to_use <- c("cosine")  
  
  for (score_type in names(processed)) { 
    
    st_base_figure_path <- sprintf("%s\\%s\\", base_figure_path, score_type)
    dir.create(st_base_figure_path)
    statistics_df <- data.frame()
    all_trial_types <- unique(processed[[score_type]]$decoder$TrialType)
    
    for (tt in all_trial_types) {
      
      df_final <- processed[[score_type]]$decoder[processed[[score_type]]$decoder$TrialType ==  tt,]
      df_pval_final <- processed[[score_type]]$pvalues[processed[[score_type]]$pvalues$TrialType ==  str_split(tt, sprintf("_%s", metrics_to_use))[[1]][1],]
      g <- 
        ggplot(df_final, 
               aes(x=factor(Class, levels=c("clusters", "activity", "clusters_shuffle",  "activity_shuffle")), 
                   y=Score), 
               group=Class)+ 
        #geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.25) + 
        scale_shape_manual(values=18) +
        geom_boxplot(fill=NA, outlier.shape=NA)  + 
        theme_light() + 
        #big_text_base_plot_theme + 
        ylab("Decoding accuracy (%)") + 
        #ylim(0,1) +
        xlab("")  +
        ggtitle(sprintf("%s_%s", score_type, tt)) +
        theme(plot.title = element_text(size=9),
              text=element_text(size=10.5, colour="black"),
              axis.text=element_text(size=10.5, colour="black"))
        
          
      
      
      gbox <- 
        ggplot(df_final, 
               aes(x=factor(Class, levels=c("clusters", "activity", "clusters_shuffle",  "activity_shuffle")), 
                   y=Score), 
               group=Class)+
        geom_boxplot(fill=NA, width=.5, outlier.shape=NA)
      
      
      for (d_idx in unique(df_final$DatasetIdx)) {
        
        row_df <- data.frame(x=c(1,2), 
                             y=c(df_final[df_final$DatasetIdx == d_idx & df_final$Class == "clusters","Score"],
                                 df_final[df_final$DatasetIdx == d_idx & df_final$Class == "activity","Score"]))
        
        gbox <- gbox + geom_line(data=row_df, aes(x=x, y=y))
      }
      
      gboxf <- 
        gbox +
        geom_point() + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Decoding accuracy (%)") + 
        ylim(0,1) +
        xlab("")  +
        theme(plot.title = element_text(size=9),
              text=element_text(size=10.5, colour="black"),
              axis.text=element_text(size=10.5, colour="black"))
      
      
      gpval <- 
        ggplot(df_pval_final, 
               aes(x=factor(Class, levels=c("clusters", "activity")), 
                   y=Pvalue, 
                   fill=factor(Class, levels=c("clusters", "activity"))))+ 
        #geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.25) +
        geom_violin(color=NA, fill="gray60") + 
        geom_jitter(color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
        #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
        stat_summary(size=1) + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Pvalue") + 
        xlab("") 
        ggtitle(sprintf("%s_%s", score_type, tt)) +
        theme(plot.title = element_text(size=9),
              text=element_text(size=10.5, colour="black"),
              axis.text=element_text(size=10.5, colour="black"))
      
      for (size_name in names(a4_sizes)) {
        

        size = a4_sizes[[size_name]]
        pdf(sprintf("%s\\%s_%s.pdf",
                    st_base_figure_path, 
                    size_name, 
                    tt), height=size, width=size)
        plot(g)
        dev.off()
        
        pdf(sprintf("%s\\%s_%s_connected.pdf",
                    st_base_figure_path, 
                    size_name, 
                    tt), height=size, width=size)
        plot(gboxf)
        dev.off()
        
        pdf(sprintf("%s\\%s_%s_pvalues.pdf",
                    st_base_figure_path, 
                    size_name, 
                    tt), height=size, width=size)
        plot(gpval)
        dev.off()
        
      }

      ##### Statistics
      clusters <- df_final[df_final$Class == "clusters",]
      activity <- df_final[df_final$Class == "activity",]
      clusters_shuffle <- df_final[df_final$Class == "clusters_shuffle",]
      activity_shuffle <- df_final[df_final$Class == "activity_shuffle",]
      
      
      
      clusters_order <-  order(clusters[,"DatasetIdx"])
      activity_order <-  order(activity[,"DatasetIdx"])
      clusters_shuffle_order <-  order(clusters_shuffle[,"DatasetIdx"])
      activity_shuffle_order <-  order(activity_shuffle[,"DatasetIdx"])
      
      
      comp_configurations <- list(list(name = "clusters_activity", 
                                       group_A = clusters[clusters_order,"Score"],
                                       group_B = activity[activity_order ,"Score"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "clusters_shuffle", 
                                       group_A = clusters[clusters_order,"Score"],
                                       group_B = clusters_shuffle[clusters_shuffle_order,"Score"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "activity_shuffle", 
                                       group_A = activity[activity_order,"Score"],
                                       group_B = activity_shuffle[activity_shuffle_order,"Score"],
                                       paired= T,
                                       alternative="greater"))
      
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
                              metric_name="Cosine",
                              score_type=score_type,
                              TrialType=tt)
        
        statistics_df <- rbind(statistics_df,
                               stat_df)
      }
    }
    
    dir.create(sprintf("%s//statistics",st_base_figure_path))
    write.csv(file=sprintf("%s//statistics//statistics.csv", st_base_figure_path),
              statistics_df)
    
  
  }
  
  

}

plot_manipulation_decoding_results <- function(prefix = "thirst_thirst", is_mode=F, metadata, chunks)
{
  if (is_mode) {
    prefix = sprintf("%s\\mode", prefix)
  }
  
  base_figure_path <- sprintf("%s\\general_decoding_results", output_path)
  dir.create(base_figure_path)
  base_figure_path <- sprintf("%s\\general_decoding_results\\%s", output_path, prefix)
  dir.create(base_figure_path)
  
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\true_vs_decoded_all.Rda", base_output_path, prefix), verbose=T)
  load(file=sprintf("%s\\data\\figure_3\\new_decoder\\%s\\activity_true_vs_decoded_all.Rda", base_output_path, prefix), verbose=T)

  
  
  Responses = c("0"="Hit", "1"="Miss", "2"="NeutCR", "3"="NeutFA", "4"="CorrectRejection", "5"="FalseAlarm")
  
  pre <- list()
  post <- list()
  
  
  for (decoded_pair_idx in 1:len(true_vs_decoded_all$cosine)) {
    
    decoded_pair <- true_vs_decoded_all$cosine[[decoded_pair_idx]]
    activity_decoded_pair <- activity_true_vs_decoded_all$cosine[[decoded_pair_idx]]
    
    relevant_metadata <- metadata[[decoded_pair$test]]
    relevant_stim_mat <- relevant_metadata$stim_master_mat[relevant_metadata$stim_master_mat[,"TrialType"] %in% c(3:5),]
    
    assert(nrow(relevant_stim_mat) == len(decoded_pair$decoded))
    assert(decoded_pair$train == activity_decoded_pair$train)
    assert(decoded_pair$test == activity_decoded_pair$test)
    
    pre_indices <- !(relevant_stim_mat[,ncol(relevant_stim_mat)] %in% chunks[[decoded_pair$test]])
    post_indices <- relevant_stim_mat[,ncol(relevant_stim_mat)] %in% chunks[[decoded_pair$test]]
    
    pre_decoded <- decoded_pair$decoded[pre_indices]
    pre_gt <- decoded_pair$gt[pre_indices]
    post_decoded <- decoded_pair$gt[post_indices]
    post_gt <- decoded_pair$decoded[post_indices]
    
    activity_pre_decoded <- activity_decoded_pair$decoded[pre_indices]
    activity_pre_gt <- activity_decoded_pair$gt[pre_indices]
    activity_post_decoded <- activity_decoded_pair$gt[post_indices]
    activity_post_gt <- activity_decoded_pair$decoded[post_indices]
    
    
    pre_f1_scores <- get_f1_scores(pre_gt, pre_decoded)
    post_f1_scores <- get_f1_scores(post_gt, post_decoded)
    activity_pre_f1_scores <- get_f1_scores(activity_pre_gt, activity_pre_decoded)
    activity_post_f1_scores <- get_f1_scores(activity_post_gt, activity_post_decoded)
    
    
    for (score_name in names(pre_f1_scores)) {
      post[[score_name]] <- rbind(post[[score_name]], c(decoded_pair$train, 
                                                        decoded_pair$test, 
                                                        post_f1_scores[[score_name]],
                                                        activity_post_f1_scores[[score_name]]))
      
      pre[[score_name]] <- rbind(pre[[score_name]], c(decoded_pair$train, 
                                                      decoded_pair$test, 
                                                      pre_f1_scores[[score_name]],
                                                      activity_pre_f1_scores[[score_name]]))
    }
  }
  
  
  processed <- 
    lapply(names(pre),
           
           function(score_name) {
             
             pre_score_df <- as.data.frame(pre[[score_name]])
             post_score_df <- as.data.frame(post[[score_name]])
             
             colnames(pre_score_df)[1] <- "train"
             colnames(post_score_df)[1] <- "train"
             
             post_score_df_avg <- 
               ddply(post_score_df,
                     .(train),
                     function(train_df) {
                       colMeans(train_df[,-2])
                     })
             
             
             pre_score_df_avg <- 
               ddply(pre_score_df,
                     .(train),
                     function(train_df) {
                       colMeans(train_df[,-2])
                     })
             
             
             pre_clusters <- pre_score_df_avg[,c(1,2:7)]
             pre_activity <- pre_score_df_avg[,c(1,8:13)]
             post_clusters <- post_score_df_avg[,c(1,2:7)]
             post_activity <- post_score_df_avg[,c(1,8:13)]
             
             colnames(pre_clusters) <- c("DatasetIdx", Responses)
             colnames(pre_activity) <- c("DatasetIdx", Responses)
             colnames(post_clusters) <- c("DatasetIdx", Responses)
             colnames(post_activity) <- c("DatasetIdx", Responses)
             
             melted_pre_clusters <- melt(pre_clusters, id.vars="DatasetIdx")
             melted_pre_activity <- melt(pre_activity, id.vars="DatasetIdx")
             melted_post_clusters <- melt(post_clusters, id.vars="DatasetIdx")
             melted_post_activity <- melt(post_activity, id.vars="DatasetIdx")
             
             colnames(melted_pre_clusters) <- c("DatasetIdx", "TrialType", "Score")
             colnames(melted_pre_activity) <- c("DatasetIdx", "TrialType", "Score")
             colnames(melted_post_clusters) <- c("DatasetIdx", "TrialType", "Score")
             colnames(melted_post_activity) <- c("DatasetIdx", "TrialType", "Score")
             
             melted_pre_clusters$Class <- "clusters"
             melted_pre_activity$Class <- "activity"
             melted_post_clusters$Class <- "clusters"
             melted_post_activity$Class <- "activity"
             
             melted_pre_clusters$Manipulation <- "pre"
             melted_pre_activity$Manipulation <- "pre"
             melted_post_clusters$Manipulation <- "post"
             melted_post_activity$Manipulation <- "post"
             
             score_df_final <- rbind(melted_pre_clusters,
                                     melted_pre_activity,
                                     melted_post_clusters,
                                     melted_post_activity)
             
             
             return(score_df_final)
             
             
             
           }
          )
  
  names(processed) <- names(pre)
  
  
  metrics_to_use <- c("cosine")  
  
  for (score_type in names(processed)) { 
    
    st_base_figure_path <- sprintf("%s\\%s\\", base_figure_path, score_type)
    dir.create(st_base_figure_path)
    statistics_df <- data.frame()
    all_trial_types <- unique(processed[[score_type]]$TrialType)
    
    for (tt in all_trial_types) {
      
      df_final <- processed[[score_type]][processed[[score_type]]$TrialType ==  tt,]
      class_levels=c("clusters.pre", "clusters.post", "activity.pre", "activity.post")
      g <- 
        ggplot(df_final, 
               aes(x=factor(interaction(Class, Manipulation), levels=class_levels), 
                   y=Score), 
               group=interaction(Class, Manipulation))+ 
        #geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=.25) + 
        scale_shape_manual(values=18) +
        geom_boxplot(fill=NA, outlier.shape=NA)  + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Decoding accuracy (%)") + 
        #ylim(0,1) +
        xlab("")  +
        ggtitle(sprintf("%s_%s", score_type, tt)) +
        theme(plot.title = element_text(size=9),
              text=element_text(size=10.5, colour="black"),
              axis.text=element_text(size=10.5, colour="black"))
      
      
      
      
      gbox <- 
        ggplot(df_final, 
               aes(x=factor(interaction(Class, Manipulation), levels=class_levels), 
                   y=Score), 
               group=interaction(Class, Manipulation))+ 
        geom_boxplot(fill=NA, width=.5, outlier.shape=NA)
      
      
      for (d_idx in unique(df_final$DatasetIdx)) {
        
        row_df <- data.frame(x=c(1,2,3,4), 
                             y=c(df_final[df_final$DatasetIdx == d_idx & df_final$Class == "clusters" & df_final$Manipulation == "pre","Score"],
                                 df_final[df_final$DatasetIdx == d_idx & df_final$Class == "clusters"  & df_final$Manipulation == "post","Score"],
                                 df_final[df_final$DatasetIdx == d_idx & df_final$Class == "activity" & df_final$Manipulation == "pre","Score"],
                                 df_final[df_final$DatasetIdx == d_idx & df_final$Class == "activity"  & df_final$Manipulation == "post","Score"]))
        
        gbox <- gbox + geom_line(data=row_df, aes(x=x, y=y))
        #gbox <- gbox + geom_line(data=activity_row_df, aes(x=x, y=y))
      }
      
      gboxf <- 
        gbox +
        geom_point() + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Decoding accuracy (%)") + 
        ylim(0,1) +
        xlab("")  +
        theme(plot.title = element_text(size=9),
              text=element_text(size=10.5, colour="black"),
              axis.text=element_text(size=10.5, colour="black"))

      
      for (size_name in names(a4_sizes)) {
        
        
        size = a4_sizes[[size_name]]
        pdf(sprintf("%s\\%s_%s.pdf",
                    st_base_figure_path, 
                    size_name, 
                    tt), height=size, width=size)
        plot(g)
        dev.off()
        
        pdf(sprintf("%s\\%s_%s_connected.pdf",
                    st_base_figure_path, 
                    size_name, 
                    tt), height=size, width=size)
        plot(gboxf)
        dev.off()
        

        
      }
      
      ##### Statistics
      pre_clusters <- df_final[df_final$Class == "clusters" & df_final$Manipulation == "pre",]
      post_clusters <- df_final[df_final$Class == "clusters" & df_final$Manipulation == "post",]
      pre_activity <- df_final[df_final$Class == "activity" & df_final$Manipulation == "pre",]
      post_activity <- df_final[df_final$Class == "activity" & df_final$Manipulation == "post",]
      
      stat_all <- 
      rbind(pre_clusters,
            post_clusters,
            pre_activity,
            post_activity)
      
      pre_clusters_order <- order(pre_clusters$DatasetIdx)
      post_clusters_order <- order(post_clusters$DatasetIdx)
      pre_activity_order <- order(pre_activity$DatasetIdx)
      post_activity_order <- order(post_activity$DatasetIdx)
      

      
      dtest <- dunn.test::dunn.test(stat_all$Score, paste(stat_all$Class, stat_all$Manipulation, sep="_"), method = "bonferroni")
      
      comp_configurations <- list(list(name = "pre_post_clusters",
                                       group_A = pre_clusters[pre_clusters_order,"Score"],
                                       group_B = post_clusters[post_clusters_order ,"Score"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "pre_post_activity",
                                       group_A = pre_activity[pre_activity_order,"Score"],
                                       group_B = post_activity[post_activity_order,"Score"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "pre_clusters_activity",
                                       group_A = pre_clusters[pre_clusters_order,"Score"],
                                       group_B = pre_activity[pre_activity_order,"Score"],
                                       paired= T,
                                       alternative="greater"),
                                  list(name = "post_clusters_activity",
                                       group_A = post_clusters[post_clusters_order,"Score"],
                                       group_B = post_activity[post_activity_order,"Score"],
                                       paired= T,
                                       alternative="greater"))
      
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
                               metric_name="Cosine",
                               score_type=score_type,
                               TrialType=tt)
      
        
         statistics_df <- rbind(statistics_df, stat_df)
         
       }
      
      # stat_df <- as.data.frame(dtest)
      # stat_df$TrialType <- tt
      # stat_df$score_type <- score_type
      # stat_df$signif_corrected <- signif.num(stat_df$P.adjusted)
      # stat_df$signif <- signif.num(stat_df$P)
      
      
    }
    
    dir.create(sprintf("%s//statistics",st_base_figure_path))
    write.csv(file=sprintf("%s//statistics//wilc_statistics.csv", st_base_figure_path),
              statistics_df)
    
    
  }
}

plot_satiation_free_vs_dynamics_examples_controled_lick_rate <- function()
{
  base_figure_path <- sprintf("%s\\satiation_free_vs_dynamics_examples_controled_lick_rate", output_path)
  dir.create(base_figure_path)
  
  
  pooled_task_similarity <- list()
  pooled_task_freely_similarity <- list()
  satiation_run_indices_all = c(5, 4, 3, 3, 3)
  satiation_run_indices_interval <- c(30, 16, 30, 16,30)
  satiation_paths <- get_satiation_paths()
  satiation_paths <- get_satiation_paths()
  metadata_satiation <- across_mice_decoding_build_metadata(satiation_paths)
  datasets_names = paste("SATIATION_", get_datasets_names(satiation_paths, sep="_"), sep="")
  
  for (i in 1:len(metadata_satiation)) {
    res = get_task_vs_freely_lick_distribution(metadata_satiation[[i]], 
                                               satiation_paths[i], 
                                               satiation_run_indices_all[i], 
                                               bouts_interval = satiation_run_indices_interval[i])
    

    task_lick_rate <- apply(res$task_licking_mat[,62:122], 1, mean)
    free_lick_rate <- apply(res$freely_licking_mat[,62:122], 1, mean)
    
    task_trials_to_use <- which(task_lick_rate <= (mean(task_lick_rate) + sd(task_lick_rate)) & task_lick_rate >= (mean(task_lick_rate) - sd(task_lick_rate)))
    free_trials_to_use <- which(free_lick_rate <= (mean(task_lick_rate) + sd(task_lick_rate)) & free_lick_rate >= (mean(task_lick_rate) - sd(task_lick_rate)))
    
    if (len(free_trials_to_use) == 0)    {
      task_trials_to_use <- 1:len(task_lick_rate)
      free_trials_to_use <- which(free_lick_rate <= max(task_lick_rate) & free_lick_rate >= min(task_lick_rate))
    }
    
    cmt <- rbind(res$task_cluster_mat[task_trials_to_use,], res$freely_cluster_mat[free_trials_to_use,])
    lmt <- rbind(res$task_licking_mat[task_trials_to_use,], res$freely_licking_mat[free_trials_to_use,])
    structured_unstructured <- data.frame(type=rep(c("Structured", "Unstructed"), times=c(len(task_trials_to_use),
                                                                                          len(free_trials_to_use))))
    

    
    
    
    
    
    
    rownames(structured_unstructured) <- 1:nrow(cmt)
    rownames(cmt) <- 1:nrow(cmt)
    colnames(cmt) <- rep("", times=ncol(cmt))
    colnames(lmt) <- rep("", times=ncol(lmt))
    colnames(cmt)[4] <- "O"
    colnames(lmt)[61] <- "O"
    
    
    cluster_color_label <- spec_cg(len(unique(metadata_satiation[[i]]$cluster_mat$labs)))
    names(cluster_color_label) <- c(-1, 1:(len(unique(metadata_satiation[[i]]$cluster_mat$labs)) - 1))
    
    phlmt <- pheatmap(lmt, cluster_rows=F, cluster_cols=F, legend=F, col=c(`0`="white",`1`="black"), border_col=NA)
    phcmt <- pheatmap(cmt, cluster_rows=F, cluster_cols=F, legend=F, border_co=NA, annotation_row = structured_unstructured, annotation_names_row = F, annotation_legend = F, show_rownames = F,
                      col=cluster_color_label)
    
    
  
    
    plt <- plot_grid(phcmt[[4]], phlmt[[4]])
    
    
    
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
    }
    
    
    sim_mat <- cosine(t(cmt))
    
    
    task_sim_mat <- sim_mat[1:len(task_trials_to_use),
                            1:len(task_trials_to_use)]
    
    free_sim_mat <- sim_mat[(len(task_trials_to_use) + 1):nrow(cmt), 1:len(task_trials_to_use)]
    
    pooled_task_similarity <- append(pooled_task_similarity,
                                     list(task_sim_mat[lower.tri(task_sim_mat)]))
    
    pooled_task_freely_similarity <- append(pooled_task_freely_similarity,
                                            list(c(free_sim_mat)))
  }
  
  
  pooled_smilarity_list <- list(task=unlist(pooled_task_similarity),
                                freely=unlist(pooled_task_freely_similarity))
  
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
    pdf(sprintf("%s\\dist_and_box\\%s_trial_similarity.pdf",base_figure_path, size_name), height=size, width=size)
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

STUB <- function() {
  
  base_figure_path <- sprintf("%s\\STUB", output_path)
  dir.create(base_figure_path)
  
  gSTUB <- ggplot()
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_STUB.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gSTUB)
    
    dev.off()
    
  }
}