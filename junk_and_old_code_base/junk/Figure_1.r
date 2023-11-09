library(latex2exp)

figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"

figure_1_illustration_traces <- function() {
  write_path <- sprintf("%s\\figure_1\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\illustration_traces", write_path)
  dir.create(write_path)  
  
  paths <- get_thirsty_quenched_paths()
  original_mat <- get_reduced_mat_full_day(paths[[1]], just_original_mat = T, window_size=15)
  
  sizes = list(big=c(width=3,
                     height=4),
               medium=c(width=2.5,
                        height=3),
               small=c(width=2,
                       height=2.75))
  
  indices <- 500:3000
  n_neurons <- 6
  
  for (trace_example_idx in 1:8) {

  max_y = (n_neurons -  1) * plot_jump_factor + 1
  
  traces_df <- c()
  random_neurons <- sample(1:nrow(original_mat), n_neurons)
  
  
  plot_jump_factor = 1.5
  for (neur_idx in 1:len(random_neurons)) {
    neur <- random_neurons[neur_idx]
    trace <- original_mat[neur,indices]
    trace <- (trace - min(trace)) / (max(trace) - min(trace))  + (neur_idx - 1) * plot_jump_factor 
  
    traces_df  <- cbind(traces_df,
                       trace)
  }
  
  
  poly_x1 <- c(200,200,400,400)
  poly_x2 <- c(1500,1500,1700,1700)
  gtrace <- ggplot() +
            geom_polygon(data=data.frame(x=poly_x1,
                                         y=c(0,max_y,max_y,0)),
                         aes(x=x,y=y),
                         color="NA",
                         fill=adjustcolor("gray60", alpha=.4)) +
  
            geom_polygon(data=data.frame(x=poly_x2,
                                         y=c(0,max_y,max_y,0)),
                         aes(x=x,y=y),
                         color="NA",
                         fill=adjustcolor("gray60", alpha=.4)) + 
            annotate(geom="text", label=TeX("$t_1$"),
                     x=mean(poly_x1),
                     y=-.35) + 
            annotate(geom="text", label=TeX("$t_n$"),
                     x=mean(poly_x2),
                     y=-.35)
    
  for (trace_idx in 1:ncol(traces_df)) {
    
    trace_tmp_df <- data.frame(dff=traces_df[,trace_idx],
                               time=1:nrow(traces_df))
    
    
    if (trace_idx == 1) {
      neur_label = "Neuron K"
    } else {
      neur_label=sprintf("Neuron %d", rev(1:ncol(traces_df))[trace_idx])
    }
    
    gtrace <- gtrace + 
              geom_line(data=trace_tmp_df, aes(x=time, y=dff),
                        size=.5,
                        color="#2281C4") +
              annotate(geom="text", 
                       label=neur_label,
                       x=-120,
                       y=min(trace_tmp_df$dff) + 0.25)
              
  }
  

  
  gtrace <- gtrace + 
            theme(line = element_blank(),
                  axis.text = element_blank(),
                  axis.title.y = element_blank(),
                  panel.background = element_blank(),
                  text=element_text(size=15)) + 
            xlim(-140, 2502) + 
            xlab("Time")
  
  for (size_name in names(sizes)) {
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    pdf(sprintf("%s\\%s\\trace_example_%d.pdf",
                write_path,
                size_name,
                trace_example_idx),
        height=sizes[[size_name]]["height"],
        width=sizes[[size_name]]["width"])
    
    plot(gtrace)
    dev.off()
    
    pdf(sprintf("%s\\%s\\transposed_trace_example_%d.pdf",
                write_path,
                size_name,
                trace_example_idx),
        height=sizes[[size_name]]["width"],
        width=sizes[[size_name]]["height"])
    
    plot(gtrace)
    dev.off()
  }
  }
}


figure_1_illustration_reduced_structure <- function() {
  
  write_path <- sprintf("%s\\figure_1\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\illustration_structure", write_path)
  dir.create(write_path)  
  
  paths <- get_thirsty_quenched_paths()
  reduced_mat <- get_mat_with_preset(paths[[1]], "dalshtaim")
  clust_mat_all <- get_clusters_mat_kmeans(reduced_mat, 
                                           20)
  
  p1 <- pair_plots_colored_clusters(reduced_mat, 
                                    cluster_labels = clust_mat_all$labs,
                                    to_label = 1,
                                    specific_dims = c(1,3))
  
  
  sizes = list(big=c(width=2,
                     height=2),
               medium=c(width=1.25,
                        height=1.25),
               medium_1=c(width=1,
                          height=1),
               small=c(width=.75,
                       height=.75))
  
  p1 <- pair_plots_colored_clusters(reduced_mat, 
                                    cluster_labels = clust_mat_all$labs,
                                    to_label = 1,
                                    specific_dims = c(1,3), 
                                    remove_outliers=T, 
                                    plot_labels = F,
                                    stroke_alpha=.25,
                                    return_grid = F)
  
  
  p2 <- pair_plots_colored_clusters(reduced_mat, 
                                    cluster_labels = clust_mat_all$labs,
                                    to_label = 1,
                                    specific_dims = c(4,6), 
                                    remove_outliers=T, 
                                    plot_labels = F,
                                    stroke_alpha=.25,
                                    return_grid = F)
  
  
  p3 <- pair_plots_colored_clusters(reduced_mat, 
                                    cluster_labels = clust_mat_all$labs,
                                    to_label = 1,
                                    specific_dims = c(1,2), 
                                    remove_outliers=T, 
                                    plot_labels = F,
                                    stroke_alpha=.25,
                                    return_grid = F)
  
  pf <- lapply(append(append(list(p1), list(p2)), list(p3)), function(p) {p[[1]] + ylab("") + xlab("") + theme(line=element_blank(),
                                                                                                               rect=element_blank(),
                                                                                                               plot.background = element_rect("white"))})
  pf$nrow <- 1
  gf <- do.call(arrangeGrob, pf)
  plot(gf)
  
  
  for (size_name in names(sizes)) {
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    pdf(sprintf("%s\\%s\\example_clusters.pdf",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]],
        width=sizes[[size_name]][["width"]]*3)
    
    plot(gf)
    dev.off()
    
    png(sprintf("%s\\%s\\example_clusters.png",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]]*2,
        width=sizes[[size_name]][["width"]]*6,
        unit="in",
        res=1500)
    
    plot(gf)
    dev.off()
    
  }
}


figure_1_example_dimensionality_estimate <- function() {
  
  write_path <- sprintf("%s\\figure_1\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\illustration_dimension", write_path)
  dir.create(write_path)  
  
  sizes = list(big=c(width=3,
                     height=3),
               medium=c(width=2.75,
                        height=2.75),
               medium_2=c(width=2.5,
                          height=2.5),
               medium_1=c(width=2,
                          height=2),
               small=c(width=1.5,
                       height=1.5))
  
  paths <- get_thirsty_quenched_paths()
  reduced_mat <- get_mat_with_preset(paths[[1]], "dalshtaim")
  
  
  dist_mat <- as.matrix(dist(reduced_mat))
  mus <- apply(dist_mat, 
               1, 
               function(r) {
                    mus_temp <- r[order(r)[2:3]]
              return(mus_temp[2] / mus_temp[1]) })
  
  mus <- sort(mus)
  
  linreg_x <- log(mus[1:(floor(len(mus)) * 0.96)])
  linreg_y <- -log(1- (0:(len(linreg_x) -1) / len(linreg_x)))
  linreg <- lm(linreg_y ~ linreg_x - 1)
  
  slope = linreg$coefficients[[1]]
  
  ind <- sort(sample(1:len(linreg_x), 1500))
  
  plot_df <- data.frame(x=linreg_x[ind],
                        y=linreg_y[ind])
  
  x_to_use <- seq(-5, max(linreg_x), by=.3)
  
  


  
  gdim_all <- ggplot() +
              theme_classic() +
              scale_x_continuous(expand=c(0,0)) +
              scale_y_continuous(expand=c(0,0))
  
  plot_slope = 3
  
  for(intercept_x in x_to_use) {
    
    b = -plot_slope * intercept_x 
    xt <- ifelse(intercept_x > 0, intercept_x, 0)
    xt <- c(xt, max(linreg_x))
    yt <- plot_slope * xt + b
    
    if (yt[1] <= max(linreg_y)) {
      if (yt[2] >= max(linreg_y)){
        
        yt[2] <- max(linreg_y)
        xt[2] <- (max(linreg_y) - b) / plot_slope
      }
    slope_lines_df <- data.frame(x=xt,
                                 y=yt)
    
    gdim_all <- gdim_all +
                geom_line(data=slope_lines_df, aes(x=x,y=y), linetype="dashed", col="gray80")
    
    }
  }
            
  
  
  true_xt <- c(0, max(linreg_x))
  true_yt <- c(0, max(linreg_x) * slope)
  
  true_slope_df <- data.frame(x=true_xt,
                              y=true_y)
  
  gdim_f <- 
    gdim_all +
    geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
    geom_point(data=plot_df, aes(x=x, y=y), alpha=.5, stroke=0) +
    xlab(TeX("$log(NN_{distance})$")) + 
    ylab(TeX("$-log(1-NN_{quantile})$")) + 
    ggtitle("Dimensionality estima1tion") +
    theme(text=element_text(size=14),
          plot.title=element_text(size=10))
  
  for (size_name in names(sizes)) {
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    pdf(sprintf("%s\\%s\\dim_estimation.pdf",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]],
        width=sizes[[size_name]][["width"]])
    
    plot(gdim_f)
    dev.off()
    
    
  }
  
  
}


figure_1_example_dimensionality_dataset <- function() {
  
  write_path <- sprintf("%s\\figure_1\\", figures_base_path)
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
  reduced_mat <- get_mat_with_preset(paths[[1]], "dalshtaim")  
  
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
    png(sprintf("%s\\%s\\lem2_example_3_rows.png",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]] * 3,
        width=sizes[[size_name]][["width"]] * 5,
        unit="in",
        res=1500)

    plot(p_lem2_3)
    dev.off()
    
    png(sprintf("%s\\%s\\lem2_example_5_rows.png",
                write_path,
                size_name),
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
    
    png(sprintf("%s\\%s\\lem2_example_3_rows_no_annot.png",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]] * 3,
        width=sizes[[size_name]][["width"]] * 5,
        unit="in",
        res=1500)
    
    plot(pf_3)
    dev.off()
    
    
    png(sprintf("%s\\%s\\lem2_example_5_rows_no_annot.png",
                write_path,
                size_name),
        height=sizes[[size_name]][["height"]] * 5,
        width=sizes[[size_name]][["width"]] * 3,
        unit="in",
        res=1500)
    
    plot(pf_5)
    dev.off()
    
  }
  
  original_mat <- get_reduced_mat_full_day(paths[[1]], window_size = 15, just_original_mat = T)
  
  for (knn in c(75,85,95,105,110,125,115,120)) {
    
    hlle_reduced <- dimRed::embed(t(original_mat), "HLLE", ndim=6, knn=knn)
    
    p_hlle_3 <- pair_plots_colored(hlle_reduced@data@data, 
                                   rep(adjustcolor("gray30", alpha=.3), times=nrow(hlle_reduced@data@data)),
                                   plot_nrow = 3)
    
    p_hlle_5 <- pair_plots_colored(hlle_reduced@data@data, 
                                   rep(adjustcolor("gray30", alpha=.3), times=nrow(hlle_reduced@data@data)),
                                   plot_nrow = 5)
    
    
    p_hlle_list <- pair_plots_colored(hlle_reduced@data@data, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(hlle_reduced@data@data)),
                                      return_grid = F,
                                      pt_size=.55)  
    
    p_hlle_list <- lapply(1:len(p_hlle_list), 
                          function(i) {return(p_hlle_list[[i]] + 
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
      png(sprintf("%s\\%s\\hlle_example_%d_3_rows.png",
                  write_path,
                  size_name,
                  knn),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      plot(p_hlle_3)
      dev.off()
      
      png(sprintf("%s\\%s\\hlle_example_%d_5_rows.png",
                  write_path,
                  size_name,
                  knn),
          height=sizes[[size_name]][["height"]] * 5,
          width=sizes[[size_name]][["width"]] * 3,
          unit="in",
          res=1500)
      
      plot(p_hlle_5)
      dev.off()
      
      p_hlle_5_no_annot <- p_hlle_list
      p_hlle_3_no_annot <- p_hlle_list
      
      p_hlle_5_no_annot$nrow <- 5
      pf_5 <- do.call(arrangeGrob, p_hlle_5_no_annot)
      
      p_hlle_3_no_annot$nrow <- 3
      pf_3 <- do.call(arrangeGrob, p_hlle_3_no_annot)
      
      png(sprintf("%s\\%s\\hlle_example_%d_3_rows_no_annot.png",
                  write_path,
                  size_name,
                  knn),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      plot(pf_3)
      dev.off()
      
      
      png(sprintf("%s\\%s\\hlle_example_%d_5_rows_no_annot.png",
                  write_path,
                  size_name,
                  knn),
          height=sizes[[size_name]][["height"]] * 5,
          width=sizes[[size_name]][["width"]] * 3,
          unit="in",
          res=1500)
      
      plot(pf_5)
      dev.off()
      
    }
  }
}


figure_1_cool_high_dim_plots <- function() {
  write_path <- sprintf("%s\\figure_1\\", figures_base_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\cool_high_dimension_plots", write_path)
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
  
  for (path_idx in 1:len(paths)) {
    
    p <- paths[[path_idx]]
    #original_mat <- get_reduced_mat_full_day(p, "da", window_size=15, just_original_mat = T)
    #reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, time_shuffled=F)
    red_mat <- get_mat_with_preset(p, "dalshtaim")
    
    for (size_name in names(sizes)) {
      
      dir.create(sprintf("%s\\%s",
                         write_path,
                         size_name))

      png(sprintf("%s\\%s\\no_shadow_cool_plots_%s.png",
                  write_path,
                  size_name,
                  datasets_names[[path_idx]]),
          height=sizes[[size_name]][["height"]] * 4,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      cool_high_dimensional_plots(red_mat)
      dev.off()
      
    }
    
  }
  
}

cool_high_dimensional_plots <- function(a, cp_from_home=NA, color_list_from_home=list()) {
  
  
  
  #combinations_mat <- combn(sample(1:ncol(a), ncol(a)), 3)
  combinations_mat <- combn(1:ncol(a), 3)
  par(mfrow=c(4,5))
  par(mar = c(0,0,0,0))
  par(oma = c(0,0,0,0))
  par(xpd=T)
  
  
  for (j in 1:20) {
    
    xdim <- combinations_mat[1,j]  
    ydim <- combinations_mat[2,j]  
    zdim <- combinations_mat[3,j]
    
    
    
    
    xlab_t = sprintf("Dim %d", xdim)
    ylab_t = sprintf("Dim %d", ydim)
    zlab_t = sprintf("Dim %d", zdim)
    mz = a[,zdim]
    mz_f <- mz + 3 * sd(mz)
    mx <- a[,xdim]
    my <- a[,ydim]
    my_f <- c(my, my)
    mx_f <- c(mx, mx)
    mz_f_2 <- c(mz_f, rep(min(mz), times=nrow(a)))
    
    nra <- nrow(a)
    # 
    # col_f <- c(1:nra)
    # col_f[which(mincol== "red")] <- (nra + 1)
    # col_f_2 <- c(col_f,
    #              rep((nra + 2), times=nra))
    # 
    # actual_colors <- c(adjustcolor(viridis(nra), alpha=0.2), "red", adjustcolor("gray50", 0.05))
    # 
    # 
    # if (len(color_list_from_home) == 2) {
    #   col_f_2 <- color_list_from_home$labels
    #   actual_colors <- color_list_from_home$colors
    # }
    # 
    # scatter3D(x=mx_f, 
    #           y=my_f, 
    #           z=mz_f_2, expand=0.85, pch=19, 
    #           colvar=col_f_2, col=actual_colors,
    #           cex=0.3, 
    #           xlab=xlab_t, ylab=ylab_t, zlab=zlab_t,
    #           theta = 25, phi = 40, colkey=F, bty="b2",
    #           
    #           zlim=c(min(mz), max(mz_f)))
    
    # 
    scatter3D(x=mx_f, 
              y=my_f, 
              z=mz_f_2, expand=0.55, pch=19, 
              col=adjustcolor(c("gray60", "white"), alpha=.15), colvar=c(rep(1, times=nrow(a)), rep(2, times=nrow(a))),
              cex=0.3, 
              xlab=xlab_t, ylab=ylab_t, zlab=zlab_t,
              theta = 25, phi = 40, colkey=F, 
              zlim=c(min(mz), max(mz_f)))
  }
  
}

