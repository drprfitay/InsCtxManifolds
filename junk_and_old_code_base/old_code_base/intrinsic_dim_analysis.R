

spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
bupu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "BuPu")))
ylgn_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "YlGn")))
blues_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "Blues")))
rdylbu_cg <- colorRampPalette(rev(brewer.pal(n = 9,  name = "RdYlBu")))
insula_hunger_paths <- c("Y:\\livneh\\itayta\\data\\IC19\\day_150911\\",
                         "Y:\\livneh\\itayta\\data\\IC17\\day_150615\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150406\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150407\\",
                         "Y:\\livneh\\itayta\\data\\IC32\\day_161214\\",
                         "Y:\\livneh\\itayta\\data\\IC42\\day_161117\\")

insula_thirst_paths <- get_thirsty_quenched_paths()

v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")

por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
               "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
               "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\")
#"Y:\\livneh\\itayta\\por_controls\\fov4\\day_150112\\")


paths_all <- c(insula_hunger_paths,
               insula_thirst_paths,
               v1_paths,
               por_paths)


path_categories <- c(rep("Hunger", times=len(insula_hunger_paths)),
                     rep("Thirst", times=len(insula_thirst_paths)),
                     rep("V1", times=len(v1_paths)),
                     rep("POR", times=len(por_paths)))

paths_all <- paths_all[-15]
path_categories <- path_categories[-15]

categories <- path_categories
categories[categories %in% c("V1", "POR")] <- "Control"

id_analysis <- function(path) {
  
  
  id_dim_df <- c()
  ext=""
  for (preset in c("skippy", "dalnatran", "jiffie", "bnd", "shufersal",   "hummus", "nitzat", "taaman")){ #,"milka", "hashahar")) {
    idims <- c()
    for (idx in 1:len(paths_all)) {
    print(idx)
    
    
    work_path <- paths_all[idx]
    control <- (categories[idx] %in% c("Control"))#
    activity_thres <- ifelse(categories[idx] %in% c("Hunger", "Control"),
                              0.25,
                              0.2)
    
    # neur_mat_func <- ifelse(control,
    #                         get_reduced_mat_full_day_control,
    #                         get_reduced_mat_full_day)
    
    mat <- get_mat_with_preset(work_path, "hummus", oldscope=T, activity_threshold=activity_thres)
    
    nr <- max(dim(mat))
    if (control) {
      stim_master_mat <- 
        get_color_palettes_control(work_path, just_mat = T, window_size = 15)
    } else {
      stim_master_mat <- 
        get_stim_mat(work_path,just_mat = T, window_size = 15)             
    }
    
    
    behav <- get_behavior_variables(stim_master_mat, mat_size = nr, window_size = 15, seconds_after=1)
    cp_base <- viridis(nr)
    reward_col <- adjustcolor(cp_base, alpha=0.1); reward_col[behav[,"Reward"] == 1] <- "red"
    tt3_col <- adjustcolor(cp_base, alpha=0.1); tt3_col[behav[,"TrialType 3"] == 1] <- "red"
    tt4_col <- adjustcolor(cp_base, alpha=0.1); tt4_col[behav[,"TrialType 4"] == 1] <- "red"
    tt5_col <- adjustcolor(cp_base, alpha=0.1); tt5_col[behav[,"TrialType 5"] == 1] <- "red"
    cumreward_col <- adjustcolor(rev(spec_cg(len(unique(behav[,"Cumulative reward"]))))[behav[,"Cumulative reward"]],
                                 0.1)
    
    cols = list(Time=cp_base,
             Reward=reward_col,
             TrialType3=tt3_col,
             TrialType4=tt4_col,
             TrialType5=tt5_col,
             CumulReward=cumreward_col)
    
    
    gen_pdf=F
    

    
    dir.create("~/IntrinsicDim")
      
      fmat <- get_mat_with_preset(work_path, preset, oldscope = control, activity_threshold=activity_thres)
      
      
      if (path_categories[idx] %in% c("V1", "POR")) {
        pattern = "fov"
      } else {
        pattern = "IC"
      }
      
      split_r <- str_split(str_split(work_path, pattern)[[1]][[2]], "\\\\")
      
      trim_005 <-  twonn(fmat, method = "linfit", c_trimmed = 0.05)$est[2]
      trim_001 <-  twonn(fmat, method = "linfit", c_trimmed = 0.01)$est[2]
      trim_none <-  twonn(fmat, method = "linfit", c_trimmed = 0)$est[2]
      
      idims <- rbind(idims, c(trim_005, trim_001, trim_none))
      
      dir.create(sprintf("~/IntrinsicDim/%s%s", pattern, split_r[[1]][1]) )
      day_path <- sprintf("~/IntrinsicDim/%s%s/%s", pattern, split_r[[1]][1], split_r[[1]][[2]])
      dir.create(day_path)
      
      if (preset != "hummus") {
        for (var_name in names(cols)) {
          final_path <- sprintf("%s/%s", day_path, var_name)
          print(sprintf("Creating %s!", final_path))
          dir.create(final_path)
          var_col <- cols[[var_name]]
          
          
          png(file=sprintf("%s/%s.png", final_path, preset),height=10,width=10,
              res=300, units="in")
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          png(file=sprintf("%s/%s_small.png", final_path, preset),height=5,width=5,
              res=300, units="in")
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          pdf(sprintf("%s/%s.pdf", final_path, preset),height=10,width=10)
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          downs <- sample(size = 3000, x=1:nrow(fmat))
          
          pdf(sprintf("%s/%s_downsampled.pdf", final_path, preset),height=6,width=6, pointsize = 0.5)
          pairs(fmat[downs,], col=var_col[downs], pch=19)
          dev.off()
        }
      }
    }
    
    
    id_dim_df <- cbind(id_dim_df, idims)
  }
}





get_clusters_mat <- function(mat, nbins=100, min_frames_for_minc=50, min_frames_for_outward=500) {
  
  cm <- colMeans(mat)
  out_n_arr <- floor(seq(min_frames_for_outward, nrow(mat), length.out=nbins))
  
  distance_mat <- apply(mat, 1, function(r) {euc_dist(r, cm)})
  ordered_distance <- order(distance_mat, decreasing = T)
  
  min_c_size_arr <- floor(seq(min_frames_for_minc, floor(nrow(mat) * 0.2), length.out=nbins))
  nc <- c()
  
  for (np in out_n_arr) {
    outmost_n <- ordered_distance[1:np]
    dbc <- dbscan(mat[outmost_n,], min(distance_mat))
    
    clusters_of_interest <- table(dbc$cluster)
    
    
    nclusters_all_sizes <- c()
    for (mc in min_c_size_arr) {
      rs <- as.numeric(names(clusters_of_interest[clusters_of_interest > mc & names(clusters_of_interest) != "0"]))
      nclusters_all_sizes <- c(nclusters_all_sizes, len(rs))
    }
    
    nc <- rbind(nc, nclusters_all_sizes)
  }
  
  colnames(nc) <- min_c_size_arr#sprintf("%.4f%%",seq(min_frames_for_minc, floor(nrow(mat) * 0.2), length.out=nbins) / nrow(mat))
  rownames(nc) <- out_n_arr #sprintf("%.2f%%",seq(min_frames_for_outward, nrow(mat), length.out=nbins) / nrow(mat))
  
  return(nc) 
}

get_clusters_mat_kmeans <- function(mat, nbins=40, min_frames=500, max_nc=20, nreps=20, new_method=F) {
  
  cm <- colMeans(mat)
  out_n_arr <- floor(seq(min_frames, nrow(mat), length.out=nbins))
  
  distance_mat <- apply(mat, 1, function(r) {euc_dist(r, cm)})
  ordered_distance <- order(distance_mat, decreasing = T)
  
  
  cluster_kmeans_mat <- c()
  
  for (np in out_n_arr) {
    outmost_n <- ordered_distance[1:np]
    working_mat <- mat[outmost_n,]
    mse <- c()
    
    for (num_of_clusters in 1:max_nc) {
      
      print(sprintf("Running for %d: %d", np, num_of_clusters))
      reps_mse <- c()  
      
      for (i in 1:nreps) {
        km <- kmeans(working_mat, centers=num_of_clusters, iter.max=100)
        
        reps_mse <- c(reps_mse,
                      sum(km$withinss))
        
      }
      
      
      mse <- cbind(mse,
                   reps_mse)
      
    }
    
    cluster_kmeans_mat <- rbind(cluster_kmeans_mat,
                                colMeans(mse))
  }
  
  colnames(cluster_kmeans_mat) <-  sprintf("# Clusters: %d", 1:max_nc)
  rownames(cluster_kmeans_mat) <- sprintf("Radius %%: %.2f", (out_n_arr) / nrow(mat))
  #colnames(nc) <- min_c_size_arr#sprintf("%.4f%%",seq(min_frames_for_minc, floor(nrow(mat) * 0.2), length.out=nbins) / nrow(mat))
  #rownames(nc) <- out_n_arr #sprintf("%.2f%%",seq(min_frames_for_outward, nrow(mat), length.out=nbins) / nrow(mat))
  
  clusters <- colMeans(cluster_kmeans_mat)
  sphere_radii <- rowMeans(cluster_kmeans_mat)
  
  clusters <- (clusters - min(clusters)) / (max(clusters) - min(clusters))
  sphere_radii <- (sphere_radii - min(sphere_radii)) / (max(sphere_radii) - min(sphere_radii))
  sphere_radii <- sphere_radii * -1  + 1
  
  sphere_radii <- sphere_radii * len(sphere_radii)
  clusters <- clusters * len(clusters)
  
  optimal_cluster <- unlist(lapply(1:len(clusters), function(i) {euc_dist(c(i, clusters[i]), c(0,0))}))
  optimal_sphere_radii <- unlist(lapply(1:len(sphere_radii), function(i) {euc_dist(c(i, sphere_radii[i]), c(0,0))}))
  energy_mat <- do.call(rbind,lapply(1:nrow(old_ckm), function(i) {
  row <- cluster_kmeans_mat[i,]
  v1 = seq(0,1, length.out=nrow(cluster_kmeans_mat))
  v2 = seq(0,16, length.out=len(row))
    
  return(unlist(lapply(1:len(row), function(j) {
      return(euc_dist(c(0,0,0), c(v1[i],v2[j],row[j])));})))
  }))
  
  if (new_method) {

    n_cluster <- which.min(colMeans(energy_mat))
  } else {  
    n_cluster <- which.min(optimal_cluster)
    }
  n_radius <- which.min(optimal_sphere_radii)
  sphere_radius <- out_n_arr[n_radius]
  ind <- ordered_distance[1:sphere_radius]
  working_mat <- mat[ind,]
  
  labels <- c()
  
  for (i in 1:500) {
    km <- kmeans(working_mat, centers=n_cluster, iter.max=100)
    reps_mse <- c(reps_mse,sum(km$withinss))
    labels <- rbind(labels, km$cluster)
  }
  
  hc <- hclust(dist(t(labels)), method = "ward.D2")
  labs <- cutree(hc, k=n_cluster)
  final_labels <- rep(-1, times=nrow(mat))
  final_labels[ind] <- labs
  
  return(list(labs=final_labels,
              clust_labs=labs,
              ind=ind,
              clust_mat=cluster_kmeans_mat,
              nc=n_cluster,
              nr=n_radius,
              energy_mat=energy_mat,
              optimal_cluster=optimal_cluster,
              optimal_radii=optimal_sphere_radii,
              clusters=clusters,
              sphere_radii=sphere_radii))
  }


get_centroid_mat <- function(mat, cluster_mat, nc, min=T) 
{
  cm <- colMeans(mat)
  out_n_arr <- floor(seq(500, nrow(mat), length.out=40))
  
  distance_mat <- apply(mat, 1, function(r) {euc_dist(r, cm)})
  ordered_distance <- order(distance_mat, decreasing = T)
  
  relevant_rows <- which(apply(cluster_mat, 1, function(r) {sum(r == nc)}) != 0)
  if (min) {
    relevant_rows <- relevant_rows[ifelse(len(relevant_rows) > 1, which.min(relevant_rows) + 1, which.min(relevant_rows))]
    relevant_cols <- min(which(cluster_mat[relevant_rows,] == nc))
  } else {
    relevant_rows <- relevant_rows[ifelse(len(relevant_rows) > 1, which.max(relevant_rows) - 1, which.max(relevant_rows))]
    relevant_cols <- max(which(cluster_mat[relevant_rows,] == nc))
  }
  
  
  cluster_sizes <- as.numeric(colnames(cluster_mat))
  outward_d <- as.numeric(rownames(cluster_mat))
  min_cluster_size <- cluster_sizes[relevant_cols]
  outward_d <- outward_d[relevant_rows]
    
  outmost_n <- ordered_distance[1:outward_d]
  dbc <- dbscan(mat[outmost_n,], min(distance_mat))

  clusters_of_interest <- table(dbc$cluster)
  clusters_of_interest <- as.numeric(names(clusters_of_interest[clusters_of_interest > min_cluster_size & names(clusters_of_interest) != "0"]))
  
  timepoint_labels <- dbc$cluster[dbc$cluster %in% clusters_of_interest]
  
  final_labels <- rep(-1, len(timepoint_labels))
  unique_tpl <- unique(timepoint_labels)
  order_tpl <- order(unique_tpl)
  
  for (idx in 1:len(unique_tpl)) {
    final_labels[which(timepoint_labels == unique_tpl[idx])] <- order_tpl[idx]
  }
  
  #assert(table(timepoint_labels) == table(final_labels))
  outmost_n_final <- outmost_n[dbc$cluster %in% clusters_of_interest]
  
  col_all <- rep(adjustcolor("gray70", alpha=0.1), times=nrow(mat))
  col_all[outmost_n_final] <-  spec_cg(max(final_labels))[final_labels]
  all_labels <- rep(-1, times=nrow(mat)) 
  all_labels[outmost_n_final] <- final_labels
  
  centroid_mat<- c()
  
  for (lbl in unique(final_labels)) {
    cluster_ind <- outmost_n_final[which(lbl == final_labels)]
    
    centroid_mat <- rbind(centroid_mat, 
                                 colMeans(mat[cluster_ind,]))
    
  }
  
  return(list(centroid_mat, col_all, all_labels))
}


# 
# pairwise_centroid_mat <- function(paths) {
#   
#   
#   mats <- lapply(insula_thirst_paths, 
#                  function(p) {get_mat_with_preset(p, "dalshtaim")})
#   
#   cluster_mat <- lapply(mats,
#                         function(mt) {get_clusters_mat(mt)})
#   
#   
#   for (nc in c(3:8)) {
#     
#     centroid_mats <- lapply(1:length(insula_thirst_paths),
#                             function())
#     
#   }
# }


cool_high_dimensional_plots <- function(a, cp_from_home=NA, color_list_from_home=list()) {
  
  if (!all(is.na(cp_from_home))) {
    mincol <- cp_from_home
  } else {
    mincol <- adjustcolor('viridis(nrow(a))', alpha=0.75)
    mincol <- adjustcolor(rep("gray70", times=nrow(a)), alpha=0.75)
  }
  
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
    scatter3D(x=mx_f, 
              y=my_f, 
              z=mz_f_2, expand=0.85, pch=19, 
              col=mincol,
              cex=0.3, 
              xlab=xlab_t, ylab=ylab_t, zlab=zlab_t,
              theta = 25, phi = 40, colkey=F, bty="b2",
              
              zlim=c(min(mz), max(mz_f)))
    
  }
  
}
cool_high_dimensional_plots_2 <- function(a, cp_from_home=NA) {
  
  if (!all(is.na(cp_from_home))) {
    mincol <- cp_from_home
  } else {
    mincol <- adjustcolor(viridis(nrow(a)), alpha=0.75)
  }
  #combinations_mat <- combn(sample(1:ncol(a), ncol(a)), 3)
  combinations_mat <- combn(1:6, 3)
  
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
    mx <- a[,xdim]
    my <- a[,ydim]
    nra <- nrow(a)
    
    sdall <- sd(c(mz, mx, my))
    sdmx <- sd(mx)
    sdmy <- sd(my)
    sdmz <- sd(mz)
    
    xlim_f <- c(min(mx) - sdmx, 
                max(mx) + sdmx)
    
    ylim_f <- c(min(my) - sdmy, 
                max(my) + sdmy)
    
    zlim_f <- c(min(mz) - sdmz, 
                max(mz) + sdmz)
    
    
    my_f <- c(my,                      my,                      rep(max(my), times=nra))
    mx_f <- c(mx,                      rep(min(mx), times=nra), mx)
    mz_f <- c(rep(min(mz), times=nra), mz,                      mz)
    
    col_f <- c(1:nra,
               1:nra,
               1:nra)
    
    mincol_f <- c(mincol, mincol, mincol)
    
    col_f[which(mincol_f == "red")] <- (nra + 1)
    
    
    scatter3D(x=mx_f, 
              y=my_f, 
              z=mz_f, expand=0.85, pch=19, colvar=(col_f), col=c(adjustcolor(viridis(nra), alpha=0.1), "red"), cex=0.5, 
              xlab=xlab_t, ylab=ylab_t, zlab=zlab_t,
              theta = 25, phi = 40, colkey=F, bty="b2")
    #xlim=xlim_f, ylim=ylim_f, zlim=zlim_f)
    
  }
  
}



get_colors_for_mat <- function(path, mat, chunk=-1) {

  stim_master_mat <- 
    get_stim_mat(path,just_mat = T, window_size = 15, chunk=chunk)    
  behav <- get_behavior_variables(stim_master_mat, mat_size = nrow(mat), window_size = 15, seconds_after=1)
  cp_base <- viridis(nrow(mat))
  reward_col <- adjustcolor(cp_base, alpha=0.1); reward_col[behav[,"Reward"] == 1] <- "red"
  tt3_col <- adjustcolor(cp_base, alpha=0.1); tt3_col[behav[,"TrialType 3"] == 1] <- "red"
  tt4_col <- adjustcolor(cp_base, alpha=0.1); tt4_col[behav[,"TrialType 4"] == 1] <- "red"
  tt5_col <- adjustcolor(cp_base, alpha=0.1); tt5_col[behav[,"TrialType 5"] == 1] <- "red"
  cumreward_col <- adjustcolor(rev(spec_cg(len(unique(behav[,"Cumulative reward"]))))[behav[,"Cumulative reward"]],
                               0.1)
  
  
  return(list(time=cp_base,
              reward=reward_col,
              tt3_col=tt3_col,
              tt4_col=tt4_col,
              tt5_col=tt5_col,
              cumreward_col=cumreward_col))
}



facial_expression_analysis <- function(path) {
  facial_expression_paths <- c("Y:\\livneh\\itayta\\data\\IC71\\day_190903",
                               "Y:\\livneh\\itayta\\data\\IC71\\day_190904",
                               "Y:\\livneh\\itayta\\data\\IC71\\day_190906",
                               "Y:\\livneh\\itayta\\data\\IC77\\day_190905",
                               "Y:\\livneh\\itayta\\data\\IC77\\day_190906",
                               "Y:\\livneh\\itayta\\data\\IC77\\day_190908")
  
  for (p in facial_expression_paths) {
        mat <- get_mat_with_preset(p, "dalshtaim", chunk=1) 
        org <- get_reduced_mat_full_day(p, window_size=15, chunk=1, just_original_mat = T)
        org <- t(org)
        color_palettes <- get_colors_for_mat(p, mat, chunk=1)
        
        png(file=sprintf("%s\\pairs_reward.png", p), units="in", height=7, width=7, res=300)
        pairs(mat, col=color_palettes$reward, pch=19)
        dev.off()
  
        png(file=sprintf("%s\\pairs_time.png", p), units="in", height=7, width=7, res=300)
        pairs(mat, col=color_palettes$time, pch=19)
        dev.off()
        
        clust_mat <- get_clusters_mat(mat, min_frames_for_minc = 5, min_frames_for_outward = 5)
        #nc=4
        centroid_mat <- get_centroid_mat(mat, clust_mat, max(c(clust_mat)))
        
        
        
        png(file=sprintf("%s\\pairs_clustered.png", p), units="in", height=7, width=7, res=300)
        pairs(mat, col=centroid_mat[[2]], pch=19)
        dev.off()
        

        
        
        timepoints <- 
        lapply(unique(centroid_mat[[3]][centroid_mat[[3]] != -1]), 
               function(clst) {as.numeric(which(centroid_mat[[3]] != -1 & centroid_mat[[3]] == clst) * 15)})
        
        max_len <- max(unlist(lapply(timepoints, length)))
        
        for(idx in 1:length(timepoints)) {
          length(timepoints[[idx]]) <-  max_len
        }
        
        timepoints_df <- do.call(cbind, timepoints)
        
        
        write.csv(file=sprintf("%s\\time_points.csv", p),
                  timepoints_df)
        
        colors <- c()
        cluster_labels <- c()
        activity_mat <- c()
        neuronal_clust_labels <- centroid_mat[[3]]
        num_of_rewards <- c()
        
        for (lbl in sort(unique(neuronal_clust_labels[neuronal_clust_labels != - 1]))) {
            indices <- which(neuronal_clust_labels == lbl)
            colors <- c(colors, unique(centroid_mat[[2]][indices]))
            cluster_labels <- c(lbl, cluster_labels)
            activity_mat <- cbind(activity_mat ,rowMeans(org[,indices]))
            
            if ("red" %in% names(table(color_palettes$reward[indices]))) {
              num_of_rewards <- c(table(color_palettes$reward[indices])["red"], num_of_rewards)
            } else {
              num_of_rewards <- c(0, num_of_rewards)
            }
        }
        
        
        
        annot_df = data.frame(cluster_labels=cluster_labels)
        annot_col = list(cluster_labels=colors)
        rownames(annot_df) <- cluster_labels
        colnames(activity_mat) <- cluster_labels
        
        reward_ratio <- data.frame(reward=num_of_rewards / sum(num_of_rewards))
        rownames(reward_ratio) <- cluster_labels
        
        png(file=sprintf("%s\\neurons_cluster_activity.png", p), units="in", height=7, width=7, res=300)
        pheatmap(activity_mat, cluster_rows=F, cluster_cols=F, col=adjustcolor(rev(inferno(20)), alpha=0.8), border_color="NA",
                 annotation_col = annot_df,
                 annotation_colors = annot_col,
                 show_colnames = T)
        dev.off()
        
        png(file=sprintf("%s\\cluster_reward_ratio.png", p), units="in", height=4, width=4, res=300)
        pheatmap(reward_ratio, 
                 cluster_rows=F, cluster_cols=F, border_color="NA",
                 annotation_row = annot_df,
                 annotation_colors = annot_col,
                 show_colnames = F)
        dev.off()
  }
  }

