get_clusters_mat_kmeans <- function(path,
                                    nbins=40,
                                    min_frames=500,
                                    max_nc=20,
                                    nreps=20,
                                    override=F,
                                    cluster_maps=F,
                                    chunk=-1,
                                    red_preset="dalshtaim",
                                    clustermaps_subpath="cluster_maps",
                                    new_method=F) {

    if (all(chunk != -1)) {
      if(len(chunk) > 1) {
        chunk_str <- paste(chunk, collapse="_")
      } else {
        chunk_str = sprintf("%d", chunk)
      }
    } else {
      chunk_str=""
    }

    cluster_map_name <- sprintf("%s_nbins%d_minframes%d_maxnc%d_nreps%d%s",
                                red_preset,
                                nbins,
                                min_frames,
                                max_nc,
                                nreps,
                                chunk_str)



  # Firstly check whether there exists a matrice for that day
  if (!override && clustermaps_subpath %in% list.dirs(path, recursive = F, full.names = F)) {

    print(path)
    print(cluster_map_name)
    if (sprintf("%s.R", cluster_map_name) %in%
        list.files(sprintf("%s\\%s\\", path, clustermaps_subpath), full.names = F)) {

      print(sprintf("Found cluster map %s, loading", cluster_map_name))
      print(sprintf("%s\\%s\\%s.R", path, clustermaps_subpath, cluster_map_name))
      load(sprintf("%s\\%s\\%s.R", path, clustermaps_subpath, cluster_map_name), verbose=F)
      return(cluster_map_final)
    }
  }


  if (red_preset == "Isomap") {
    mat <- get_reduced_mat_full_day(path, ndim = 6, type = "isomap", knn1=100, knn2=0, window_size=15)
  } else {
    mat <- get_mat_with_preset(path, preset_name =  red_preset, chunk=chunk)
  }

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

  energy_mat <- NA

  if (new_method) {
    energy_mat <- do.call(rbind,lapply(1:nrow(old_ckm), function(i) {
      row <- cluster_kmeans_mat[i,]
      v1 = seq(0,1, length.out=nrow(cluster_kmeans_mat))
      v2 = seq(0,16, length.out=len(row))

      return(unlist(lapply(1:len(row), function(j) {
        return(euc_dist(c(0,0,0), c(v1[i],v2[j],row[j])));})))
    }))

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


  cluster_map_final <- list(labs=final_labels,
                            clust_labs=labs,
                            ind=ind,
                            clust_mat=cluster_kmeans_mat,
                            nc=n_cluster,
                            nr=n_radius,
                            energy_mat=energy_mat,
                            optimal_cluster=optimal_cluster,
                            optimal_radii=optimal_sphere_radii,
                            clusters=clusters,
                            sphere_radii=sphere_radii)

  dir.create(sprintf("%s\\%s", path, clustermaps_subpath))
  save(cluster_map_final, file=sprintf("%s\\%s\\%s.R", path, clustermaps_subpath, cluster_map_name))
  print(sprintf("Saving cluster map! %s", cluster_map_name))

  return(cluster_map_final)
}

