get_pupil_vec <- function(path, window_size=15, raw_vec=F) {
  pupil_files <- list.files(sprintf("%s\\pupil_files", path), full.names = T)
  
  pupil_mat <- lapply(pupil_files,
                      function(p) {
                        mt <- readMat(p)
                        return(mt)
                      })
  
  
  pupil_vec <- unlist(lapply(pupil_mat,
                       function(mt) {
                         return(mt$pupil.mat[,3])
                       }))
  
  if (raw_vec) {
    return(pupil_vec)
  }
  
  pupil_vec_d <- rep(pupil_vec, each=2)
  binned_vec <- time_bin_average_vec(pupil_vec_d, window_size)
  
  return(binned_vec)
}


convert_pupil_vec_into_color <- function(pupil_vec, nbins = 100) {
  spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
  return(spec_cg(nbins)[as.numeric(cut(pupil_vec, breaks = nbins))])
}


output_path <- sprintf("%s//pupil_analysis//", output_path_f)
dir.create(output_path)

path_list = insula_thirst_paths[1:3]



mice_name_indices <- unlist(gregexpr("IC[0-9]{2}", path_list))
mice_names <- unlist(lapply(1:len(mice_name_indices), 
                            function(i) {substr(path_list[[i]], mice_name_indices[[i]], mice_name_indices[[i]] + 3)}))
days_indices <- unlist(gregexpr("day_[0-9]{6}", path_list))
days <- unlist(lapply(1:len(days_indices),
                      function(i) {substr(path_list[[i]], days_indices[[i]] + 4, days_indices[[i]] + 9)}))
datasets_names <- paste(mice_names, days, sep = " ")
annotation_df <- data.frame(Mice=mice_names)


for (path_idx in 1:len(path_list)) {

  
  path = path_list[path_idx]

  pupil_vec <- get_pupil_vec(path)
  raw_vec <- get_pupil_vec(path, raw_vec=T)
  rolled_vc <- rollmean(pupil_vec, 20)
  rolled_vc <- c(rolled_vc, rep(rolled_vc[len(rolled_vc)], 19))
  org <- get_reduced_mat_full_day(path, window_size = 15, just_original_mat = T)
  
  pup_col <- convert_pupil_vec_into_color(rolled_vc)
  pup_col_raw <- convert_pupil_vec_into_color(raw_vec)  
  if (ncol(org) > nrow(org)) {
    org <- t(org)
  }
  mt <- get_mat_with_preset(path, "dalshtaim")
  
  cc <- cancor(org,  as.matrix(rolled_vc))
  projected_mode <- org %*% cc$xcoef[,1]
  
  
  cor_value <- cor(projected_mode, rolled_vc)
  coeff <- (cor_value / abs(cor_value))[1]
  df <- data.frame(mode=scale(coeff * projected_mode),
                   variable=scale(rolled_vc),
                   time=1:len(rolled_vc),
                   col=pup_col)
  
  df_raw <- data.frame(variable=raw_vec,
                       time=1:len(raw_vec),
                       col=pup_col_raw)
  
  gf <- 
  ggplot(df, aes(x=time, y=mode)) + geom_line(col="red", alpha=0.8) + 
    geom_line(aes(y=variable), col="black", alpha=0.8) +
    theme_classic() + 
    theme(plot.title=element_text(size=9)) + 
    xlab("Time") + 
    ylab("Mode") + 
    ggtitle(sprintf("Variable (black), Neural mode (red), %f",
                    abs(cor_value)))
  
  
  gf2 <-  
    ggplot(df, aes(x=time, y=variable)) + geom_line(color=df$col) + 
    theme_classic() + 
    theme(plot.title=element_text(size=9)) +
    xlab("Time") + 
    ylab("Pupil size") + 
    ggtitle("Pupil- smoothed")
  
  
  gf_raw <-  
    ggplot(df_raw, aes(x=time, y=variable)) + geom_line(color=df_raw$col) + 
    theme_classic() + 
    theme(plot.title=element_text(size=9)) + 
    xlab("Time") + 
    ylab("Pupil size") +
    ggtitle("Pupil - raw")
  
  
  gfinal <- grid.arrange(gf, gf2, gf_raw, nrow=1)
    
  g <- pair_plots_colored(mt, adjustcolor(df$col, alpha=.5))
  
  dir.create(sprintf("%s//%s", 
                     output_path,
                     str_replace(datasets_names[path_idx], " ", "_")))


  png(sprintf("%s\\%s\\structure_by_pupil.png", 
              output_path,
              str_replace(datasets_names[path_idx], " ", "_")),
      units="in",
      res=1000,
      height=8,
      width=6)
  
  plot(g)
  dev.off()
  
  
  
  png(sprintf("%s\\%s\\pupil_mode_cca.png", 
              output_path,
              str_replace(datasets_names[path_idx], " ", "_")),
      units="in",
      res=1000,
      height=3,
      width=9)
  
  plot(gfinal)
  dev.off()
}
