figures_base_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\figures\\"
base_output_path = "Y:\\livneh\\itayta\\Itay_group_meeting_dfs\\code_base_for_paper\\final_figures\\"
output_path <- sprintf("%s\\figures_structure_comparision\\", base_output_path) 

dir.create(base_output_path)
dir.create(output_path)


pairwise_topological_permutation_generic <- function(metadata_1, 
                                                     metadata_2, 
                                                     remove_diagonal = F, 
                                                     niter = 20, 
                                                     single_example=F)
{
  library(TDAstats)
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
  
  save(file=sprintf("%s\\data\\figure_1\\%s",
                    base_output_path,
                    output_name),
       final_res)
  
  return(final_res)
}

TT_shuffle_topo_similarity_statistics <- function(df_final, metric)
{
  
  
  statistics_df <- data.frame()
  
  df_tt_across <- df_final[df_final$group == "T-T_Across",]
  df_tshuff_across <- df_final[df_final$group == "T-Shuff_Across",]
  df_tt_within <- df_final[df_final$group == "T-T_Within",]
  df_tshuff_within <- df_final[df_final$group == "T-Shuff_Within",]
  
  # Sort dataframes for pairing
  df_tt_across <- df_tt_across[order(paste(df_tt_across$path1, df_tt_across$path2)),]
  df_tshuff_across <- df_tshuff_across[order(paste(df_tshuff_across$path1, df_tshuff_across$path2)),]
  df_tt_within <- df_tt_within[order(paste(df_tt_within$path1, df_tt_within$path2)),]
  df_tshuff_within <- df_tshuff_within[order(paste(df_tshuff_within$path1, df_tshuff_within$path2)),]
  
  
  N_comparisions = 3      
  # Statistical tests
  wilcox_across_b0 <- 
    wilcox.test(df_tt_across$betti0,
                df_tshuff_across$betti0,
                paired=T,
                correct = F,
                alternative = "less")
  
  
  wilcox_across_b1 <- 
    wilcox.test(df_tt_across$betti1,
                df_tshuff_across$betti1,
                paired=T,
                correct = F,
                alternative = "less")
  
  
  ks_across_b0 <- 
    ks.test(df_tt_across$betti0,
            df_tshuff_across$betti0,
            exact = T,
            alternative = "greater")
  
  
  ks_across_b1 <- 
    ks.test(df_tt_across$betti1,
            df_tshuff_across$betti1,
            exact = T,
            alternative = "greater")
  
  
  wilcox_within_b0 <- 
    wilcox.test(df_tt_within$betti0,
                df_tshuff_within$betti0,
                paired=T,
                correct = F,
                alternative = "less")
  
  
  wilcox_within_b1 <- 
    wilcox.test(df_tt_within$betti1,
                df_tshuff_within$betti1,
                paired=T,
                correct = F,
                alternative = "less")
  
  
  ks_within_b0 <- 
    ks.test(df_tt_within$betti0,
            df_tshuff_within$betti0,
            exact = T,
            alternative = "greater")
  
  
  ks_within_b1 <- 
    ks.test(df_tt_within$betti1,
            df_tshuff_within$betti1,
            exact = T,
            alternative = "greater")
  
  wilcox_within_across_b0 <- 
    wilcox.test(df_tt_within$betti0,
                df_tt_across$betti0,
                correct = F,
                alternative = "two.sided")
  
  
  wilcox_within_across_b1 <- 
    wilcox.test(df_tt_within$betti1,
                df_tt_across$betti1,
                correct = F,
                alternative = "two.sided")
  
  
  ks_within_across_b0 <- 
    ks.test(df_tt_within$betti0,
            df_tt_across$betti0,
            exact = T,
            alternative = "two.sided")
  
  
  ks_within_across_b1 <- 
    ks.test(df_tt_within$betti1,
            df_tt_across$betti1,
            exact = T,
            alternative = "two.sided")
  
  
  
  stat_df <-
    data.frame(comparision="wilcox_across_b0",
               method=wilcox_across_b0$method,
               alternative=wilcox_across_b0$alternative,
               pval=wilcox_across_b0$p.value,
               adjusted.pval=wilcox_across_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_across$betti0),
               mean_group_b=mean(df_tshuff_across$betti0),
               sd_group_a=sd(df_tt_across$betti0),
               sd_group_b=sd(df_tshuff_across$betti0),
               sem_group_a=sem(df_tt_across$betti0),
               sem_group_b=sem(df_tshuff_across$betti0),
               N_group_a=len(df_tt_across$betti0),
               N_group_b=len(df_tshuff_across$betti0))
  
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_across_b1",
               method=wilcox_across_b1$method,
               alternative=wilcox_across_b1$alternative,
               pval=wilcox_across_b1$p.value,
               adjusted.pval=wilcox_across_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_across$betti1),
               mean_group_b=mean(df_tshuff_across$betti1),
               sd_group_a=sd(df_tt_across$betti1),
               sd_group_b=sd(df_tshuff_across$betti1),
               sem_group_a=sem(df_tt_across$betti1),
               sem_group_b=sem(df_tshuff_across$betti1),
               N_group_a=len(df_tt_across$betti1),
               N_group_b=len(df_tshuff_across$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_across_b0",
               method=ks_across_b0$method,
               alternative=ks_across_b0$alternative,
               pval=ks_across_b0$p.value,
               adjusted.pval=ks_across_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_across$betti0),
               mean_group_b=mean(df_tshuff_across$betti0),
               sd_group_a=sd(df_tt_across$betti0),
               sd_group_b=sd(df_tshuff_across$betti0),
               sem_group_a=sem(df_tt_across$betti0),
               sem_group_b=sem(df_tshuff_across$betti0),
               N_group_a=len(df_tt_across$betti0),
               N_group_b=len(df_tshuff_across$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_across_b1",
               method=ks_across_b1$method,
               alternative=ks_across_b1$alternative,
               pval=ks_across_b1$p.value,
               adjusted.pval=ks_across_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_across$betti1),
               mean_group_b=mean(df_tshuff_across$betti1),
               sd_group_a=sd(df_tt_across$betti1),
               sd_group_b=sd(df_tshuff_across$betti1),
               sem_group_a=sem(df_tt_across$betti1),
               sem_group_b=sem(df_tshuff_across$betti1),
               N_group_a=len(df_tt_across$betti1),
               N_group_b=len(df_tshuff_across$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_within_b0",
               method=wilcox_within_b0$method,
               alternative=wilcox_within_b0$alternative,
               pval=wilcox_within_b0$p.value,
               adjusted.pval=wilcox_within_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti0),
               mean_group_b=mean(df_tshuff_within$betti0),
               sd_group_a=sd(df_tt_within$betti0),
               sd_group_b=sd(df_tshuff_within$betti0),
               sem_group_a=sem(df_tt_within$betti0),
               sem_group_b=sem(df_tshuff_within$betti0),
               N_group_a=len(df_tt_within$betti0),
               N_group_b=len(df_tshuff_within$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_within_b1",
               method=wilcox_within_b1$method,
               alternative=wilcox_within_b1$alternative,
               pval=wilcox_within_b1$p.value,
               adjusted.pval=wilcox_within_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti1),
               mean_group_b=mean(df_tshuff_within$betti1),
               sd_group_a=sd(df_tt_within$betti1),
               sd_group_b=sd(df_tshuff_within$betti1),
               sem_group_a=sem(df_tt_within$betti1),
               sem_group_b=sem(df_tshuff_within$betti1),
               N_group_a=len(df_tt_within$betti1),
               N_group_b=len(df_tshuff_within$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_within_b0",
               method=ks_within_b0$method,
               alternative=ks_within_b0$alternative,
               pval=ks_within_b0$p.value,
               adjusted.pval=ks_within_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti0),
               mean_group_b=mean(df_tshuff_within$betti0),
               sd_group_a=sd(df_tt_within$betti0),
               sd_group_b=sd(df_tshuff_within$betti0),
               sem_group_a=sem(df_tt_within$betti0),
               sem_group_b=sem(df_tshuff_within$betti0),
               N_group_a=len(df_tt_within$betti0),
               N_group_b=len(df_tshuff_within$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  stat_df <-
    data.frame(comparision="ks_within_b1",
               method=ks_within_b1$method,
               alternative=ks_within_b1$alternative,
               pval=ks_within_b1$p.value,
               adjusted.pval=ks_within_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti1),
               mean_group_b=mean(df_tshuff_within$betti1),
               sd_group_a=sd(df_tt_within$betti1),
               sd_group_b=sd(df_tshuff_within$betti1),
               sem_group_a=sem(df_tt_within$betti1),
               sem_group_b=sem(df_tshuff_within$betti1),
               N_group_a=len(df_tt_within$betti1),
               N_group_b=len(df_tshuff_within$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_within_across_b0",
               method=wilcox_within_across_b0$method,
               alternative=wilcox_within_across_b0$alternative,
               pval=wilcox_within_across_b0$p.value,
               adjusted.pval=wilcox_within_across_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti0),
               mean_group_b=mean(df_tt_across$betti0),
               sd_group_a=sd(df_tt_within$betti0),
               sd_group_b=sd(df_tt_across$betti0),
               sem_group_a=sem(df_tt_within$betti0),
               sem_group_b=sem(df_tt_across$betti0),
               N_group_a=len(df_tt_within$betti0),
               N_group_b=len(df_tt_across$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_within_across_b1",
               method=wilcox_within_across_b1$method,
               alternative=wilcox_within_across_b1$alternative,
               pval=wilcox_within_across_b1$p.value,
               adjusted.pval=wilcox_within_across_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti1),
               mean_group_b=mean(df_tt_across$betti1),
               sd_group_a=sd(df_tt_within$betti1),
               sd_group_b=sd(df_tt_across$betti1),
               sem_group_a=sem(df_tt_within$betti1),
               sem_group_b=sem(df_tt_across$betti1),
               N_group_a=len(df_tt_within$betti1),
               N_group_b=len(df_tt_across$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_within_across_b0",
               method=ks_within_across_b0$method,
               alternative=ks_within_across_b0$alternative,
               pval=ks_within_across_b0$p.value,
               adjusted.pval=ks_within_across_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti0),
               mean_group_b=mean(df_tt_across$betti0),
               sd_group_a=sd(df_tt_within$betti0),
               sd_group_b=sd(df_tt_across$betti0),
               sem_group_a=sem(df_tt_within$betti0),
               sem_group_b=sem(df_tt_across$betti0),
               N_group_a=len(df_tt_within$betti0),
               N_group_b=len(df_tt_across$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)  
  
  stat_df <-
    data.frame(comparision="ks_within_across_b1",
               method=ks_within_across_b1$method,
               alternative=ks_within_across_b1$alternative,
               pval=ks_within_across_b1$p.value,
               adjusted.pval=ks_within_across_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_tt_within$betti1),
               mean_group_b=mean(df_tt_across$betti1),
               sd_group_a=sd(df_tt_within$betti1),
               sd_group_b=sd(df_tt_across$betti1),
               sem_group_a=sem(df_tt_within$betti1),
               sem_group_b=sem(df_tt_across$betti1),
               N_group_a=len(df_tt_within$betti1),
               N_group_b=len(df_tt_across$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  
  return(statistics_df)
}

plot_TT_shuffle_topo_similarity <- function(average_across_datastes=T)
{
  
  base_figure_path <- sprintf("%s\\TT_shuffle_topo_similarity_boxplots", output_path)
  dir.create(base_figure_path)
  
  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path)
  
  
  output_name = "thirst_thirst_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  tt <- final_res
  
  output_name = "thirst_shuffle_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  tshuff <- final_res
  
  
  statistics_df <- data.frame()
  
  
  
  for (metric in c("JSD", "KL", "wasserstein")) {
      if (average_across_datastes) {
        figure_path <- sprintf("%s\\averaged_over_datasets\\", base_figure_path)
        dir.create(figure_path)
        figure_path <- sprintf("%s\\%s", figure_path, metric)
      } else {
        figure_path <- sprintf("%s\\%s", base_figure_path, metric)
      }
      dir.create(figure_path)

      df_all <- rbind(tt$topological_similarity[tt$topological_similarity$metric == metric,],
                      tshuff$topological_similarity[tshuff$topological_similarity$metric == metric,])
      
      cmps <- c(nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]),
                nrow(tshuff$topological_similarity[tshuff$topological_similarity$metric == metric,]))
      
      df_all$group = rep(c("T-T", "T-Shuff"), cmps)
      df_all <- as.data.frame(df_all)
      df_all$betti0 <- as.numeric(df_all$betti0)
      df_all$betti1 <- as.numeric(df_all$betti1)
      df_final <- df_all[!(df_all$path1 == df_all$path2),]
      
      mice_names_for_path <- unlist(lapply(str_split(get_datasets_names(get_thirsty_quenched_paths()), " "), function(lst) {lst[[1]]}))
      names(mice_names_for_path) <- 1:14
      
      df_final$mice_path1 <- mice_names_for_path[df_final$path1]
      df_final$mice_path2 <- mice_names_for_path[df_final$path2]
      df_final$original_group <- df_final$group
      df_final$group <- paste(df_final$group, ifelse(as.numeric(df_final$mice_path1 == df_final$mice_path2) > 0, "Within", "Across"), sep="_")
      
      
      if (average_across_datastes) {
        
      df_final <- 
      ddply(df_final,
            
            .(original_group),
            function(group_df) {
              res <- 
              ddply(group_df,
                    .(path1),
                    function(mice_path_df) {
                      
                      #print(mice_path_df)
                      
                      across <- mice_path_df[mice_path_df$mice_path1 != mice_path_df$mice_path2,]
                      within  <- mice_path_df[mice_path_df$mice_path1 == mice_path_df$mice_path2,]
                      
                      across_df <- data.frame(betti0=mean(as.numeric(across$betti0)),
                                    betti1=mean(as.numeric(across$betti1)),
                                    group=across$group[1])
                      
                      #print(nrow(within))
                      
                      if (nrow(within) > 0) {
                        within_df <- data.frame(betti0=mean(as.numeric(within$betti0)),
                                       betti1=mean(as.numeric(within$betti1)),
                                       group=within$group[1])
                        
                        across_df <- rbind(across_df,
                                           within_df)
                      }
                      
                      
                      
                      
                      return(across_df)
                    })
              
            })
      }
      
      plot_groups <- c("T-T_Within", "T-Shuff_Within", "T-T_Across", "T-Shuff_Across")
      
      gb0 <- 
      ggplot(df_final, 
             aes(x=factor(group, levels=plot_groups), 
                 y=log2(as.numeric(betti0))), 
                group=group)
        #geom_violin(aes(fill=group),color=NA, fill="gray65")
        #geom_boxplot(fill=NA)
      
        if (average_across_datastes) {
           #gb0 <- 
          gb0 +
            # stat_summary(
            #   fun.min = function(z) { quantile(z,0.25) },
            #   fun.max = function(z) { quantile(z,0.75) },
            #   fun = mean,
            #   size=1) +
            geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
             scale_shape_manual(values=18) +
            geom_boxplot(fill=NA, width=.5, ) 
        } else {
          gb0 <- gb0 + 
          geom_violin(aes(fill=group),color=NA, fill="gray65")
          geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.25)) + 
          stat_summary(size=1) 
        }
        
        gb0 <- gb0 +
        theme_light() + 
        big_text_base_plot_theme + 
        ylab(TeX(sprintf("%s (%s)", "Topologic divergence", "$\\beta_{0}$"))) + 
        xlab("")  +
        ggtitle(metric) +
        theme(plot.title = element_text(size=9))
        
      
        
        gb1 <- 
          ggplot(df_final, 
                 aes(x=factor(group, levels=plot_groups), 
                     y=log2(as.numeric(betti1))), 
                 group=group) 
          #geom_violin(aes(fill=group),color=NA, fill="gray65")
          #geom_boxplot(fill=NA)
        
        if (average_across_datastes) {
          gb1 <- gb1 +
            # stat_summary(
            #   fun.min = function(z) { quantile(z,0.25) },
            #   fun.max = function(z) { quantile(z,0.75) },
            #   fun = mean,
            #   size=1) +
            geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
            geom_boxplot(fill=NA, width=.5, ) 
        } else {
          gb1 <- gb1 + 
            geom_violin(aes(fill=group),color=NA, fill="gray65")
            geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.25)) + 
            stat_summary(size=1) 
        }
        
        gb1 <- gb1 +
          theme_light() + 
          big_text_base_plot_theme + 
          ylab(TeX(sprintf("%s (%s)", "Topologic divergence", "$\\beta_{1}$"))) + 
          xlab("")  +
          ggtitle(metric) +
          theme(plot.title = element_text(size=9))

      stats <- TT_shuffle_topo_similarity_statistics(df_final, metric)
      stats$signif = signif.num(stats$pval)
      stats$adjusted.signif = signif.num(stats$adjusted.pval)
      statistics_df <- rbind(statistics_df,
                             stats)
      
      colnames(tt$pvalue_b0[[metric]]) <- 1:ncol(tt$pvalue_b0[[metric]])
      rownames(tt$pvalue_b0[[metric]]) <- 1:ncol(tt$pvalue_b0[[metric]])
      colnames(tt$pvalue_b1[[metric]]) <- 1:ncol(tt$pvalue_b1[[metric]])
      rownames(tt$pvalue_b1[[metric]]) <- 1:ncol(tt$pvalue_b1[[metric]])
      meltedb0 <- melt(tt$pvalue_b0[[metric]])
      meltedb1 <- melt(tt$pvalue_b1[[metric]])
      
      
      pvalues_df <- data.frame(pval_b0=meltedb0$value[meltedb0$Var1 != meltedb0$Var2],
                               pval_b1=meltedb1$value[meltedb0$Var1 != meltedb1$Var2])
      pvalues_df$group = "T-T"
      
      #group=paste("T-T", ifelse(as.numeric(melted$Var1 == melted$Var2) > 0, "Within", "Across"), sep="_"))
      
      
      colnames(tshuff$pvalue_b0[[metric]]) <- 1:ncol(tshuff$pvalue_b0[[metric]])
      rownames(tshuff$pvalue_b0[[metric]]) <- 1:ncol(tshuff$pvalue_b0[[metric]])
      colnames(tshuff$pvalue_b1[[metric]]) <- 1:ncol(tshuff$pvalue_b1[[metric]])
      rownames(tshuff$pvalue_b1[[metric]]) <- 1:ncol(tshuff$pvalue_b1[[metric]])
      meltedb0 <- melt(tshuff$pvalue_b0[[metric]])
      meltedb1 <- melt(tshuff$pvalue_b1[[metric]])
      
      tmp_df <- data.frame(pval_b0=meltedb0$value[meltedb0$Var1 != meltedb0$Var2],
                           pval_b1=meltedb1$value[meltedb0$Var1 != meltedb1$Var2])
      #group=paste("T-Shuff", ifelse(as.numeric(melted$Var1 == melted$Var2) > 0, "Within", "Across"), sep="_")))
      tmp_df$group = "T-Shuffle"
      
      pvalues_df <- rbind(pvalues_df,
                          tmp_df)
      
      pvalues_fraction_hist_df <- 
        ddply(pvalues_df, .(group), 
              function(group_df) {
                hb0 <-  hist(group_df$pval_b0, breaks=seq(0, 1,length.out=21), plot=F)
                hb1 <-  hist(group_df$pval_b1, breaks=seq(0, 1,length.out=21), plot=F)
                
                
                hist_df <- 
                  data.frame(frac_pval_b0=hb0$counts / sum(hb0$counts), 
                             breaks_pval_b0=hb0$breaks[-1], 
                             frac_pval_b1=hb1$counts / sum(hb1$counts), 
                             breaks_pval_b1=hb1$breaks[-1], 
                             group=rep(group_df$group[1],len(hb0$counts)))
                
                return(hist_df)
              })
      
      
      gpvb0 <- 
      ggplot(pvalues_fraction_hist_df) + 
        geom_bar(aes(y=frac_pval_b0, x=breaks_pval_b0, group=group, fill=group), color="black", stat="identity", position=position_nudge(.05), width=.04, alpha=.55) + 
        scale_y_continuous(expand=c(0,0)) +
        ylab("Fraction (%)") +
        xlab("P-value") + 
        theme_light() + 
        #geom_vline(xintercept = .1) +
        big_text_base_plot_theme + 
        ggtitle(TeX(sprintf("%s (%s)", metric, "$\\beta_{0}$"))) +
        theme(plot.title = element_text(size=9))
      
      
      gpvb1 <- 
        ggplot(pvalues_fraction_hist_df) + 
        geom_bar(aes(y=frac_pval_b1, x=breaks_pval_b1, group=group, fill=group), color="black", stat="identity", position=position_nudge(.05), width=.04, alpha=.55) + 
        scale_y_continuous(expand=c(0,0)) +
        ylab("Fraction (%)") +
        xlab("P-value") + 
        theme_light() + 
        #geom_vline(xintercept = .1) +
        big_text_base_plot_theme + 
        ggtitle(TeX(sprintf("%s (%s)", metric, "$\\beta_{1}$"))) +
        theme(plot.title = element_text(size=9))
      
      
      
      gpval_violin_b0 <- 
      ggplot(pvalues_df, aes(x=group, y=pval_b0, group=group)) +
        geom_violin(aes(fill=group),color=NA, fill="gray60") + 
        geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
        #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
        stat_summary(size=1) + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Pvalue") + 
        xlab("") 
      
      gpval_violin_b1 <- 
        ggplot(pvalues_df, aes(x=group, y=pval_b1, group=group)) +
        geom_violin(aes(fill=group),color=NA, fill="gray60") + 
        geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
        #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
        stat_summary(size=1) + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Pvalue") + 
        xlab("") 
      
      for (size_name in names(a4_sizes)) {
        size = a4_sizes[[size_name]]
        pdf(sprintf("%s\\%s_b0_violin.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gb0)
        
        dev.off()
        
        pdf(sprintf("%s\\%s_b1_violin.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gb1)
        
        dev.off()
        
        
        pdf(sprintf("%s\\%s_b0_pvalue.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gpvb0)
        
        dev.off()
        
        pdf(sprintf("%s\\%s_b1_pvalue.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gpvb1)
        
        dev.off()
        
        pdf(sprintf("%s\\%s_b0_violin_pvalue.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gpval_violin_b0)
        
        dev.off()
        
        pdf(sprintf("%s\\%s_b1_violin_pvalue.pdf",
                    figure_path,
                    size_name),
            height=size,
            width=size)
        
        plot(gpval_violin_b1)
        
        dev.off()
      }
  }
  
  write.csv(file=sprintf("%s\\statistics.csv", statistics_figure_path),
            statistics_df)
}

plot_TT_structure_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\TT_structure_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
    
  sizes = list(medium=c(width=1,
                          height=1),
               small=c(width=.75,
                       height=.75))
  
  for (dataset_idx in 1:len(paths)) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[dataset_idx]])
    
    dir.create(figure_path)
    
    reduced_mat <- get_mat_with_preset(paths[[dataset_idx]], "dalshtaim")  
    
    
    
    
    p_lem2_list <- pair_plots_colored(reduced_mat, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                      return_grid = F,
                                      pt_size=.55, 
                                      just_annotation = F)  
    
    
    p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                            rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                            return_grid = F,
                                            pt_size=.55, 
                                            just_annotation = T)  
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    
    
    for (size_name in names(sizes)) {
      

      
      p_lem2_no_annot <- p_lem2_list
      p_lem2_no_annot$nrow <- 3
      pf_structure <- do.call(plot_grid, p_lem2_no_annot)
      
      
      p_lem2_annot_only <- p_lem2_list_annot
      p_lem2_annot_only$nrow <- 3
      pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
      
      
      
      png(sprintf("%s\\%s_lem2_example_no_annot.png",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      plot(pf_structure)
      dev.off()
      
      pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5)
      
      plot(pf_annot_only)
      dev.off()
      
    }
  }
}

plot_HS_structure_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\HS_structure_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_hungry_sated_paths()
  datasets_names <- get_datasets_names(paths, sep="_")
  
  sizes = list(medium=c(width=1,
                        height=1),
               small=c(width=.75,
                       height=.75))
  
  for (dataset_idx in 1:len(paths)) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[dataset_idx]])
    
    dir.create(figure_path)
    
    reduced_mat <- get_mat_with_preset(paths[[dataset_idx]], "dalshtaim")  
    
    
    
    
    p_lem2_list <- pair_plots_colored(reduced_mat, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                      return_grid = F,
                                      pt_size=.55, 
                                      just_annotation = F)  
    
    
    p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                            rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                            return_grid = F,
                                            pt_size=.55, 
                                            just_annotation = T)  
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    
    
    for (size_name in names(sizes)) {
      
      
      
      p_lem2_no_annot <- p_lem2_list
      p_lem2_no_annot$nrow <- 3
      pf_structure <- do.call(plot_grid, p_lem2_no_annot)
      
      
      p_lem2_annot_only <- p_lem2_list_annot
      p_lem2_annot_only$nrow <- 3
      pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
      
      
      
      png(sprintf("%s\\%s_lem2_example_no_annot.png",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      plot(pf_structure)
      dev.off()
      
      pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5)
      
      plot(pf_annot_only)
      dev.off()
      
    }
  }
}

plot_control_structure_examples <- function()
{
  
  
  base_figure_path <- sprintf("%s\\control_structure_examples", output_path)
  dir.create(base_figure_path)
  
  paths <- get_control_paths()
  datasets_names <- get_datasets_names(paths, sep="_", control = T)
  
  sizes = list(medium=c(width=1,
                        height=1),
               small=c(width=.75,
                       height=.75))
  
  for (dataset_idx in 1:len(paths)) {
    
    figure_path <- sprintf("%s\\%s",
                           base_figure_path,
                           datasets_names[[dataset_idx]])
    
    dir.create(figure_path)
    
    reduced_mat <- get_mat_with_preset(paths[[dataset_idx]], "dalshtaim")  
    
    
    
    
    p_lem2_list <- pair_plots_colored(reduced_mat, 
                                      rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                      return_grid = F,
                                      pt_size=.55, 
                                      just_annotation = F)  
    
    
    p_lem2_list_annot <- pair_plots_colored(reduced_mat, 
                                            rep(adjustcolor("gray30", alpha=.3), times=nrow(reduced_mat)),
                                            return_grid = F,
                                            pt_size=.55, 
                                            just_annotation = T)  
    
    p_lem2_list <- lapply(1:len(p_lem2_list), 
                          function(i) {return(p_lem2_list[[i]] + 
                                                theme(line=element_blank(),
                                                      #rect=element_rect(color="white"),
                                                      plot.margin=margin(t = 8, r = 8, b = 8, l = 8, unit = "pt"),
                                                      panel.border=element_rect(color="white"),
                                                      plot.background = element_rect("white"),
                                                      axis.title=element_text(color="white")))})
    
    
    
    
    for (size_name in names(sizes)) {
      
      
      
      p_lem2_no_annot <- p_lem2_list
      p_lem2_no_annot$nrow <- 3
      pf_structure <- do.call(plot_grid, p_lem2_no_annot)
      
      
      p_lem2_annot_only <- p_lem2_list_annot
      p_lem2_annot_only$nrow <- 3
      pf_annot_only <- do.call(plot_grid, p_lem2_annot_only)
      
      
      
      png(sprintf("%s\\%s_lem2_example_no_annot.png",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5,
          unit="in",
          res=1500)
      
      plot(pf_structure)
      dev.off()
      
      pdf(sprintf("%s\\%s_lem2_example_no_annot.pdf",
                  figure_path,
                  size_name),
          height=sizes[[size_name]][["height"]] * 3,
          width=sizes[[size_name]][["width"]] * 5)
      
      plot(pf_annot_only)
      dev.off()
      
    }
  }
}


plot_TT_dimensionality_examples <- function(create=F, 
                                            epsilon_neighborhood = T)
{
  
  base_figure_path <- sprintf("%s\\TT_dimensionality_examples", output_path)
  dir.create(base_figure_path)
  
  
  paths <- get_thirsty_quenched_paths()
  
  
  
  if (epsilon_neighborhood) {
    output_name = "TT_calculated_dimensionality_linreg_eps_neighbor.Rda"
  } else {
    output_name = "TT_calculated_dimensionality_linreg.Rda"
  }
  
  max_x  <- c()
  plot_df <- data.frame() 
  if (create) {
    
  
  for (idx in 1:len(paths)) {
  
  #reduced_mat <- get_mat_with_preset(, "dalshtaim")
    reduced_mat <- get_reduced_mat_full_day(paths[[idx]], "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, time_shuffled=F)
  
  
  
  if (epsilon_neighborhood) { 
    dm <- as.matrix(dist(reduced_mat[sample(1:nrow(reduced_mat), 5000),]))
    eps <- seq(min(dm), max(dm), length.out=200)
    nn <- unlist(lapply(eps, function(ep) {sum(c(dm) <= ep)}))
    
    df <- data.frame(x=log10(eps),
                     y=log10(nn),
                     i=rep(idx, 200))
    max_x <- c(max_x, max(df$x))
  } else {
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
    
    ind <- sort(sample(1:len(linreg_x), 300))
    
    df <- data.frame(x=linreg_x[ind],
                  y=linreg_y[ind],
                  i=rep(idx, 300))
    max_x <- c(max_x, max(linreg_x))
  }
  
  
  plot_df <- rbind(plot_df, df)
  
  
}
    
    res_list <- list(plot_df=plot_df,
                     max_x=max_x)
    
    save(file=sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), res_list)
  } else {
    
    #load(sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), verbose=T)
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
    max_x = res_list$max_x
    plot_df = res_list$plot_df
  }
  
  #x_to_use <- seq(-5, max(max_x), by=.015)
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
  
  plot_slope = 6
  
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
      #next
      xt <- c(intercept_x, max_x)
      yt <- c(min_y, get_y_from_x(max_x))
      
      if (yt[2] > max_y) {
        yt[2] <- max_y
        xt[2]  <- get_x_from_y(max_y)
      }
      print(xt)
      print(yt)
    }
    

    # xt <- ifelse(intercept_x > 0, intercept_x, min(x_to_use))
    # xt <- c(xt, max(x_to_use))
    # yt <- plot_slope * xt + b
    # 
    # if (yt[1] <= max(plot_df$y)) {
    #   if (yt[2] >= max(plot_df$y)){
    #     
    #     yt[2] <- max(plot_df$y)
    #     xt[2] <- (max(plot_df$y) - b) / plot_slope
    #   }
    #   
      slope_lines_df <- data.frame(x=xt,
                                   y=yt)
      
      #print(slope_lines_df)
      gdim_all <- gdim_all +
        geom_line(data=slope_lines_df, aes(x=x,y=y), linetype="dashed", col="gray80")
      
    # } else 
  }

  #true_xt <- c(0, max_x)
  #true_yt <- c(0, max_x * slope)
  
  #true_slope_df <- data.frame(x=true_xt,
                              #y=true_yt)
  
  if (epsilon_neighborhood) {
    gdim_f <- 
      gdim_all +
      #geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
      geom_point(data=plot_df, aes(x=x, y=y, color=i), alpha=.2, stroke=0, size=2) +
      xlab(TeX("$log(neighbours)$")) + 
      ylab(TeX("$log(radius)$")) + 
      ggtitle("Dimensionality estimation") +
      big_text_base_plot_theme +
      theme(plot.title=element_text(size=10)) +
      scale_color_distiller(palette="Spectral")    
  } else {
  gdim_f <- 
    gdim_all +
    #geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
    geom_point(data=plot_df, aes(x=x, y=y, color=i), alpha=.5, stroke=0, size=2) +
    xlab(TeX("$log(NN_{distance})$")) + 
    ylab(TeX("$-log(1-NN_{quantile})$")) + 
    ggtitle("Dimensionality estimation") +
    big_text_base_plot_theme +
    theme(plot.title=element_text(size=10)) +
    scale_color_distiller(palette="Spectral")
  }
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s%s_dimensionality_plot.pdf",
                base_figure_path,
                size_name,
                ifelse(epsilon_neighborhood, "_eps", "")),
        height=size,
        width=size)
    
    plot(gdim_f)
    
    dev.off()
    
  }
  
}


plot_HS_dimensionality_examples <- function(create=F, 
                                            epsilon_neighborhood = T)
{
  
  base_figure_path <- sprintf("%s\\HS_dimensionality_examples", output_path)
  dir.create(base_figure_path)
  
  
  paths <- get_hungry_sated_paths()
  
  
  plot_slope = 6
  if (epsilon_neighborhood) {
    output_name = "HS_calculated_dimensionality_linreg_eps_neighbor.Rda"
  } else {
    output_name = "HS_calculated_dimensionality_linreg.Rda"
  }
  
  max_x  <- c()
  plot_df <- data.frame() 
  if (create) {
    
    
    for (idx in 1:len(paths)) {
      
      #reduced_mat <- get_mat_with_preset(, "dalshtaim")
      reduced_mat <- get_reduced_mat_full_day(paths[[idx]], "lem", ndim=20, activity_threshold = .5, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, time_shuffled=F)
      
      
      
      if (epsilon_neighborhood) { 
        dm <- as.matrix(dist(reduced_mat[sample(1:nrow(reduced_mat), 5000),]))
        eps <- seq(min(dm), max(dm), length.out=200)
        nn <- unlist(lapply(eps, function(ep) {sum(c(dm) <= ep)}))
        
        df <- data.frame(x=log10(eps),
                         y=log10(nn),
                         i=rep(idx, 200))
        max_x <- c(max_x, max(df$x))
      } else {
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
        
        ind <- sort(sample(1:len(linreg_x), 300))
        
        df <- data.frame(x=linreg_x[ind],
                         y=linreg_y[ind],
                         i=rep(idx, 300))
        max_x <- c(max_x, max(linreg_x))
      }
      
      
      plot_df <- rbind(plot_df, df)
      
      
    }
    
    res_list <- list(plot_df=plot_df,
                     max_x=max_x)
    
    save(file=sprintf("%s\\data\\figure_2\\%s", base_output_path, output_name), res_list)
  } else {
    
    #load(sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), verbose=T)
    load(sprintf("%s\\data\\figure_2\\%s", base_output_path, output_name), verbose=T)
    max_x = res_list$max_x
    plot_df = res_list$plot_df
  }
  
  #x_to_use <- seq(-5, max(max_x), by=.015)
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
  
  plot_slope = 6
  
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
      #next
      xt <- c(intercept_x, max_x)
      yt <- c(min_y, get_y_from_x(max_x))
      
      if (yt[2] > max_y) {
        yt[2] <- max_y
        xt[2]  <- get_x_from_y(max_y)
      }
      print(xt)
      print(yt)
    }
    
    
    # xt <- ifelse(intercept_x > 0, intercept_x, min(x_to_use))
    # xt <- c(xt, max(x_to_use))
    # yt <- plot_slope * xt + b
    # 
    # if (yt[1] <= max(plot_df$y)) {
    #   if (yt[2] >= max(plot_df$y)){
    #     
    #     yt[2] <- max(plot_df$y)
    #     xt[2] <- (max(plot_df$y) - b) / plot_slope
    #   }
    #   
    slope_lines_df <- data.frame(x=xt,
                                 y=yt)
    
    #print(slope_lines_df)
    gdim_all <- gdim_all +
      geom_line(data=slope_lines_df, aes(x=x,y=y), linetype="dashed", col="gray80")
    
    # } else 
  }
  
  #true_xt <- c(0, max_x)
  #true_yt <- c(0, max_x * slope)
  
  #true_slope_df <- data.frame(x=true_xt,
  #y=true_yt)
  
  if (epsilon_neighborhood) {
    gdim_f <- 
      gdim_all +
      #geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
      geom_point(data=plot_df, aes(x=x, y=y, color=i), alpha=.2, stroke=0, size=2) +
      xlab(TeX("$log(neighbours)$")) + 
      ylab(TeX("$log(radius)$")) + 
      ggtitle("Dimensionality estimation") +
      big_text_base_plot_theme +
      theme(plot.title=element_text(size=10)) +
      scale_color_distiller(palette="Spectral")    
  } else {
    gdim_f <- 
      gdim_all +
      #geom_line(data=true_slope_df, aes(x=x,y=y), col="#F05A28", size=1) +
      geom_point(data=plot_df, aes(x=x, y=y, color=i), alpha=.5, stroke=0, size=2) +
      xlab(TeX("$log(NN_{distance})$")) + 
      ylab(TeX("$-log(1-NN_{quantile})$")) + 
      ggtitle("Dimensionality estimation") +
      big_text_base_plot_theme +
      theme(plot.title=element_text(size=10)) +
      scale_color_distiller(palette="Spectral")
  }
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s%s_dimensionality_plot.pdf",
                base_figure_path,
                size_name,
                ifelse(epsilon_neighborhood, "_eps", "")),
        height=size,
        width=size)
    
    plot(gdim_f)
    
    dev.off()
    
  }
  
}

plot_TT_dimension_boxplots <- function(create=F) 
{
  
  base_figure_path <- sprintf("%s\\TT_dimension_boxplots", output_path)
  dir.create(base_figure_path)

  statistics_df <- data.frame()
  
  output_name = "TT_calculated_dimensionality_estimation_full.Rda"
  paths <- get_thirsty_quenched_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  if (create) {
    dim_df <- data.frame()
    for (path_idx in 1:len(paths)) {
      
      p <- paths[[path_idx]]
      original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, just_original_mat = T, control = F, override=F)
      shuffled_original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, just_original_mat = T, control = F, shuffled = T, time_shuffled = F, override=F)
      reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, control=T, shuffled=F, time_shuffled=F)
      shuffled_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, control=T, shuffled=T, time_shuffled=F)
      #reg_id <- get_intrinsic_dim(reg_mat)
      
      shuff_org_id <- twonn(shuffled_original_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      org_id <- twonn(original_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      reg_id <- twonn(reg_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      shuffled_id <- twonn(shuffled_mat, c_trimmed = 0.04, method="linfit")
    
      print(sprintf("Original %.3f, Shuffled %.3f, Red %.3f, Shuffled %.3f",
               org_id$est[[2]],
               shuff_org_id$est[[2]],
               reg_id$est[[2]],
               shuffled_id$est[[2]]))
      
      dim_df <- rbind(dim_df,
                      data.frame(name=datasets_names[path_idx],
                                 n_neurons=min(dim(original_mat)),
                                 n_timepoints=max(dim(original_mat)),
                                 original_dim=org_id$est[[2]],
                                 shuffled_original_dim=shuff_org_id$est[[2]],
                                 reg_dim=reg_id$est[[2]],
                                 shuffled_dim=shuffled_id$est[[2]]))
    }
    
    save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), dim_df)
  } else {
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  }
  
  
  tmp_df <- dim_df[,c("original_dim", "shuffled_original_dim", "reg_dim","shuffled_dim")]
  colnames(tmp_df) <- c("Data", "Shuffle", "Data (denoised)", "Shuffle (denoised")
  
  melted_df <- melt(tmp_df)
  colnames(melted_df) <- c("x", "y")
  
  wilc <- wilcox.test(tmp_df[,"Data"], tmp_df[,"Shuffle"],
                      correct = F,
                      alternative = "less", 
                      paired = T)
  
  statistics_df <- rbind(statistics_df,
                         data.frame(statistic=wilc$statistic,
                                    pval=wilc$p.value,
                                    method=wilc$method,
                                    alternative=wilc$alternative,
                                    mean_group_reg=mean(tmp_df[,"Data"]),
                                    mean_group_shuffle=mean(tmp_df[,"Shuffle"]),
                                    sd_group_reg=sd(tmp_df[,"Data"]),
                                    sd_group_shuffle=sd(tmp_df[,"Shuffle"]),
                                    sem_group_reg=sem(tmp_df[,"Data"]),
                                    sem_group_shuffle=sem(tmp_df[,"Shuffle"]),
                                    R=NA,
                                    signif=signif.num(wilc$p.value)))
                         
  crtest <- cor.test(dim_df[,"reg_dim"], dim_df[,"n_neurons"])
  
  statistics_df <- rbind(statistics_df,
                         data.frame(statistic=crtest$statistic,
                                    pval=crtest$p.value,
                                    method=crtest$method,
                                    alternative=crtest$alternative,
                                    mean_group_reg=mean(tmp_df[,"Data"]),
                                    mean_group_shuffle=mean(dim_df[,"n_neurons"]),
                                    sd_group_reg=sd(tmp_df[,"Data"]),
                                    sd_group_shuffle=sd(dim_df[,"n_neurons"]),
                                    sem_group_reg=sem(tmp_df[,"Data"]),
                                    sem_group_shuffle=sem(dim_df[,"n_neurons"]),
                                    R=cor(dim_df[,"reg_dim"], dim_df[,"n_neurons"]),
                                    signif=signif.num(crtest$p.value)))
  
  gbar <- 
    ggplot(melted_df, aes(x=x,y=y)) +
    # geom_bar(stat="summary", width=.45,
    #          fill="gray70") +
    # geom_errorbar(stat="summary", width = .3) +
    # geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
    #             size=1.5) + 
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75, binwidth = 1/2)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape=NA) +
    xlab("") +
    ylab("Estimated dimension") +
    big_text_base_plot_theme + 
    #ylim(0,30) +
    # geom_line(data=data.frame(x=c(1,1,2,2),
    #                           y=c(24.5,25,25,24.5)),
    #           aes(x=x,y=y)) +
    # 
    # geom_point(data=data.frame(x=1.5,y=26.5),
    #            color="white",fill="white",stroke=0) + 
    # geom_text(label=signif.num(wilcox.test(dim_df[,"reg_dim"], dim_df[,"shuffled_dim"], alternative="less", paired=T)$p.value),
    #           x=1.5,
    #           y=25.5) +
    scale_color_manual(values=c("#405EAB", "#CE3736"))
    #scale_y_continuous(expand=c(0,0)) 
  
  tmp_df_2 <- dim_df[,c("n_neurons", "reg_dim", "shuffled_dim")]
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
    big_text_base_plot_theme + 
    ylim(0,30)  +
    xlim(150,430) +
    scale_color_manual(values=c("#405EAB", "#CE3736")) + 
    geom_hline(yintercept=mean(dim_df$reg_dim), linetype="dashed")
  
  
  gf <- plot_grid(gbar, gdim, nrow=1)
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_dimensionality_example.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=2 * size)
    
    plot(gf)
    
    dev.off()
    
  }
  
  write.csv(file=sprintf("%s\\statistics.csv", base_figure_path),
            statistics_df)
}

plot_HS_TT_dimension_boxplots <- function(create=F) 
{
  
  base_figure_path <- sprintf("%s\\HS_dimension_boxplots", output_path)
  dir.create(base_figure_path)
  
  statistics_df <- data.frame()
  statistics_df_cor <- data.frame()
  
  output_name = "HS_calculated_dimensionality_estimation.Rda"
  tt_output_name = "TT_calculated_dimensionality_estimation.Rda"
  
  paths <- get_hungry_sated_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  if (create) {
    dim_df <- data.frame()
    for (path_idx in 1:len(paths)) {
      
      activity_thresh <- ifelse(path_idx >=9, .25, .5)
      p <- paths[[path_idx]]
      original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, activity_threshold = activity_thresh, just_original_mat = T)
      reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, activity_threshold = activity_thresh, time_shuffled=F)
      shuffled_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, activity_threshold = activity_thresh, time_shuffled=F)
      #reg_id <- get_intrinsic_dim(reg_mat)
      
      reg_id <- twonn(reg_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      shuffled_id <- twonn(shuffled_mat, c_trimmed = 0.04, method="linfit")
      
      print(reg_id$est[[2]])
      
      dim_df <- rbind(dim_df,
                      data.frame(name=datasets_names[path_idx],
                                 n_neurons=min(dim(original_mat)),
                                 n_timepoints=max(dim(original_mat)),
                                 reg_dim=reg_id$est[[2]],
                                 shuffled_dim=shuffled_id$est[[2]]))
    }
    
    dir.create(sprintf("%s\\data\\figure_1\\", base_output_path))
    
    save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), dim_df)
  } else {
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
    hs_dim_df <- dim_df
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, tt_output_name), verbose=T)
    tt_dim_df <- dim_df
  }
  
  
  tmp_df <- data.frame(dim=c(hs_dim_df[,"reg_dim"], tt_dim_df[,"reg_dim"]),
                       group=rep(c("HS", "TT"), times=c(nrow(hs_dim_df), nrow(tt_dim_df))))
  
  colnames(tmp_df) <- c("y", 
                        "x")
  

  
  wilc <- wilcox.test(hs_dim_df[,"reg_dim"], tt_dim_df[,"reg_dim"],
                      correct = F,
                      alternative = "two.sided", 
                      paired = F)
  
  statistics_df <- rbind(statistics_df,
                         data.frame(statistic=wilc$statistic,
                                    pval=wilc$p.value,
                                    method=wilc$method,
                                    alternative=wilc$alternative,
                                    mean_group_HS=mean(hs_dim_df[,"reg_dim"]),
                                    mean_group_TT=mean(tt_dim_df[,"reg_dim"]),
                                    sd_group_HS=sd(hs_dim_df[,"reg_dim"]),
                                    sd_group_TT=sd(tt_dim_df[,"reg_dim"]),
                                    sem_group_HS=sem(hs_dim_df[,"reg_dim"]),
                                    sem_group_TT=sem(tt_dim_df[,"reg_dim"]),
                                    R=NA,
                                    signif=signif.num(wilc$p.value)))
  
  crtest <- cor.test(hs_dim_df[,"reg_dim"], hs_dim_df[,"n_neurons"])
  
  statistics_df_cor <- rbind(statistics_df_cor,
                         data.frame(statistic=crtest$statistic,
                                    pval=crtest$p.value,
                                    method=crtest$method,
                                    alternative=crtest$alternative,
                                    mean_group_hs=mean(hs_dim_df[,"reg_dim"]),
                                    mean_group_hsneur=mean(hs_dim_df[,"n_neurons"]),
                                    sd_group_hs=sd(hs_dim_df[,"reg_dim"]),
                                    sd_group_hsneur=sd(hs_dim_df[,"n_neurons"]),
                                    sem_group_hs=sem(hs_dim_df[,"reg_dim"]),
                                    sem_group_hsneur=sem(hs_dim_df[,"n_neurons"]),
                                    R=cor(hs_dim_df[,"reg_dim"], hs_dim_df[,"n_neurons"]),
                                    signif=signif.num(crtest$p.value)))
  
  gbar <- 
    ggplot(tmp_df, aes(x=x,y=y)) +
    # geom_bar(stat="summary", width=.45,
    #          fill="gray70") +
    # geom_errorbar(stat="summary", width = .3) +
    # geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
    #             size=1.5) + 
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75, binwidth = 1/2)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape=NA) +
    xlab("") +
    ylab("Estimated dimension") +
    big_text_base_plot_theme + 
    ylim(0,13) +
    geom_hline(yintercept=mean(hs_dim_df$shuffled_dim), col="purple", linetype="dashed") +
    geom_hline(yintercept=mean(tt_dim_df$shuffled_dim), col="orange", linetype="dashed")
    # geom_line(data=data.frame(x=c(1,1,2,2),
    #                           y=c(24.5,25,25,24.5)),
    #           aes(x=x,y=y)) +
    # 
    # geom_point(data=data.frame(x=1.5,y=26.5),
    #            color="white",fill="white",stroke=0) + 
    # geom_text(label=signif.num(wilcox.test(dim_df[,"reg_dim"], dim_df[,"shuffled_dim"], alternative="less", paired=T)$p.value),
    #           x=1.5,
    #           y=25.5) +
  #   scale_color_manual(values=c("#405EAB", "#CE3736"))
  # #scale_y_continuous(expand=c(0,0)) 
  # 
  # tmp_df_2 <- dim_df[,c("n_neurons", "reg_dim", "shuffled_dim")]
  # colnames(tmp_df_2) <- c("Neurons", "Regular", "Shuffle")
  # 
  # melted_df_2 <- melt(tmp_df_2, id.vars = "Neurons")
  # colnames(melted_df_2)[2:3] <- c("group", "Dim")
  # 
  # gdim <- 
  #   ggplot(melted_df_2) +
  #   geom_point(aes(y=Dim, x=Neurons, color=group), alpha=.75,
  #              size=3.5,
  #              stroke=0,
  #              fill=NA) + 
  #   theme_classic() +
  #   xlab("Number of neurons") +
  #   ylab("Estimated dimension") +
  #   big_text_base_plot_theme + 
  #   ylim(0,30)  +
  #   xlim(150,430) +
  #   scale_color_manual(values=c("#405EAB", "#CE3736")) + 
  #   geom_hline(yintercept=mean(dim_df$reg_dim), linetype="dashed")
  # 
  
  #gf <- plot_grid(gbar, gdim, nrow=1)
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%s_hs_vs_tt_dimensionality_example_no_shuffle.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gbar)
    
    dev.off()
    
  }
  
  write.csv(file=sprintf("%s\\hs_tt_statistics.csv", base_figure_path),
            statistics_df)
  write.csv(file=sprintf("%s\\hs_tt_statistics_cor.csv", base_figure_path),
            statistics_df_cor)
  
}

plot_HS_v1_dimension_boxplots <- function(create=F) 
{
  
  base_figure_path <- sprintf("%s\\HS_control_dimension_boxplots", output_path)
  dir.create(base_figure_path)
  
  statistics_df <- data.frame()
  statistics_df_cor <- data.frame()
  
  
  output_name = "v1_calculated_dimensionality_estimation.Rda"
  hs_output_name = "HS_calculated_dimensionality_estimation.Rda"
  
  
  if (create) {
    dim_df <- data.frame()
    for (path_idx in 1:len(get_control_paths())) {
      
      
      p <- paths[[path_idx]]
      original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, activity_threshold = .25, just_original_mat = T, control=T)
      reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, activity_threshold = .25, time_shuffled=F, control = T)
      shuffled_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, activity_threshold = .25, time_shuffled=F, control=T)
      #reg_id <- get_intrinsic_dim(reg_mat)
      
      reg_id <- twonn(reg_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      shuffled_id <- twonn(shuffled_mat, c_trimmed = 0.04, method="linfit")
      
      print(reg_id$est[[2]])
      
      dim_df <- rbind(dim_df,
                      data.frame(name=datasets_names[path_idx],
                                 n_neurons=min(dim(original_mat)),
                                 n_timepoints=max(dim(original_mat)),
                                 reg_dim=reg_id$est[[2]],
                                 shuffled_dim=shuffled_id$est[[2]]))
    }
    
    dir.create(sprintf("%s\\data\\figure_1\\", base_output_path))
    save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), dim_df)
  } else {
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, hs_output_name), verbose=T)
    hs_dim_df <- dim_df
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
    v1_dim_df <- dim_df
  }
  
  
  tmp_df <- data.frame(dim=c(hs_dim_df[,"reg_dim"], v1_dim_df[,"reg_dim"]),
                       group=rep(c("InsCtx", "Visual areas"), times=c(nrow(hs_dim_df), nrow(v1_dim_df))))
  
  colnames(tmp_df) <- c("y", 
                        "x")
  
  
  
  wilc <- wilcox.test(hs_dim_df[,"reg_dim"], v1_dim_df[,"reg_dim"],
                      correct = F,
                      alternative = "less", 
                      paired = F)
  
  statistics_df <- rbind(statistics_df,
                         data.frame(statistic=wilc$statistic,
                                    pval=wilc$p.value,
                                    method=wilc$method,
                                    alternative=wilc$alternative,
                                    N_group_HS=len(hs_dim_df[,"reg_dim"]),
                                    N_group_V1=len(v1_dim_df[,"reg_dim"]),
                                    mean_group_HS=mean(hs_dim_df[,"reg_dim"]),
                                    mean_group_V1=mean(v1_dim_df[,"reg_dim"]),
                                    sd_group_HS=sd(hs_dim_df[,"reg_dim"]),
                                    sd_group_V1=sd(v1_dim_df[,"reg_dim"]),
                                    sem_group_HS=sem(hs_dim_df[,"reg_dim"]),
                                    sem_group_V1=sem(v1_dim_df[,"reg_dim"]),
                                    R=NA,
                                    signif=signif.num(wilc$p.value)))

  
  gbar <- 
    ggplot(tmp_df, aes(x=x,y=y)) +
    # geom_bar(stat="summary", width=.45,
    #          fill="gray70") +
    # geom_errorbar(stat="summary", width = .3) +
    # geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
    #             size=1.5) + 
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75, binwidth = 1/2)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape=NA) +
    xlab("") +
    ylab("Estimated dimension") +
    big_text_base_plot_theme + 
    ylim(0,13) +
    geom_hline(yintercept=mean(hs_dim_df$shuffled_dim), col="purple", linetype="dashed") +
    geom_hline(yintercept=mean(v1_dim_df$shuffled_dim), col="orange", linetype="dashed")
  # geom_line(data=data.frame(x=c(1,1,2,2),
  #                           y=c(24.5,25,25,24.5)),
  #           aes(x=x,y=y)) +
  # 
  # geom_point(data=data.frame(x=1.5,y=26.5),
  #            color="white",fill="white",stroke=0) + 
  # geom_text(label=signif.num(wilcox.test(dim_df[,"reg_dim"], dim_df[,"shuffled_dim"], alternative="less", paired=T)$p.value),
  #           x=1.5,
  #           y=25.5) +
  #   scale_color_manual(values=c("#405EAB", "#CE3736"))
  # #scale_y_continuous(expand=c(0,0)) 
  # 
   tmp_df_2 <- hs_dim_df[,c("n_neurons", "reg_dim", "shuffled_dim")]
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
    big_text_base_plot_theme +
    ylim(0,30)  +
    xlim(min(melted_df_2$Neurons),max(melted_df_2$Neurons)) +
    scale_color_manual(values=c("#405EAB", "#CE3736")) +
    geom_hline(yintercept=mean(hs_dim_df$reg_dim), linetype="dashed")

  
  #gf <- plot_grid(gbar, gdim, nrow=1)
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    # pdf(sprintf("%s\\%s_hs_vs_v1_dimensionality_example_no_shuffle.pdf",
    #             base_figure_path,
    #             size_name),
    #     height=size,
    #     width=size)
    # 
    # plot(gbar)
    # 
    # dev.off()
    
    pdf(sprintf("%s\\%s_hs_dim.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gdim)
    
    dev.off()
    
    
  }
  
  write.csv(file=sprintf("%s\\hs_v1_statistics.csv", base_figure_path),
            statistics_df)
}

plot_HS_TT_structure <- function(create=F) 
{
  
  base_figure_path <- sprintf("%s\\HS_dimension_boxplots", output_path)
  dir.create(base_figure_path)
  
  statistics_df <- data.frame()
  statistics_df_cor <- data.frame()
  
  output_name = "HS_calculated_dimensionality_estimation.Rda"
  tt_output_name = "TT_calculated_dimensionality_estimation.Rda"
  
  paths <- get_hungry_sated_paths()
  datasets_names <- get_datasets_names(paths, sep = "_")
  
  if (create) {
    dim_df <- data.frame()
    for (path_idx in 1:len(paths)) {
      
      activity_thresh <- ifelse(path_idx >=9, .25, .5)
      p <- paths[[path_idx]]
      original_mat <- get_reduced_mat_full_day(p, "lem", window_size=15, activity_threshold = activity_thresh, just_original_mat = T)
      reg_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=F, activity_threshold = activity_thresh, time_shuffled=F)
      shuffled_mat <- get_reduced_mat_full_day(p, "lem", ndim=20, window_size=15, knn1=0.075, knn2 = 0, shuffled=T, activity_threshold = activity_thresh, time_shuffled=F)
      #reg_id <- get_intrinsic_dim(reg_mat)
      
      reg_id <- twonn(reg_mat, c_trimmed = 0.04, method="linfit") # for some reason their method is faster!!! probably all the cpp stuff :[
      shuffled_id <- twonn(shuffled_mat, c_trimmed = 0.05, method="linfit")
      
      print(reg_id$est[[2]])
      
      dim_df <- rbind(dim_df,
                      data.frame(name=datasets_names[path_idx],
                                 n_neurons=min(dim(original_mat)),
                                 n_timepoints=max(dim(original_mat)),
                                 reg_dim=reg_id$est[[2]],
                                 shuffled_dim=shuffled_id$est[[2]]))
    }
    
    dir.create(sprintf("%s\\data\\figure_1\\", base_output_path))
    
    save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), dim_df)
  } else {
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
    hs_dim_df <- dim_df
    load(sprintf("%s\\data\\figure_1\\%s", base_output_path, tt_output_name), verbose=T)
    tt_dim_df <- dim_df
  }
  
  
  tmp_df <- data.frame(dim=c(hs_dim_df[,"reg_dim"], tt_dim_df[,"reg_dim"]),
                       group=rep(c("HS", "TT"), times=c(nrow(hs_dim_df), nrow(tt_dim_df))))
  
  colnames(tmp_df) <- c("y", 
                        "x")
  
  
  
  wilc <- wilcox.test(hs_dim_df[,"reg_dim"], tt_dim_df[,"reg_dim"],
                      correct = F,
                      alternative = "two.sided", 
                      paired = F)
  
  statistics_df <- rbind(statistics_df,
                         data.frame(statistic=wilc$statistic,
                                    pval=wilc$p.value,
                                    method=wilc$method,
                                    alternative=wilc$alternative,
                                    mean_group_HS=mean(hs_dim_df[,"reg_dim"]),
                                    mean_group_TT=mean(tt_dim_df[,"reg_dim"]),
                                    sd_group_HS=sd(hs_dim_df[,"reg_dim"]),
                                    sd_group_TT=sd(tt_dim_df[,"reg_dim"]),
                                    sem_group_HS=sem(hs_dim_df[,"reg_dim"]),
                                    sem_group_TT=sem(tt_dim_df[,"reg_dim"]),
                                    R=NA,
                                    signif=signif.num(wilc$p.value)))
  
  crtest <- cor.test(hs_dim_df[,"reg_dim"], hs_dim_df[,"n_neurons"])
  
  statistics_df_cor <- rbind(statistics_df_cor,
                             data.frame(statistic=crtest$statistic,
                                        pval=crtest$p.value,
                                        method=crtest$method,
                                        alternative=crtest$alternative,
                                        mean_group_hs=mean(hs_dim_df[,"reg_dim"]),
                                        mean_group_hsneur=mean(hs_dim_df[,"n_neurons"]),
                                        sd_group_hs=sd(hs_dim_df[,"reg_dim"]),
                                        sd_group_hsneur=sd(hs_dim_df[,"n_neurons"]),
                                        sem_group_hs=sem(hs_dim_df[,"reg_dim"]),
                                        sem_group_hsneur=sem(hs_dim_df[,"n_neurons"]),
                                        R=cor(hs_dim_df[,"reg_dim"], hs_dim_df[,"n_neurons"]),
                                        signif=signif.num(crtest$p.value)))
  
  gbar <- 
    ggplot(tmp_df, aes(x=x,y=y)) +
    # geom_bar(stat="summary", width=.45,
    #          fill="gray70") +
    # geom_errorbar(stat="summary", width = .3) +
    # geom_jitter(position=position_jitterdodge(.5), aes(fill=x),
    #             size=1.5) + 
    geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75, binwidth = 1/2)+
    scale_shape_manual(values=18) +
    geom_boxplot(fill=NA, width=.5, outlier.shape=NA) +
    xlab("") +
    ylab("Estimated dimension") +
    big_text_base_plot_theme + 
    ylim(0,30) +
    geom_hline(yintercept=mean(hs_dim_df$shuffled_dim), col="purple", linetype="dashed") +
    geom_hline(yintercept=mean(tt_dim_df$shuffled_dim), col="orange", linetype="dashed")
  # geom_line(data=data.frame(x=c(1,1,2,2),
  #                           y=c(24.5,25,25,24.5)),
  #           aes(x=x,y=y)) +
  # 
  # geom_point(data=data.frame(x=1.5,y=26.5),
  #            color="white",fill="white",stroke=0) + 
  # geom_text(label=signif.num(wilcox.test(dim_df[,"reg_dim"], dim_df[,"shuffled_dim"], alternative="less", paired=T)$p.value),
  #           x=1.5,
  #           y=25.5) +
  #   scale_color_manual(values=c("#405EAB", "#CE3736"))
  # #scale_y_continuous(expand=c(0,0)) 
  # 
  # tmp_df_2 <- dim_df[,c("n_neurons", "reg_dim", "shuffled_dim")]
  # colnames(tmp_df_2) <- c("Neurons", "Regular", "Shuffle")
  # 
  # melted_df_2 <- melt(tmp_df_2, id.vars = "Neurons")
  # colnames(melted_df_2)[2:3] <- c("group", "Dim")
  # 
  # gdim <- 
  #   ggplot(melted_df_2) +
  #   geom_point(aes(y=Dim, x=Neurons, color=group), alpha=.75,
  #              size=3.5,
  #              stroke=0,
  #              fill=NA) + 
  #   theme_classic() +
  #   xlab("Number of neurons") +
  #   ylab("Estimated dimension") +
  #   big_text_base_plot_theme + 
  #   ylim(0,30)  +
  #   xlim(150,430) +
  #   scale_color_manual(values=c("#405EAB", "#CE3736")) + 
  #   geom_hline(yintercept=mean(dim_df$reg_dim), linetype="dashed")
  # 
  
  #gf <- plot_grid(gbar, gdim, nrow=1)
  
  for (size_name in names(a4_sizes)) {
    size = a4_sizes[[size_name]]
    pdf(sprintf("%s\\%shs_vs_tt_dimensionality_example.pdf",
                base_figure_path,
                size_name),
        height=size,
        width=size)
    
    plot(gbar)
    
    dev.off()
    
  }
  
  write.csv(file=sprintf("%s\\hs_tt_statistics.csv", base_figure_path),
            statistics_df)
}

average_across_datastes_across_within  <- function(df) 
{
  
  return(
    ddply(df,
          .(path1),
          function (path_df) {
            
            across <- path_df[path_df$path1 != path_df$path2]
            within <- path_df[path_df$path1 == path_df$path2]
            across_df <- 
              ddply(across, .(metric),
                    function(metric_df) {
                      return(
                        c(mean(as.numeric(metric_df$betti0)),
                          mean(as.numeric(metric_df$betti1))))
                    })
            
            within_ <- 
              ddply(across, .(metric),
                    function(metric_df) {
                      return(
                        c(mean(as.numeric(metric_df$betti0)),
                          mean(as.numeric(metric_df$betti1))))
                    })
          }))
}

average_similarity_across_datastes  <- function(df, path1=T)
{
  
  if (path1) {
  return(
  ddply(df,
        .(path1),
        function (path_df) {
          return(
            ddply(path_df, .(metric),
                  function(metric_df) {
                    return(
                      c(betti0=mean(as.numeric(metric_df$betti0)),
                        betti1=mean(as.numeric(metric_df$betti1))))
                  }))
        }))
  } 
  
  return(
    ddply(df,
          .(path2),
          function (path_df) {
            return(
              ddply(path_df, .(metric),
                    function(metric_df) {
                      return(
                        c(betti0=mean(as.numeric(metric_df$betti0)),
                          betti1=mean(as.numeric(metric_df$betti1))))
                    }))
          }))
}

generic_topo_similarity_statistics <- function(df_final, metric, isPaired=T, wilc_alt_to_use="less", ks_alt_to_use="greater")
{
  
  
  statistics_df <- data.frame()
  
  df_group_A <- df_final[df_final$group == "Group_A",]
  df_group_B <- df_final[df_final$group == "Group_B",]

  N_comparisions = 1     
  # Statistical tests
  wilcox_b0 <- 
    wilcox.test(df_group_A$betti0,
                df_group_B$betti0,
                paired=isPaired,
                correct = F,
                alternative = wilc_alt_to_use)
  
  
  wilcox_b1 <- 
    wilcox.test(df_group_A$betti1,
                df_group_B$betti1,
                paired=isPaired,
                correct = F,
                alternative = wilc_alt_to_use)
  
  
  ks_b0 <- 
    ks.test(df_group_A$betti0,
            df_group_B$betti0,
            exact = T,
            alternative = ks_alt_to_use)
  
  
  ks_b1 <- 
    ks.test(df_group_A$betti1,
            df_group_B$betti1,
            exact = T,
            alternative = ks_alt_to_use)
  
  stat_df <-
    data.frame(comparision="wilcox_b0",
               method=wilcox_b0$method,
               alternative=wilcox_b0$alternative,
               pval=wilcox_b0$p.value,
               adjusted.pval=wilcox_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_group_A$betti0),
               mean_group_b=mean(df_group_B$betti0),
               sd_group_a=sd(df_group_A$betti0),
               sd_group_b=sd(df_group_B$betti0),
               sem_group_a=sem(df_group_A$betti0),
               sem_group_b=sem(df_group_B$betti0),
               N_group_a=len(df_group_A$betti0),
               N_group_b=len(df_group_B$betti0))
  
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="wilcox_b1",
               method=wilcox_b1$method,
               alternative=wilcox_b1$alternative,
               pval=wilcox_b1$p.value,
               adjusted.pval=wilcox_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_group_A$betti1),
               mean_group_b=mean(df_group_B$betti1),
               sd_group_a=sd(df_group_A$betti1),
               sd_group_b=sd(df_group_B$betti1),
               sem_group_a=sem(df_group_A$betti1),
               sem_group_b=sem(df_group_B$betti1),
               N_group_a=len(df_group_A$betti1),
               N_group_b=len(df_group_B$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_b0",
               method=ks_b0$method,
               alternative=ks_b0$alternative,
               pval=ks_b0$p.value,
               adjusted.pval=ks_b0$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_group_A$betti0),
               mean_group_b=mean(df_group_B$betti0),
               sd_group_a=sd(df_group_A$betti0),
               sd_group_b=sd(df_group_B$betti0),
               sem_group_a=sem(df_group_A$betti0),
               sem_group_b=sem(df_group_B$betti0),
               N_group_a=len(df_group_A$betti0),
               N_group_b=len(df_group_B$betti0))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  stat_df <-
    data.frame(comparision="ks_b1",
               method=ks_b1$method,
               alternative=ks_b1$alternative,
               pval=ks_b1$p.value,
               adjusted.pval=ks_b1$p.value * N_comparisions,
               metric=metric,
               mean_group_a=mean(df_group_A$betti1),
               mean_group_b=mean(df_group_B$betti1),
               sd_group_a=sd(df_group_A$betti1),
               sd_group_b=sd(df_group_B$betti1),
               sem_group_a=sem(df_group_A$betti1),
               sem_group_b=sem(df_group_B$betti1),
               N_group_a=len(df_group_A$betti1),
               N_group_b=len(df_group_B$betti1))
  
  statistics_df <- rbind(statistics_df, 
                         stat_df)
  
  
  return(statistics_df)
}

plot_TT_TH_HH_HC_shuffle_topo_similarity <- function()
{
  
  base_figure_path <- sprintf("%s\\TT_TH_HH_HC_shuffle_topo_similarity", output_path)
  dir.create(base_figure_path)
  
  statistics_figure_path <- sprintf("%s\\statistics", base_figure_path)
  dir.create(statistics_figure_path)
  
  
  output_name = "thirst_thirst_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  tt <- final_res
  
  output_name = "hunger_thirst_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  th <- final_res
  

  output_name = "hunger_control_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose = T)
  hc <- final_res
  

  output_name = "hunger_hunger_topological_comparision.Rda"
  load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  hh <- final_res
  
  
  tt$topological_similarity <-  average_similarity_across_datastes(tt$topological_similarity)
  th$topological_similarity <-  average_similarity_across_datastes(th$topological_similarity)
  hc$topological_similarity <-  average_similarity_across_datastes(hc$topological_similarity, path1=F)
  colnames(hc$topological_similarity)[1] <- "path1"
  hh$topological_similarity <-  average_similarity_across_datastes(hh$topological_similarity)
  
  statistics_df <- data.frame()
  values_df_all <- data.frame()
  
  
  
  for (metric in c("JSD", "KL", "wasserstein")) {

    
    figure_path <- sprintf("%s\\%s", base_figure_path, metric)
  
    dir.create(figure_path)
    
    df_tt_th <- rbind(tt$topological_similarity[tt$topological_similarity$metric == metric,],
                      th$topological_similarity[th$topological_similarity$metric == metric,])
    
    df_hh_hc <- rbind(hh$topological_similarity[hh$topological_similarity$metric == metric,],
                      hc$topological_similarity[hc$topological_similarity$metric == metric,])
    
    cmps_tt_th <- c(nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]),
                    nrow(th$topological_similarity[th$topological_similarity$metric == metric,]))
    
    cmps_hh_hc <- c(nrow(hh$topological_similarity[hh$topological_similarity$metric == metric,]),
                    nrow(hc$topological_similarity[hc$topological_similarity$metric == metric,]))
    
    df_tt_th$group = rep(c("Group_A", "Group_B"), cmps_tt_th)
    df_tt_th <- as.data.frame(df_tt_th)
    df_tt_th$betti0 <- as.numeric(df_tt_th$betti0)
    df_tt_th$betti1 <- as.numeric(df_tt_th$betti1)
    
    df_hh_hc$group = rep(c("Group_A", "Group_B"), cmps_hh_hc)
    df_hh_hc <- as.data.frame(df_hh_hc)
    df_hh_hc$betti0 <- as.numeric(df_hh_hc$betti0)
    
    gb0_tt_th <- 
      ggplot(df_tt_th, 
             aes(x=group, 
                 y=log2(as.numeric(betti0))), 
             group=group) + 
      geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
      scale_shape_manual(values=18) +
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Topological features\ndivergence (log2)") + 
      xlab("")  +
      ggtitle(metric) +
      theme(plot.title = element_text(size=9))
    
    
    gb1_tt_th <- 
      ggplot(df_tt_th, 
             aes(x=group, 
                 y=log2(as.numeric(betti1))), 
             group=group) + 
      geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
      scale_shape_manual(values=18) +
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Topological features\ndivergence (log2)") + 
      xlab("")  +
      ggtitle(metric) +
      theme(plot.title = element_text(size=9))
    

    gb0_hh_hc <- 
      ggplot(df_hh_hc, 
             aes(x=group, 
                 y=log2(as.numeric(betti0))), 
             group=group) + 
      geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
      scale_shape_manual(values=18) +
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Topological features\ndivergence (log2)") + 
      xlab("")  +
      ggtitle(metric) +
      theme(plot.title = element_text(size=9))
    
    
    gb1_hh_hc <- 
      ggplot(df_hh_hc, 
             aes(x=group, 
                 y=log2(as.numeric(betti1))), 
             group=group) + 
      geom_dotplot(binaxis='y', stackdir='center', alpha=.4, stroke=0, dotsize=1.75)+
      scale_shape_manual(values=18) +
      geom_boxplot(fill=NA, width=.5, outlier.shape=NA) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Topological features\ndivergence (log2)") + 
      xlab("")  +
      ggtitle(metric) +
      theme(plot.title = element_text(size=9))
    
    
    stats <- generic_topo_similarity_statistics(df_tt_th, metric, isPaired = F, wilc_alt_to_use = "two.sided")
    stats$signif = signif.num(stats$pval)
    stats$adjusted.signif = signif.num(stats$adjusted.pval)
    stats$comp <- c("TT vs TH")
    statistics_df <- rbind(statistics_df,
                           stats)
    
    stats <- generic_topo_similarity_statistics(df_hh_hc, metric, isPaired = F)
    stats$signif = signif.num(stats$pval)
    stats$adjusted.signif = signif.num(stats$adjusted.pval)
    stats$comp <- c("HH vs HC")
    statistics_df <- rbind(statistics_df,
                           stats)
    
    all_df <- data.frame()
    colnames(hh$pvalue_b0[[metric]]) <- 1:ncol(hh$pvalue_b0[[metric]])
    rownames(hh$pvalue_b0[[metric]]) <- 1:ncol(hh$pvalue_b0[[metric]])
    colnames(hh$pvalue_b1[[metric]]) <- 1:ncol(hh$pvalue_b1[[metric]])
    rownames(hh$pvalue_b1[[metric]]) <- 1:ncol(hh$pvalue_b1[[metric]])
    meltedb0 <- melt(hh$pvalue_b0[[metric]])
    meltedb1 <- melt(hh$pvalue_b1[[metric]])
    pvalues_df <- data.frame(pval_b0=meltedb0$value[meltedb0$Var1 != meltedb0$Var2],
                             pval_b1=meltedb1$value[meltedb0$Var1 != meltedb1$Var2])
    pvalues_df$group = "HH"
    all_df <- rbind(all_df,
                    pvalues_df)
    
    colnames(tt$pvalue_b0[[metric]]) <- 1:ncol(tt$pvalue_b0[[metric]])
    rownames(tt$pvalue_b0[[metric]]) <- 1:ncol(tt$pvalue_b0[[metric]])
    colnames(tt$pvalue_b1[[metric]]) <- 1:ncol(tt$pvalue_b1[[metric]])
    rownames(tt$pvalue_b1[[metric]]) <- 1:ncol(tt$pvalue_b1[[metric]])
    meltedb0 <- melt(tt$pvalue_b0[[metric]])
    meltedb1 <- melt(tt$pvalue_b1[[metric]])
    pvalues_df <- data.frame(pval_b0=meltedb0$value[meltedb0$Var1 != meltedb0$Var2],
                             pval_b1=meltedb1$value[meltedb0$Var1 != meltedb1$Var2])
    pvalues_df$group = "TT"
    all_df <- rbind(all_df,
                    pvalues_df)
    
    colnames(th$pvalue_b0[[metric]]) <- 1:ncol(th$pvalue_b0[[metric]])
    rownames(th$pvalue_b0[[metric]]) <- 1:nrow(th$pvalue_b0[[metric]])
    colnames(th$pvalue_b1[[metric]]) <- 1:ncol(th$pvalue_b1[[metric]])
    rownames(th$pvalue_b1[[metric]]) <- 1:nrow(th$pvalue_b1[[metric]])
    meltedb0 <- melt(th$pvalue_b0[[metric]])
    meltedb1 <- melt(th$pvalue_b1[[metric]])
    pvalues_df <- data.frame(pval_b0=meltedb0$value,
                             pval_b1=meltedb1$value)
    pvalues_df$group = "TH"
    all_df <- rbind(all_df,
                    pvalues_df)
    
    colnames(hc$pvalue_b0[[metric]]) <- 1:ncol(hc$pvalue_b0[[metric]])
    rownames(hc$pvalue_b0[[metric]]) <- 1:nrow(hc$pvalue_b0[[metric]])
    colnames(hc$pvalue_b1[[metric]]) <- 1:ncol(hc$pvalue_b1[[metric]])
    rownames(hc$pvalue_b1[[metric]]) <- 1:nrow(hc$pvalue_b1[[metric]])
    meltedb0 <- melt(hc$pvalue_b0[[metric]])
    meltedb1 <- melt(hc$pvalue_b1[[metric]])
    pvalues_df <- data.frame(pval_b0=meltedb0$value,
                             pval_b1=meltedb1$value)
    pvalues_df$group = "HC"
    all_df <- rbind(all_df,
                    pvalues_df)
    

    values_df <- ddply(all_df, .(group), function(group_df) {c(mean(group_df[,"pval_b0"]), sd(group_df[,"pval_b0"]), sem(group_df[,"pval_b0"]),
                                                               mean(group_df[,"pval_b1"]), sd(group_df[,"pval_b1"]), sem(group_df[,"pval_b1"]),
                                                               nrow(group_df))})
    
    colnames(values_df) <- c("Group", paste(c("mean", "sd", "sem"), rep(c("b0", "b1"), each=3), sep="_"), "N_comp")
    values_df$metric = metric
    
    values_df_all <- rbind(values_df_all,
                           values_df)
     gviolinb0 <- 
      ggplot(all_df, aes(x=group, y=pval_b0, group=group)) +
      geom_violin(aes(fill=group),color=NA, fill="gray60") + 
      geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
      #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
      stat_summary(size=1) + 
      theme_light() + 
      big_text_base_plot_theme + 
      ylab("Pvalue") + 
      xlab("") 
      
      
      gviolinb1 <- 
      ggplot(all_df, aes(x=group, y=pval_b1, group=group)) +
        geom_violin(aes(fill=group),color=NA, fill="gray60") + 
        geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.5)) + 
        #geom_dotplot(binaxis='y', stackdir='center', binwidth = 1/150) +
        stat_summary(size=1) + 
        theme_light() + 
        big_text_base_plot_theme + 
        ylab("Pvalue") + 
        xlab("") 
    

    
    for (size_name in names(a4_sizes)) {
      size = a4_sizes[[size_name]]
      pdf(sprintf("%s\\%s_b0_violin_tt_th.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=size)
      plot(gb0_tt_th)
      dev.off()
      
      pdf(sprintf("%s\\%s_b1_violin_tt_th.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=size)
      plot(gb1_tt_th)
      dev.off()
      
      pdf(sprintf("%s\\%s_b0_violin_hh_hc.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=size)
      plot(gb0_hh_hc)
      dev.off()
      
      pdf(sprintf("%s\\%s_b1_violin_hh_hc.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=size)
      plot(gb1_hh_hc)
      dev.off()
      
      
      pdf(sprintf("%s\\%s_b0_pvalue.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=1.5 * size)
      
      plot(gviolinb0)
      
      dev.off()
      
      pdf(sprintf("%s\\%s_b1_pvalue.pdf",
                  figure_path,
                  size_name),
          height=size,
          width=1.5 * size)
      
      plot(gviolinb1)
      
      dev.off()
      
    }
  }
  
  write.csv(file=sprintf("%s\\statistics.csv", statistics_figure_path),
            statistics_df)
  
  write.csv(file=sprintf("%s\\values.csv", statistics_figure_path),
            values_df_all)
}

commented_junk <- function()
{
  # output_name = "thirst_SFO_new_topological_comparision.Rda"
  # save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), final_res)
  output_name = "hunger_AGRP_new_topological_comparision.Rda"
  save(file=sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), final_res)
  
  # 
  # output_name = "thirst_thirst_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # tt <- final_res
  # 
  # output_name = "hunger_thirst_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # th <- final_res
  # 
  # output_name = "hunger_v1_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # hv1 <- final_res
  # 
  # output_name = "hunger_control_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose = T)
  # hc <- final_res
  # 
  # output_name = "thirst_control_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose = T)
  # tc <- final_res
  # 
  # 
  # output_name = "thirst_v1_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # tv1 <- final_res
  # 
  # 
  # output_name = "hunger_hunger_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # hh <- final_res
  # 
  # 
  # output_name = "thirst_shuffle_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # tshuff <- final_res
  # 
  # output_name = "v1_v1_topological_comparision.Rda"
  # load(sprintf("%s\\data\\figure_1\\%s", base_output_path, output_name), verbose=T)
  # v1v1 <- final_res
  # 
  # 
  # 
  # tt$topological_similarity <-  average_similarity_across_datastes(tt$topological_similarity)
  # th$topological_similarity <-  average_similarity_across_datastes(th$topological_similarity, path1 = F)
  # colnames(th$topological_similarity)[1] <- "path1"
  # hv1$topological_similarity <-  average_similarity_across_datastes(hv1$topological_similarity)
  # hc$topological_similarity <-  average_similarity_across_datastes(hc$topological_similarity)
  # #tc$topological_similarity <-  average_similarity_across_datastes(tc$topological_similarity)
  # tv1$topological_similarity <-  average_similarity_across_datastes(tv1$topological_similarity)
  # hh$topological_similarity <-  average_similarity_across_datastes(hh$topological_similarity)
  # #v1v1$topological_similarity <-  average_similarity_across_datastes(v1v1$topological_similarity)
  # tshuff$topological_similarity <-  average_similarity_across_datastes(tshuff$topological_similarity)
  # 
  # metric = "wasserstein"
  # # pvalues_df <- data.frame(pval=c(c(tt$pvalue_b0[[metric]]),
  # #                                 c(tshuff$pvalue_b0[[metric]])),
  # #                          group=rep(c("T-T", "T-Shuff"), c(len(tt$pvalue_b0[[metric]]),
  # #                                                           len(tshuff$pvalue_b0[[metric]]))))
  # 
  # 
  # 
  # 
  # pvalues_df <- data.frame(pval=melted$value,
  #                          group=paste("T-T", ifelse(as.numeric(melted$Var1 == melted$Var2) > 0, "Within", "Across"), sep="_"))
  # 
  # 
  # 
  # colnames(tt$pvalue_b0[[metric]]) <- 1:14
  # rownames(tt$pvalue_b0[[metric]]) <- 1:14
  # melted <- melt(tt$pvalue_b0[[metric]])
  # melted <- melted[melted$Var1 != melted$Var2,]
  # 
  # pvalues_df <- data.frame(pval=c(melted$value, 
  #                                 c(tv1$pvalue_b0[[metric]]),
  #                                 c(th$pvalue_b0[[metric]]),
  #                                 c(h$pvalue_b0[[metric]])),
  #                          group=rep(c("T-T", "T-V1", "T-H", "H-C"), c(len(melted$value), 
  #                                                                      len(tv1$pvalue_b0[[metric]]),
  #                                                                      len(th$pvalue_b0[[metric]]),
  #                                                                      len(hv1$pvalue_b0[[metric]]))))
  # 
  # 
  # 
  # 
  # pvalues_fraction_hist_df <- 
  #   ddply(pvalues_df, .(group), 
  #         function(group_df) {
  #           h <-  hist(group_df$pval, breaks=seq(0, 1,length.out=51), plot=F)
  #           
  #           
  #           hist_df <- 
  #             data.frame(frac_pval=h$counts / sum(h$counts), 
  #                        breaks_pval=h$breaks[-len(h$counts)], 
  #                        group=rep(group_df$group[1],len(h$counts)))
  #           
  #           return(hist_df)
  #         })
  # 
  # 
  # 
  # ggplot(pvalues_fraction_hist_df) + 
  #   geom_line(aes(y=frac_pval, x=breaks_pval, group=factor(group), fill=group), stat="identity", position=position_dodge(), alpha=.75) + 
  #   scale_y_continuous(expand=c(0,0)) +
  #   ylab("Fraction") +
  #   xlab("P-value") + 
  #   theme_light() + 
  #   base_plot_theme
  # 
  # 
  # 
  # 
  # 
  # 
  # metric="wasserstein"
  # dfa <- rbind(th$topological_similarity[th$topological_similarity$metric == metric,],
  #              #tv1$topological_similarity[tv1$topological_similarity$metric == metric,],
  #              tt$topological_similarity[tt$topological_similarity$metric == metric,],
  #              hc$topological_similarity[hc$topological_similarity$metric == metric,],
  #              hh$topological_similarity[hh$topological_similarity$metric == metric,])
  # 
  # cmps <- c(nrow(th$topological_similarity[th$topological_similarity$metric == metric,]),
  #           #nrow(tv1$topological_similarity[tv1$topological_similarity$metric == metric,]),
  #           nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]),
  #           nrow(hc$topological_similarity[hc$topological_similarity$metric == metric,]),
  #           nrow(hh$topological_similarity[hh$topological_similarity$metric == metric,]))
  # 
  # dfa$group = rep(c("T-H", 
  #                   #"T-C", 
  #                   "T-T",
  #                   "H-C", 
  #                   "H-H"), cmps)
  # dfa <- as.data.frame(dfa)
  # dfa$betti0 <- as.numeric(dfa$betti0)
  # dfa$betti1 <- as.numeric(dfa$betti1)
  # dff <- dfa
  # #dff <- dfa[!(dfa$path1 == dfa$path2 & (dfa$group == "T-T" | dfa$group == "H-H")),]
  # 
  # ggplot(dff, aes(x=factor(group, levels=c("T-T", "T-H", "T-C", "H-H", "H-C")), y=log2(as.numeric(betti0))), group=group) +
  #   #geom_violin(aes(fill=group),color=NA, fill="gray60") + 
  #   #geom_jitter(aes(fill=group), color="black", alpha=.25, stroke=0, position=position_jitterdodge(.25)) + 
  #   geom_dotplot(binaxis='y', stackdir='center', dotsize=.25) +
  #   geom_boxplot(fill=NA, width=.5, outlier.shape = NA)  +
  #   #stat_summary(size=1) + 
  #   theme_light() + 
  #   base_plot_theme + 
  #   #ylab(TeX(sprintf("JSD (%s)", "$\\beta_{0}$"))) + 
  #   ylab("log2(top)") + 
  #   xlab("") 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # geom_line(data=data.frame(x=c(1,3), y=c(max(dff$betti0) - .035, max(dff$betti0) - .035)),
  # #           aes(x=x,y=y)) + 
  # # geom_line(data=data.frame(x=c(1,3), y=c(max(dff$betti0) + 0.035, max(dff$betti0) + 0.035)),
  # #           aes(x=x,y=y)) +
  # # geom_line(data=data.frame(x=c(4,5), y=c(max(dff$betti0) + 0.025, max(dff$betti0) + 0.025)),
  # #           aes(x=x,y=y)) +
  # # geom_line(data=data.frame(x=c(4,5), y=c(max(dff$betti0) + 0.095, max(dff$betti0) + 0.095)),
  # #           aes(x=x,y=y))
  # dt <- dunn.test:::dunn.test(dfa$betti0[!grepl("V1", dfa$group)], dfa$group[!grepl("V1", dfa$group)], method = "bonferroni")
  # dt <- dunn.test:::dunn.test(dfa$betti0, dfa$group, method = "bonferroni")
  # names(dt$P.adjusted) <- dt$comparisons
  # names(dt$P) <- dt$comparisons
  # barplot(dt$P.adjusted)
  # barplot(dt$P)
  # abline(h=.05)
  # 
  # 
  # fraction_hist_df <- 
  #   ddply(dff, .(group), 
  #         function(group_df) {
  #           hb0 <-  hist(group_df$betti0, breaks=seq(min(dff$betti0),
  #                                                    max(dff$betti0),
  #                                                    length.out=31), plot=F)
  #           
  #           hb1 <-  hist(group_df$betti1, breaks=seq(min(dff$betti1),
  #                                                    max(dff$betti1),
  #                                                    length.out=31), plot=F)
  #           
  #           
  #           hist_df <- 
  #             data.frame(frac_betti0=hb0$counts / sum(hb0$counts), 
  #                        breaks_betti0=hb0$breaks[-1], 
  #                        frac_betti1=hb1$counts / sum(hb1$counts), 
  #                        breaks_betti1=hb1$breaks[-1], 
  #                        group=rep(group_df$group[1],len(hb0$counts)))
  #           print(sum(hist_df$frac_betti0))
  #           
  #           return(hist_df)
  #         })
  # 
  # # 
  # ggplot(fraction_hist_df, aes(y=frac_betti0, x=breaks_betti0, group=group)) + 
  #   geom_bar(aes(fill=group), stat="identity", size=2, alpha=.5,  position=position_nudge()) + 
  #   #geom_line(aes(color=group), size=.2, alpha=1, position=position_dodge(.0001)) + 
  #   theme_light() +
  #   scale_y_continuous(expand=c(0,0)) +
  #   base_plot_theme
  # 
  # 
  # 
  # 
  # 
  # 
  # output_name = "hunger_v1_topological_comparision.Rda"
  # load(sprintf("%s\\figure_2\\data\\%s", figures_base_path, output_name), verbose=T)
  # hv1 <- final_res
  # 
  # metric="wasserstein"
  # dfa <- rbind(th$topological_similarity[th$topological_similarity$metric == metric,],
  #              tv1$topological_similarity[tv1$topological_similarity$metric == metric,],
  #              tt$topological_similarity[tt$topological_similarity$metric == metric,],
  #              hv1$topological_similarity[hv1$topological_similarity$metric == metric,])
  # 
  # cmps <- c(nrow(th$topological_similarity[th$topological_similarity$metric == metric,]),
  #           nrow(tv1$topological_similarity[tv1$topological_similarity$metric == metric,]),
  #           nrow(tt$topological_similarity[tt$topological_similarity$metric == metric,]),
  #           nrow(hv1$topological_similarity[hv1$topological_similarity$metric == metric,]))
  # 
  # dfa$group = rep(c("T-H", "T-V1", "T-T", "H-V1"), cmps)
  # dfa <- as.data.frame(dfa)
  # dfa$betti0 <- as.numeric(dfa$betti0)
  # dfa$betti1 <- as.numeric(dfa$betti1)
  # dff <- dfa[!(dfa$path1 == dfa$path2 & dfa$group == "T-T"),]
}
