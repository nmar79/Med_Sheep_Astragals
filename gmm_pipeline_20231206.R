#' gmm_pipeline_voles
#' Coded by N Marom <nimrod.arch@gmail.com> 
#' @param spec_info_file specimen info file with "specID" and "group" column names.
#' @param sliders_file a csv file with column names "before", "slide", "after" defining semisliding landmarks ("slide") and the landmarks they slide between. Using the BE default for sliding.
#' @param tps_file tps file with 2D landmark data. I'm using here the landmark/semisliding landmark configuration from Piras et al. (2012) for water voles, with a single change: landmark 11 on the posterior of the tooth coded as a semisliding landmark between landmarks 10 and 12. 
#' @param order_vect a vector of the order in which the landmarks are to be connected in concensus shape plots. I'm using the order in Piras et al. as the default (c(1,19,20,2,21,3,22,4,23,5,24,6,25,7,26,8,27,9,28,10,11,12,29,13,30,14,31,15,32,16,33,17,34,18,35) )
#' @param library_check whether to install/update the required libraries automatically when the function is run. The default is FALSE. The required packages can be downloaded from CRAN 
#'
#' @return 
#' @export
#'
#' @examples gmm_pipeline("redig_specinfo_20230222.csv", "sliders_11.csv", "lmrks_20230222.tps", c(1,19,20,2,21,3,22,4,23,5,24,6,25,7,26,8,27,9,28,10,11,12,29,13,30,14,31,15,32,16,33,17,34,18,35), library_check = TRUE)
#' 
gmm_pipeline <- function(spec_info_file, sliders_file, tps_file, order_vect = c(1,19,20,2,21,3,22,4,23,5,24,6,25,7,26,8,27,9,28,10,11,12,29,13,30,14,31,15,32,16,33,17,34,18,35), library_check = FALSE) { # nolint: line_length_linter.
  
  #library checkup
if (library_check == TRUE) {print("checking geomorph...")
  if (isTRUE((grep("geomorph", installed.packages())) > 0) == FALSE) {install.packages("geomorph")} else {print("geomorph is installed")}
  
  print("checking Morpho...")
  if (isTRUE((grep("Morpho", installed.packages())) > 0) == FALSE) {install.packages("devtools")} else {print("Morpho is installed")}
  
  print("checking ggsci...")
  if (isTRUE((grep("ggsci", installed.packages())) > 0) == FALSE) {install.packages("ggsci")} else {print("ggsci is installed")}
  
  print("checking ggplot2...")
  if (isTRUE((grep("ggplot2", installed.packages())) > 0) == FALSE) {install.packages("ggplot2")} else {print("ggplot is installed")}
  
}
  
  
  require(geomorph)
  require(Morpho)
  require(ggplot2)
  require(ggsci)
  require(ape)
  
  #set directories
  
  working_dir <- getwd()
  print(paste0("Working directory is ", working_dir))
  dir.create(paste0("results_", format(Sys.time(), "%H_%M_%S_%Y_%M")))
  output_dir <- max(list.dirs())
  print(paste0("Output directory is ", output_dir))
 
  
  #read data files
  print("performing Generalized Procrustes Analysis in geomorph...")
  spec_info <- read.csv(spec_info_file)
  landmarks <- readland.tps(tps_file, readcurves = TRUE, specID = 'ID')
  sliders <- read.csv(sliders_file)
  #GPA
  gpa_out <<- gpagen(landmarks, curves = sliders)
  #size
  print("plotting centroid sizes...")
  size <- gpa_out$Csize
  size_plot_df <- cbind.data.frame(factor(spec_info$group), size)
  colnames(size_plot_df) <- c("group", "size")
  size_plot <- ggplot(size_plot_df, aes(x = group, y = size, color = group)) + geom_boxplot(aes(fill = NA)) + geom_jitter() + scale_color_npg() + scale_fill_npg() + theme_bw()
  plot(size_plot)
  size_foo <- aov(gpa_out$Csize ~ spec_info$group)
  summary(size_foo)
  print("Saving results...")
  setwd(output_dir)
  sink("size_ANOVA.txt")
  print(summary(size_foo))
  print(TukeyHSD(size_foo))
  sink()
  setwd(working_dir)
  
  #Shape ANOVA
  print("calculating Procrustes distance linear model in geomorph, shape ~ group * centroid size...")
  model_df <- geomorph.data.frame(gpa_out, group = spec_info$group)
  mod_1 <<- procD.lm(coords ~ group * Csize, data = model_df, RRPP = TRUE)
  print(summary(mod_1))
 print("created 'procLM.txt with the results.")
 
  setwd(output_dir) 
  sink(file = "procLM.txt")
  print(summary(mod_1))
  sink()
  setwd(working_dir)
 
  #Disparity
  print("calculating morphological disparity with geomorph, call: shape ~ group * size")
  disparity_1 <<- morphol.disparity(mod_1, groups = model_df$group)
  print(disparity_1)
  plot(nj(dist(disparity_1$PV.dist)), "u")
  print("calculating morphological disparity with geomorph, call: shape ~ group")
  disparity_2 <<- morphol.disparity(gpa_out$coords ~ 1, groups = model_df$group)
  print(disparity_2)
  plot(nj(dist(disparity_2$PV.dist)), "u")
  setwd(output_dir)
  
  sink("disparity_group_size.txt")
  print(disparity_1)
  sink()
  
  sink("disparity_group.txt")
  print(disparity_2)
  sink()
  
  setwd(working_dir)
  
  
  #PCA
  print("doing pca in geomorph...")
  pca_out <<- gm.prcomp(gpa_out$coords)
  #Prepare 10 PC data frame
  ten_pcs <- pca_out$x[,1:10]
  ten_pcs <- cbind.data.frame(spec_info$specID, factor(spec_info$group),
   ten_pcs)
  colnames(ten_pcs)[1:2] <- c("ID", "group")
  #Plot first PCs
  pc_variances <- pca_out$sdev^2/sum(pca_out$sdev^2)
  #get relative proportion of variance explained
  pc_variances_axes <- as.numeric((round(pc_variances[1:2], digits = 2))*100) 
  #round and convert to percentages
  pc_plot <- ggplot(ten_pcs, aes(x = Comp1, y = Comp2, colour = group, 
  label = ID)) + geom_point(size = 8, alpha = 0.3) + theme_classic() + scale_fill_aaas() + scale_colour_aaas() + xlab(paste0("PC1"," ",pc_variances_axes[1], "%")) + ylab(paste0("PC2"," ",pc_variances_axes[2], "%")) + geom_text()
  
  plot(pc_plot)
  
  print("Plotting PC1-PC2 shape extremes...")
  
  plot_opposites <- function(gpa_out, pca_out, ord){
    
    nlmks <- length(gpa_out$coords[ ,1 ,1])
    
    #find extreme specimens  
    min_pc1 <- grep(min(pca_out$x[, 1]), pca_out$x[, 1])
    max_pc1 <- grep(max(pca_out$x[, 1]), pca_out$x[, 1])
    min_pc2 <- grep(min(pca_out$x[, 2]), pca_out$x[, 2])
    max_pc2 <- grep(max(pca_out$x[, 2]), pca_out$x[, 2])
    
    min_pc1_coords <- gpa_out$coords[ , , min_pc1]
    max_pc1_coords <- gpa_out$coords[ , , max_pc1]
    min_pc2_coords <- gpa_out$coords[ , , min_pc2]
    max_pc2_coords <- gpa_out$coords[ , , max_pc2]
    
    extreme_coords_pc1 <- list(min_pc1_coords, max_pc1_coords)
    names(extreme_coords_pc1) <- c("min PC1", "max PC1")
    extreme_coords_pc2 <- list(min_pc2_coords, max_pc2_coords)
    names(extreme_coords_pc2) <- c("min PC2", "max PC2")
    extreme_coords <- list(extreme_coords_pc1, extreme_coords_pc2)
    
    for (h in 1:2){
      
      
      for (i in 1:2){
        temp_coords <- extreme_coords[[h]][[i]]
        new_coords <- extreme_coords[[h]][[i]]
        
        for (j in 1:nlmks){new_coords[j,] <- temp_coords[ord[j],]}
        
        plot(new_coords, main = paste0(names(extreme_coords[[h]][i])))
        
        newer_coords <- array(1:2*(length(ord)+1), dim = (c(length(ord)+1, 2)))
        newer_coords[1:length(ord), ] <- new_coords
        newer_coords[length(ord) + 1, ] <- new_coords[1, ]
        if (i == 1){
          lines(newer_coords, col = "black")
        } else {
          lines(newer_coords, col = "red")
        }
      }
    }
  }
  
  plot_opposites(gpa_out, pca_out, order_vect)
  
  
  #run CVA
  print("running cva in Morpho...")
  cva_input_mx <- as.matrix(pca_out$x[, 1:10])
  groups <- ten_pcs$group
  cva_out <<- CVA(cva_input_mx, groups, rounds = 10000, cv = TRUE)
  #cross-validated classification results
  classifications <- typprobClass(cva_out$CVscores, groups = groups)
  specimen_probs <- cbind.data.frame(spec_info[, 1], classifications$probsCV)
  colnames(specimen_probs)[1] <- "ID"
  classify_out <<- classify(cva_out)
  confusion_df <- data.frame(classify_out$groups, classify_out$class)
  print("writing classification probabilities to file 'spec_probs.csv'")
  setwd(output_dir)
  write.csv(specimen_probs, file = "spec_probs.csv")
  write.csv(confusion_df, file = "classification.csv")
  setwd(working_dir)
  
  #Plot CVA
  print("cva plot...")
  axis_variances <- as.numeric(c(round(cva_out$Var[1, 2], digits = 2),
  round(cva_out$Var[2,2],digits = 2)))
  cva_plot_df <- cbind.data.frame(spec_info$specID, groups,
  cva_out$CVscores[, 1:2])
  colnames(cva_plot_df) <- c("ID", "group", "CV1", "CV2")
cva_plot <- ggplot(cva_plot_df, aes(x = CV1, y = CV2,
colour = group, label = ID)) + geom_point(size = 8, alpha = 0.2) + theme_classic() + scale_fill_ucscgb() + scale_colour_ucscgb() + xlab(paste0("CV1", " ", axis_variances[1], "%")) + ylab(paste0("CV2", " ", axis_variances[2], "%")) + geom_text()
print(cva_plot)


  #dendrogram of between group Mahalanobis distances from CVA
print("clustering the Mahalanobis distance between group centroids...")
  dendroS = hclust(cva_out$Dist$GroupdistMaha)
  dendroS$labels = levels(groups)
  par(mar = c(4, 4.5, 1, 1))
  dendroS = as.dendrogram(dendroS)
  mahalanobis_dist <- plot(dendroS, main = '', sub = '',
                           ylab = 'Mahalahobis distance')
  
 mahalanobis_dist
  
 
 print("clustering the Euclidean distance between group centroids...")
 dendroS2 = hclust(cva_out$Dist$GroupdistEuclid)
 dendroS$labels = levels(groups)
 par(mar = c(4, 4.5, 1, 1))
 dendroS2 = as.dendrogram(dendroS2)
 euclid_dist <- plot(dendroS2, main = '', sub = '',
                          ylab = 'Euclidean distance')
 
 euclid_dist
  
  #plot concensus shapes for each group
  
  print("plotting mean shape for each group...")
  
  plot_concensus_shapes <- function (spec_info_file, gpa_object, order_vect){
    spec_info <- read.csv(spec_info_file)
    spec_info$group <- factor(spec_info$group)
    gp_vole <- gpa_object
    ngroups <- length(unique(spec_info$group))
    group_names <- unique(spec_info$group)
    group_members <- list()
    for (i in 1:ngroups) {
      group_members[[i]] <- which(spec_info$group == group_names[i])
    }
    names(group_members) <- group_names
    group_coords <- list()
    for (i in 1:length(group_members)) {
      group_coords[[i]] <- gp_vole$coords[, , c(group_members[[i]])]
    }
    group_concensus <- list()
    lmks <- 1:length(gp_vole$coords[, 1, 1])
    XY <- 1:2
    spec <- 1:length(group_members)
    concensus <- array(c(lmks, XY, spec), dim = c(length(lmks), 
                                                  2, length(group_members)))
    for (i in 1:length(group_members)) {
      for (j in 1:length(lmks)) {
        concensus[j, 1:2, i] <- c(mean(group_coords[[i]][j, 
                                                         1, ]), mean(group_coords[[i]][j, 2, ]))
      }
    }
    for (i in 1:length(group_members)) {
      plot(concensus[, , i], main = paste("concensus shape", 
                                          group_names[i], sep = " "))
      temp_coords <- concensus[, , i]
      new_coords <- concensus[, , i]
      for (j in 1:length(lmks)) {
        new_coords[j, ] <- temp_coords[order_vect[j], ]
      }
      newer_coords <- array(1:2 * (length(order_vect) + 1), 
                            dim = (c(length(order_vect) + 1, 2)))
      newer_coords[1:length(order_vect), ] <- new_coords
      newer_coords[length(order_vect) + 1, ] <- new_coords[1, 
      ]
      lines(newer_coords)
    }
  }
  
  spread_foo <- plot_concensus_shapes(spec_info_file, gpa_out, order_vect)
  spread_foo
  print("DONE!")
}


