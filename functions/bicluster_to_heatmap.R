bicluster_to_heatmap <- function(bicluster_list, df=NULL, eset=NULL, sample_info=NULL, sample_ID, multiple=NULL, values, color_legend, gradient, title_components,text_size=5){
  require(ComplexHeatmap)
  require(grid)
  myplots <- list()
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  mosbiConvert <- function(bicluster_list){
    require(mosbi)
    outList <- list()
    
    for(k in 1:length(bicluster_list)){
      outList[[k]] <- list(bicluster_list[[k]]@row, 
                           bicluster_list[[k]]@column,
                           bicluster_list[[k]]@rowname,
                           bicluster_list[[k]]@colname,
                           bicluster_list[[k]]@algorithm)
      names(outList[[k]]) <- c("row", "column", "rowname", "colname", "algorithm")
      names(outList)[k] <- paste0("Bicluster ", k)
    }
    return(outList)
  }
  if(all(lapply(bicluster_list, function(x) class(x)=="bicluster")==TRUE)){
    bicluster_list <- mosbiConvert(bicluster_list)
  }
  # extracting info from eset
  if(hasArg(eset) & hasArg(df)) stop(print("Please specify only one of `df` or `eset`."))
  if(hasArg(eset)){
    if(all(eset$samples$norm.factors==1)){eset <- calcNormFactors(eset, method="TMM")}
    data_matrix = cpm(eset, log=TRUE)
    data_matrix <- t(apply(data_matrix, 1, cal_z_score))
    sample_info = eset$samples
  }
  if(hasArg(df)){
    data_matrix <- t(apply(df, 1, cal_z_score))
    sample_info <- sample_info
  }
  if(!hasArg(sample_ID)) stop(print("Please specify the column name of `sample_info` containing sample IDs"))

  for (i in 1:length(bicluster_list)) {
    #message(i)
    myplots[[i]] <- local({
      #i <- i
      df_cluster <- data_matrix[(bicluster_list[[i]]["row"] %>% unlist()),] %>% as.data.frame() %>%
        dplyr::select(as.vector(unlist(bicluster_list[[i]][["colname"]])), everything()) %>% as.matrix()
      # create a cluster legend
      sample_col <- data.frame(row.names=colnames(as.data.frame(df_cluster)),
                               cluster=ifelse(colnames(as.data.frame(df_cluster)) %in% unlist(bicluster_list[[i]][["colname"]]), 
                                              "BiCluster", 
                                              "NonCluster"))
      # adding groups of interest to cluster legend
      sample_col <- sample_col %>% rownames_to_column(all_of(sample_ID)) %>% 
        left_join(sample_info %>% dplyr::select(all_of(sample_ID), all_of(values)), by=all_of(sample_ID)) %>% column_to_rownames(all_of(sample_ID))
      
      title <-paste0(names(bicluster_list)[1], " (", length(bicluster_list[[i]][["colname"]]), " samples)", 
                     "\nalgorithm: ",bicluster_list[[i]][["algorithm"]],
                     if(any("DEG" %in% title_components)){
                       paste0("\n", length(bicluster_list[[i]][["shared_genes"]]),"/", length(bicluster_list[[i]][["rowname"]]), " DE bicluster genes")
                       },
                     if(any("overlap" %in% title_components)){
                       paste0(
                         "\n overlap coefficient: ", format(signif(bicluster_list[[i]][["overlap"]], digits=3), nsmall=3)
                       )
                     },
                     if(any("jaccard" %in% title_components)){
                       paste0(
                         "\n jaccard similarity: ", format(signif(bicluster_list[[i]][["jaccard"]], digits=3), nsmall=3)
                       )
                     },
                     if(any("SNR" %in% title_components)){
                       paste0(
                         "\n Signal-to-Noise Ratio: ", format(signif(bicluster_list[[i]][["SNR"]], digits=3), nsmall=3)
                       )
                     },
                     if(any("K-W_p" %in% title_components)){
                       paste0(
                         "\n Kruskal-Wallis p-value = ", formatC(bicluster_list[[i]][[grep("KW",names(bicluster_list[[i]]))]]$p.value, format = "e", digits = 2)
                       )
                     },
                     if(any("K-W_fdr" %in% title_components)){
                       paste0(
                         "\n Kruskal-Wallis FDR = ", formatC(bicluster_list[[i]][[grep("KW",names(bicluster_list[[i]]))]]$fdr, format = "e", digits = 2)
                       )
                     },
                     if(any("chi_p" %in% title_components)){
                       paste0(
                         "\n ChiSq p-value = ", formatC(bicluster_list[[i]][[grep("chisq",names(bicluster_list[[i]]))]]$p.value, format = "e", digits = 2)
                       )
                     },
                     if(any("chi_fdr" %in% title_components)){
                       paste0(
                         "\n ChiSq FDR = ", formatC(bicluster_list[[i]][[grep("chisq",names(bicluster_list[[i]]))]]$fdr, format = "e", digits = 2)
                       )
                     },
                     if(any("fisher_p" %in% title_components)){
                       paste0(
                         "\n Fisher p-value = ", formatC(bicluster_list[[i]][[grep("fisher",names(bicluster_list[[i]]))]]$p.value, format = "e", digits = 2)
                       )
                     },
                     if(any("fisher_fdr" %in% title_components)){
                       paste0(
                         "\n Fisher FDR = ", formatC(bicluster_list[[i]][[grep("fisher",names(bicluster_list[[i]]))]]$fdr, format = "e", digits = 2)
                       )
                     }
                               
      )
      
      if(i %% multiple == 0){
        p1 <- suppressMessages(ComplexHeatmap::pheatmap(mat = df_cluster,
                                       border_color=NA,
                                       col=gradient,
                                       cluster_rows = TRUE,
                                       cluster_cols = FALSE,
                                       show_rownames = FALSE,
                                       show_colnames = FALSE,
                                       annotation_legend = TRUE,
                                       annotation_names_col = TRUE,
                                       annotation_colors = color_legend,
                                       annotation_col = sample_col,
                                       main = title,
                                       fontsize=text_size)
        )
        return(grid.grabExpr(draw(p1)))
      }
      else{
        p1 <- suppressMessages(ComplexHeatmap::pheatmap(mat = df_cluster,
                                       border_color=NA,
                                       col=gradient,
                                       cluster_rows = TRUE,
                                       cluster_cols = FALSE,
                                       show_rownames = FALSE,
                                       show_colnames = FALSE,
                                       annotation_legend = FALSE,
                                       annotation_names_col = FALSE,
                                       annotation_colors = color_legend,
                                       annotation_col = sample_col,
                                       main = title,
                                       fontsize=text_size)
        )
        return(grid.grabExpr(draw(p1)))
      }
    })
  }
  names(myplots) <- names(bicluster_list)
  return(myplots)
}