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
    names(outList)[k] <- paste0("Community ", k)
  }
  return(outList)
}