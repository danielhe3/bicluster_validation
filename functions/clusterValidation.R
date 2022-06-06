clusterValidation <- function(biclusters, eset=NULL, df=NULL, type, fdr=0.05, fold=1.0){
  require(mosbi)
  require(limma)
  require(dplyr)
  require(edgeR)
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
  # Jaccard index function
  jaccard <- function(s1, s2, overlap) {
    return(((overlap) / (s1 + s2 - overlap)))
  }
  # overlap coefficient function
  overlap <- function(s1, s2, overlap) {
    return((overlap / min(s1, s2)))
  }
  if(missing(type)) stop(print("Error: Please specify whether your data is `rnaseq` or `microarray` via the `type` argument."))
  # sanity checks
  if(type=="microarray" & is.null(df)) stop(print("Error: Please use `df` to specify a log-normalized matrix of microarray data, not eset."))
  #if(type=="rnaseq" & is.null(eset)) stop(print("Error: Please use `eset` to specify a DGEList containing counts with normalization factors, not df.")) 
  if(type=="rnaseq" & hasArg(eset) & class(eset)!="DGEList") stop(print("Error: eset object is not of class DGEList. To create a DGEList object, please use the `DGEList` function from the edgeR package and run the function `calcNormFactors` to obtain norm.factors."))
  if(type=="rnaseq" & hasArg(eset) & all(eset$samples$norm.factors==1)) {
    eset <- calcNormFactors(eset, method="TMM")
  }
  #if(!all(lapply(biclusters, function(x) c("colname", "rowname", "column", "row") %in% names(x) %>% unlist() ))) stop('Bicluster list does not contain bicluster sample names under `colname`')
  if(type=="rnaseq" & hasArg(eset)){
    data_matrix <- as.matrix(cpm(eset, log=TRUE))
  }
  if(type=="rnaseq" & hasArg(df)){
    data_matrix <- df
  }
  if(type=="microarray"){
    data_matrix <- df
  }
  
  clusterList <- list()
  if(all(lapply(biclusters, function(x) class(x)=="bicluster")==TRUE)){
    biclusters <- mosbiConvert(biclusters)
  }
  for(i in 1:length(biclusters)){
    
    skip_to_next <- FALSE
    
    design <- data.frame(row.names = colnames(data_matrix),
                         cluster = ifelse(colnames(data_matrix) %in% biclusters[[i]][["colname"]], "bicluster","noncluster"))
    cluster <- ifelse(colnames(data_matrix) %in% biclusters[[i]][["colname"]], "bicluster","noncluster")
    tryCatch(model.matrix(~0+cluster), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }     
    design <- model.matrix(~0+cluster)
    
    if(type=="rnaseq" & hasArg(eset)){
      y <- voom(eset, design, plot = F)
      fit <- lmFit(y, design)
      contr.matrix <- makeContrasts(cluster = clusterbicluster - clusternoncluster, levels=colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contr.matrix)
      tmp <- eBayes(tmp)
    }
    
    if(type=="microarray" | (type=="rnaseq" & hasArg(df))){
      fit <- lmFit(data_matrix, design)
      contr.matrix <- makeContrasts(cluster = clusterbicluster - clusternoncluster, levels=colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contr.matrix)
      tmp <- eBayes(tmp)
    }
    
    shared_genes <- intersect(biclusters[[i]][["rowname"]],
                              topTable(tmp, coef="cluster", n=Inf) %>% dplyr::filter(adj.P.Val<=fdr) %>% dplyr::filter(logFC>=abs(fold)) %>% rownames())
    absent_genes <- setdiff(biclusters[[i]][["rowname"]],
                            topTable(tmp, coef="cluster", n=Inf) %>% dplyr::filter(adj.P.Val<=fdr) %>% dplyr::filter(logFC>=abs(fold)) %>% rownames())
    num_deg <- length(topTable(tmp, coef="cluster", n=Inf) %>% dplyr::filter(adj.P.Val<=fdr) %>% dplyr::filter(logFC>=abs(fold)) %>% rownames())
    num_bicluster_deg <- length(shared_genes) 
    num_bicluster_gene <- length(biclusters[[i]][["rowname"]])
    num_samples <- length(biclusters[[i]][["colname"]])
    cluster_jaccard <- jaccard(length(biclusters[[i]][["rowname"]]),
                               length(topTable(tmp, coef="cluster", n=Inf) %>% dplyr::filter(adj.P.Val<=fdr) %>% dplyr::filter(logFC>=abs(fold)) %>% rownames()),
                               overlap=length(shared_genes)
    )
    cluster_overlap <- overlap(length(biclusters[[i]][["rowname"]]),
                               length(topTable(tmp, coef="cluster", n=Inf) %>% dplyr::filter(adj.P.Val<=fdr) %>% dplyr::filter(logFC>=abs(fold)) %>% rownames()),
                               overlap=length(shared_genes)
    )
    SNR <- local({
      z_score <- function(x, margin = 2) {
        z_fun <- function(y) {
          (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
        }
        
        if (margin == 2) {
          return(apply(x, margin, z_fun))
        } else if (margin == 1) {
          return(t(apply(x, margin, z_fun)))
        }
      }
      #data_matrix <- t(apply(data_matrix, 1, cal_z_score))
      data_matrix = z_score(data_matrix, 1)
      return(mean(data_matrix[biclusters[[i]][["row"]],biclusters[[i]][["column"]]])/sd(data_matrix[biclusters[[i]][["row"]],biclusters[[i]][["column"]]]))
    })
    
    clusterList[[i]] <- list(paste0("Bicluster ", i),
                             biclusters[[i]][["row"]], 
                             biclusters[[i]][["column"]],
                             biclusters[[i]][["rowname"]],
                             biclusters[[i]][["colname"]],
                             biclusters[[i]][["algorithm"]],
                             num_deg, # total num of DEG
                             num_bicluster_deg, # num of bicluster genes that are DE
                             num_bicluster_gene, # num of bicluster genes
                             num_samples,
                             cluster_jaccard,
                             cluster_overlap,
                             SNR,
                             shared_genes, # genes in bicluster that appear in DEG list
                             absent_genes # genes in bicluster that do not appear in DEG list
    )
    names(clusterList[[i]]) <- c("bicluster","row", "column", "rowname", "colname", "algorithm","num_deg", "num_bicluster_deg","num_bicluster_gene","num_samples","jaccard","overlap","SNR", "shared_genes", "absent_genes")
    names(clusterList)[i] <- paste0("Bicluster ", i)
  }
  return(clusterList)
}