# also need to specify test... will implement chisquare, kruskal wallis, and fisher's
bicluster_stat_test <- function(bicluster_list, sample_info, sample_id, outcome, outcome_cat = "factor", order, test){
  require(dplyr)
  require(stats)
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
  # bicluster_list checks
  if(missing(bicluster_list)) stop('Please provide a list of biclusters, with each list element containing bicluster information such as bicluster sample names (`colname`), bicluster gene names (`rowname`). For mosbi objects, use the `convertBicluster` function to change the output list of mosbi bicluster run_ methods.')
  if(!is.list(bicluster_list)) stop('`bicluster_list` must be a list of biclusters, with each list element containing bicluster information such as bicluster sample names (`colname`), bicluster gene names (`rowname`). Please use the `convertBicluster` function to change output list of mosbi bicluster run_ methods.')
  if(!all(lapply(bicluster_list, function(x) "colname" %in% names(x)) %>% unlist())) stop('Bicluster list does not contain bicluster sample names under `colname`')
  
  # sample_info checks
  if(missing(sample_info)) stop('Please provide a dataframe with columns corresponding to sample IDs and outcome of interest.')
  if(!is.data.frame(sample_info)) stop('`sample_info` must be a dataframe with columns of sample names and outcome of interest')
  
  # sample_id checks
  if(missing(sample_id)) stop('Please provide the name of the column in your sample_info dataframe that contains your sample names.')
  if(!sample_id %in% colnames(sample_info)) stop('`sample_id` column does not exist in provided sample_info file, please check for typos!')
  if(!any( (lapply( bicluster_list, function(x) x$colname) %>% unlist() ) %in% sample_info[[all_of(sample_id)]])) stop('Sample names in `bicluster_list` and `sample_info` do not match.')
  
  # outcome checks
  if(missing(outcome)) stop('Please provide the name of the column in your sample_info dataframe for which you would like to perform a statistical test.')
  if(length(outcome)>1) stop('Only calculate one outcome at a time, please.')
  if(!outcome %in% colnames(sample_info)) stop('`Outcome` column does not exist in provided sample_info file, please check for typos!')
  
  # test checks
  if(missing(test)) stop('Please specify a statistical test (kruskal-wallis, chisq, and/or fisher) under the `test` argument.')
  if(!test %in% c("kruskal-wallis", "chisq", "fisher")) stop('Please specify "kruskal-wallis", "chisq", and/or "fisher" under the `test` argument.')
  if(test=="kruskal-wallis" & outcome_cat == "factor" & missing(order)) stop('Please provide a vector with the order of outcome categories from lowest to highest using the `order` argument.')
  if(test=="fisher" & length(levels(sample_info[[all_of(outcome)]])) > 2) stop("Fisher's exact test is not recommended for more than two outcomes.")
  
  if(all(lapply(bicluster_list, function(x) class(x)=="bicluster")==TRUE)){
    bicluster_list <- mosbiConvert(bicluster_list)
  }
  
  for(i in 1:length(bicluster_list)){
    if("kruskal-wallis" %in% test){
      bicluster_list[[i]][[paste0("KW_",outcome)]] <- local({
        # identifying samples in bicluster
        bicluster_samples <- bicluster_list[[i]][["colname"]]
        # creating column indicating bicluster status
        sample_info$cluster <- ifelse(sample_info[[all_of(sample_id)]] %in% bicluster_samples, "bicluster","noncluster")
        # selecting only outcome and cluster info
        df <- sample_info %>% dplyr::select(outcome, cluster)
        if(outcome_cat=="factor"){
          # setting the outcome as a factor
          df$outcome.f <- factor(df[[all_of(outcome)]], levels = order, ordered=TRUE)
        } else{
          df$outcome.f <- df[[all_of(outcome)]] %>% as.integer()
        }
        return(kruskal.test(outcome.f ~ cluster, data = df))
      })
    }
    if("chisq" %in% test){
      bicluster_list[[i]][[paste0("chisq_",outcome)]] <- local({
        bicluster_samples <- bicluster_list[[i]][["colname"]]
        sample_info$cluster <- ifelse(sample_info[[all_of(sample_id)]] %in% bicluster_samples, "bicluster","noncluster")
        contingency <- table(sample_info[[all_of(outcome)]],sample_info[["cluster"]])
        chisq <- suppressWarnings(chisq.test(contingency))
        return(chisq)
      })
    }
    if("fisher" %in% test){
      bicluster_list[[i]][[paste0("fisher_",outcome)]] <- local({
        bicluster_samples <- bicluster_list[[i]][["colname"]]
        sample_info$cluster <- ifelse(sample_info[[all_of(sample_id)]] %in% bicluster_samples, "bicluster","noncluster")
        contingency <- table(sample_info[[all_of(outcome)]],sample_info$cluster)
        fisher <- fisher.test(contingency)
        return(fisher)
      })
    }
  }
  
  # applying FDR correction to p.values
  if("kruskal-wallis" %in% test){
    fdr_values_kw <- p.adjust(lapply(bicluster_list, function(x) x[[paste0("KW_",outcome)]][["p.value"]] ), method = "fdr")
    for(i in 1:length(bicluster_list)){
      bicluster_list[[i]][[paste0("KW_",outcome)]][["fdr"]] <- fdr_values_kw[i] %>% unname()
    }
  }
  if("chisq" %in% test){
    fdr_values_chi <- p.adjust(lapply(bicluster_list, function(x) x[[paste0("chisq_",outcome)]][["p.value"]] ), method = "fdr")
    for(i in 1:length(bicluster_list)){
      bicluster_list[[i]][[paste0("chisq_",outcome)]][["fdr"]] <- fdr_values_chi[i] %>% unname()
    }
  }
  if("fisher" %in% test){
    fdr_values_fisher <- p.adjust(lapply(bicluster_list, function(x) x[[paste0("fisher_",outcome)]][["p.value"]] ), method = "fdr")
    for(i in 1:length(bicluster_list)){
      bicluster_list[[i]][[paste0("fisher_",outcome)]][["fdr"]] <- fdr_values_fisher[i] %>% unname()
    }
  }
  
  return(bicluster_list)
}