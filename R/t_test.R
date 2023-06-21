#' t-test 
#' 
#' @param data FA substructure profile. Output of "lipid_sub_matrix".
#' @param ctrl An integer vector specifying samples in the control group in 
#'  lipid expression table. Note that the first column containing lipid names 
#'  is not counted.
#' @param exp An integer vector specifying samples in the experimental group in 
#'  lipid expression table. Note that the first column containing lipid names is
#'   not counted.
#' @param method A character string indicating which statistical method to be. 
#'  Currently, iLipidome supports two-sample t-tests ("t.test"), Wilcoxon tests 
#'  ("wilcox.test"), or moderated t-tests ("mod.t.test").
#' @param significant "p_value" or "adj_p_value" can be used for selecting 
#'  significant lipid features.
#'  
#' @return statistical results table
#' 
#' @importFrom stats p.adjust
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' 
#' @export
t_test <- function(data, ctrl, exp, method='t.test', significant='p_value'){
  lipid <- character()
  mean_ctrl <- numeric()
  mean_exp <- numeric()
  mean_all <- numeric()
  
  FC <- numeric()
  
  sd_ctrl <- numeric()
  sd_exp <- numeric()
  p_value <- numeric()
  statistics <- numeric()
  
  # if(method=='mod.t.test'){
  #   
  #   g2 <- rep("group 2", length(c(ctrl, exp)))
  #   g2[exp] <- "group 1"
  #   g2 <- factor(g2)
  #   data_m <- as.matrix(data)
  #   mod_t_test <- tryCatch({mod.t.test(data_m, group = g2)},
  #                          error=function(e){return(NULL)})
  #   
  #   if(!is.null(mod_t_test)){
  #     p_value <- mod_t_test$p.value
  #     statistics <- mod_t_test$t
  #   }else{
  #     p_value <- NA
  #     statistics <- NA
  #   }
  #   
  # }
  
  
  for (a  in 1:nrow(data)) {
    lipid[a] <- rownames(data)[a]
    mean_ctrl[a] <- data[a, ctrl] %>% unlist() %>% mean(na.rm=T)
    mean_exp[a] <- data[a, exp] %>% unlist() %>% mean(na.rm=T)
    mean_all[a] <- data[a, c(ctrl, exp)] %>% unlist() %>% mean(na.rm=T)
    FC[a] <- mean_exp[a]/mean_ctrl[a]
    sd_ctrl[a] <- data[a, ctrl] %>% unlist() %>% sd(na.rm=T)
    sd_exp[a] <- data[a, exp] %>% unlist() %>% sd(na.rm=T)
    if(method=='t.test'){
      p_value[a] <- tryCatch({t.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({t.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    else if(method=='wilcox.test'){
      p_value[a] <- tryCatch({wilcox.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({wilcox.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    
  }
  adj_p_value <- p.adjust(p_value, method = 'fdr', n = length(p_value))
  result <- data.frame(lipid=lipid, mean_all=mean_all, mean_ctrl=mean_ctrl, mean_exp=mean_exp, sd_ctrl=sd_ctrl,
                       sd_exp=sd_exp, FC=FC, log2FC=log2(FC), statistics=statistics, p_value=p_value, mlog10p=-log10(p_value),
                       adj_p_value=adj_p_value, mlog10padj=-log10(adj_p_value))
  
  if(significant=='p_value'){
    result <- result %>% mutate(sig=ifelse(p_value<0.05, 'yes', 'no'))
  }
  else{
    result <- result %>% mutate(sig=ifelse(adj_p_value<0.05, 'yes', 'no'))
  }
  
  result$lipid <- str_replace(result$lipid, '_FA\\d+','')
  
  return(result)
}