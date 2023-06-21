#' Uses output of "build_char_table" to perform differential expression for 
#'  three types of data: 1.) lipid species, 2.) fatty acids, and 3.) lipid 
#'  classes
#' 
#' @param exp_data Lipid expression table. Output of "build_char_table".
#' @param char_table Lipid characteristics table. Output of "build_char_table".
#' @param method A character string indicating which statistical method to be. 
#'  Currently, iLipidome supports two-sample t-tests ("t.test"), Wilcoxon tests 
#'  ("wilcox.test"), or moderated t-tests ("mod.t.test").
#' @param significant "p_value" or "adj_p_value" can be used for selecting 
#'  significant lipid features.
#' @param ctrl_group An integer vector specifying samples in the control group 
#'  in lipid expression table. Note that the first column containing lipid names
#'  is not counted.
#' @param exp_group An integer vector specifying samples in the experimental 
#'  group in lipid expression table. Note that the first column containing lipid
#'  names is not counted.
#'  
#' @return list of expression tables
#' 
#' @importFrom dplyr everything
#' @importFrom dplyr group_by
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' 
#' @export
unprocessed_data_test <- function(exp_data, char_table, method='t.test',
                                  significant=p_value, ctrl_group, 
                                  exp_group){
  
  ctrl_group <- ctrl_group+1
  exp_group <- exp_group+1
  
  #FA expression
  
  char_sel <- char_table %>% filter(feature %in% exp_data$feature)
  
  
  each_FA <- str_extract_all(char_sel$FA_split, '\\d+:\\d+;\\d+')
  
  FA_mapping <- unique(unlist(each_FA)) %>% sort()
  
  Lipid_FA_freq_list <- each_FA  %>% 
    map(.f = function(x){factor(sort(x[x!='0:0;0']), levels = FA_mapping) %>% 
        table() %>% as.integer()})
  
  Lipid_FA_freq_matrix <- matrix(unlist(Lipid_FA_freq_list), ncol = length(Lipid_FA_freq_list))
  
  exp_matrix <- exp_data[-1] %>% as.matrix()
  
  
  FA_exp <- Lipid_FA_freq_matrix %*% exp_matrix %>% as.data.frame()
  
  FA_exp <- FA_exp %>% mutate(feature=FA_mapping, type='FA') %>% 
    dplyr::select(feature, type, everything())
  
  #lipid class expression
  
  exp_data[-1][is.na(exp_data[-1])] <- 0
  
  lipid_class_exp <- exp_data %>% left_join(char_table[c('feature', 'class')], by='feature') %>% 
    dplyr::select(-1)
  
  lipid_class_exp <- lipid_class_exp[!is.na(lipid_class_exp[['class']]),]
  if(nrow(lipid_class_exp)==0){
    lipid_class_exp <- data.frame()
  }else{
    lipid_class_exp <- lipid_class_exp %>% 
      aggregate(as.formula(str_c('. ~ class')), ., sum)
  }
  
  colnames(lipid_class_exp)[1] <- 'feature'
  
  lipid_class_exp <- lipid_class_exp %>% mutate(type='class') %>% 
    dplyr::select(feature, type, everything())
  
  lipid_species_exp <- exp_data %>% mutate(type='species') %>% 
    dplyr::select(feature, type, everything())
  
  all_exp_data <- rbind(lipid_species_exp, lipid_class_exp, FA_exp)
  
  #statistics test
  lipid <- character()
  type <- character()
  
  mean_ctrl <- numeric()
  mean_exp <- numeric()
  mean_all <- numeric()
  
  FC <- numeric()
  
  sd_ctrl <- numeric()
  sd_exp <- numeric()
  p_value <- numeric()
  statistics <- numeric()
  
  ctrl_group <- ctrl_group+1
  exp_group <- exp_group+1
  # if(method=='mod.t.test'){
  #   
  #   g2 <- rep("group 2", length(c(ctrl_group, exp_group)))
  #   g2[(exp_group-2)] <- "group 1"
  #   g2 <- factor(g2)
  #   print(g2)
  #   mod_t_test <- tryCatch({mod.t.test(as.matrix(all_exp_data[-c(1:2)]), group = g2)},
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
  for (a  in 1:nrow(all_exp_data)) {
    lipid[a] <- all_exp_data$feature[a]
    type[a] <- all_exp_data$type[a]
    
    mean_ctrl[a] <- all_exp_data[a, ctrl_group] %>% unlist() %>% mean(na.rm=T)
    mean_exp[a] <- all_exp_data[a, exp_group] %>% unlist() %>% mean(na.rm=T)
    mean_all[a] <- all_exp_data[a, c(ctrl_group, exp_group)] %>% unlist() %>% mean(na.rm=T)
    
    FC[a] <- mean_exp[a]/mean_ctrl[a]
    sd_ctrl[a] <- all_exp_data[a, ctrl_group] %>% unlist() %>% sd(na.rm=T)
    sd_exp[a] <- all_exp_data[a, exp_group] %>% unlist() %>% sd(na.rm=T)
    
    if(method=='t.test'){
      p_value[a] <- tryCatch({t.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({t.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }else if(method=='wilcox.test'){
      p_value[a] <- tryCatch({wilcox.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({wilcox.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    
  }
  
  result <- data.frame(lipid=lipid, type=type, mean_all=mean_all, mean_ctrl=mean_ctrl, mean_exp=mean_exp, sd_ctrl=sd_ctrl,
                       sd_exp=sd_exp, FC=FC, log2FC=log2(FC), statistics=statistics, p_value=p_value, mlog10p=-log10(p_value))
  
  result <- result %>% group_by(type) %>% 
    mutate(adj_p_value = stats::p.adjust(p_value, method = 'fdr')) %>% 
    mutate(mlog10padj=-log10(adj_p_value))
  
  
  if(significant=='p_value'){
    result <- result %>% mutate(sig=ifelse(p_value<0.05, 'yes', 'no'))
  }else{
    result <- result %>% mutate(sig=ifelse(adj_p_value<0.05, 'yes', 'no'))
  }
  
  return(list(all_exp_data, result))
}