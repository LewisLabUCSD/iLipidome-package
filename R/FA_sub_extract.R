#' maps FA substructures in each pathway with fold changes from the 
#'  “unprocessed_data_test” result and extracts them
#'  
#' @param char_table Lipid characteristics table. Output of "build_char_table".
#' @param FA_substructure FA substructure table. Output of "FA_sub_transform".
#' @param unprocessed_data_result Differential expression for unprocessed 
#'  lipidomics data. Output of "unprocessed_data_test".
#' @param exact_FA "yes" or "no" to decide if the exact FA identity 
#'  (e.g., w3 or w6) is known or not.
#' @param exo_lipid An character vector specifying exogenous lipid addition. The
#'  exogenous lipids and their adjacent nodes in FA network will not be 
#'  decomposed.
#'  
#' @return extracted FA substructures
#' @export
FA_sub_extract <- function(char_table, FA_substructure,
                           unprocessed_data_result,
                           exact_FA = "no", exo_lipid = NULL) {
  if (class(unprocessed_data_result)[1] == "list") {
    unprocessed_data_result <- unprocessed_data_result[[2]] %>% filter(type == "FA")
  }
  non_processed_data_result <- unprocessed_data_result %>% filter(lipid %in% FA_substructure$FA)
  if (exact_FA == "no") {
    FA_t_test <- FA_substructure[-1] %>%
      t() %>%
      as.data.frame() %>%
      apply(MARGIN = 2, FUN = function(x) {
        as.data.frame(x) %>%
          filter(x != "") %>%
          mutate(Lipid = str_extract(x, "\\d+:\\d+;\\d+")) %>% # modify
          left_join(non_processed_data_result, by = c("Lipid" = "lipid"))
      })
  } else {
    FA_t_test <- FA_substructure[-1] %>%
      t() %>%
      as.data.frame() %>%
      apply(MARGIN = 2, FUN = function(x) {
        as.data.frame(x) %>%
          filter(x != "") %>%
          mutate(Lipid = x) %>% # modify
          mutate(Lipid = str_replace_all(Lipid, "endo: ", "")) %>%
          mutate(Lipid = str_replace_all(Lipid, "exo: ", "")) %>%
          left_join(non_processed_data_result, by = c("Lipid" = "lipid"))
      })
  }
  
  if (!is.null(exo_lipid)) {
    exo_lipid_with_neighbor <- unlist(filter(FA_network, S1 == exo_lipid | P1 == exo_lipid)[c("S1", "P1")]) %>% unique()
  } else {
    exo_lipid_with_neighbor <- ""
  }
  FA_sub_stop <- list()
  FA_name <- character()
  for (num in 1:length(FA_t_test)) {
    FA_sub <- FA_t_test[[num]]
    FA_name[num] <- last(FA_sub$x)
    
    if (is.na(last(FA_sub$log2FC))) {
      FA_sub_stop[[num]] <- ""
    } else if (last(FA_sub$x) == "16:0;0") {
      FA_sub_stop[[num]] <- FA_sub$x
    } else if (exo_lipid_with_neighbor != "" && last(FA_sub$x) %in% exo_lipid_with_neighbor) {
      FA_sub_stop[[num]] <- last(FA_sub$x)
    } else {
      if (last(FA_sub$log2FC > 0)) {
        stop_point <- which(FA_sub$log2FC < 0)
        if (length(stop_point) == 0) {
          FA_sub_stop[[num]] <- FA_sub$x
        } else if (max(stop_point) == (nrow(FA_sub) - 1)) {
          FA_sub_stop[[num]] <- ""
        } else {
          stop_loc <- max(stop_point) + 1
          FA_sub_stop[[num]] <- FA_sub$x[stop_loc:length(FA_sub$x)]
        }
      } else {
        stop_point <- which(FA_sub$log2FC > 0)
        if (length(stop_point) == 0) {
          FA_sub_stop[[num]] <- FA_sub$x
        } else if (max(stop_point) == (nrow(FA_sub) - 1)) {
          FA_sub_stop[[num]] <- ""
        } else {
          stop_loc <- max(stop_point) + 1
          FA_sub_stop[[num]] <- FA_sub$x[stop_loc:length(FA_sub$x)]
        }
      }
    }
    names(FA_sub_stop) <- FA_name # modify
  }
  FA_sub_stop <- plyr::ldply(FA_sub_stop, rbind)
  FA_sub_stop[is.na(FA_sub_stop)] <- ""
  
  colnames(FA_sub_stop) <- c("FA", str_c("Unit", 1:(ncol(FA_sub_stop) - 1)))
  FA_sub_stop <- unique(FA_sub_stop) %>% filter(Unit1 != "")
  
  exist_FA <- FA_substructure %>%
    filter(FA %in% non_processed_data_result$lipid) %>%
    apply(MARGIN = 1, FUN = function(x) {
      last(x[x != ""])
    }) %>%
    unique()
  FA_sub_exist_FA <- FA_sub_stop %>%
    apply(MARGIN = 1, FUN = function(x) {
      last(x[x != ""])
    }) %>%
    unique()
  lost_FA <- exist_FA[!exist_FA %in% FA_sub_exist_FA]
  
  if (length(lost_FA) != 0) {
    FA_sub_stop <- FA_sub_stop %>% add_row(FA = lost_FA, Unit1 = lost_FA)
    FA_sub_stop[is.na(FA_sub_stop)] <- ""
  }
  if (exact_FA == "no") {
    FA_sub_stop <- FA_sub_stop %>%
      mutate(FA = str_extract(FA, "\\d+:\\d+;\\d+"))
  } else {
    FA_sub_stop <- FA_sub_stop %>%
      mutate(FA = ifelse(str_detect(FA, "(exo)|(endo)"), str_extract(FA, "\\d+:\\d+;\\d+"), FA))
  }
  
  FA_substructure_transform <- function(char_table, FA_substructure) {
    # Filter lipid and FA with substructure
    char <- char_table %>% filter(FA_split != "")
    
    lipid_name <- list()
    lipid_sub_list <- list()
    
    add <- 1
    
    for (var1 in 1:nrow(char)) {
      char[var1, "feature"]
      FA <- str_split(char[var1, "FA_split"], "_") %>% unlist()
      FA_sub <- FA_sub_trans(char[var1, "feature"], FA, FA_substructure)
      
      
      lipid_name[[add]] <- FA_sub[[1]]
      lipid_sub_list[[add]] <- FA_sub[[2]]
      add <- add + 1
    }
    lipid_name <- unlist(lipid_name, recursive = F)
    lipid_sub_list <- unlist(lipid_sub_list, recursive = F)
    
    return(list(lipid_name, lipid_sub_list))
  }
  FA_sub_trans <- function(lipid, FA_list, FA_substructure) {
    add <- 1
    # FA substructure: all possible path method
    all_FA_sub <- list()
    for (FA_num in 1:length(FA_list)) {
      FA_sub_all <- FA_substructure %>% filter(FA == FA_list[FA_num])
      
      if (nrow(FA_sub_all) == 0) {
        FA_sub_all <- data.frame(FA = FA_list[FA_num], Unit1 = FA_list[FA_num])
      }
      all_FA_sub[[FA_num]] <- FA_sub_all
    }
    
    lipid_name <- character()
    FA_sub_list <- list()
    
    for (var1 in 1:nrow(all_FA_sub[[1]])) {
      if (length(all_FA_sub) > 1) {
        for (var2 in 1:nrow(all_FA_sub[[2]])) {
          if (length(all_FA_sub) > 2) {
            for (var3 in 1:nrow(all_FA_sub[[3]])) {
              if (length(all_FA_sub) > 3) {
                for (var4 in 1:nrow(all_FA_sub[[4]])) {
                  FA1_sub <- all_FA_sub[[1]][var1, ] %>%
                    dplyr::select(-FA) %>%
                    unlist()
                  FA1_sub <- FA1_sub[FA1_sub != ""]
                  
                  FA2_sub <- all_FA_sub[[2]][var2, ] %>%
                    dplyr::select(-FA) %>%
                    unlist()
                  FA2_sub <- FA2_sub[FA2_sub != ""]
                  
                  FA3_sub <- all_FA_sub[[3]][var3, ] %>%
                    dplyr::select(-FA) %>%
                    unlist()
                  FA3_sub <- FA3_sub[FA3_sub != ""]
                  
                  FA4_sub <- all_FA_sub[[4]][var4, ] %>%
                    dplyr::select(-FA) %>%
                    unlist()
                  FA4_sub <- FA4_sub[FA4_sub != ""]
                  
                  FA_sub <- c(FA1_sub, FA2_sub, FA3_sub, FA4_sub)
                  
                  names(FA_sub) <- c(
                    rep("FA1", length(FA1_sub)),
                    rep("FA2", length(FA2_sub)),
                    rep("FA3", length(FA3_sub)),
                    rep("FA4", length(FA4_sub))
                  )
                  
                  lipid_name[add] <- lipid
                  FA_sub_list[[add]] <- FA_sub
                  add <- add + 1
                }
              } else {
                FA1_sub <- all_FA_sub[[1]][var1, ] %>%
                  dplyr::select(-FA) %>%
                  unlist()
                FA1_sub <- FA1_sub[FA1_sub != ""]
                
                FA2_sub <- all_FA_sub[[2]][var2, ] %>%
                  dplyr::select(-FA) %>%
                  unlist()
                FA2_sub <- FA2_sub[FA2_sub != ""]
                
                FA3_sub <- all_FA_sub[[3]][var3, ] %>%
                  dplyr::select(-FA) %>%
                  unlist()
                FA3_sub <- FA3_sub[FA3_sub != ""]
                
                
                FA_sub <- c(FA1_sub, FA2_sub, FA3_sub)
                
                names(FA_sub) <- c(
                  rep("FA1", length(FA1_sub)),
                  rep("FA2", length(FA2_sub)),
                  rep("FA3", length(FA3_sub))
                )
                
                lipid_name[add] <- lipid
                FA_sub_list[[add]] <- FA_sub
                add <- add + 1
              }
            }
          } else {
            FA1_sub <- all_FA_sub[[1]][var1, ] %>%
              dplyr::select(-FA) %>%
              unlist()
            FA1_sub <- FA1_sub[FA1_sub != ""]
            FA2_sub <- all_FA_sub[[2]][var2, ] %>%
              dplyr::select(-FA) %>%
              unlist()
            FA2_sub <- FA2_sub[FA2_sub != ""]
            
            FA_sub <- c(FA1_sub, FA2_sub)
            
            names(FA_sub) <- c(
              rep("FA1", length(FA1_sub)),
              rep("FA2", length(FA2_sub))
            )
            
            lipid_name[add] <- lipid
            FA_sub_list[[add]] <- FA_sub
            add <- add + 1
          }
        }
      } else {
        FA_sub <- all_FA_sub[[1]][var1, ] %>%
          dplyr::select(-FA) %>%
          unlist()
        # 20:2;0: 18:1;0,9, 18:2;0,6,9, 20:2;0,8,11
        FA_sub <- FA_sub[FA_sub != ""]
        
        names(FA_sub) <- c(rep("FA1", length(FA_sub)))
        
        lipid_name[add] <- lipid
        FA_sub_list[[add]] <- FA_sub
        add <- add + 1
      }
    }
    return(list(lipid_name, FA_sub_list))
  }
  
  FA_sub_stop <- FA_substructure_transform(char_table, FA_sub_stop)
  
  return(FA_sub_stop)
}
