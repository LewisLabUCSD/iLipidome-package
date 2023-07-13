#' maps species substructures in each pathway with fold changes from the 
#'  “unprocessed_data_test” result and extracts them
#' 
#' @param lipid_substructure Lipid species substructure table. Output of 
#'  "species_sub_transform".
#' @param unprocessed_data_result Differential expression for unprocessed 
#'  lipidomics data. Output of "unprocessed_data_test".
#' @param type "species" or "class" to indicate which type of substructure.
#' @param pct_limit A threshold of non-NA percent (0~1) in one biosynthetic 
#'  route can be set to control decompostion process and reduce artifacts. If 
#'  the proportion of non-missing substructures exceeds the threshold in one 
#'  biosynthetic route, the target lipid species will not be decomposed to 
#'  substructures.
#' @param exo_lipid An character vector specifying exogenous lipid addition. 
#'  The exogenous lipids and their adjacent nodes in FA network will not be 
#'  decomposed.
#'  
#' @return lipid species
#' 
#' @importFrom dplyr arrange
#' @importFrom dplyr everything
#' @importFrom dplyr filter
#' @importFrom dplyr last
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom plyr ldply
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' 
#' @export
species_sub_extract <- function(lipid_substructure, unprocessed_data_result,
                                type = "species", pct_limit = 0.3,
                                exo_lipid = NULL) {
  if (class(unprocessed_data_result)[1] == "list") {
    unprocessed_data_result <- unprocessed_data_result[[2]] %>% filter(type == type)
  }
  
  non_pro <- unprocessed_data_result
  # non_pro <- non_processed_data_result %>% filter(lipid %in% lipid_substructure$Lipid)
  lipid_substructure <- lipid_substructure %>% filter(Lipid %in% non_pro$lipid)
  
  lipid_test <- lipid_substructure %>%
    mutate(num = 1:nrow(.)) %>%
    gather(-Lipid, -num, key = "Unit", value = "sub") %>%
    mutate(Unit = factor(Unit, levels = unique(Unit))) %>%
    mutate(sub = str_replace(sub, "_FA\\d", "")) %>%
    left_join(non_pro[c("lipid", "log2FC")], by = c("sub" = "lipid")) %>%
    .[-4] %>%
    spread(key = "Unit", value = "log2FC") %>%
    arrange(num) %>%
    dplyr::select(-num)
  
  lipid_sub_trans <- list(nrow(lipid_test))
  lipid_name <- character(nrow(lipid_test))
  for (num in 1:nrow(lipid_test)) {
    # print(num)
    lipid_sub <- lipid_substructure[-1][num, ] %>%
      unlist() %>%
      str_replace("_FA\\d", "")
    lipid_sub_FC <- lipid_test[-1][num, ] %>% unlist()
    
    lipid_sub_FC <- lipid_sub_FC[lipid_sub != ""]
    lipid_sub <- lipid_sub[lipid_sub != ""]
    
    lipid_name[num] <- last(lipid_sub)
    
    Not_NA_pct <- sum(!is.na(lipid_sub_FC)) / length(lipid_sub_FC)
    if (Not_NA_pct < pct_limit) {
      lipid_sub_trans[[num]] <- ""
    } else {
      if (last(lipid_sub_FC) > 0) {
        stop_point <- which(lipid_sub_FC < 0)
        if (length(stop_point) == 0) {
          lipid_sub_trans[[num]] <- lipid_sub
        } else if (max(stop_point) == (length(lipid_sub) - 1)) {
          lipid_sub_trans[[num]] <- ""
        } else {
          stop_loc <- max(which(lipid_sub_FC < 0)) + 1
          lipid_sub_trans[[num]] <- lipid_sub[stop_loc:length(lipid_sub)]
        }
      } else {
        stop_point <- which(lipid_sub_FC > 0)
        if (length(stop_point) == 0) {
          lipid_sub_trans[[num]] <- lipid_sub
        } else if (max(stop_point) == (length(lipid_sub) - 1)) {
          lipid_sub_trans[[num]] <- ""
        } else {
          stop_loc <- max(which(lipid_sub_FC > 0)) + 1
          lipid_sub_trans[[num]] <- lipid_sub[stop_loc:length(lipid_sub)]
        }
      }
    }
  }
  names(lipid_sub_trans) <- lipid_name
  
  lipid_sub_trans <- plyr::ldply(lipid_sub_trans, rbind)
  lipid_sub_trans[is.na(lipid_sub_trans)] <- ""
  
  
  colnames(lipid_sub_trans) <- c("Lipid", str_c("Unit", 1:(ncol(lipid_sub_trans) - 1)))
  lipid_sub_trans <- unique(lipid_sub_trans) %>% filter(Unit1 != "")
  
  
  if (type == "class") {
    exist_lipid <- non_pro %>%
      filter(type == "class") %>%
      .$lipid
    lost_lipid <- exist_lipid[!exist_lipid %in% lipid_sub_trans$Lipid]
    lipid_sub_trans <- lipid_sub_trans %>% add_row(
      Lipid = lost_lipid,
      Unit1 = lost_lipid
    )
    lipid_sub_trans[is.na(lipid_sub_trans)] <- ""
  } else if (type == "species") {
    exist_lipid <- non_pro %>%
      filter(type == "species") %>%
      .$lipid
    lost_lipid <- exist_lipid[!exist_lipid %in% lipid_sub_trans$Lipid]
    lipid_sub_trans <- lipid_sub_trans %>% add_row(
      Lipid = lost_lipid,
      Unit1 = lost_lipid
    )
    lipid_sub_trans[is.na(lipid_sub_trans)] <- ""
  }
  
  species_sub_stop <- list(lipid_sub_trans[[1]], apply(lipid_sub_trans[-1], MARGIN = 1, FUN = function(x) {
    x[x != ""]
  }))
  if (!is.null(exo_lipid)) {
    exo_lipid_with_neighbor <- c(
      which(map_chr(species_sub_stop[[2]], ~ last(.)) %in% exo_lipid),
      which(map_chr(species_sub_stop[[2]], ~ last(., 2)[1]) %in% exo_lipid)
    ) %>% unique()
    
    species_sub_stop[[2]][exo_lipid_with_neighbor] <-
      species_sub_stop[[2]][exo_lipid_with_neighbor] %>% map(.f = function(x) {
        y <- last(x)
        names(y) <- "Unit1"
        y
      })
  }
  return(species_sub_stop)
}
