#' Fatty acid biosynthetic network transformation
#' 
#' @param FA_network A data.frame describing the reference FA biosynthetic 
#'  network in iLipdiome.
#' @param unprocessed_data_result FA_network: A data.frame describing the 
#'  reference FA biosynthetic network in iLipdiome.
#'  
#' @return transformed fatty acid network
#' @export
build_FA_net <- function(FA_network, unprocessed_data_result) {
  non_processed_data_result <- unprocessed_data_result
  all_FA <- non_processed_data_result[[2]] %>%
    filter(type == "FA") %>%
    .$lipid %>%
    unique()
  
  FA_below_16 <- min(which(FA_network$S1[1:7] %in% all_FA))
  
  if (is.infinite(FA_below_16)) {
    FA_network <- FA_network[-c(1:7), ]
    endo_start <- data.frame(
      S1 = str_c("endo: 2:0;0--14:0;0"),
      P1 = "16:0;0",
      S1_detail = "endo: 2:0;0--14:0;0",
      P1_detail = "16:0;0", pathway = "Non_essential_FA_synthesis"
    )
    FA_network <- FA_network %>%
      rbind(endo_start)
    exo <- ""
  } else {
    exo <- FA_network[FA_below_16:7, ]
    exo$S1 <- str_c("exo: ", exo$S1)
    exo$P1 <- str_c("exo: ", exo$P1)
    exo$S1_detail <- str_c("exo: ", exo$S1_detail)
    exo$P1_detail <- str_c("exo: ", exo$P1_detail)
    endo <- data.frame(lapply(exo, function(x) {
      gsub("exo", "endo", x)
    }))
    endo_start <- data.frame(
      S1 = str_c("endo: 2:0;0--", 2 * (FA_below_16 - 1), ":0;0"),
      P1 = endo[1, 1],
      S1_detail = str_c("endo: 2:0;0--", 2 * (FA_below_16 - 1), ":0;0"),
      P1_detail = endo[1, 1], pathway = "Non_essential_FA_synthesis"
    )
    
    FA_network <- FA_network[-c(1:7), ] %>%
      rbind(exo) %>%
      rbind(endo) %>%
      rbind(endo_start)
    FA_network <- data.frame(lapply(FA_network, function(x) {
      gsub("exo: 16:0;0|endo: 16:0;0", "16:0;0", x)
    }))
  }
  
  return(FA_network)
}
