#' decomposes FAs into FA substructures based on the FA biosynthetic network.
#' 
#' @param FA_network Trimmed FA biosynthetic network. Output of "build_FA_net".
#' @param unprocessed_data_result Differential expression for unprocessed 
#'  lipidomics data. Output of "unprocessed_data_test".
#' @param unmapped_FA An character vector allowing users to decide which FAs 
#'  should be ignored since some FAs can be mapped to more than one node in FA 
#'  network (e.g., 20:4;0 for w6-20:4;0 and w3-20:4;0). 
#'  
#' @return FA substructure table
#' @export
FA_sub_transform <- function(FA_network, unprocessed_data_result,
                             unmapped_FA = NULL) {
  non_processed_data_result <- unprocessed_data_result
  all_FA <- non_processed_data_result[[2]] %>%
    filter(type == "FA") %>%
    .$lipid %>%
    unique()
  
  endo_start <- FA_network %>% filter(str_detect(S1, "endo: 2:0;0"))
  exo <- FA_network %>%
    filter(str_detect(S1, "exo: ")) %>%
    .[1, ]
  
  graph <- graph_from_data_frame(FA_network[c("S1", "P1")],
                                 directed = T,
                                 vertices = unique(c(FA_network$S1, FA_network$P1))
  )
  
  
  De_novo <- FA_network %>%
    filter(pathway == "Non_essential_FA_synthesis") %>%
    unlist()
  Omega6 <- FA_network %>%
    filter(pathway == "Omega_6_FA_synthesis") %>%
    unlist()
  Omega3 <- FA_network %>%
    filter(pathway == "Omega_3_FA_synthesis") %>%
    unlist()
  
  lipid <- character()
  G3P_start_path <- list()
  lipid_list <- unique(c(FA_network$S1, FA_network$P1))
  num_class <- 1
  num_path <- 1
  
  while (!is.na(lipid_list[num_class])) {
    if (str_detect(lipid_list[num_class], "exo|w6-18:2;0|w3-18:3;0")) {
      lipid[num_path] <- lipid_list[num_class]
      G3P_start_path[[num_path]] <- lipid_list[num_class]
      num_path <- num_path + 1
      num_class <- num_class + 1
      next
    }
    if (lipid_list[num_class] %in% De_novo) {
      start <- c(exo[1, 1], endo_start[1, 1])
    } else if (lipid_list[num_class] %in% Omega6) {
      start <- "w6-18:2;0"
    } else if (lipid_list[num_class] %in% Omega3) {
      start <- "w3-18:3;0"
    }
    for (starts in start) {
      shortest_path <- tryCatch(
        {
          all_simple_paths(graph, starts,
                           lipid_list[num_class],
                           mode = "out"
          )
        },
        error = function(e) {
          NULL
        }
      )
      
      
      if (!is.null(shortest_path)) {
        if (length(shortest_path) != 0) {
          for (path_num in 1:length(shortest_path)) {
            lipid[num_path] <- lipid_list[num_class]
            G3P_start_path[[num_path]] <- shortest_path[[path_num]] %>%
              attributes() %>%
              .$names
            num_path <- num_path + 1
          }
        } else {
          lipid[num_path] <- lipid_list[num_class]
          G3P_start_path[[num_path]] <- ""
          num_path <- num_path + 1
        }
      } else {
        lipid[num_path] <- lipid_list[num_class]
        G3P_start_path[[num_path]] <- ""
        num_path <- num_path + 1
      }
    }
    num_class <- num_class + 1
  }
  
  names(G3P_start_path) <- lipid
  FA_substructure <- plyr::ldply(G3P_start_path, rbind) %>% unique()
  
  FA_substructure[is.na(FA_substructure)] <- ""
  
  colnames(FA_substructure) <- c("FA", str_c("Unit", 1:(ncol(FA_substructure) - 1)))
  
  
  FA_substructure <- FA_substructure %>% filter(Unit1 != "")
  FA_substructure$FA <- str_extract(FA_substructure$FA, "\\d+:.+")
  
  special_FA <- all_FA[all_FA != ""]
  special_FA <- special_FA[!special_FA %in% FA_substructure$FA]
  
  FA_substructure <- FA_substructure %>% bind_rows(data.frame(FA = special_FA, Unit1 = special_FA))
  
  FA_substructure[is.na(FA_substructure)] <- ""
  
  if (!is.null(unmapped_FA)) {
    rm_FA <- which(apply(FA_substructure, MARGIN = 1, FUN = function(x) {
      last(x[x != ""])
    }) %in% unmapped_FA)
    
    FA_substructure <- FA_substructure[-rm_FA, ]
  }
  
  return(FA_substructure)
}
