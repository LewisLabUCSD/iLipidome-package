#' build biosynthetic network using substructures, pathway and reaction scoring 
#'  results
#' 
#' @param network_data Trimmed FA biosynthetic network. Output of "build_FA_net"
#' @param DE_data Statistical result table for substructures. Output of "t_test"
#' @param if_species "True" or "False" to indicate data type is lipid species or
#'  the others (FA or class).
#' @param significant "p_value" or "adj_p_value" can be used for selecting 
#'  significant lipid features.
#' @param path_scoring_result Pathway scoring result table. Output of 
#'  "path_scoring".
#' @param reaction_scoring_result Reaction scoring result table. Output of 
#'  "reaction_scoring".
#' @param top_n 1 to 5 to label top N pathways and reactions.
#' @param path_type "active", "suppressed", or "both" is used to select specific
#'  types of top N pathways and reactions.
#'  
#' @return biosynthetic network node
#' @export
draw_network <- function(network_data, DE_data, if_species = F, 
                         significant = "p_value", path_scoring_result, 
                         reaction_scoring_result, top_n = 3, 
                         path_type = "both") {
  if (path_type == "both") {
    active <- T
    suppressed <- T
  } else if (path_type == "active") {
    active <- T
    suppressed <- F
  } else {
    active <- F
    suppressed <- T
  }
  FA_node <- data.frame(lipid = unique(c(network_data$S1, network_data$P1)), stringsAsFactors = F) %>%
    left_join(DE_data, by = c("lipid"))
  
  
  replace_inf <- max(abs(FA_node$log2FC[!is.infinite(FA_node$log2FC)]), na.rm = T) + 0.5
  
  if (significant == "p_value") {
    FA_node <- FA_node %>% mutate(sig = ifelse(p_value < 0.05, "yes", "no"))
  } else {
    FA_node <- FA_node %>%
      mutate(sig = ifelse(adj_p_value < 0.05, "yes", "no")) %>%
      mutate(mlog10p = -log10(adj_p_value))
  }
  
  FA_node <- FA_node %>%
    mutate(in_ref = ifelse(!is.na(FC), "yes", "no")) %>%
    mutate(sig = ifelse(is.na(sig), "na", sig)) %>%
    mutate(log2FC = ifelse(is.infinite(log2FC) & log2FC > 0, replace_inf, log2FC)) %>%
    mutate(log2FC = ifelse(is.infinite(log2FC) & log2FC < 0, -replace_inf, log2FC)) %>%
    mutate(mlog10p = ifelse(is.na(mlog10p), 0, mlog10p)) %>%
    mutate(mean_all = ifelse(is.na(mean_all), 0, mean_all))
  
  FC_color <- cut(FA_node$log2FC, breaks = seq(-replace_inf - 0.1, replace_inf + 0.1, length.out = 100), labels = bluered(99)) %>%
    as.character()
  
  
  
  FC_color[is.na(FC_color)] <- "black"
  
  
  node <- data.frame(
    id = FA_node$lipid,
    label = FA_node$lipid,
    color.background = FC_color,
    color.border = recode(FA_node$sig,
                          "yes" = "purple",
                          "no" = "black",
                          "na" = "transparent"
    ),
    borderWidth = recode(FA_node$sig,
                         "yes" = 2.5,
                         "no" = 1,
                         "na" = 0
    ),
    # value=FA_node$mlog10p,
    value = (FA_node$mlog10p),
    font.size = 35
  )
  edge <- data.frame(
    from = network_data$S1, to = network_data$P1,
    # label = paste("Edge", 1:8),
    color = "gray",
    arrows = c("to"),
    length = 100
  )
  
  path_color <- function(network_edge, path_scoring_result,
                         reaction_scoring_result, top_n = 3,
                         active = T, suppressed = T) {
    if (active == T) {
      # top paths
      
      sig_active_path <- path_scoring_result %>%
        filter(Type == "Active") %>%
        .[!duplicated(.$rep_sub_path), ] %>%
        filter(Significant == "yes")
      if (nrow(sig_active_path) == 0) {
        top_n_new <- nrow(sig_active_path)
        topN_active_path <- data.frame()
      } else if (nrow(sig_active_path) < top_n) {
        top_n_new <- nrow(sig_active_path)
        topN_active_path <- sig_active_path[1:top_n_new, ] %>%
          mutate(rank = 1:top_n_new)
      } else {
        top_n_new <- top_n
        topN_active_path <- sig_active_path[1:top_n_new, ] %>%
          mutate(rank = 1:top_n_new)
      }
      
      
      # top edges
      
      sig_active_edge <- reaction_scoring_result %>%
        filter(Mode == "Increase") %>%
        filter(p_value < 0.05)
      
      if (nrow(sig_active_edge) == 0) {
        topN_active_edge <- data.frame()
      } else if (nrow(sig_active_edge) < top_n) {
        top_n_new <- nrow(sig_active_edge)
        topN_active_edge <- sig_active_edge[1:top_n_new, ] %>%
          mutate(rank = 1:top_n_new)
      } else {
        top_n_new <- top_n
        topN_active_edge <- sig_active_edge[1:top_n_new, ] %>%
          mutate(rank = 1:top_n_new)
      }
    } else {
      topN_active_path <- data.frame()
      topN_active_edge <- data.frame()
    }
    
    if (suppressed == T) {
      # top paths
      
      sig_suppressed_path <- path_scoring_result %>%
        filter(Type == "Suppressed") %>%
        arrange(cal_score) %>%
        .[!duplicated(.$rep_sub_path), ] %>%
        filter(Significant == "yes")
      if (nrow(sig_suppressed_path) == 0) {
        topN_suppressed_path <- data.frame()
      } else if (nrow(sig_suppressed_path) < top_n) {
        top_n_new <- nrow(sig_suppressed_path)
        topN_suppressed_path <- sig_suppressed_path[1:top_n_new, ] %>%
          mutate(rank = c(rep(5, top_n_new) + 1:top_n_new))
      } else {
        top_n_new <- top_n
        topN_suppressed_path <- sig_suppressed_path[1:top_n_new, ] %>%
          mutate(rank = c(rep(5, top_n_new) + 1:top_n_new))
      }
      
      
      # top edges
      
      sig_suppressed_edge <- reaction_scoring_result %>%
        filter(Mode == "Decrease") %>%
        arrange(perturbation_score) %>%
        filter(p_value < 0.05)
      if (nrow(sig_suppressed_edge) == 0) {
        topN_suppressed_edge <- data.frame()
      } else if (nrow(sig_suppressed_edge) < top_n) {
        top_n_new <- nrow(sig_suppressed_edge)
        topN_suppressed_edge <- sig_suppressed_edge[1:top_n_new, ] %>%
          mutate(rank = c(rep(5, top_n_new) + 1:top_n_new))
      } else {
        top_n_new <- top_n
        topN_suppressed_edge <- sig_suppressed_edge[1:top_n_new, ] %>%
          mutate(rank = c(rep(5, top_n_new) + 1:top_n_new))
      }
    } else {
      topN_suppressed_path <- data.frame()
      topN_suppressed_edge <- data.frame()
    }
    
    
    topN_path <- rbind(topN_active_path, topN_suppressed_path)
    
    topN_edge <- rbind(topN_active_edge, topN_suppressed_edge)
    
    path <- data.frame(
      rank = 1:10, edge_color = bluered(100)[c(seq(100, 56, length.out = 5), seq(1, 45, length.out = 5))],
      width = rep(seq(12, 4, length.out = 5), 2)
    )
    
    #
    edge <- data.frame(
      rank = 1:10, perturb_type = c(rep("Increase", 5), rep("Decrease", 5)),
      fontsize = rep(rep(40, 5), 2),
      fontcolor = c(rep("red", 5), rep("blue", 5))
    )
    if (nrow(topN_path) != 0) {
      path_score <- topN_path %>% left_join(path, by = "rank")
    } else {
      path_score <- data.frame()
    }
    if (nrow(topN_edge) != 0) {
      perturbation_score <- topN_edge %>% left_join(edge, by = "rank")
    } else {
      perturbation_score <- data.frame()
    }
    
    edge_color <- character()
    width <- numeric()
    perturb_type <- character()
    fontsize <- numeric()
    fontcolor <- character()
    
    for (edge_num in 1:nrow(network_edge)) {
      from <- network_edge$from[edge_num]
      to <- network_edge$to[edge_num]
      
      edge_in_path <- str_detect(path_score$path, str_c(from, " --> ", to))
      
      if (sum(edge_in_path) == 0) {
        edge_color[edge_num] <- "black"
        width[edge_num] <- 1
      } else {
        names(edge_in_path) <- path_score$rank
        path_rank <- names(edge_in_path)[which.max(edge_in_path)]
        edge_color[edge_num] <- filter(path_score, rank == path_rank)$edge_color
        width[edge_num] <- filter(path_score, rank == path_rank)$width
      }
      
      edge_in_perturb <- perturbation_score$edge_name == str_c(from, " --> ", to)
      
      if (sum(edge_in_perturb) == 0) {
        perturb_type[edge_num] <- ""
        fontsize[edge_num] <- 0
        fontcolor[edge_num] <- "black"
      } else {
        names(edge_in_perturb) <- perturbation_score$rank
        perturb_rank <- names(edge_in_perturb)[which.max(edge_in_perturb)]
        perturb_type[edge_num] <- filter(perturbation_score, rank == perturb_rank)$perturb_type
        fontsize[edge_num] <- filter(perturbation_score, rank == perturb_rank)$fontsize
        fontcolor[edge_num] <- filter(perturbation_score, rank == perturb_rank)$fontcolor
      }
    }
    
    network_edge <- network_edge %>% mutate(
      color = edge_color,
      width = width,
      label = perturb_type,
      font.size = fontsize,
      font.color = fontcolor
    )
    return(network_edge)
  }
  
  if (if_species == T) {
    if (active == T && suppressed == T) {
      topN_rep_path <- c(
        path_scoring_result %>% filter(Type == "Active") %>%
          .[!duplicated(.$rep_sub_path), ] %>% .[1:top_n, ] %>% .$rep_sub_path,
        path_scoring_result %>% filter(Type == "Suppressed") %>%
          arrange(cal_score) %>% .[!duplicated(.$rep_sub_path), ] %>%
          .[1:top_n, ] %>% .$rep_sub_path
      )
    } else if (active == T && suppressed == F) {
      topN_rep_path <- path_scoring_result %>%
        filter(Type == "Active") %>%
        .[!duplicated(.$rep_sub_path), ] %>%
        .[1:top_n, ] %>%
        .$rep_sub_path
    } else if (active == F && suppressed == T) {
      topN_rep_path <- path_scoring_result %>%
        filter(Type == "Suppressed") %>%
        arrange(cal_score) %>%
        .[!duplicated(.$rep_sub_path), ] %>%
        .[1:top_n, ] %>%
        .$rep_sub_path
    }
    
    topN_rep_path <- path_scoring_result %>%
      filter(rep_sub_path %in% topN_rep_path, Significant == "yes") %>%
      .$path %>%
      str_split(" --> ") %>%
      unlist() %>%
      unique()
    
    edge_in_topN_path <- reaction_scoring_result$edge_name %>%
      str_split(" --> ") %>%
      map_lgl(.f = function(x) {
        x[1] %in% topN_rep_path && x[2] %in% topN_rep_path
      })
    edge_in_topN_path <- reaction_scoring_result[edge_in_topN_path, ]$edge_name
    
    topN_net_edge <- network_data %>% filter(S1 %in% topN_rep_path, P1 %in% topN_rep_path)
    
    topN_net_edge <- topN_net_edge %>% mutate(color = "gray", arrows = "to", length = 100)
    
    colnames(topN_net_edge)[1:2] <- c("from", "to")
    
    edge <- topN_net_edge
    
    topN_net_node <- node %>%
      filter(id %in% unlist(topN_net_edge))
    
    
    reaction_scoring_result <- reaction_scoring_result %>%
      filter(edge_name %in% edge_in_topN_path)
  }
  
  edge <- path_color(
    edge, path_scoring_result,
    reaction_scoring_result, top_n,
    active, suppressed
  )
  
  if (if_species == T) {
    node <- topN_net_node %>% filter(id %in% c(edge$from, edge$to))
  }
  return(list(node, edge))
}
