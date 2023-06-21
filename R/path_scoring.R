#' uses substructures to score pathways in biosynthetic network
#' 
#' @param network Trimmed FA biosynthetic network. Output of "build_FA_net".
#' @param sub_t Statistical result table for substructures. Output of "t_test".
#' @param calibrate "True" or "False" to calibrate pathway scores or not.
#' @param data_type "FA", "Class", or "Species" to indicate which type of 
#'  substructure.
#'  
#' @return pathway scoring result table
#' 
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr everything
#' @importFrom gtools permutations
#' @importFrom purrr map_dbl
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom utils tail
#' 
#' 
#' @export
path_scoring <- function(network, sub_t, calibrate = T, data_type = "FA") {
  options(scipen = 999)
  
  sub_t <- sub_t %>% dplyr::select(-sig)
  sub_t <- sub_t %>%
    mutate(zscore = qnorm(1 - sub_t$p_value / 2)) %>%
    mutate(zscore = ifelse(is.infinite(zscore), 8, zscore)) %>%
    mutate(zscore = ifelse(statistics > 0, zscore, -zscore))
  
  Network_w <- data.frame(node = unique(unlist(network[c("S1", "P1")]))) %>%
    left_join(sub_t, by = c("node" = "lipid")) %>%
    mutate(weight = zscore) %>%
    dplyr::select(node, weight)
  
  graph <- graph_from_data_frame(network[c("S1", "P1")], directed = T, vertices = unique(c(network$S1, network$P1)))
  
  all_node <- network[c("S1", "P1")] %>%
    unlist() %>%
    unique()
  
  combo <- permutations(length(all_node), 2, all_node)
  path <- character()
  path_full <- character()
  score <- numeric()
  type <- character()
  from <- character()
  to <- character()
  
  num_path <- 1
  for (num in 1:nrow(combo)) {
    if ((combo[num, 1] %in% sub_t$lipid) && (combo[num, 2] %in% sub_t$lipid)) {
      all_path <- all_simple_paths(graph, combo[num, 1], combo[num, 2], mode = "out")
      
      
      if (length(all_path) != 0) {
        for (path_num in 1:length(all_path)) {
          path_node <- all_path[[path_num]] %>%
            attributes() %>%
            .$names
          # c(2:0;0, 4:0;0, 6:0;0, 8:0;0, 10:0;0, 12:0;0)
          
          path[num_path] <- str_c(path_node, collapse = " --> ")
          # c(2:0;0_4:0;0_6:0;0_8:0;0_10:0;0_12:0;0)
          
          from[num_path] <- path_node[1]
          to[num_path] <- tail(path_node, 1)
          
          S1_P1_w <- Network_w %>%
            filter(node %in% path_node) %>%
            .$weight
          score[num_path] <- sum(S1_P1_w) / sqrt(length(path_node))
          
          num_path <- num_path + 1
        }
      }
    }
  }
  path_score <- data.frame(path = path, from = from, to = to, score = score) %>%
    filter(!is.na(score)) %>%
    mutate(Type = ifelse(score > 0, "Active", "Suppressed")) %>%
    arrange(score)
  
  
  if (calibrate == T) {
    all_z <- sub_t$zscore
    max_length <- path_score$path %>%
      str_split(" --> ") %>%
      map(length) %>%
      unlist() %>%
      max()
    
    path <- numeric()
    mean_each_length <- numeric()
    sd_each_length <- numeric()
    
    for (path_length in 2:max_length) {
      path[path_length - 1] <- path_length
      
      cal <- replicate(10000, sample(all_z, path_length, replace = F)) %>%
        colSums()
      cal <- cal / sqrt(path_length)
      
      mean_each_length[path_length - 1] <- mean(cal)
      sd_each_length[path_length - 1] <- sd(cal)
    }
    
    cal_score <- numeric()
    for (each_path in 1:nrow(path_score)) {
      path_length <- path_score$path[each_path] %>%
        str_split(" --> ") %>%
        unlist() %>%
        length()
      cal_score[each_path] <- (path_score$score[each_path] - mean_each_length[path_length - 1]) / sd_each_length[path_length - 1]
    }
    path_score <- path_score %>%
      mutate(cal_score = cal_score) %>%
      dplyr::select(path, path, from, to, score, cal_score, everything())
    
    path_score <- path_score %>%
      mutate(Significant = ifelse(abs(cal_score) > 1.96, "yes", "no")) %>%
      arrange(desc(cal_score)) %>%
      mutate(Type = ifelse(cal_score > 0, "Active", "Suppressed"))
  } else {
    path_score <- path_score %>%
      mutate(Significant = ifelse(abs(score) > 1.96, "yes", "no")) %>%
      arrange(desc(score))
  }
  
  if (data_type == "FA") {
    rep_sub <- function(path_score, type) {
      path_score <- path_score %>% filter(Type == type)
      
      if (nrow(path_score) == 0) {
        return(NULL)
      }
      score_list <- path_score$path
      score_list_name <- path_score$path
      
      rep_sub_path_list <- list()
      rep_sub_path_num <- 1
      rep_sub_path_name <- character()
      # path_name <- character()
      
      
      for (num in 1:length(score_list)) {
        sub_path <- unlist(str_split(score_list[num], " --> "))
        sub_path <- str_c(sub_path[-length(sub_path)], sub_path[-1], sep = "_")
        overlap_prop <- rep_sub_path_list %>% map_dbl(.f = function(x) {
          sum(sub_path %in% x) / length(sub_path)
        })
        
        if (sum(overlap_prop > 0.5) == 0) {
          # path_name[num] <- score_list_name[num]
          
          rep_sub_path_list[[rep_sub_path_num]] <- sub_path
          names(rep_sub_path_list)[rep_sub_path_num] <- score_list_name[num]
          
          rep_sub_path_name[num] <- score_list_name[num]
          rep_sub_path_num <- rep_sub_path_num + 1
        } else {
          rep_sub_path_name[num] <- names(rep_sub_path_list)[which.max(overlap_prop)]
          # path_name[num] <- score_list_name[num]
        }
      }
      return(rep_sub_path_name)
    }
    
    rep_sub_path <- c(rep_sub(path_score, "Active"), rev(rep_sub(path_score[nrow(path_score):1, ], "Suppressed")))
    path_score <- path_score %>% mutate(rep_sub_path = rep_sub_path)
  } else if (data_type == "Class") {
    rep_sub <- function(path_score, type) {
      path_score <- path_score %>% filter(Type == type)
      if (nrow(path_score) == 0) {
        return(NULL)
      }
      score_list <- path_score$path
      score_list_name <- path_score$path
      
      rep_sub_path_list <- list()
      rep_sub_path_num <- 1
      rep_sub_path_name <- character()
      # path_name <- character()
      
      
      for (num in 1:length(score_list)) {
        sub_path <- unlist(str_split(score_list[num], " --> "))
        
        sub_path <- str_c(sub_path[-length(sub_path)], sub_path[-1], sep = "_")
        
        overlap_prop <- rep_sub_path_list %>% map_dbl(.f = function(x) {
          sum(sub_path %in% x) / length(sub_path)
        })
        
        if (sum(overlap_prop >= 0.5) == 0) {
          # path_name[num] <- score_list_name[num]
          
          rep_sub_path_list[[rep_sub_path_num]] <- sub_path
          names(rep_sub_path_list)[rep_sub_path_num] <- score_list_name[num]
          
          rep_sub_path_name[num] <- score_list_name[num]
          rep_sub_path_num <- rep_sub_path_num + 1
        } else {
          rep_sub_path_name[num] <- names(rep_sub_path_list)[which.max(overlap_prop)]
          # path_name[num] <- score_list_name[num]
        }
      }
      return(rep_sub_path_name)
    }
    
    rep_sub_path <- c(rep_sub(path_score, "Active"), rev(rep_sub(path_score[nrow(path_score):1, ], "Suppressed")))
    path_score <- path_score %>% mutate(rep_sub_path = rep_sub_path)
  } else {
    if (sum(str_detect(path_score$path, "\\d+:\\d+;0_\\d+:\\d+;0")) != 0) {
      rep_sub_path <- path_score$path %>%
        str_extract_all("\\d+:\\d+;0_\\d+:\\d+;0") %>%
        map(.f = function(x) {
          k <- table(x) %>% sort(decreasing = T)
          if (length(k) == 0) {
            ""
          } else if (sum(k != 1) == 0) {
            sort(names(k), decreasing = T)[1]
          } else if (sum(k != 1) != 0) {
            names(k)[1]
          }
        }) %>%
        unlist()
    } else {
      rep_sub_path <- path_score$path %>%
        str_extract_all("\\d+:\\d+;\\d+") %>%
        map(.f = function(x) {
          k <- table(x) %>% sort(decreasing = T)
          if (length(k) == 0) {
            ""
          } else if (sum(k != 1) == 0) {
            sort(names(k), decreasing = T)[1]
          } else if (sum(k != 1) != 0) {
            sort(names(k[k == max(k)]), decreasing = T)[1]
          }
        }) %>%
        unlist()
    }
    
    
    
    path_score <- path_score %>%
      mutate(rep_sub_path = rep_sub_path)
  }
  return(path_score)
}
