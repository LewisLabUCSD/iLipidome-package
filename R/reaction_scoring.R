#' evaluates each reaction in biosynthetic network using substructures
#' 
#' @param network Trimmed FA biosynthetic network. Output of "build_FA_net".
#' @param sub_exp Substructure profile. Output of "lipid_sub_matrix".
#' @param sub_t Statistical result table for substructures. Output of "t_test".
#' @param ctrl An integer vector specifying samples in the control group in 
#'  lipid expression table. Note that the first column containing lipid names is
#'  not counted.
#' @param exp An integer vector specifying samples in the experimental group in 
#'  lipid expression table. Note that the first column containing lipid names is
#'  not counted.
#' @param Species "human", "mouse", or "rat" can be used to label 
#'  species-specific genes for lipid reactions.
#'  
#' @return reaction scoring result table
#' 
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom purrr map2_chr
#' @importFrom stats t.test
#' @importFrom stats p.adjust
#' 
#' @export
reaction_scoring <- function(network, sub_exp, sub_t, ctrl = 1:7,
                             exp = 8:13, Species = "human") {
  sub_exp <- sub_exp %>% apply(MARGIN = 1, FUN = function(x) {
    y <- x
    y[y == 0] <- min(y[y != 0]) * 0.5
    y
  }, simplify = T)
  sub_exp <- as.data.frame(t(sub_exp))
  
  
  edge_name <- character(nrow(network))
  
  
  edge_type <- character(nrow(network))
  node1_log2FC <- numeric(nrow(network))
  node2_log2FC <- numeric(nrow(network))
  p_value <- numeric(nrow(network))
  statistics <- numeric(nrow(network))
  FC_ctrl <- numeric(nrow(network))
  FC_exp <- numeric(nrow(network))
  FC_exp_ctrl <- numeric(nrow(network))
  
  num <- 1
  for (edge in 1:nrow(network)) {
    if ((network[edge, 1] %in% rownames(sub_exp)) && (network[edge, 2] %in% rownames(sub_exp))) {
      edge_name[num] <- str_c(network[edge, 1], " --> ", network[edge, 2])
      
      a <- sub_t %>%
        filter(lipid == network[edge, 1]) %>%
        .$log2FC
      b <- sub_t %>%
        filter(lipid == network[edge, 2]) %>%
        .$log2FC
      
      node1_log2FC[num] <- a
      node2_log2FC[num] <- b
      
      if (a < 0 && b > 0) {
        edge_type[num] <- "Increase"
      } else if (a > 0 && b < 0) {
        edge_type[num] <- "Decrease"
      } else {
        edge_type[num] <- "Non-change"
      }
      
      from <- which(rownames(sub_exp) == network[edge, 1])
      to <- which(rownames(sub_exp) == network[edge, 2])
      
      edge_FC <- unlist(sub_exp[to, ]) / unlist(sub_exp[from, ])
      
      FC_ctrl[num] <- mean(edge_FC[ctrl])
      FC_exp[num] <- mean(edge_FC[exp])
      FC_exp_ctrl[num] <- FC_exp[num] / FC_ctrl[num]
      
      p_value[num] <- tryCatch(
        {
          t.test(edge_FC[exp], edge_FC[ctrl], var.equal = T)$p.value
        },
        error = function(e) {
          return(NA)
        }
      )
      statistics[num] <- tryCatch(
        {
          t.test(edge_FC[exp], edge_FC[ctrl], var.equal = T)$statistic
        },
        error = function(e) {
          return(NA)
        }
      )
      
      num <- num + 1
    } else {}
  }
  
  
  result <- data.frame(
    edge_name = edge_name, FC_ctrl = FC_ctrl, FC_exp = FC_exp,
    FC_exp_ctrl = FC_exp_ctrl, statistics = statistics,
    p_value = p_value, edge_type = edge_type,
    node1_log2FC = node1_log2FC, node2_log2FC = node2_log2FC
  ) %>%
    filter(edge_name != "") %>%
    mutate(
      adj_p_value = p.adjust(.$p_value, method = "fdr"),
      log2FC_exp_ctrl = log2(FC_exp_ctrl),
      mlog10p = -log10(p_value), mlog10padj = -log10(adj_p_value)
    )
  
  result <- result %>%
    mutate(perturbation_score = log2FC_exp_ctrl * mlog10p) %>%
    mutate(Mode = ifelse(perturbation_score > 0, "Increase", "Decrease"))
  
  result <- result %>%
    dplyr::select(
      edge_name, FC_ctrl, FC_exp,
      FC_exp_ctrl, statistics,
      p_value, mlog10p,
      adj_p_value, mlog10padj,
      perturbation_score, Mode, edge_type,
      node1_log2FC, node2_log2FC
    ) %>%
    mutate(perturbation_score = ifelse(is.na(perturbation_score), 0, perturbation_score)) %>%
    arrange(desc(perturbation_score))
  
  
  edge_lipid1 <- str_split(result$edge_name, " --> ") %>% map_chr(~ .[1])
  edge_lipid2 <- str_split(result$edge_name, " --> ") %>% map_chr(~ .[2])
  
  
  FA_change <- map2_chr(edge_lipid1, edge_lipid2, .f = function(x, y) {
    a <- str_extract_all(x, "\\d+:\\d+;\\d+") %>% unlist()
    if (length(a) == 0) {
      a <- "0:0;0"
    }
    a1 <- str_extract_all(a, "\\d+:") %>%
      str_sub(end = -2) %>%
      as.integer() %>%
      sum()
    a2 <- str_extract_all(a, "\\d+;") %>%
      str_sub(end = -2) %>%
      as.integer() %>%
      sum()
    a3 <- str_extract_all(a, "\\d+$") %>%
      as.integer() %>%
      sum()
    b <- str_extract_all(y, "\\d+:\\d+;\\d+") %>% unlist()
    if (length(b) == 0) {
      b <- "0:0;0"
    }
    b1 <- str_extract_all(b, "\\d+:") %>%
      str_sub(end = -2) %>%
      as.integer() %>%
      sum()
    b2 <- str_extract_all(b, "\\d+;") %>%
      str_sub(end = -2) %>%
      as.integer() %>%
      sum()
    b3 <- str_extract_all(b, "\\d+$") %>%
      as.integer() %>%
      sum()
    return(str_c(" (", abs(a1 - b1), ":", abs(a2 - b2), ";", abs(a3 - b3), ")"))
  })
  
  FA_change[FA_change == " (0:0;0)"] <- ""
  
  result <- result %>% mutate(FA_change = str_replace(FA_change, ";0", ""))
  
  
  lipid_reaction <- str_split(result$edge_name, " --> ") %>%
    map(~ str_replace_all(.x, "endo: ", "") %>%
          str_replace_all("exo: ", "") %>%
          str_replace_all("2:0;0--", "")) %>%
    map_chr(.f = function(x) {
      a <-
        if (str_detect(x[1], "_")) {
          str_extract(x[1], "[A-Za-z -]+")
        } else {
          x[1]
        }
      b <- if (str_detect(x[2], "_")) {
        str_extract(x[2], "[A-Za-z -]+")
      } else {
        x[2]
      }
      str_c(a, b, sep = "_")
    })
  
  genes <- data.frame(Reaction = lipid_reaction) %>%
    left_join(filter(reaction_gene_mapping, species == Species), by = "Reaction") %>%
    .$gene
  
  result <- result %>% mutate(genes = genes)
  
  
  return(result)
}
