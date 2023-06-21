#' Converts expression of FAs to expression of FA substructures.
#' 
#' @param exp_data Lipid expression table. Output of "build_char_table".
#' @param sub_data Extracted substructure data. Output of "FA_sub_extract".
#' @param sub_type "FA", "Class", or "Species" to indicate which type of 
#'  substructure.
#'  
#' @return FA substructure matrix encoding the frequency of each substructure 
#' @export
lipid_sub_matrix <- function(exp_data, sub_data,
                             sub_type = "species") {
  sub_mat <- list()
  
  if (sub_type == "FA") {
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x) {
      x[str_detect(names(x), "FA")]
    })
  } else if (sub_type == "Class") {
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x) {
      x <- x[!str_detect(names(x), "FA")]
      x <- unique(str_replace(x, "_.+", ""))
    })
  } else {
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x) {
      x <- x[!str_detect(names(x), "FA")]
    })
  }
  
  all_sub <- sub_data[[2]] %>%
    unlist() %>%
    unique() %>%
    sort()
  
  for (num in 1:nrow(exp_data)) {
    loc <- which(sub_data[[1]] == exp_data[num, 1])
    if (length(loc) != 0) {
      # Substructure weight: all possible path divided by path number
      sub_mat[[num]] <- as.integer(table(factor(unlist(sub_data[[2]][loc]), levels = all_sub))) / length(loc)
    } else {
      sub_mat[[num]] <- integer(length(all_sub))
    }
  }
  
  sub_mat <- sub_mat %>% as.data.frame()
  
  colnames(sub_mat) <- exp_data$feature
  
  rownames(sub_mat) <- sub_data[[2]] %>%
    unlist() %>%
    unique() %>%
    sort()
  
  sub_mat <- sub_mat[rowSums(sub_mat) != 0, ]
  
  sub_mat <- as.matrix(sub_mat)
  
  exp_mat <- exp_data[-1] %>% as.matrix()
  
  sub_exp <- sub_mat %*% exp_mat
  
  return(list(sub_mat, exp_mat, sub_exp))
}
