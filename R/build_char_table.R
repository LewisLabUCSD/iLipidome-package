#' Transforms lipid expression table
#' 
#' Transforms lipid expression table (“exp”) into two iLipidome inputs: 
#'  selected lipid expression table (“exp_sel”) and 
#'  selected lipid characteristics table (“char_sel”).
#'
#' @param raw_data A data.frame where rows are lipid species and columns are 
#'  samples. Lipid names should be provided in the first column. Lipid values 
#'  should be numeric.
#' @param network_node A data.frame recording iLipidome-supported lipid classes 
#'  and related information.
#'  
#' @return list containing selected lipid expression table and selected lipid 
#'  characteristics table
#'  
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr map_int
#' @importFrom stringr str_c
#' @importFrom stringr str_extract
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_sub
#' 
#' @export
build_char_table <- function(raw_data, network_node) {
  class <- raw_data$feature %>% str_extract("[A-Za-z]+( O-)*")
  class_included <- class %in% network_node$Abbreviation
  
  
  totallength <- raw_data$feature %>%
    str_extract_all("\\d+:") %>%
    map_int(~ str_sub(.x, end = -2) %>%
              as.integer() %>%
              sum())
  
  totaldb <- raw_data$feature %>%
    str_extract_all("\\d+;") %>%
    map_int(~ str_sub(.x, end = -2) %>%
              as.integer() %>%
              sum())
  
  totaloh <- raw_data$feature %>%
    str_extract_all(";\\d+") %>%
    map_int(~ str_sub(.x, start = 2) %>%
              as.integer() %>%
              sum())
  
  FA_sum <- str_c(totallength, ":", totaldb, ";", totaloh)
  
  FA_split <- raw_data$feature %>%
    str_extract_all("\\d+:\\d+;\\d+") %>%
    map_chr(~ str_c(.x, collapse = "_"))
  
  char_table <- data.frame(
    feature = raw_data$feature, class = class, totallength = totallength,
    totaldb = totaldb, totaloh = totaloh, FA_sum = FA_sum, FA_split = FA_split
  ) %>%
    .[class_included, ] %>%
    left_join(network_node[c("Abbreviation", "FA")], by = c("class" = "Abbreviation"))
  
  colnames(char_table)[8] <- "FA_num"
  FA_exact <- map_int(str_split(char_table$FA_split, "_"), ~ length(.x)) == char_table$FA_num
  char_table$FA_split[!FA_exact] <- ""
  
  each_FA <- str_extract_all(char_table$FA_split, "\\d+:\\d+;\\d+") %>% map(.f = function(x) {
    x[x != "0:0;0"]
  })
  
  char_table <- char_table %>% mutate(each_FA = each_FA)
  
  raw_data <- raw_data[class_included, ]
  
  return(list(raw_data, char_table))
}
