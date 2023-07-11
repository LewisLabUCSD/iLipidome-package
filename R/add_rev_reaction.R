#' Completes all reversible reactions in lipid species biosynthetic network
#' 
#' @param network_edge A data.frame recording edge information of the reference 
#'  lipid biosynthetic network in iLipdiome.
#' @param species_net Lipid species biosynthetic network. Output of 
#'  "build_species_net".
#'  
#' @return Lipid species biosynthetic network with complete reversible reactions
#' 
#' @export
add_rev_reaction <- function(network_edge, species_net) {
  rev_reaction <- network_edge %>%
    filter(reverse == 1) %>%
    apply(MARGIN = 1, FUN = function(x) {
      c(x[1], x[2])
    })
  
  species_net_pair <- map2(species_net$S1, species_net$P1,
                           .f = function(x, y) {
                             c(str_split(x, "_")[[1]][1], str_split(y, "_")[[1]][1])
                           }
  )
  
  add_reaction <- species_net_pair %>% map_lgl(.f = function(x) {
    max(colSums(matrix(rev_reaction %in% x, nrow = 2))) == 2
  })
  
  species_net_w_rev <- data.frame(S1 = species_net$P1[add_reaction], P1 = species_net$S1[add_reaction]) %>%
    rbind(species_net) %>%
    unique()
  
  return(species_net_w_rev)
}

