#' Construct lipid biosynthetic network
#' 
#' Uses species substructures to construct lipid biosynthetic network
#' 
#' @param species_substructure Output of "species_sub_transform"
#' 
#' @return lipid species biosynthetic network
#' 
#' @export
build_species_net <- function(species_substructure) {
  species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
  species_network <- species_network %>% map(.f = function(x) {
    x[x != ""]
  })
  
  species_network <- species_network %>% map(.f = function(x) {
    head(rep(x, each = 2)[-1], -1)
  })
  species_network <- matrix(unlist(species_network), ncol = 2, byrow = T, dimnames = list(NULL, c("S1", "P1"))) %>%
    as.data.frame() %>%
    unique()
  
  species_network <- species_network %>% mutate(
    S1 = str_replace(S1, "_FA\\d", ""),
    P1 = str_replace(P1, "_FA\\d", "")
  )
  return(species_network)
}
