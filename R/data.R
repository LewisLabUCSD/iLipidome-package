#' Dataset Reference fatty acid network
#' 
#' Contains the reference fatty acid  network data
#' 
#'@format A dataframe with 51 rows and 5 variables
#' \describe{
#'    \item{S1}{fatty acid substrate}
#'    \item{P1}{fatty acid product}
#'    \item{S1_detail}{fatty acid substrate with double bond specification}
#'    \item{P1_detail}{fatty acid product with double bond specification}
#'    \item{pathway}{fatty acid pathway}
#'    }
#'  
#' @source {Created in-house}
#' 
#' @examples
#' data(FA_network)
"FA_network"


#' Dataset Lipid substructure
#' 
#' Contains the substructure decomposition rules based on reference lipid class network 
#' 
#'@format dataframe with 108 rows and 13 variables
#' \describe{
#'    \item{Lipid}{reference lipid class}
#'    \item{Unit1}{substructure unit 1}
#'    \item{Unit2}{substructure unit 2}
#'    \item{Unit3}{substructure unit 3}
#'    \item{Unit4}{substructure unit 4}
#'    \item{Unit5}{substructure unit 5}
#'    \item{Unit6}{substructure unit 6}
#'    \item{Unit7}{substructure unit 7}
#'    \item{Unit8}{substructure unit 8}
#'    \item{Unit9}{substructure unit 9}
#'    \item{Unit10}{substructure unit 10}
#'    \item{Unit11}{substructure unit 11}
#'    \item{Unit12}{substructure unit 12}  
#'    }
#'  
#' @source {Created in house}
#' 
#' @examples
#' data(lipid_substructure)
"lipid_substructure"


#' Dataset Reference lipid class network node 
#' 
#' Contains the reference lipid class network node data
#' 
#'@format dataframe with 64 rows and 6 variables
#' \describe{
#'    \item{Abbreviation}{abbreviation of lipid class}
#'    \item{Full_name}{Full name of lipid class}
#'    \item{Class}{category of lipid class}
#'    \item{LIPIDMAPS_subclass}{LipidMaps classification of lipid class}
#'    \item{LIPIDMAPS_subclass_id}{LipidMaps classification ID of lipid class}
#'    \item{FA}{Fatty acid number of lipid class}
#'    }
#'  
#' @source {Created in house}
#' 
#' @examples
#' data(network_node)
"network_node"

#' Dataset Reference lipid class network edge 
#' 
#' Contains the reference lipid class network edge  data
#' 
#'@format dataframe with 99 rows and 3 variables
#' \describe{
#'    \item{S1}{lipid class substract}
#'    \item{P1}{lipid class product}
#'    \item{reverse}{A binary indicator to specify if the reaction is reversible or not}
#'    }
#'  
#' @source {Created in house}
#' 
#' @examples
#' data(network_edge)
"network_edge"

#' Dataset LipidMaps mapping data
#' 
#' Contains LipidMaps mapping data for lipids without an exact fatty acid identity
#' 
#'@format dataframe with 16334 rows and 4 variables
#' \describe{
#'    \item{name}{lipid name}
#'    \item{sub_class_id}{LipidMaps classification ID of lipid}
#'    \item{FA_sum}{summed fatty acid composition of lipid}
#'    \item{FA_split}{individual fatty acid composition of lipid}
#'    }
#'  
#' @source {https://www.lipidmaps.org/databases/lmsd/download}
#' 
#' @examples
#' data(lipidmaps_sum)
"lipidmaps_sum"

#' Dataset LipidMaps mapping data for cardiolipin
#' 
#' Contains LipidMaps mapping data for cardiolipin without an exact fatty acid identity
#'
#'@format dataframe with 10102 rows and 8 variables
#' \describe{
#'    \item{name}{lipid name}
#'    \item{sub_class_id}{LipidMaps classification ID of lipid}
#'    \item{FA_sum}{summed fatty acid of lipid}
#'    \item{FA_split}{individual fatty acid of lipid}
#'    \item{possible_MLCL_sum}{summed fatty acid composition of matched monocardiolipin}
#'    \item{possible_MLCL_split}{individual fatty acid composition of matched monocardiolipin}
#'    \item{possible_PG_sum}{summed fatty acid composition of matched phosphatidylglycerol}
#'    \item{possible_PG_split}{individual fatty acid composition of matched phosphatidylglycerol}
#'    }
#'  
#' @source {https://www.lipidmaps.org/databases/lmsd/download}
#' 
#' @examples
#' data(MLCL_mapping)
"MLCL_mapping"

#' Dataset Reaction gene
#' 
#' Contains lipid reactions and their involved genes
#' 
#'@format dataframe with 390 rows and 3 variables
#' \describe{
#'    \item{Reaction}{lipid reaction}
#'    \item{species}{species name}
#'    \item{gene}{involved gene name}
#'    }
#'  
#' @source {https://www.lipidmaps.org/resources/tools/reactions}
#' 
#' @examples
#' data(reaction_gene_mapping)
"reaction_gene_mapping"
