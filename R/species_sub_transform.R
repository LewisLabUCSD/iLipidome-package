#' decomposes lipids into species substructures based on the lipid biosynthetic 
#'  network.
#'
#' @param char Lipid characteristics table. Output of Output of 
#'  "build_char_table".
#' @param lipid_substructure A data.frame recording lipid class substructures 
#'  computed from the reference lipid biosynthetic network in iLipdiome.
#' @param network_node A data.frame recording iLipidome-supported lipid classes 
#'  and related information.
#'  
#' @return lipid species substructure table
#' @export
species_sub_transform <- function(char, lipid_substructure, network_node) {
  char <- char %>% filter(class %in% lipid_substructure$Lipid)
  lipid_name <- character()
  lipid_sub_list <- list()
  
  add <- 1
  add2 <- 1
  
  FA_pool <- char$FA_split[char$FA_split != ""] %>%
    str_split("_") %>%
    unlist() %>%
    unique()
  
  for (var1 in 1:nrow(char)) {
    # PC, PE
    # print(var1)
    
    LMID <- filter(network_node, Abbreviation == char$class[var1])$LIPIDMAPS_subclass_id
    
    # lipid substructure: all possible path method
    lipid_sub_all <- lipid_substructure %>% filter(Lipid %in% char$class[var1])
    # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
    # PE: G3P, LPA_FA1, PA_FA2, DAG_FA2, PE_FA2
    if (char[var1, "FA_num"] == 2) {
      if (char[var1, "FA_split"] == "") {
        FA2 <- char$FA_sum[var1] # c('34:1;0')
        
        # chain_db_oh <- str_extract_all(char[var1,'FA_sum'], '\\d+') %>% unlist()
        # c('34','1','0')
        if_sphingolipid <- filter(network_node, Abbreviation == char$class[var1])$Class
        
        if (if_sphingolipid %in% c("Sphingolipid", "Deoxysphingolipid")) {
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[str_detect(lipid_sub, "FA2")] <- str_c(lipid_sub[str_detect(lipid_sub, "FA2")], FA2, sep = "_")
            
            lipid_sub_FA1 <- lipid_sub
            
            if (!char$class[var1] %in% c("dhCer", "dhSM", "doxdhCer", "doxdhSM")) {
              length_db_oh <- str_extract_all(FA2, "\\d+") %>%
                unlist() %>%
                as.integer()
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "^dhCer")] <- str_c(
                "dhCer_FA2_", length_db_oh[1], ":",
                length_db_oh[2] - 1, ";2"
              )
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "doxdhCer")] <- str_c(
                "doxdhCer_FA2_", length_db_oh[1], ":",
                length_db_oh[2] - 1, ";1"
              )
            }
            if (str_detect(char$class[var1], "dox")) {
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")], "18:0;0", sep = "_")
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")], "18:0;1", sep = "_")
            } else {
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")], "18:0;1", sep = "_")
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")], "18:0;2", sep = "_")
            }
            
            
            names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
            
            lipid_sub_list[[add]] <- lipid_sub_FA1
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
        } else {
          LM_mapping <- lipidmaps_sum %>% filter(sub_class_id == LMID, FA_sum == str_extract(FA2, "\\d+:\\d+"))
          
          if (nrow(LM_mapping) != 0) {
            possible_FA <- LM_mapping$FA_split %>%
              str_split("_") %>%
              unlist() %>%
              unique()
            
            possible_FA <- str_c(possible_FA, ";0")
            
            if (length(FA_pool) != 0) {
              possible_FA <- possible_FA[possible_FA %in% FA_pool]
            }
          } else {
            for (var2 in 1:nrow(lipid_sub_all)) {
              lipid_sub <- lipid_sub_all[var2, ] %>%
                dplyr::select(-Lipid) %>%
                unlist()
              # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
              
              lipid_sub <- lipid_sub[lipid_sub != ""]
              
              lipid_sub[str_detect(lipid_sub, "FA2")] <- str_c(lipid_sub[str_detect(lipid_sub, "FA2")], FA2, sep = "_")
              
              names(lipid_sub) <- rep("Lipid", length(lipid_sub))
              
              lipid_sub_list[[add]] <- lipid_sub
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
            next
          }
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[str_detect(lipid_sub, "FA2")] <- str_c(lipid_sub[str_detect(lipid_sub, "FA2")], FA2, sep = "_")
            
            last_FA1_lipid <- last(lipid_sub[str_detect(lipid_sub, "FA1")]) %>%
              str_replace("_FA1", "")
            
            possible_FA1_lipid <- char %>% filter(class == last_FA1_lipid, FA_split %in% possible_FA)
            
            if (nrow(possible_FA1_lipid) != 0) {
              possible_FA1 <- possible_FA1_lipid$FA_split
            } else {
              possible_FA1 <- possible_FA
            }
            for (var3 in 1:length(possible_FA1)) {
              lipid_sub_FA1 <- lipid_sub
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")],
                                                                       possible_FA1[var3],
                                                                       sep = "_"
              )
              
              names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
          }
        }
        
        # FA_pool_chain_db_oh <- str_extract_all(FA_pool, '\\d+')
        # c('16','0','0'), c('18','1','0')...
        
        # possible_FA <- FA_pool_chain_db_oh %>% map(.f = function(x){str_c(as.integer(chain_db_oh)-as.integer(x), collapse = ':')}) %>%
        #  unlist()
        # c('18:1:0', '16:0:0')
        
        # str_sub(possible_FA, -2,-2) <- ';'
        # c('18:1;0', '16:0;0')
      } else if (char[var1, "FA_split"] != "") {
        FA2 <- unlist(char[var1, "each_FA"]) # c('16:0;0', '18:0;0')
        
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[str_detect(lipid_sub, "FA2")] <- str_c(lipid_sub[str_detect(lipid_sub, "FA2")], str_c(FA2, collapse = "_"), sep = "_")
          
          if_sphingolipid <- filter(network_node, Abbreviation == char$class[var1])$Class
          
          if (if_sphingolipid %in% c("Sphingolipid", "Deoxysphingolipid")) {
            lipid_sub_FA1 <- lipid_sub
            
            if (str_detect(char$class[var1], "dox")) {
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "doxdhCer")] <- str_c("doxdhCer_FA2_18:0;1_", FA2[2])
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")], "18:0;0", sep = "_")
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")], "18:0;1", sep = "_")
            } else {
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhCer")] <- str_c("dhCer_FA2_18:0;2_", FA2[2])
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "ketoSPB")], "18:0;1", sep = "_")
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "dhSPB")], "18:0;2", sep = "_")
            }
            
            
            names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
            
            lipid_sub_list[[add]] <- lipid_sub_FA1
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
            
            FA2 <- FA2[2]
          } else {
            if (filter(network_node, Abbreviation == char$class[var1])$Class == "Ether lipid") {
              lipid_sub_FA1 <- lipid_sub
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")], FA2[1], sep = "_")
              
              names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
              next
            } else {
              # first FA
              possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% FA2)
              
              if (nrow(possible_FA1_lipid) != 0) {
                FA1 <- possible_FA1_lipid$FA_split
              } else {
                FA1 <- FA2
              }
              for (var3 in 1:length(FA1)) {
                lipid_sub_FA1 <- lipid_sub
                
                lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA1")], FA1[var3], sep = "_")
                
                names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
                
                lipid_sub_list[[add]] <- lipid_sub_FA1
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          }
        }
      }
    } else if (char[var1, "FA_num"] == 1) {
      FA1 <- unlist(char[var1, "each_FA"]) # c('16:0;0')
      
      if_sphingolipid <- filter(network_node, Abbreviation == char$class[var1])$Class
      
      
      for (var2 in 1:nrow(lipid_sub_all)) {
        lipid_sub <- lipid_sub_all[var2, ] %>%
          dplyr::select(-Lipid) %>%
          unlist()
        # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
        
        lipid_sub <- lipid_sub[lipid_sub != ""]
        lipid_sub[str_detect(lipid_sub, "FA1")] <- str_c(lipid_sub[str_detect(lipid_sub, "FA1")], FA1, sep = "_")
        
        if (if_sphingolipid %in% c("Sphingolipid", "Deoxysphingolipid")) {
          if (!str_detect(char$class[var1], "dox")) {
            lipid_sub <- str_replace(lipid_sub, "ketoSPB_FA1_.+", "ketoSPB_FA1_18:0;1")
            lipid_sub <- str_replace(lipid_sub, "dhSPB_FA1_.+", "dhSPB_FA1_18:0;2")
          } else {
            lipid_sub <- str_replace(lipid_sub, "doxketoSPB_FA1_.+", "doxketoSPB_FA1_18:0;0")
            lipid_sub <- str_replace(lipid_sub, "doxdhSPB_FA1_.+", "doxdhSPB_FA1_18:0;1")
          }
        }
        
        # LPE: LPA_FA1_16:0;0
        
        
        if (sum(str_detect(lipid_sub, "FA2")) != 0) {
          last_FA2_lipid <- last(lipid_sub[str_detect(lipid_sub, "FA2")]) %>%
            str_replace("_FA2", "")
          # PC for LPC
          if (filter(network_node, Abbreviation == last_FA2_lipid)$Class == "Ether lipid") {
            possible_FA2_lipid <- char %>%
              filter(class == last_FA2_lipid) %>%
              filter(str_detect(FA_split, str_c(FA1, "_")))
          } else {
            possible_FA2_lipid <- char %>%
              filter(class == last_FA2_lipid) %>%
              filter(str_detect(FA_split, FA1))
          }
          
          
          if (nrow(possible_FA2_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA2_lipid)) {
              lipid_sub_FA1 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA2")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA2")], possible_FA2, sep = "_")
              
              if (if_sphingolipid %in% c("Sphingolipid", "Deoxysphingolipid")) {
                lipid_sub_FA1 <- str_replace(lipid_sub_FA1, "^dhCer_FA2_18:1;2", "dhCer_FA2_18:0;2")
                lipid_sub_FA1 <- str_replace(lipid_sub_FA1, "doxdhCer_FA2_18:1;1", "dhCer_FA2_18:0;1")
              }
              
              names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
          } else if (nrow(possible_FA2_lipid) == 0) {
            LMID <- filter(network_node, Abbreviation == last_FA2_lipid)$LIPIDMAPS_subclass_id
            
            LM_mapping <- lipidmaps_sum %>%
              filter(sub_class_id == LMID) %>%
              filter(str_detect(FA_split, str_extract(FA1, "\\d+:\\d+")))
            
            if (nrow(LM_mapping) != 0) {
              possible_FA2_sum <- LM_mapping$FA_sum %>% unique()
              possible_FA2_sum <- str_c(possible_FA2_sum, ";0")
            } else {
              lipid_sub_FA1 <- lipid_sub
              
              names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
              
              next
            }
            
            possible_FA2_lipid <- char %>%
              filter(FA_split == "") %>%
              filter(class == last_FA2_lipid, FA_sum %in% possible_FA2_sum)
            
            if (nrow(possible_FA2_lipid) != 0) {
              for (var3 in 1:nrow(possible_FA2_lipid)) {
                lipid_sub_FA1 <- lipid_sub
                
                possible_FA2 <- possible_FA2_lipid$FA_sum[var3]
                
                lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA2")] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, "FA2")], possible_FA2, sep = "_")
                
                if (if_sphingolipid %in% c("Sphingolipid", "Deoxysphingolipid")) {
                  length_db_oh <- str_extract_all(possible_FA2, "\\d+") %>% as.integer()
                  
                  lipid_sub_FA1[str_detect(lipid_sub_FA1, "^dhCer")] <- str_c("dhCer_FA2_", length_db_oh[1], ":",
                                                                              length_db_oh[2] - 1, ";2",
                                                                              sep = "_"
                  )
                  lipid_sub_FA1[str_detect(lipid_sub_FA1, "doxdhCer")] <- str_c("doxdhCer_FA2_", length_db_oh[1], ":",
                                                                                length_db_oh[2] - 1, ";1",
                                                                                sep = "_"
                  )
                }
                
                
                names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
                
                lipid_sub_list[[add]] <- lipid_sub_FA1
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            } else {
              lipid_sub_FA1 <- lipid_sub
              
              names(lipid_sub_FA1) <- rep("Lipid", length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
          }
        } else {
          names(lipid_sub) <- rep("Lipid", length(lipid_sub))
          
          lipid_sub_list[[add]] <- lipid_sub
          
          lipid_name[add] <- char$feature[var1]
          
          add <- add + 1
        }
      }
    } else if (char$class[var1] == "TAG") {
      if (char[var1, "FA_split"] != "") {
        FA3 <- unlist(char[var1, "each_FA"])
        # c('16:0;0', '18:0;0','18:0;0')
        
        # combo1 <- str_c(FA3[1], FA3[2],sep = '_')
        # combo2 <- str_c(FA3[2], FA3[3],sep = '_')
        # combo3 <- str_c(FA3[1], FA3[1],sep = '_')
        # all_combo <- c(combo1, combo2,combo3)
        
        possible_FA2_lipid <- char %>%
          filter(class == "DAG") %>%
          filter(sum(!unlist(str_split(FA_split, "_")) %in% FA3) == 0) %>%
          filter(sum(FA3 %in% unlist(str_split(FA_split, "_"))) >= 2)
        
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[lipid_sub == "TAG_FA3"] <- str_c("TAG_FA3", str_c(FA3, collapse = "_"), sep = "_")
          
          if (nrow(possible_FA2_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA2_lipid)) {
              lipid_sub_FA3 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")],
                                                                       possible_FA2,
                                                                       sep = "_"
              )
              possible_FA1 <- unlist(str_split(possible_FA2))
              
              possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
              
              if (nrow(possible_FA1_lipid) != 0) {
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for (var4 in 1:length(possible_FA1)) {
                lipid_sub_FA3_2 <- lipid_sub_FA3
                
                lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")],
                                                                             possible_FA1[var4],
                                                                             sep = "_"
                )
                names(lipid_sub_FA3_2) <- rep("Lipid", length(lipid_sub_FA3_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          } else if (nrow(possible_FA2_lipid) != 0) {
            lipid_sub_FA3 <- lipid_sub
            
            names(lipid_sub_FA3) <- rep("Lipid", length(lipid_sub_FA3))
            
            lipid_sub_list[[add]] <- lipid_sub_FA3
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
        }
      } else if (char[var1, "FA_split"] == "") {
        FA3 <- char$FA_sum[var1] # c('52:2;0')
        LM_mapping <- lipidmaps_sum %>% filter(sub_class_id == LMID, FA_sum == str_extract(FA3, "\\d+:\\d+"))
        
        if (nrow(LM_mapping) != 0) {
          possible_DAG <- LM_mapping$FA_split %>%
            str_split("_") %>%
            map(.f = function(x) {
              a <- str_c(x[1], x[2], sep = "_")
              b <- str_c(x[1], x[3], sep = "_")
              c <- str_c(x[2], x[3], sep = "_")
              return(c(a, b, c))
            }) %>%
            unlist() %>%
            unique()
          possible_DAG <- str_replace_all(possible_DAG, "_", ";0_")
          possible_DAG <- str_c(possible_DAG, ";0")
          
          if (length(FA_pool) != 0) {
            FA_pool_filter <- str_split(possible_DAG, "_") %>%
              map(.f = function(x) {
                sum(!x %in% FA_pool) == 0
              }) %>%
              unlist()
            possible_DAG <- possible_DAG[FA_pool_filter]
          }
        }
        if (nrow(LM_mapping) == 0 || length(possible_DAG) == 0) {
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[lipid_sub == "TAG_FA3"] <- str_c("TAG_FA3", FA3, sep = "_")
            
            names(lipid_sub) <- rep("Lipid", length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
          next
        }
        
        possible_FA2_lipid <- char %>% filter(class == "DAG", FA_split %in% possible_DAG)
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[lipid_sub == "TAG_FA3"] <- str_c("TAG_FA3", FA3, sep = "_")
          
          
          if (nrow(possible_FA2_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA2_lipid)) {
              lipid_sub_FA3 <- lipid_sub
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")], possible_FA2, sep = "_")
              
              possible_FA1 <- unlist(str_split(possible_FA2, "_"))
              possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
              
              if (nrow(possible_FA1_lipid) != 0) {
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              
              for (var4 in 1:length(possible_FA1)) {
                lipid_sub_FA3_2 <- lipid_sub_FA3
                
                
                lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")],
                                                                             possible_FA1[var4],
                                                                             sep = "_"
                )
                names(lipid_sub_FA3_2) <- rep("Lipid", length(lipid_sub_FA3_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          } else if (nrow(possible_FA2_lipid) == 0) {
            possible_FA_sum <- possible_DAG %>%
              str_extract_all("\\d+") %>%
              map(.f = function(x) {
                a <- as.integer(x[1]) + as.integer(x[4])
                b <- as.integer(x[2]) + as.integer(x[5])
                c <- str_c(a, b, sep = ":")
                return(c)
              }) %>%
              unlist()
            
            possible_FA_sum <- str_c(possible_FA_sum, ";0")
            
            
            possible_FA2_lipid_2 <- char %>%
              filter(FA_split == "") %>%
              filter(class == "DAG", FA_sum %in% possible_FA_sum)
            
            
            if (nrow(possible_FA2_lipid_2) != 0) {
              for (var3 in 1:nrow(possible_FA2_lipid_2)) {
                lipid_sub_FA3 <- lipid_sub
                possible_FA2 <- possible_FA2_lipid_2$FA_sum[var3]
                lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, "FA2")], possible_FA2, sep = "_")
                
                possible_FA1 <- possible_DAG[which(possible_FA_sum == possible_FA2)] %>%
                  str_split("_") %>%
                  unlist()
                possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
                
                if (nrow(possible_FA1_lipid) != 0) {
                  possible_FA1 <- possible_FA1_lipid$FA_split
                }
                
                for (var4 in 1:length(possible_FA1)) {
                  lipid_sub_FA3_2 <- lipid_sub_FA3
                  
                  lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA1")],
                                                                               possible_FA1[var4],
                                                                               sep = "_"
                  )
                  names(lipid_sub_FA3_2) <- rep("Lipid", length(lipid_sub_FA3_2))
                  
                  lipid_sub_list[[add]] <- lipid_sub_FA3_2
                  
                  lipid_name[add] <- char$feature[var1]
                  
                  add <- add + 1
                }
              }
            } else if (nrow(possible_FA2_lipid_2) == 0) {
              lipid_sub_FA3 <- lipid_sub
              
              names(lipid_sub_FA3) <- rep("Lipid", length(lipid_sub_FA3))
              
              lipid_sub_list[[add]] <- lipid_sub_FA3
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
          }
        }
      }
    } else if (char$class[var1] == "CL") {
      if (char[var1, "FA_split"] != "") {
        FA4 <- unlist(char[var1, "each_FA"])
        
        possible_FA2_lipid <- char %>%
          filter(class == "PG") %>%
          filter(sum(!unlist(str_split(FA_split, "_")) %in% FA4) == 0) %>%
          filter(sum(FA4 %in% unlist(str_split(FA_split, "_"))) >= 2)
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[lipid_sub == "CL_FA4"] <- str_c("CL_FA4", str_c(FA4, collapse = "_"), sep = "_")
          
          if (nrow(possible_FA2_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA2_lipid)) {
              lipid_sub_FA4 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")],
                                                                       possible_FA2,
                                                                       sep = "_"
              )
              
              possible_FA1 <- unlist(str_split(possible_FA2, "_"))
              
              possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
              
              if (nrow(possible_FA1_lipid) != 0) {
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for (var4 in 1:length(possible_FA1)) {
                lipid_sub_FA4_2 <- lipid_sub_FA4
                
                lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")],
                                                                             possible_FA1[var4],
                                                                             sep = "_"
                )
                names(lipid_sub_FA4_2) <- rep("Lipid", length(lipid_sub_FA4_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA4_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          } else if (nrow(possible_FA2_lipid) == 0) {
            lipid_sub_FA4 <- lipid_sub
            
            names(lipid_sub_FA4) <- rep("Lipid", length(lipid_sub_FA4))
            
            lipid_sub_list[[add]] <- lipid_sub_FA4
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
        }
      } else if (char[var1, "FA_split"] == "") {
        FA4 <- char$FA_sum[var1] # c('52:2;0')
        
        LM_mapping <- lipidmaps_sum %>% filter(sub_class_id == LMID, FA_sum == str_extract(FA4, "\\d+:\\d+"))
        
        if (nrow(LM_mapping) != 0) {
          possible_PG <- MLCL_mapping %>%
            filter(FA_sum == str_extract(FA4, "\\d+:\\d+")) %>%
            .$possible_PG_split %>%
            unique()
          
          possible_PG <- str_replace_all(possible_PG, "_", ";0_")
          possible_PG <- str_c(possible_PG, ";0")
          
          if (length(FA_pool) != 0) {
            FA_pool_filter <- str_split(possible_PG, "_") %>%
              map(.f = function(x) {
                sum(!x %in% FA_pool) == 0
              }) %>%
              unlist()
            possible_PG <- possible_PG[FA_pool_filter]
          }
        }
        if (nrow(LM_mapping) == 0 || length(possible_PG) == 0) {
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[lipid_sub == "CL_FA4"] <- str_c("CL_FA4", FA4, sep = "_")
            
            names(lipid_sub) <- rep("Lipid", length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
          next
        }
        
        possible_FA2_lipid <- char %>% filter(class == "PG", FA_split %in% possible_PG)
        
        if (nrow(possible_FA2_lipid) != 0) {
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[lipid_sub == "CL_FA4"] <- str_c("CL_FA4", FA4, sep = "_")
            for (var3 in 1:nrow(possible_FA2_lipid)) {
              lipid_sub_FA4 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")], possible_FA2, sep = "_")
              
              possible_FA1 <- unlist(str_split(possible_FA2, "_"))
              
              possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
              
              if (nrow(possible_FA1_lipid) != 0) {
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for (var4 in 1:length(possible_FA1)) {
                lipid_sub_FA4_2 <- lipid_sub_FA4
                
                lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")],
                                                                             possible_FA1[var4],
                                                                             sep = "_"
                )
                names(lipid_sub_FA4_2) <- rep("Lipid", length(lipid_sub_FA4_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA4_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          }
        } else if (nrow(possible_FA2_lipid) == 0) {
          possible_PG_sum <- possible_PG %>%
            str_extract_all("\\d+") %>%
            map(.f = function(x) {
              a <- as.integer(x[1]) + as.integer(x[4])
              b <- as.integer(x[2]) + as.integer(x[5])
              c <- str_c(a, b, sep = ":")
              return(c)
            }) %>%
            unlist()
          
          possible_PG_sum <- str_c(possible_PG_sum, ";0")
          
          
          possible_FA2_lipid_2 <- char %>%
            filter(FA_split == "") %>%
            filter(class == "PG", FA_sum %in% possible_PG_sum)
          
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[lipid_sub == "CL_FA4"] <- str_c("CL_FA4", FA4, sep = "_")
            
            if (nrow(possible_FA2_lipid_2) != 0) {
              for (var3 in 1:nrow(possible_FA2_lipid_2)) {
                lipid_sub_FA4 <- lipid_sub
                
                possible_FA2 <- possible_FA2_lipid_2$FA_sum[var3]
                
                lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, "FA2")], possible_FA2, sep = "_")
                
                possible_FA1 <- possible_PG[which(possible_PG_sum == possible_FA2)] %>%
                  str_split("_") %>%
                  unlist()
                
                possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
                
                if (nrow(possible_FA1_lipid) != 0) {
                  possible_FA1 <- possible_FA1_lipid$FA_split
                }
                for (var4 in 1:length(possible_FA1)) {
                  lipid_sub_FA4_2 <- lipid_sub_FA4
                  
                  lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, "FA1")],
                                                                               possible_FA1[var4],
                                                                               sep = "_"
                  )
                  names(lipid_sub_FA4_2) <- rep("Lipid", length(lipid_sub_FA4_2))
                  
                  lipid_sub_list[[add]] <- lipid_sub_FA4_2
                  
                  lipid_name[add] <- char$feature[var1]
                  
                  add <- add + 1
                }
              }
            } else if (nrow(possible_FA2_lipid_2) == 0) {
              lipid_sub_FA4 <- lipid_sub
              
              names(lipid_sub_FA4) <- rep("Lipid", length(lipid_sub_FA4))
              
              lipid_sub_list[[add]] <- lipid_sub_FA4
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add + 1
            }
          }
        }
      }
    } else if (char$class[var1] == "MLCL") {
      if (char[var1, "FA_split"] != "") {
        FA3 <- unlist(char[var1, "each_FA"])
        # c('16:0;0', '18:0;0','18:0;0')
        
        
        possible_FA4_lipid <- char %>%
          filter(class == "CL") %>%
          filter(sum(!unlist(str_split(FA_split, "_")) %in% FA3) == 0) %>%
          filter(sum(!FA3 %in% unlist(str_split(FA_split, "_"))) == 0)
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[lipid_sub == "MLCL_FA3"] <- str_c("MLCL_FA3", str_c(FA3, collapse = "_"), sep = "_")
          
          if (nrow(possible_FA4_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA4_lipid)) {
              lipid_sub_FA3 <- lipid_sub
              
              possible_FA4 <- possible_FA4_lipid$FA_split[var3]
              
              lipid_sub_FA3[lipid_sub == "CL_FA4"] <- str_c("CL_FA4", possible_FA4, sep = "_")
              
              FA4 <- str_split(possible_FA4, "_") %>% unlist()
              
              possible_FA2_lipid <- char %>%
                filter(class == "PG") %>%
                filter(sum(!unlist(str_split(FA_split, "_")) %in% FA4) == 0) %>%
                filter(sum(FA4 %in% unlist(str_split(FA_split, "_"))) >= 2)
              
              if (nrow(possible_FA2_lipid) != 0) {
                for (var4 in 1:nrow(possible_FA2_lipid)) {
                  lipid_sub_FA3_2 <- lipid_sub_FA3
                  
                  possible_FA2 <- possible_FA2_lipid$FA_split[var4]
                  
                  lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA2")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA2")],
                                                                               possible_FA2,
                                                                               sep = "_"
                  )
                  
                  possible_FA1 <- unlist(str_split(possible_FA2, "_"))
                  
                  possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
                  
                  if (nrow(possible_FA1_lipid) != 0) {
                    possible_FA1 <- possible_FA1_lipid$FA_split
                  }
                  for (var4 in 1:length(possible_FA1)) {
                    lipid_sub_FA3_3 <- lipid_sub_FA3_2
                    
                    
                    lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")],
                                                                                 possible_FA1[var4],
                                                                                 sep = "_"
                    )
                    names(lipid_sub_FA3_3) <- rep("Lipid", length(lipid_sub_FA3_3))
                    
                    lipid_sub_list[[add]] <- lipid_sub_FA3_3
                    
                    lipid_name[add] <- char$feature[var1]
                    
                    add <- add + 1
                  }
                }
              } else if (nrow(possible_FA2_lipid) == 0) {
                names(lipid_sub_FA3) <- rep("Lipid", length(lipid_sub_FA3))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add + 1
              }
            }
          } else if (nrow(possible_FA4_lipid) == 0) {
            lipid_sub_FA3 <- lipid_sub
            
            names(lipid_sub_FA3) <- rep("Lipid", length(lipid_sub_FA3))
            
            lipid_sub_list[[add]] <- lipid_sub_FA3
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
        }
      } else if (char[var1, "FA_split"] == "") {
        FA3 <- char$FA_sum[var1] # c('52:2;0')
        
        LM_mapping <- MLCL_mapping %>%
          filter(possible_MLCL_sum == FA3) %>%
          .[3:5] %>%
          unique()
        
        if (nrow(LM_mapping) != 0) {
          possible_CL <- LM_mapping$FA_split
          possible_CL_sum <- LM_mapping$FA_sum
          
          possible_CL <- str_replace_all(possible_CL, "_", ";0_")
          possible_CL <- str_c(possible_CL, ";0")
          possible_CL_sum <- str_c(possible_CL_sum, ";0")
          
          if (length(FA_pool) != 0) {
            FA_pool_filter <- str_split(possible_CL, "_") %>%
              map(.f = function(x) {
                sum(!x %in% FA_pool) == 0
              }) %>%
              unlist()
            
            possible_CL <- possible_CL[FA_pool_filter]
            possible_CL_sum <- possible_CL_sum[FA_pool_filter]
          }
        }
        
        if (nrow(LM_mapping) == 0 || length(possible_CL) == 0) {
          for (var2 in 1:nrow(lipid_sub_all)) {
            lipid_sub <- lipid_sub_all[var2, ] %>%
              dplyr::select(-Lipid) %>%
              unlist()
            # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub != ""]
            
            lipid_sub[lipid_sub == "MLCL_FA3"] <- str_c("MLCL_FA3", FA3, sep = "_")
            
            names(lipid_sub) <- rep("Lipid", length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
          next
        }
        
        possible_FA4_lipid <- char %>% filter(class == "CL", FA_sum %in% possible_CL_sum)
        
        
        for (var2 in 1:nrow(lipid_sub_all)) {
          lipid_sub <- lipid_sub_all[var2, ] %>%
            dplyr::select(-Lipid) %>%
            unlist()
          # PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub != ""]
          
          lipid_sub[lipid_sub == "MLCL_FA3"] <- str_c("MLCL_FA3", FA3, sep = "_")
          
          
          if (nrow(possible_FA4_lipid) != 0) {
            for (var3 in 1:nrow(possible_FA4_lipid)) {
              lipid_sub_FA3 <- lipid_sub
              FA4 <- possible_FA4_lipid$FA_sum[var3]
              lipid_sub_FA3[lipid_sub_FA3 == "CL_FA4"] <- str_c("CL_FA4", FA4, sep = "_")
              
              possible_FA4 <- possible_CL[which(possible_CL_sum == FA4)]
              
              FA4s <- str_split(possible_FA4, "_")
              
              
              for (var4 in 1:length(FA4s)) {
                FA4 <- FA4s[[var4]]
                
                possible_FA2_lipid <- char %>%
                  filter(class == "PG") %>%
                  filter(sum(!unlist(str_split(FA_split, "_")) %in% FA4) == 0) %>%
                  filter(sum(FA4 %in% unlist(str_split(FA_split, "_"))) >= 2)
                
                
                if (nrow(possible_FA2_lipid) != 0) {
                  for (var5 in 1:nrow(possible_FA2_lipid)) {
                    lipid_sub_FA3_2 <- lipid_sub_FA3
                    
                    possible_FA2 <- possible_FA2_lipid$FA_split[var5]
                    
                    lipid_sub_FA4[str_detect(lipid_sub_FA3_2, "FA2")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA2")],
                                                                               possible_FA2,
                                                                               sep = "_"
                    )
                    
                    possible_FA1 <- unlist(str_split(possible_FA2, "_"))
                    
                    possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
                    
                    
                    if (nrow(possible_FA1_lipid) != 0) {
                      possible_FA1 <- possible_FA1_lipid$FA_split
                    }
                    for (var6 in 1:length(possible_FA1)) {
                      lipid_sub_FA3_3 <- lipid_sub_FA3_2
                      
                      lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")],
                                                                                   possible_FA1[var6],
                                                                                   sep = "_"
                      )
                      names(lipid_sub_FA3_3) <- rep("Lipid", length(lipid_sub_FA3_3))
                      
                      lipid_sub_list[[add]] <- lipid_sub_FA3_3
                      
                      lipid_name[add] <- char$feature[var1]
                      
                      add <- add + 1
                    }
                  }
                } else if (nrow(possible_FA2_lipid) == 0) {
                  FA4 <- str_c(FA4, collapse = "_") %>%
                    str_replace_all(pattern = ";0", "")
                  
                  possible_FA_sum_2 <- MLCL_mapping %>%
                    filter(FA_split == FA4) %>%
                    .$possible_PG_sum %>%
                    unique()
                  
                  possible_FA2_lipid_2 <- char %>%
                    filter(FA_split == "") %>%
                    filter(class == "PG", FA_sum %in% possible_FA_sum_2)
                  
                  
                  if (nrow(possible_FA2_lipid_2) != 0) {
                    for (var5 in 1:nrow(possible_FA2_lipid_2)) {
                      lipid_sub_FA3_2 <- lipid_sub_FA3
                      
                      possible_FA2 <- possible_FA2_lipid_2$FA_sum[var5]
                      lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA2")] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, "FA2")], possible_FA2, sep = "_")
                      
                      
                      LM_mapping_2 <- lipidmaps_sum %>% filter(sub_class_id == "[GP0401]", FA_sum == str_extract(possible_FA2, "\\d+:\\d+"))
                      
                      
                      if (nrow(LM_mapping_2) != 0) {
                        possible_FA1 <- LM_mapping_2$FA_split %>%
                          str_split("_") %>%
                          unlist() %>%
                          unique()
                        
                        possible_FA1 <- str_c(possible_FA1, ";0")
                        
                        if (length(FA_pool) != 0) {
                          possible_FA1 <- possible_FA1[possible_FA1 %in% FA_pool]
                        }
                      } else {
                        names(lipid_sub_FA3_2) <- rep("Lipid", length(lipid_sub_FA3_2))
                        
                        lipid_sub_list[[add]] <- lipid_sub_FA3_2
                        
                        lipid_name[add] <- char$feature[var1]
                        
                        add <- add + 1
                        next
                      }
                      
                      
                      
                      possible_FA1_lipid <- char %>% filter(class == "LPA", FA_split %in% possible_FA1)
                      
                      if (nrow(possible_FA1_lipid) != 0) {
                        possible_FA1 <- possible_FA1_lipid$FA_split
                      }
                      for (var6 in 1:length(possible_FA1)) {
                        lipid_sub_FA3_3 <- lipid_sub_FA3_2
                        
                        lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, "FA1")],
                                                                                     possible_FA1[var6],
                                                                                     sep = "_"
                        )
                        names(lipid_sub_FA3_3) <- rep("Lipid", length(lipid_sub_FA3_3))
                        
                        lipid_sub_list[[add]] <- lipid_sub_FA3_3
                        
                        lipid_name[add] <- char$feature[var1]
                        
                        add <- add + 1
                      }
                    }
                  } else if (nrow(possible_FA2_lipid) == 0) {
                    names(lipid_sub_FA3) <- rep("Lipid", length(lipid_sub_FA3))
                    
                    lipid_sub_list[[add]] <- lipid_sub_FA3
                    
                    lipid_name[add] <- char$feature[var1]
                    
                    add <- add + 1
                  }
                }
              }
            }
          } else if (nrow(possible_FA4_lipid) == 0) {
            names(lipid_sub) <- rep("Lipid", length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add + 1
          }
        }
      }
    }
  }
  sub_species <- list(lipid_name, lipid_sub_list)
  species_substructure <- sub_species[[2]] %>%
    map(.f = function(x) {
      names(x) <- 1:length(x)
      return(x)
    }) %>%
    plyr::ldply(rbind) %>%
    mutate(Lipid = sub_species[[1]]) %>%
    dplyr::select(Lipid, everything())
  
  species_substructure[is.na(species_substructure)] <- ""
  
  colnames(species_substructure) <- c("Lipid", str_c("Unit", 1:(ncol(species_substructure) - 1)))
  
  
  return(species_substructure)
}
