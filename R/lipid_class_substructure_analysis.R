#' Website function for lipid class substructure analysis
#' 
#' @param exp_raw raw expression input file 
#' @param method analysis method
#' @param ctrl range of values for control (1:7)
#' @param exp range of values for exp (8:13)
#' @param exo_lipid vector of fatty acids selected to prevent substructure decomposition
#' @param species the analysis species, (human, mouse, rat)
#' @param add_reaction add a reaction
#' @param delete_reaction delete a reaction
#' 
#' @return list of tables and plots relevant to the analysis
#' 
#' @import ggplot2
#' 
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom ggtext element_markdown
#' @importFrom igraph all_simple_paths
#' @importFrom magrittr %>%
#' @importFrom plyr ldply
#' @importFrom stats reorder 
#' @importFrom stringr str_c
#' @importFrom stringr str_split
#' @importFrom visNetwork visNetwork
#' 
#' @export

lipid_class_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                              exo_lipid=NULL,species='rat',
                                              add_reaction=NULL, delete_reaction=NULL){
  
  if(!is.null(add_reaction)){
    
    colnames(add_reaction) <- c('S1','P1')
    add_reaction <- add_reaction %>% filter(S1 %in% network_node$Abbreviation,
                                            P1 %in% network_node$Abbreviation)
    new_network_edge <- network_edge[,1:2] %>% rbind(add_reaction[,1:2]) %>% 
      unique()
    
  }
  if(!is.null(delete_reaction)){
    
    colnames(delete_reaction) <- c('S1','P1')
    delete_reaction <- delete_reaction %>% filter(S1 %in% network_node$Abbreviation,
                                                  P1 %in% network_node$Abbreviation)
    delete_reaction <- str_c(delete_reaction[[1]],'_', delete_reaction[[2]])
    
    new_network_edge <- new_network_edge %>% mutate(reaction=str_c(S1, '_', P1)) %>% 
      filter(!reaction %in% delete_reaction) %>% 
      dplyr::select(-reaction) %>% 
      unique()
  }
  
  if(exists('new_network_edge')){
    new_network_edge_sort <- new_network_edge %>% 
      apply(MARGIN = 1, FUN = function(x){str_c(sort(x), collapse = '_')})
    
    
    
    new_network_edge <- new_network_edge %>% mutate(reverse=0)
    
    new_network_edge$reverse[which(new_network_edge_sort %in% (names(table(new_network_edge_sort))[table(new_network_edge_sort)==2]))] <- 1
    
    network_edge <- new_network_edge
    
    graph <- graph_from_data_frame(new_network_edge[c('S1','P1')], directed = T, vertices = network_node$Abbreviation)
    
    find_path_from_G3P <- function(network_node, network_edge){
      lipid <- character()
      G3P_start_path <- list()
      lipid_list <- unique(c(network_edge$S1,network_edge$P1))
      num_class <- 1
      num_path <- 1
      
      while(!is.na(lipid_list[num_class])){
        class <- network_node %>% filter(Abbreviation==lipid_list[num_class]) %>% .$Class
        
        if(class=='Sphingolipid'){
          shortest_path <- tryCatch({all_simple_paths(graph,'Serine+Palmitoyl CoA',lipid_list[num_class], mode='out')},
                                    error=function(e){NULL})
        }
        else if(class=='Deoxysphingolipid'){
          shortest_path <- tryCatch({all_simple_paths(graph,'Alanine+Palmitoyl CoA',lipid_list[num_class], mode='out')},
                                    error=function(e){NULL})
        }
        else{
          shortest_path <- tryCatch({all_simple_paths(graph,'G3P',lipid_list[num_class], mode='out')},
                                    error=function(e){NULL})
        }
        
        if(!is.null(shortest_path)){
          if(length(shortest_path)!=0){
            
            for (path_num in 1:length(shortest_path)) {
              lipid[num_path] <- lipid_list[num_class]
              G3P_start_path[[num_path]] <- shortest_path[[path_num]] %>% 
                attributes() %>% .$names
              num_path <- num_path+1
              
            }
          }
          else{
            lipid[num_path] <- lipid_list[num_class]
            G3P_start_path[[num_path]] <- ''
            num_path <- num_path+1
          }
          
        }
        else{
          lipid[num_path] <- lipid_list[num_class]
          G3P_start_path[[num_path]] <- ''
          num_path <- num_path+1
        }
        num_class <- num_class+1
      }
      names(G3P_start_path) <- lipid
      
      
      lipid_substructure <- plyr::ldply(G3P_start_path, rbind) %>% unique()
      
      lipid_substructure[is.na(lipid_substructure)] <- ''
      
      colnames(lipid_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(lipid_substructure)-1)))
      
      
      for(a in 1:nrow(lipid_substructure)){
        FA_num <- network_node %>% filter(Abbreviation==lipid_substructure$Lipid[a]) %>% .$FA
        if(FA_num!=0){
          lipid_substructure[-1][lipid_substructure[-1]==lipid_substructure$Lipid[a]] <- str_c(lipid_substructure$Lipid[a], '_FA', as.character(FA_num))
        }
      }
      
      return(lipid_substructure)
      
    }
    
    lipid_substructure <- find_path_from_G3P(network_node, new_network_edge)
    
  }
  
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]
  
  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]
  
  no_sub_t <- unprocessed_data_test(exp_data = exp_data,
                                    char_table = char_data,
                                    method = method,
                                    significant='adj_p_value',
                                    ctrl_group = ctrl, exp_group = exp)
  
  #-------------------Lipid class substructure analysis-------------------
  
  #Extract class substructures using fold changes
  
  class_sub_stop <- species_sub_extract(lipid_substructure =lipid_substructure,
                                        unprocessed_data_result = no_sub_t,
                                        type = 'class', pct_limit = 0.01,
                                        exo_lipid=exo_lipid)
  
  #Transform lipid exp into substructures exp
  
  class_exp <- no_sub_t[[1]] %>% filter(type=='class') %>% 
    dplyr::select(-type)
  
  
  class_sub_exp <- lipid_sub_matrix(exp_data = class_exp, 
                                    sub_data = class_sub_stop,
                                    sub_type = 'Class')
  
  
  
  #Differential expression analysis for substructures
  
  class_sub_exp_t <- t_test(data = class_sub_exp[[3]], ctrl = ctrl, exp = exp,
                            method = method, significant = 'adj_p_value')
  
  print('Substructure transformation complete.')
  
  #Class biosynthetic network data transformation
  
  class_network <- network_edge[c('S1','P1')] %>% 
    filter(S1 %in% class_sub_exp_t$lipid, P1 %in% class_sub_exp_t$lipid)
  
  #Essential pathway analysis for class substructures
  
  set.seed(1)
  path_score <-  path_scoring(network = class_network,
                              sub_t = class_sub_exp_t,
                              calibrate = T, data_type = 'Class')
  
  
  path_score_sel <- path_score %>% filter(Significant=='yes') %>% 
    mutate(Type=ifelse(Type=='Active','Increase','Decrease'))
  
  if(nrow(path_score_sel)==0){
    path_data_fig <- NA
    #path_score_sel <- NA
  }
  else{
    path_data <- rbind(path_score_sel %>% 
                         filter(Type=='Decrease') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                       path_score_sel %>% 
                         filter(Type=='Increase') %>% 
                         .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
      filter(Significant=='yes')  %>% 
      dplyr::select(-score)
    
    
    add_suffix <- function(strings) {
      counts <- table(strings)
      duplicated_indices <- which(duplicated(strings))
      count <- 2
      for (index in duplicated_indices) {
        strings[index] <- paste(strings[index], paste0("(", count, ")"), sep = "")
        count <- count + 1
      }
      
      return(strings)
    }
    path_data$path <- str_split(path_data$path, ' --> ') %>% map_chr(~str_c(.x[[1]], ' --> ', last(.x)))
    path_data$path <- add_suffix(path_data$path)
    
    
    path_data_fig <- path_data %>% 
      mutate(path=factor(.$path, levels = .$path)) %>% 
      ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=Type))+
      geom_bar(stat='identity')+
      geom_hline(yintercept = 0)+
      geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
      coord_flip()+
      theme_bw()+
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue"))+
      theme(legend.position='top',
            plot.title = element_text(hjust = 0.5))+
      labs(x='', y='Pathway score',
           title='Top 5 representative pathways')
  }
  
  #path_score_sel <- path_score_sel[,c('path', 'from', 'to', 'cal_score', 'Type', 'rep_sub_path')]
  path_score_sel <- path_score[,c('path', 'from', 'to', 'cal_score','Significant', 'Type', 'rep_sub_path')]
  
  colnames(path_score_sel) <- c('Pathway', 'From', 'To', 'Score','Significant', 'Type', 'Representative pathway')
  
  
  
  print('Pathway analysis complete.')
  
  
  #Essential edges (reactions) analysis for class substructures
  
  reaction_score <- reaction_scoring(network = class_network,
                                     sub_exp = class_sub_exp[[3]],
                                     sub_t = class_sub_exp_t,
                                     ctrl=ctrl, exp=exp,
                                     Species = species)
  
  
  reaction_score_sel <- 
    reaction_score[,c("edge_name", "p_value", "mlog10p",
                      'perturbation_score' ,'Mode', 'genes')] %>% 
    filter(p_value<0.05)
  
  if(nrow(reaction_score_sel)==0){
    reaction_data_fig <- NA
    #reaction_score_sel <- NA
  }
  else{
    reaction_data <- rbind(reaction_score %>% filter(perturbation_score>0) %>% .[1:5,],
                           reaction_score %>% filter(perturbation_score<0) %>% 
                             arrange(perturbation_score) %>%.[1:5,]) %>% 
      filter(p_value<0.05)
    
    
    reaction_data <- reaction_data %>% 
      mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
             node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
      mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                                paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
      mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                                paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
      mutate(edge_color=paste0(node1_color,node2_color))
    
    
    
    reaction_data_fig <- reaction_data %>%
      ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                 fill=Mode))+
      geom_bar(stat='identity', size=0.8)+
      scale_y_discrete(
        labels=rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0)+
      theme_bw()+
      theme(legend.position = 'top', 
            axis.text.y = element_markdown(),
            plot.title = element_text(hjust = 0.5))+
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue"))+
      labs(y='', fill='Reaction', x='Perturbation score', 
           title='Top 5 significant reactions')
    
    #reaction_score_sel <- reaction_score_sel[,c('edge_name', 'p_value', 'perturbation_score', 'Mode', 'genes')]
    
  }
  reaction_score_sel <- reaction_score %>% mutate(Significant=ifelse(p_value<0.05, 'yes','no'))
  reaction_score_sel <- reaction_score_sel[,c('edge_name', 'p_value', 'perturbation_score','Significant', 'Mode', 'genes')]
  
  colnames(reaction_score_sel) <- c('Reaction', 'P-value', 'Perturbation score', 'Significant','Type', 'Gene')
  
  print('Reaction analysis complete.')
  
  
  #Lipid class biosynthetic network construction
  
  class_network_data <- draw_network(network_data = class_network,
                                     DE_data = class_sub_exp_t,
                                     if_species = F,significant = 'adj_p_value',
                                     path_scoring_result = path_score,
                                     reaction_scoring_result = reaction_score,
                                     top_n = 5, path_type = 'both')
  
  network <- visNetwork(class_network_data[[1]], class_network_data[[2]])
  print('Network analysis complete.')
  
  
  #---------------------new code--------------------------------
  
  sub_result <- class_sub_exp_t[,c("lipid", "mean_ctrl", "mean_exp", "FC",
                                   "log2FC","p_value", "adj_p_value", "sig")] %>% 
    `colnames<-`(c('Substructure', 'Mean(ctrl)','Mean(exp)','FC', 'Log2(FC)',
                   'P-value','Adjusted p-value', 'Significance')) %>% 
    arrange(`Adjusted p-value`,desc(Substructure))
  
  network_node <- class_network_data[[1]] %>% left_join(sub_result, by=c('id'='Substructure'))
  network_edge <- class_network_data[[2]] %>% 
    mutate(Reaction=str_c(from, ' --> ', to)) %>% 
    left_join(reaction_score[,c("edge_name", "p_value",
                                'perturbation_score' ,'Mode', 'genes')] %>% 
                `colnames<-`(c('Reaction', 'P-value', 'Perturbation score', 'Type', 'Gene')),
              by='Reaction')  %>% arrange(desc(label))
  
  DE_volcano_plot <- class_sub_exp_t %>%
    mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
    mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
    mutate(Significance=ifelse(sig=='yes', Significance, 'No change')) %>% 
    ggplot(aes(x=log2FC, y=mlog10padj, col=Significance)) +
    geom_jitter(width=0.3) +
    #scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
    scale_color_manual(values = c("Increase" = "red", "Decrease" = "blue", "No change" = "gray"))+
    geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
    theme_classic()+
    theme(legend.position = 'top')+
    labs(y='-Log10 (Adjusted p-value)', x='Log2 (Fold change)')
  
  return(list(sub_result, DE_volcano_plot, path_score_sel, path_data_fig, reaction_score_sel,
              reaction_data_fig, network_node, network_edge, network))
}