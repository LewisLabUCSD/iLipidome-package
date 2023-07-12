#' Website function for fatty acid substructure analysis
#' 
#' @param exp_raw raw expression file
#' @param method analysis method
#' @param ctrl range of values for control (1:7)
#' @param exp range of values for exp (8:13)
#' @param unmapped_FA vector of lowly expressed fatty acid isomers to exclude
#' @param exo_lipid vector of fatty acids selected to prevent substructure decomposition
#' @param species the analysis species, (human, mouse, rat)
#' @param add_reaction add a reaction
#' @param delete_reaction delete a reaction
#' 
#' @import ggplot2
#' @importFrom ggtext element_markdown
#' @importFrom stats reorder
#' @importFrom visNetwork visIgraphLayout
#' @importFrom visNetwork visNetwork
#' 
#' @return list of tables and plots relevant to the analysis
#' 
#' @export

FA_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                     unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                     exo_lipid='w3-22:6;0',species='rat',
                                     add_reaction=NULL, delete_reaction=NULL){
  
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]
  
  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]
  
  no_sub_t <- unprocessed_data_test(exp_data = exp_data,
                                    char_table = char_data,
                                    method = method,
                                    significant='adj_p_value',
                                    ctrl_group = ctrl, exp_group = exp)
  
  #-------------------FA substructure analysis-------------------
  
  #FA biosynthetic network data transformation
  if(!is.null(add_reaction)){
    colnames(add_reaction) <- c('S1','P1','pathway')
    add_reaction <- add_reaction %>% mutate(S1_detail=NA, P1_detail=NA) %>% 
      dplyr::select(everything(), pathway)
    FA_network <- rbind(FA_network, add_reaction) %>% unique()
  }
  if(!is.null(delete_reaction)){
    delete_reaction <- str_c(delete_reaction[[1]],'_', delete_reaction[[2]])
    FA_network <- FA_network %>% mutate(reaction=str_c(S1, '_', P1)) %>% 
      filter(!reaction %in% delete_reaction) %>% 
      dplyr::select(-reaction) %>% unique()
  }
  FA_network_new <- build_FA_net(FA_network = FA_network,
                                 unprocessed_data_result = no_sub_t)
  
  #Decompose lipids into FA substructures
  #18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them
  
  FA_substructure <- FA_sub_transform(FA_network = FA_network_new,
                                      unprocessed_data_result = no_sub_t,
                                      unmapped_FA = unmapped_FA)
  
  #Extract FA substructures using fold changes
  
  FA_sub_stop <- FA_sub_extract(char_table = char_data,
                                FA_substructure = FA_substructure,
                                unprocessed_data_result = no_sub_t,
                                exact_FA='no', exo_lipid=exo_lipid)
  
  #Transform FA exp into substructure exp
  
  FA_sub_exp <- lipid_sub_matrix(exp_data = exp_data, sub_data = FA_sub_stop,
                                 sub_type = 'FA')
  
  
  #Differential expression analysis for FA substructures
  
  FA_sub_exp_t <- t_test(data = FA_sub_exp[[3]], ctrl = ctrl, exp = exp,
                         method = method, significant = 'adj_p_value')
  
  print('Substructure transformation complete.')
  #Essential pathway analysis for FA substructures
  
  set.seed(1)
  
  path_score_FA <- path_scoring(network = FA_network_new, sub_t = FA_sub_exp_t, 
                                calibrate = T, data_type = 'FA')
  path_score_FA_sel <- 
    path_score_FA %>% filter(Significant=='yes') %>% 
    mutate(Type=ifelse(Type=='Active','Increase','Decrease'))
  
  if(nrow(path_score_FA_sel)!=0){
    path_data <- rbind(path_score_FA_sel %>% 
                         filter(Type=='Decrease') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                       path_score_FA_sel %>% 
                         filter(Type=='Increase') %>% 
                         .[!duplicated(.$rep_sub_path),]%>% .[1:5,])%>% 
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
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'top')+
      labs(x='', y='Pathway score',
           title='Top 5 representative pathways')
    
    #path_score_FA_sel <- path_score_FA_sel[,c('path', 'from', 'to', 'cal_score', 'Type', 'rep_sub_path')]
    
  }
  else{
    path_data_fig <- NA
    #path_score_FA_sel <- NA
  }
  
  path_score_FA_sel <- path_score_FA[,c('path', 'from', 'to', 'cal_score','Significant', 'Type', 'rep_sub_path')]
  
  colnames(path_score_FA_sel) <- c('Pathway', 'From', 'To', 'Score','Significant', 'Type', 'Representative pathway')
  
  
  print('Pathway analysis complete.')
  
  #Essential edges (reactions) analysis for FA substructures
  
  reaction_score_FA <- reaction_scoring(network = FA_network_new, 
                                        sub_exp = FA_sub_exp[[3]],
                                        sub_t = FA_sub_exp_t, 
                                        ctrl = ctrl, exp = exp, 
                                        Species = species)
  
  reaction_score_FA_sel <- 
    reaction_score_FA[,c("edge_name", "p_value", "mlog10p",
                         'perturbation_score' ,'Mode', 'genes')] %>% 
    filter(p_value<0.05)
  
  if(nrow(reaction_score_FA_sel)!=0){
    reaction_data <- rbind(reaction_score_FA %>% filter(perturbation_score>0) %>% .[1:5,],
                           reaction_score_FA %>% filter(perturbation_score<0) %>% 
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
           title='Top 5 significant reations')
    
    #reaction_score_FA_sel <- reaction_score_FA_sel[,c('edge_name', 'p_value', 'perturbation_score', 'Mode', 'genes')]
    
  }
  else{
    #reaction_score_FA_sel <- NA
    reaction_data_fig <- NA
  }
  
  reaction_score_FA_sel <- reaction_score_FA %>% mutate(Significant=ifelse(p_value<0.05, 'yes','no'))
  reaction_score_FA_sel <- reaction_score_FA_sel[,c('edge_name', 'p_value', 'perturbation_score','Significant', 'Mode', 'genes')]
  
  colnames(reaction_score_FA_sel) <- c('Reaction', 'P-value', 'Perturbation score', 'Significant','Type', 'Gene')
  
  print('Reaction analysis complete.')
  
  #FA biosynthetic network construction
  FA_network_data <- draw_network(network_data = FA_network_new,
                                  DE_data = FA_sub_exp_t,
                                  if_species = F, significant = 'adj_p_value',
                                  path_scoring_result = path_score_FA,
                                  reaction_scoring_result = reaction_score_FA,
                                  top_n = 5, path_type = 'both')
  
  
  network <- visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>% 
    visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                    physics = F, smooth = TRUE, randomSeed =5) 
  
  print('Network analysis complete.')
  
  #---------------------new code--------------------------------
  
  sub_result <- FA_sub_exp_t[,c("lipid", "mean_ctrl", "mean_exp", "FC",
                                "log2FC","p_value", "adj_p_value", "sig")] %>% 
    `colnames<-`(c('Substructure', 'Mean(ctrl)','Mean(exp)','FC', 'Log2(FC)',
                   'P-value','Adjusted p-value', 'Significance')) %>% 
    arrange(`Adjusted p-value`,desc(Substructure))
  
  network_node <- FA_network_data[[1]] %>% left_join(sub_result, by=c('id'='Substructure'))
  network_edge <- FA_network_data[[2]] %>% 
    mutate(Reaction=str_c(from, ' --> ', to)) %>% 
    left_join(reaction_score_FA[,c("edge_name", "p_value",
                                   'perturbation_score' ,'Mode', 'genes')] %>% 
                `colnames<-`(c('Reaction', 'P-value', 'Perturbation score', 'Type', 'Gene')),
              by='Reaction')  %>% arrange(desc(label))
  
  DE_volcano_plot <- FA_sub_exp_t %>%
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
  
  return(list(sub_result,DE_volcano_plot, path_score_FA_sel, path_data_fig, reaction_score_FA_sel,
              reaction_data_fig, network_node,network_edge, network))
}