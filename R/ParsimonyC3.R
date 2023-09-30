#' ParsimonyC3 methods
#'
#' @param cci The dataframe with LR interactions
#' @param lr_combination The ligand-receptor(LR) combination list
#' @param start the starting celltype name
#' @param b  A parameter that adjusts the relative strength in the shortest path mode
#' @param r  A parameter that adjusts the relative strength in the diffusion-kernel based mode
#'
#' @return The priortized LR pair combinations. c_lr: LR combination ID; z_shortest: z score (can be found in the Paper Methods) in shortest path mode;
#' z_diffusion: z score in diffusion-kernel based mode; n_lr: number of LRs in the specific combination
#' @examples
#' \donttest{
#' result <- ParsimonyC3(cci,lr_combination)
#' }
#'
#' @import igraph
#' @import dplyr
#' @export

ParsimonyC3 = function(
    cci,
    lr_combination,
    start = 'Vn_ACKR1',
    b = 0.02,
    r = 1
) {
  c_lr=c();c_shortest=c();c_diffusion=c();f_shortest=c();f_diffusion=c();n_lr=c()
  for (m in lr_combination) {
    print(m)
    num_pr <- str_count(m, "\\|")
    if (num_pr == 0) {
      lr_pairs <- m
    } else if (num_pr>0) {
      lr_pairs <- as.vector(str_split_fixed(m, '\\|', num_pr+1))
    }
    sub_cci <- cci[cci$inter_lr %in% lr_pairs,]
    coln <- c("from","to","inter_cc","inter_lr","means_e","FC_ec","weight","weight_d")
    if ( ! all( coln %in% colnames(cci) ) ) {
      print('error: columns "from","to","inter_cc","inter_lr","means_e","FC_ec","weight","weight_d" not in cci')
    } else {
      df_sub <- sub_cci[,c("from","to","inter_cc","inter_lr","means_e","FC_ec","weight","weight_d")]
      check_dup <- data.frame()
      for (mc in unique(df_sub$inter_cc)) {
        df_dup = as.data.frame(df_sub[df_sub$inter_cc==mc,])
        if (dim(df_dup)[1]>=1) {
          tmp <- df_dup[1,]
          tmp[,c("means_e","FC_ec","weight","weight_d")] <- colSums(df_dup[,c("means_e","FC_ec","weight","weight_d")])
          check_dup <- rbind(check_dup,tmp)
        } else {
          check_dup <- rbind(check_dup,df_sub[df_sub$inter_cc==mc,])
        }
      }
      if ( sum(duplicated(check_dup$inter_cc)) != 0 ) {
        print('error: data still have duplicates')
      } else {
        df_sub <- check_dup
        node_n <- setdiff(unique(cci$from),c(unique(df_sub$from),unique(df_sub$to)))
        df_sub <- as.data.frame(df_sub)
        if ( length(node_n)>=1 ) {
          for ( i in 1:length(node_n) ) {
            df_sub[nrow(df_sub) + 1,] <- c(node_n[i], node_n[i],paste0(node_n[i],'|',node_n[i]),unique(df_sub$inter_lr[1]),0,0,0,0)
          }
        }
        df_sub$weight <- as.numeric(df_sub$weight)
        df_sub$weight_d <- as.numeric(df_sub$weight_d)
        graph <- igraph::graph_from_data_frame(df_sub,directed = TRUE)
        shortest_dis <- DistGraph(graph,method = 'shortest.paths')
        if (! start %in% colnames(shortest_dis)) {
          print('error: start cell type was not found!')
        } else {
          shortest_dis <- shortest_dis[,!grepl(start,colnames(shortest_dis))]
          diffusion_dis <- DistGraph(graph,method='diffusion',weight_name = 'weight_d',a=0.5)
          diffusion_dis <- diffusion_dis[,!grepl(start,colnames(diffusion_dis))]
          if ( sum( !is.finite( shortest_dis[start,] ) )>0 ) {
            next; } else {
              shortest_sum <- sum( as.vector( shortest_dis[start,] ) )
            }
          if ( sum( diffusion_dis[start,]==0 )>0 ){
            next; } else {
              diffusion_sum <- sum( as.vector( diffusion_dis[start,] ) )
            }
          z_shortest <- length(lr_pairs) + b*shortest_sum
          z_diffusion <- length(lr_pairs) + r*(-1*log10(diffusion_sum))
          c_lr <- c(c_lr, paste(lr_pairs, collapse = "|"))
          n_lr <- c(n_lr,length(lr_pairs))
          c_shortest <- c( c_shortest, shortest_sum )
          c_diffusion <- c( c_diffusion, (-1*log10(diffusion_sum)) )
          f_shortest <- c( f_shortest, z_shortest)
          f_diffusion <- c(f_diffusion, z_diffusion)
        }
      }
    }
  }
  res <- data.frame(pairs=c_lr,
                    z_shortest=f_shortest,
                    z_diffusion=f_diffusion,
                    n_lr=n_lr)
  return(res)
}
