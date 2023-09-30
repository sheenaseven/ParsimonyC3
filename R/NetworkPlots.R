#' Plotting ParsimonyC3 results with Network Plot
#'
#' @param cci The dataframe with LR interactions
#' @param start the starting celltype name
#' @param lr_select The selected LR combination
#'
#' @return Plotting the selected LR combination
#' @examples
#' \donttest{
#' NetworkPlot(cci,lr_select='PECAM1_CD38|integrin_a4b7_complex_FN1')
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import igraph
#' @import viridis
#' @import ggh4x
#' @import RColorBrewer
#' @import patchwork
#' @export

NetworkPlot = function(
    cci,
    start='Vn_ACKR1',
    lr_select='GRN_SORT1|PECAM1_CD38'
) {
  df_sub=cci[,c("from","to","inter_cc","inter_lr","means_e","FC_ec","weight","weight_d")]
  num_pr <- str_count(lr_select, "\\|")
  if (num_pr == 0) {
    lr_pairs <- lr_select
  } else if (num_pr>0) {
    lr_pairs <- as.vector(str_split_fixed(lr_select, '\\|', num_pr+1))
  }
  df_sub=df_sub[df_sub$inter_lr %in% lr_pairs,]
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
    mode <- "out"
    node_to <- setdiff(unique( c(df_sub$from,df_sub$to) ), start)
    plot_list <- list()
    for ( ct in 1:length(node_to) ) {
      plot_list[[ct]] <- shortest_paths(graph, from = start,to=node_to[ct], mode = mode,algorithm ="dijkstra",output ='both')
    }
    ShortestpathPlot(graph,df_sub,plot_list)
  }
}

#' Eqarrow Plot
#'
#' @param graph An graph created by igraph package
#' @param layout layout
#' @param edge.width edge.width
#' @param edge.lty edge.lty
#' @param edge.arrow.size edge.arrow.size
#' @param vertex.shape vertex.shape
#' @param edge.curved edge.curved
#' @param ... other parameters
#'
#' @return draw Eqarrow Plot
#' @import ggplot2
#' @import dplyr
#' @import igraph
#' @import viridis
#' @import ggh4x
#' @import RColorBrewer
#' @import patchwork
#' @export
EqarrowPlot <- function(graph, layout, edge.width, edge.lty=rep(1, ecount(graph)),
                        edge.arrow.size=rep(1, ecount(graph)),
                        vertex.shape="circle",
                        edge.curved=autocurve.edges(graph), ...) {
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape="none",vertex.label=NA)
  for (e in seq_len(ecount(graph))) {
    graph2 <- delete.edges(graph, E(graph)[(1:ecount(graph))[-e]])
    plot(graph2, edge.lty=edge.lty[e], edge.arrow.size=edge.arrow.size[e],
         edge.curved=edge.curved[e],edge.width=edge.width[e], layout=layout, vertex.shape="none",
         vertex.label=NA, add=TRUE)
  }
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape=vertex.shape, vertex.label.color= "black",
       vertex.label.family="Helvetica",vertex.label.font=6,add=TRUE, ...)
  invisible(NULL)
}

#' ShortestpathPlot
#'
#' @param graph An graph created by igraph package
#' @param df_sub The dataframe with the LR combinations
#' @param plot_list The paths between the starting node to the remaining nodes
#'
#' @return A Shortestpath Plot
#' @import ggplot2
#' @import dplyr
#' @import igraph
#' @import viridis
#' @import ggh4x
#' @import RColorBrewer
#' @import patchwork
#' @export
ShortestpathPlot=function(
    graph,
    df_sub,
    plot_list) {
  E(graph)$weight=df_sub$weight
  ew <- rep(0.5, ecount(graph))
  for (p in 1:length(plot_list)) {
    ew[unlist(plot_list[[p]]$epath)] <- 4
  }
  E(graph)$size <- ew
  ea <- rep(0.3, ecount(graph))
  for (p in 1:length(plot_list)) {
    ea[unlist(plot_list[[p]]$epath)] <- 0.9
  }
  ENDS = ends(graph, E(graph))
  Curv = rep(F, nrow(ENDS))
  for(i in 1:nrow(ENDS)) {
    Curv[i] = are.connected(graph, ENDS[i,2], ENDS[i,1]) }
  color.range <- brewer.pal(nlevels(as.factor(df_sub$inter_lr)), name = "Accent")
  col_l = color.range[as.numeric(as.factor(levels(as.factor(df_sub$inter_lr))))]
  col_e = color.range[as.numeric(as.factor(df_sub$inter_lr))]
  E(graph)$color <- col_e
  EqarrowPlot(graph, layout.auto(graph), edge.arrow.size=ea,
              edge.width=E(graph)$size, vertex.color="#d6d6d6",vertex.frame.width=0.5,
              cex=0.7,edge.curved=Curv,vertex.size=20)
  df_sub$inter_lr = factor(df_sub$inter_lr)
  legend("bottom",bty = "n",
         legend=levels(df_sub$inter_lr),
         fill=col_l, border=NA, xpd=TRUE, inset=c(0, -.1),x.intersp = 0.5)
}
