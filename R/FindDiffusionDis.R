#' calculate diffusion-kernel based distance
#'
#' @param graph An graph created by igraph package
#' @param weight_name Either NULL or a character string giving an edge attribute name
#' @param a The parameter controls the likelihood that the diffusion process follows a connection originating from a node (default is 0.5)
#'
#' @return distance matrix of the nodes in the graph
#'
#' @examples
#' \donttest{
#'
#' }
#'
#' @import igraph
#' @import dplyr
#' @export

FindDiffusionDis = function(
    graph,
    weight_name = 'weight',
    a = 0.5
) {
  if (a<=0 | a>1) {
    print('error: wrong value for a, the range for a is 0 < a <= 1')
  }
  W <- as_adjacency_matrix(graph,attr = weight_name,sparse = F)
  Dis <- a*W+a**2*(W%*%W)+a**3*(W%*%W%*%W)
  return(Dis)
}
