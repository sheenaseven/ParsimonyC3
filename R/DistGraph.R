#' choose the different strategies to calculate the distance
#'
#' @param graph An graph created by igraph package
#' @param method 'shortest.paths' or 'diffusion'
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

DistGraph = function(
    graph,
    method = 'shortest.paths',
    weight_name = 'weight',
    a = 0.5
) {
  if (method == 'shortest.paths') {
    distance <- FindSPDis(graph)
  }else if (method == 'diffusion') {
    distance <- FindDiffusionDis(graph,weight_name = weight_name,a = a)
  }else {
    print('error: wrong parameter, only "shortest.paths" and "diffusion" are supported!')
  }
  return(distance)
}
