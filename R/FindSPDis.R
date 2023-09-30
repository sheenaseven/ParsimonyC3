#' calculate shortest path distance based on Dijkstra's algorithm
#'
#' @param graph An graph created by igraph package
#'
#' @return distance matrix of the nodes in the graph
#' @examples
#' \donttest{
#'
#' }
#'
#' @import igraph
#' @export

FindSPDis = function(graph){
  Dis <- distances(graph,algorithm ="dijkstra", mode = "out")
  return(Dis)
}
