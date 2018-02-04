
#' Exports a population as an igraph graph.
#' 
#' The resulting graph is directed. 
#' Depending to parameters its vertices might be attributed.
# 
#' @param pop the population to convert
#' @param with.attribute if TRUE, the attributes of entities are copied as vertex attributes 
#'
#' @return an igraph graph
#'
#' @export
#'
#' @examples
#' # TODO loading of pop
#' # load the igraph library
#' library(igraph)
#' # generate a population based on sample case 1
#' data(cas1)
#' case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#' disc <- matching.arbitrate(case.prepared, 
#'                     nA=500,nB=400, 
#'                     nu.A=0, phi.A=0, delta.A=1, 
#'                     gamma=1, 
#'                     delta.B=0, phi.B=0, nu.B=0)
#' case <- matching.generate(disc, cas1$sample.A, cas1$sample.B)
#' 
#' # convert the population as a population
#' g <- as.igraph(case$pop, with.attributes=TRUE)
#'
#' # can export the graph using one of the igraph formats
#' write.graph(g,file="mygraph.graphml", format="graphml")
#'
#' # can visualize the population as a graph
#' #lay <- layout_nicely(g)
#' #tkplot(g,layout=lay)
#'
#' # test structural properties
#' is.connected(g)
#' average.path.length(g)
#'
#' # check attributes 
#' # ... list all the attributes defined for verticies 
#' vertex_attr_names(g)
#' # ... view all attributes
#' vertex_attr(g)
#'
#' # detect communities 
#' # TODO !
#'
#' @describeIn as.igraph exports a generation result as an igraph graph
#'
#' @importFrom igraph graph_from_edgelist 
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph as.igraph
#'
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
as.igraph.dpp_population <- function(pop, with.attributes=FALSE, ...) {
    
    # pop$links is a data frame that requires convertion to a matrix

    g <- graph_from_edgelist(as.matrix(pop$links), directed=TRUE)

    if (as.logical(with.attributes)) {

        # first define values to NA 
        # ... attributes of A which are not defined
        for (name in colnames(pop$B)) {
            g <- set_vertex_attr(g, name, index=pop$A$id.A, value=rep(NA, times=nrow(pop$A)))
        }
        # ... attributes of B which are not defined
        for (name in colnames(pop$A)) {
            g <- set_vertex_attr(g, name, index=pop$B$id.B, value=rep(NA, times=nrow(pop$B)))   
        } 

        # copy the content of attributes of A
        for (name in colnames(pop$A)) {
            g <- set_vertex_attr(g, name, index=pop$A$id.A, value=pop$A[,name])   
        }
        
        # attributes of B
        for (name in colnames(pop$B)) {
            g <- set_vertex_attr(g, name, index=pop$B$id.B, value=pop$B[,name])    
        }

    }

    g

}

#' Exports a generation result as an igraph graph. 
#' 
#' @describeIn as.igraph exports a generation result as an igraph graph
#'
#'
#' @export
#'
#' @importFrom igraph as.igraph
#'
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
as.igraph.dpp_result <- function(generated, with.attributes=FALSE, ...) {

    as.igraph(generated$pop, with.attributes=with.attributes)

}