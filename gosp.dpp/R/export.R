
# ' Takes a population and returns it as an igraph graph
# ' 
# ' @param pop the population to convert
# ' @return an igraph result
# '
# ' @export
# as_igraph <- function (x, ...) {
#    UseMethod("as_igraph", x)
# }

#' @importFrom igraph graph_from_edgelist 
#' @importFrom igraph set_vertex_attr
as.igraph.dpp_population <- function(pop, with.attributes=FALSE, ...) {
    
    # pop$links is a data frame

    g <- graph_from_edgelist(as.matrix(pop$links), directed=TRUE)

    # TODO set_vertex_attr
    if (as.logical(with.attributes)) {

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