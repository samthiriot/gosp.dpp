
#' Measures the contigencies from a sample 
#'
#' measures, for a given variable, the cumulated weights for 
#' each modality of this variable in the sample
#' 
#' @param sample the sample to measure
#' @param group.colname the name of the column to check
#' 
#' @return a vector containing the contigencies 
#'  
measure_contigencies <- function(sample, group.colname) {
 
    res <- list()
    #total <- 0
    for(i in 1:nrow(sample$sample)) {
        v = toString(sample$sample[i,group.colname])

        if (is.null(res[[v]])) {
            res[v] <- sample$sample[i,sample$dictionary$colname.weight]
        } else {
            res[[v]] <- res[[v]] + sample$sample[i,sample$dictionary$colname.weight]
        }
        #total <- total + 1
    }

    # normalize
    #for (k in names(res)) {
    #    res[[k]] <- res[[k]]/total
    #}
    
    # order into a vector
    ci <- c() 
    for (code in names(sample$dictionary$encoding[[group.colname]])) {
        idx <- sample$dictionary$encoding[[group.colname]][[code]]
        ci <- c(ci, res[[toString(idx)]])
    }

    ci
}


#' Computes the average degree from a distribution of degrees
#' 
#' @param pd the probability distribution
#' @return a vector containing the average degree
#' 
compute_average_degree <- function(pd) {

    colSums(pd$data * (0:(nrow(pd$data)-1)))

}

#' Computes the maximum degree which might be obtained using a given 
#' distribution of degrees. It corresponds to the maximal average degre
#' which can be obtained by still enforcing the zeros in the initial distribution.
#'
#' @param orig.pd the probability distribution
#' @return a vector containing the maximum possible degree 
#'
assess.max.pd <- function(orig.pd) {
    max.pd <- rep(0,ncol(orig.pd$data))
    for (i in 1:ncol(orig.pd$data)) {
        max.pd[i] <- max(which(orig.pd$data[,i]!=0))-1 # because of indexing starts with 1 but degree 0
    }
    max.pd
}
assess.min.pd <- function(orig.pd) {
    min.pd <- rep(0,ncol(orig.pd$data))
    for (i in 1:ncol(orig.pd$data)) {
        min.pd[i] <- min(which(orig.pd$data[,i]!=0))-1 # because of indexing starts with 1 but degree 0
    }
    min.pd
}


#' Prepares a case for Direct Probabilistic Peering
#'
#' Prepares a case based on two samples and two sets of probabilities.
#' It takes two samples, probabilities distributions for degrees 
#' (count of links per entity) and matching probabilities. 
#' It ensures all the parameters are consistent.
#'
#' @param sample.A the sample to use for population A created using \code{\link{create_sample}}
#' @param sample.B the sample to use for population B created using \code{\link{create_sample}}
#' @param pdi the distribution of degrees for population A created with \code{\link{create_degree_probabilities_table}}
#' @param pdj the distribution of degrees for population B created with \code{\link{create_degree_probabilities_table}}
#' @param pij the matching probabilities created with \code{\link{create_matching_probabilities_table}}
#' @return a case ready to be arbitrated with \code{\link{matching.arbitrate}}
#' 
#' @examples
#' data(cas1)
#' case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#'
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.prepare <- function(sample.A, sample.B, pdi, pdj, pij) {

    inputs <- list(
                pdi=pdi, 
                pdj=pdj,
                pij=pij,
                sample.A=list(dictionary=sample.A$dictionary),
                sample.B=list(dictionary=sample.B$dictionary)
                )

  
    # store as inputs
    
    inputs$pij <- pij
    
    # analyze sample frequencies
    # ... compute frequencies fi    
    inputs$ci <- measure_contigencies(sample.A, pdi$attributes)
    inputs$cj <- measure_contigencies(sample.B, pdj$attributes)
    
    # analyze degrees
    # ... compute average degrees
    inputs$di <- compute_average_degree(pdi)
    inputs$dj <- compute_average_degree(pdj)

    inputs$max.di <- assess.max.pd(inputs$pdi)
    inputs$max.dj <- assess.max.pd(inputs$pdj)
    inputs$min.di <- assess.min.pd(inputs$pdi)
    inputs$min.dj <- assess.min.pd(inputs$pdj)

    # compute basic statistics 
    stats <- list()
    stats$fi <- inputs$ci / sum(inputs$ci)
    stats$fj <- inputs$cj / sum(inputs$cj)
        
    res <- list(    
            inputs=inputs,
            stats=stats
            )
    class(res) <- "dpp_prepared"
    res
}

#' Display a prepared case
#' 
#' @param x the case to print
#' @param ... ignored
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_prepared <- function(x,...) {
    cat("prepared case:n\n")
    print(x$inputs$pdi)
    print(x$inputs$pdj)
    print(x$inputs$pij)
    print(x$stats)
    cat("di:\n")
    print(x$inputs$di)
    # TODO
}
