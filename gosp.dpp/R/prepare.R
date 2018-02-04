
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

# for each class, compute the degree expected 
# so starting with a probability table which described p(d=n|A1)
# it returns the expected probability d(A1)
compute_average_degree <- function(pd) {
 
    # TODO move into data preparation?

    # ensure the given colnames(cas1$pdi)
    for (name in colnames(pd$data)) {
        if (abs(sum(pd$data[name]) - 1.0) >= 0.0000000001) {
            stop(paste(c("invalid probabilities for ", name, ": should sum to 1 but sums to ", sum(pd$data[name]), sep="")))
        }
    }

    # precomputation: sum the expected degree for each class
    # based on the input table
    #expected.degree = list()
    tilded <- c()
    for (name in colnames(pd$data)) {
        sum <- 0;
        for (i in 1:nrow(pd$data)) {
            sum <- sum + (i-1) * pd$data[i,name]
        }
        #expected.degree[name] = sum
        tilded <- c(tilded, sum)
    }

    #expected.degree
    tilded
}

#' Computes the maximum degree which might be obtained using a given 
#' distribution of degrees. It corresponds to the maximal average degre
#' which can be obtained by still enforcing the zeros in the initial distribution.
#'
#' @param orig.pd the probability distribution
#' @return a vector containing the maximum possible degree 
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
