
#' Measures the contigencies from a sample 
#'
#' measures, for a given variable, the cumulated weights for 
#' each modality of this variable in the sample
#' 
#' @param sample the sample to measure
#' @param strs the vector of strings containing the variable names and modalities to be measured
#' 
#' @return a vector containing the contigencies 
#'  
measure_contigencies <- function(sample, strs) {

    ci <- list()
    
    for (str in strs) {

        # parse the key and value (variable and modality) to be measured

        k2v <- extract_attributes_values(str)

        # identify which lines are relevant 
        selected_ids <- 1:nrow(sample$sample) 
        for (k in names(k2v)) {
            selected_ids <- intersect(selected_ids, which(sample$sample[k] == k2v[[k]]))
        }

        # sum the weights there
        sum_weights <- sum(sample$sample[selected_ids,sample$dictionary$colname.weight])
        ci[[str]] <- sum_weights

    }

    unlist(ci)
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


#' Extracts the list of modalities per variable
#'
#' From a list of strings, extracts the modalities possible 
#' for each variable.
#' 
#' @param strs a list of strings
#' 
#' @return a named list with for each name (variable) contains the vector of modalities
#' 
list_variabes_modalities <- function(strs) {

    idx2k2v <- lapply(strs,extract_attributes_values)

    res <- list()
    for (name in unique(unlist(lapply(idx2k2v, names)))) {

        res[[name]] <- unlist(lapply(idx2k2v, "[[", name))
    }

    res
    
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
#' @return a case ready to be solved with \code{\link{matching.solve}}
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

    # do we have the right types ? 
    if ( (class(sample.A)!="dpp_sample") || (class(sample.B)!="dpp_sample")) {
        stop("sample.A and sample.B should be created using the function create_sample. See ?create_sample .")
    }
    if ( (class(pdi)!="dpp_degree_cpt") || (class(pdi)!="dpp_degree_cpt")) {
        stop("pdi and pdj should be created using the function create_degree_probabilities_table. See ?create_degree_probabilities_table .")
    }
    if (class(pij) != "dpp_matching_probas") {
        stop("pij should be created using the function create_matching_probabilities_table. See ?create_matching_probabilities_table")
    }

    # ensure everything is consistent 

    # ... the list of attributes tackled on sample.A, pdi and colname(pij) should be consistent 
    if (all.equal(pdi$attributes, pij$Ai) != TRUE) {
        stop("the attributes and values defined in pdi (",paste(pdi$attributes,collapse=","),") ",
            "should be the same as the columns of pij (", paste(pij$Ai,collapse=","), ").")
    }
    if (!all(pij$Ai %in% colnames(sample.A$sample))) {
         stop("the Ai attributes defined pdi and pij ",
            "(",paste(pij$Ai,collapse=","),")",
            " should be present in the sample of A (",
            paste(colnames(sample.A$sample),collapse=","),").")
    }
    if (!all.equal(pdj$attributes, pij$Bi)) {
        stop("the attributes and values defined in pdj (",paste(pdj$attributes,collapse=","),") ",
            "should be the same as the rows of pij (", paste(pij$Bi,collapse=","), ").")
    }
    if (!all(pij$Bi %in% colnames(sample.B$sample))) {
        stop("the Bj attributes defined pdj and pij ",
            "(",paste(pij$Bi,collapse=","),")",
            " should be present in the sample of B (",
            paste(colnames(sample.B$sample),collapse=","),").")
    }

    # ... for sure the list of indices should be the very same 
    if (!identical(colnames(pdi$data),colnames(pij$data))) {
        stop("the variables and modalities for pdi and pij should be exactly the same")
    }
    if (!identical(colnames(pdj$data),rownames(pij$data))) {
        stop("the variables and modalities for pdj and pij should be exactly the same")
    }

    # ... the list of values for Ai should be covered everywhere. 
    list_variabes_modalities(colnames(pdi$data))

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
    inputs$ci <- measure_contigencies(sample.A, colnames(pdi$data))
    inputs$cj <- measure_contigencies(sample.B, colnames(pdj$data))
    
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
