



# import("dplyr")

# exportPattern("^[[:alpha:]]+")

# S3method(print, dpp_measure)
# S3method(print, dpp_population)
# S3method(print, dpp_result)

# S3method(print, dpp_degree_cpt)
# S3method(print, dpp_matching_probas)
# S3method(print, dpp_prepared)
# S3method(print, dpp_resolved)
# S3method(print, dpp_sample)



# creates a sample. 
# note a sample is just a specific type of list
create_sample <- function(data, encoding, weight.colname) {

    # ensure the weight column does exist
    # TODO 

    # ensure the inputs are consistent

    encoding.table <- list()

    # build the inverse dictionnary 
    decoding <- list()
    for (attname in names(encoding)) {
        values <- list()
        for (label in names(encoding[[attname]])) {
            code <- encoding[[attname]][[label]]
            
            # add the correspondance with the code used in dataframes
            encoding.table[[attname]][make.names(label)] <- code
            
            values[code] <- label
        }
        decoding[[attname]] <- values
    }

    res <- list(
            sample=data,
            dictionnary=list(
                encoding=encoding,
                encoding.table=encoding.table,
                decoding=decoding,
                colname.weight=weight.colname
                )
    )
    class(res) <- "dpp_sample"

    return(res)

}

print.dpp_sample <- function(x, ...) {
    cat("Sample containing", nrow(x$sample), "elements ")

    cat("having ",length(x$dictionnary$encoding),"columns:",names(x$dictionnary$encoding))
    
    cat(" (weight column:",x$dictionnary$colname.weight,")\n")
}

# creates a table storing probabilities for degrees
create_degree_probabilities_table <- function(probabilities, attributes.names) {
    # TODO check
    res <- list(
            data=probabilities,
            attributes=attributes.names
        )
    class(res) <- "dpp_degree_cpt"

    return(res)
}
print.dpp_degree_cpt <- function(x, ...) {
    cat("distribution of degrees depending to attributes '", x$attributes,"':\n")
    print(x$data)
   
}

create_matching_probabilities_table <- function(data, Ai, Bi) {

    # TODO check input
    res <- list(data=data, Ai=Ai,Bi=Bi)
    class(res) <- "dpp_matching_probas"
    res
}
print.dpp_matching_probas <- function(x, ...) {
    cat("matching probabilities given '", x$Ai,"' and '",x$Bi,"':\n")
    print(x$data)
}
