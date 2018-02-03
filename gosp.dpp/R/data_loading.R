#
# All the methods related to loading data and making it into the structures
# expected by the rest of the library.
#

#' Creates a sample organized for dpp manipulation. 
#' 
#' The resulting object contains the sample and a dictionary of data. 
#' If no weight column is provided, one will be created with value 1 (uniform weights).
#' 
#' @param data a data frame containing a sample (weighted list) of entities
#' @param encoding a dictionary containing information about the variables of the sample
#' @param weight.colname the name of the column containing the weights
#' @return a sample 
#' 
#' @examples 
#' # to read a CSV file as a sample
#' f <- system.file("data-raw", "logements.csv", package = "gosp.dpp")
#' m <- read.csv(f, sep=";", dec=",")
#' df <- as.data.frame(m)
#' dictionary <- list('surface'=list('small'=1, 'medium'=2, 'large'=3))
#' create_sample(data=df, encoding=dictionary, weight.colname="weight")
#' 
#' # to create a sample from random data
#' # ... create 100 entities being either male of female
#' df <- data.frame(gender=sample(1:2, size=100, replace=TRUE))
#' # ... describe the encoding of data
#' dictionary <- list("gender"=list("male"=1,"female"=2))
#' create_sample(data=df, encoding=dictionary)
#' # ... the weights columns was created automatically
#' 
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
create_sample <- function(data, encoding=NULL, weight.colname=NULL) {

    # test types if parameters
    if (!is.data.frame(data)) {
        stop("data is expected to be a dataframe")
    }

    if (is.null(weight.colname)) {
        # if no weight column is given, then create a weight column filled by 1.0 
        weight.colname <- "_auto_weight"
        data[ , weight.colname] <- rep(1, times=nrow(data))
    } else {
        # ensure the weight column does exist in the dataset
        if (!weight.colname %in% colnames(data)) {
            stop(paste("There is no column weight.colname='",weight.colname,"' in the data",sep=""))
        }
    }

    # if no encoding is provided, then create one with a direct mapping
    if (is.null(encoding)) {
        encoding <- list()
        for (name in colnames(data)) {

            uniques <- unique(data[,name])
            if (
                # skip when columns have a huge cardinality,
                # as they might be ids, weights... which don't require a dictionnary
                (length(uniques) < nrow(data))
                # don't create a dictionary for weights
                && (name != weight.colname)
                ) {
                
                    values <- sort(uniques)
                    keys <- as.character(values)
                    encoding[[name]] <- setNames(values, keys)        
            }

        }
    }

    # ensure the inputs are consistent

    encoding.table <- list()

    # build the inverse dictionary 
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
            dictionary=list(
                encoding=encoding,
                encoding.table=encoding.table,
                decoding=decoding,
                colname.weight=weight.colname
                )
    )
    class(res) <- "dpp_sample"

    return(res)

}

#' Display a sample for Direct Probabilistic Peering
#' 
#' @describeIn print prints a sample prepared for Direct Probabilistic Peering
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_sample <- function(x, ...) {
    cat("Sample containing", nrow(x$sample), "elements ")

    cat("having ",length(x$dictionary$encoding),"columns:",names(x$dictionary$encoding))
    
    cat(" (weight column:",x$dictionary$colname.weight,")\n")
}

# TODO manage multiple attributes
# 
#' Creates a table storing probabilities for degrees
#' 
#' @param probabilities a data frame containing the probabilities 
#' @param attributes.names the vector of the attributes refered to in the probabilities
#' 
#' @examples
#' 
#' p <- data.frame(
#'                  'small'=c(0.2, 0.8, 0, 0, 0),
#'                  'medium'=c(0.15, 0.8, 0.05, 0, 0),
#'                  'large'=c(0.05, 0.8, 0.1, 0.05, 0)
#'                  ),
#' create_degree_probabilities_table(probabilities=p, attributes.names="size")
#' 
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
create_degree_probabilities_table <- function(probabilities, attributes.names) {
    # TODO check
    res <- list(
            data=probabilities,
            attributes=attributes.names
        )
    class(res) <- "dpp_degree_cpt"

    return(res)
}

#' Display a sample for Direct Probabilistic Peering
#' 
#' @describeIn print prints a degrees distribution of probabilities for Direct Probabilistic Peering
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
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
