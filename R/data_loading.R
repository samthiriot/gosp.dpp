#
# All the methods related to loading data and making it into the structures
# expected by the rest of the library.
#


#' Extracts attribute values from a string
#'
#' From a string link "att1=v1,att2=v2",
#' returns a named list containing names "att1" and "att2"
#' and values "v1" and "v2". An error (\link{stop}) is raised
#' if the string is malformed.
#'
#' @param name the string to decode. 
#' 
#' @return a named list with attribute as key and value for each  
#'
#' @keywords internal
#'
extract_attributes_values <- function(name) { 

    kv <- unlist(strsplit(name,","))
    
    if (length(kv) < 1) {
        stop("the name should not be empty but was'",name,"'")
    }

    k2v <- strsplit(kv,"=")

    # check the validity (we expect pairs !)
    if (!all(lapply(k2v,length)==2)) {
        stop("invalid name ",name,"; we expect a scheme like att1=v1,att2=vx")
    }

    k <- unlist(lapply(k2v,'[',1))
    v <- unlist(lapply(k2v,'[',2))

    # return as a named list 
    setNames(v,k)
}

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
#' f <- system.file("extdata", "logements.csv", package = "gosp.dpp")
#' m <- read.csv(f, sep=";", dec=",", check.names=FALSE)
#' df <- as.data.frame(m, check.names=FALSE)
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
#'
#' @importFrom stats setNames
#'
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
    }
    for (name in colnames(data)) {

        uniques <- unique(data[,name])
        if (
            # skip the columns already defined
            is.null(encoding[[name]])
            # skip when columns have a huge cardinality,
            # as they might be ids, weights... which don't require a dictionnary
            && (length(uniques) < nrow(data))
            # don't create a dictionary for weights
            && (name != weight.colname)
            ) {
            
                values <- sort(uniques)
                keys <- as.character(values)
                encoding[[name]] <- setNames(values, keys)        

                message("no dictionary provided for the column ", name, " of the sample; creating a dictionary ",paste(encoding[[name]],collapse=","))
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
#' @param x the matching probabilities to print
#' @param ... ignored
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_sample <- function(x, ...) {
    cat("Sample containing ", nrow(x$sample), " elements ",sep="")

    cat("having ",length(x$dictionary$encoding)," columns:",paste(names(x$dictionary$encoding),collapse=","),sep="")
    
    cat(" (weight column:",x$dictionary$colname.weight,")\n",sep="")
}


#' Coerce the a dpp sample into a data frame
#' 
#' Extracts the data frame of a sample prepare by function \code{\link{create_sample}}
#' 
#' @param x the sample to convert
#' @param ... further parameters are ignored quietly
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
as.data.frame.dpp_sample <- function(x, ...) {
    x$sample
}

# TODO manage multiple attributes
# 
#' Creates a table storing probabilities for degrees
#' 
#' @param probabilities a data frame containing the probabilities 
#' @param norm if TRUE, will normalize the table so the columns sum up to 1 (defaults to TRUE)
#' 
#' @examples
#' # create a table describing degrees depending to size; the bigger the size, the highest the degree
#' p <- data.frame(
#'                  'size=0'=c(0.2, 0.8, 0, 0, 0),
#'                  'size=1'=c(0.15, 0.8, 0.05, 0, 0),
#'                  'size=2'=c(0.05, 0.8, 0.1, 0.05, 0),
#'                  check.names=FALSE
#'                  )
#' create_degree_probabilities_table(probabilities=p)
#' 
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
create_degree_probabilities_table <- function(probabilities, norm=TRUE) {

    # check inputs
    if (!is.data.frame(probabilities)) {
        stop("probabilities should be a data frame")
    }

    # ensure the given colnames(cas1$pdi)
    for (name in colnames(probabilities)) {

        # ensure probabilities sum to 1
        if (norm) {
            probabilities[,name] <- normalise(probabilities[,name])
        } else if (abs(sum(probabilities[,name]) - 1.0) >= 0.0000000001) {
            stop(paste("invalid probabilities for ", name, 
                ": should sum to 1 but sums to ", 
                sum(probabilities[,name]), sep=""))
        }

    }

    idx2k2v <- lapply(colnames(probabilities),extract_attributes_values)

    # do all the parameters have the same length of attributes ?
    if (length(unique(lapply(idx2k2v, length))) != 1) {
        stop("all the column names should contain the same count of attributes")
    }
    
    # TODO check attributes !
    
    # list the attribute names concerned by the table
    attributes.names <- unique(unlist(lapply(idx2k2v, names)))

    # forge the result 
    res <- list(
            data=probabilities,
            attributes=attributes.names
        )
    class(res) <- "dpp_degree_cpt"

    return(res)
}

#' Display a sample for Direct Probabilistic Peering
#'
#' @param x the degrees distribution of probabilities to display
#' @param ... ignored 
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_degree_cpt <- function(x, ...) {
    cat("distribution of degrees depending to attributes '", paste(x$attributes, collapse=","),"':\n",sep="")
    print(x$data)
}

#' Coerce a distribution of degrees into a data frame
#' 
#' Extracts the data frame of the distribution of degrees. 
#' 
#' @param x the distribution of degrees created by \code{\link{create_degree_probabilities_table}}
#' @param ... further parameters are ignored quietly
#' 
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
as.data.frame.dpp_degree_cpt <- function(x, ...) {
    x$data
}

# TODO manage multiple attributes
# 
#' Creates a table storing matching probabilities 
#' 
#' 
#' 
#' @param data a data frame containing matching probabilities
#' @param norm if TRUE, will normalize the table so the totals sum up to one (defaults to TRUE)
#' @return a matching probability table ready to be used for usage with \code{\link{matching.prepare}}
#' 
#' @examples
#' 
#' cas1.pij <- create_matching_probabilities_table(
#'                data.frame(
#'                    'surface=1'=c(0.2, 0.1, 0.05, 0.025),
#'                    'surface=2'=c(0.0375, 0.125, 0.1, 0.05),
#'                    'surface=3'=c(0.0125, 0.025, 0.1, 0.175),
#'                    row.names=c("size=1", "size=2", "size=3", "size=4"),
#'                    check.names=FALSE
#'                    )
#'                )
#' print(cas1.pij)
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
create_matching_probabilities_table <- function(data, norm=TRUE) {

    # check inputs
    if (!is.data.frame(data)) {
        stop("data should be a data frame")
    }

    data_normalized <- if (norm) normalise(data) else data

    if (abs(sum(data_normalized) - 1.0) >= 1e-6) {
        stop(paste("the pairing probabilities should sum up to 1.0 but sum up to", sum(data_normalized)))
    }

    # create list of Ai
    Ai.idx2k2v <- lapply(colnames(data_normalized),extract_attributes_values)
    # do all the parameters have the same length of attributes ?
    if (length(unique(lapply(Ai.idx2k2v, length))) != 1) {
        stop("all the column names should contain the same count of attributes")
    }
    # list the attribute names concerned by the table
    Ai <- unique(unlist(lapply(Ai.idx2k2v, names)))


    # create list of Bj
    Bj.idx2k2v <- lapply(rownames(data_normalized),extract_attributes_values)
    # do all the parameters have the same length of attributes ?
    if (length(unique(lapply(Bj.idx2k2v, length))) != 1) {
        stop("all the row names should contain the same count of attributes")
    }
    # list the attribute names concerned by the table
    Bj <- unique(unlist(lapply(Bj.idx2k2v, names)))

    # TODO check input

    res <- list(data=data_normalized, Ai=Ai,Bi=Bj)
    class(res) <- "dpp_matching_probas"

    res
}

#' Display matching probabilities
#' 
#' @param x the matching probabilities to print
#' @param ... ignored
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_matching_probas <- function(x, ...) {
    cat("matching probabilities given '", x$Ai,"' and '",x$Bi,"':\n",sep="")
    print(x$data)
}


#' Coerce pairing probabilities into a data frame
#'
#' Extracts the data frame of pairing probabilities prepared by function \code{\link{create_matching_probabilities_table}}
#'
#' @param x the pairing probabilities to convert
#' @param ... further parameters are ignored quietly
#'
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
as.data.frame.dpp_matching_probas <- function(x, ...) {
    x$data
}
