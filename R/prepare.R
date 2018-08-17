
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
#'
list_variabes_modalities <- function(strs) {

    idx2k2v <- lapply(strs,extract_attributes_values)

    res <- list()
    for (name in unique(unlist(lapply(idx2k2v, names)))) {

        res[[name]] <- unlist(lapply(idx2k2v, "[[", name))
    }

    res
    
}

#' Expands a data frame by adding combinations with a novel variable
#'
#' Takes the data frame df passed as a parameter, and 
#' adds to columns the combinations of existing values and the novel
#' variable values. 
#' For instance if the columns where a=1 and a=2, and you add name=b
#' and values=list(3,4), the columns will be a=1,b=3 a=1,b=4 a=2,b=3
#' a2,b=4
#' 
#' @param df a data frame
#' @param name a string
#' @param values a list of values
#' @return the same dataframe with additional columns minus the original columns
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variable_col.data.frame <- function(df, name, values) {

    to_pair <- function(v) {
        paste(",", name, "=", v, sep="")
    } 
    names_to_add <- sapply(values, to_pair)

    df2 <- df
    #df2 <- data.frame(check.names=FALSE, check.rows=FALSE,row.names=row.names(df))

    # duplicate existing columns with new col names
    for (q in colnames(df)) {
        for (n in names_to_add) {
            coln <- paste(q,n,sep="")
            #cat(q, " -> ", coln, "\n")
            df2[coln] <- df[[q]]
        }
        # delete former column
        df2[q] <- NULL
    }

    df2
}

#' Expands a data frame by adding combinations with a novel variable
#'
#' Takes the data frame df passed as a parameter, and 
#' adds to rows the combinations of existing values and the novel
#' variable values. 
#' For instance if the rows where a=1 and a=2, and you add name=b
#' and values=list(3,4), the rows will be a=1,b=3 a=1,b=4 a=2,b=3
#' a2,b=4
#' 
#' @param df a data frame
#' @param name a string
#' @param values a list of values
#' @return the same dataframe with additional rows minus the original rows
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variable_row.data.frame <- function(df, name, values) {

    to_pair <- function(v) {
        paste(",", name, "=", v, sep="")
    } 
    names_to_add <- sapply(values, to_pair)

    df2 <- df
    #df2 <- data.frame(check.names=FALSE, check.rows=FALSE,row.names=row.names(df))

    to_delete <- rownames(df2)

    # duplicate existing columns with new col names
    for (q in rownames(df)) {
        for (n in names_to_add) {
            rown <- paste(q,n,sep="")
            #cat("line: ", q, " -> ", rown, "\n")
            #print(df2[q,])
            df2[rown,] <- df2[q,]
        }
    }
    
    # delete former columns
    df2[!(rownames(df2) %in% to_delete),]

}

#' Expands a conditional probability table for degrees with an additional variable
#'
#' Takes the conditional probability table for degrees passed as a parameter, and 
#' adds to columns the combinations of existing values and the novel
#' variable values. 
#' For instance if the columns where a=1 and a=2, and you add name=b
#' and values=list(3,4), the columns will be a=1,b=3 a=1,b=4 a=2,b=3
#' a2,b=4
#' 
#' @param cpt a conditional probability table for degrees
#' @param name a string
#' @param values a list of values
#' @return the CPT with additional columns minus the original columns
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variable.dpp_degree_cpt <- function(cpt, name, values) {

    if (class(cpt) != "dpp_degree_cpt")
        stop("cpt should be a degree Conditional Probability Table.")

    cpt$data <- expand_with_variable_col.data.frame(cpt$data, name, values)
    cpt$attributes[length(cpt$attributes)+1] <- name

    cpt
}

#' Expands a conditional probability table for degrees with additional variables
#'
#' Takes the conditional probability table for degrees passed as a parameter, 
#' and extends it with various variables
#' 
#' @param cpt a conditional probability table for degrees
#' @param names a list of strings
#' @param dictionary a dictionary in which the domains of each variable of names is described
#' @return the CPT with additional columns minus the original columns
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variables.dpp_degree_cpt <- function(cpt, names, dictionary) {

    if (class(cpt) != "dpp_degree_cpt")
        stop("cpt should be a degree Conditional Probability Table.")

    cpt2 <- cpt
    for (n in names) {
        values <- dictionary$encoding[[n]] #levels
        if (is.null(values)) {
            stop("there is no variable ",n," in the dictionary")
        }
        cat("adding for ",n," values ",values,"\n")
        cpt2 <- expand_with_variable.dpp_degree_cpt(cpt2, n, values)
    }
    cpt2
}

#' Expands the columns of matching probabilities with additional variables
#'
#' Takes the probability table for matching passed as a parameter, 
#' and extends its columns with a novel variable
#' 
#' @param pt a probability table for degrees
#' @param name a string
#' @param values a list of values
#' @return the CPT with additional columns minus the original columns
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variable.dpp_matching_probas.A <- function(pt, name, values) {

    if (class(pt) != "dpp_matching_probas")
        stop("pt should be a matching Probability Table.")

    pt$data <- expand_with_variable_col.data.frame(pt$data, name, values)
    pt$Ai[length(pt$Ai)+1] <- name

    pt
}

#' Expands the columns of matching probabilities with additional variables
#'
#' Takes the probability table for degrees passed as a parameter, 
#' and extends it with various variables
#' 
#' @param pt a probability table for degrees
#' @param names a list of strings
#' @param dictionary a dictionary in which the domains of each variable of names is described
#' @return the CPT with additional columns minus the original columns
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variables.dpp_matching_probas.A <- function(pt, names, dictionary) {

    if (class(pt) != "dpp_matching_probas")
        stop("cpt should be a matching Conditional Probability Table.")

    pt2 <- pt
    for (n in names) {
        values <- dictionary$encoding[[n]] #levels
        if (is.null(values)) {
            stop("there is no variable ",n," in the dictionary")
        }
        cat("adding for", n, "values", values, "\n")
        pt2 <- expand_with_variable.dpp_matching_probas.A(pt2, n, values)
    }
    pt2
}

#' Expands the rows of matching probabilitis with additional variables
#'
#' Takes the probability table for matching passed as a parameter, 
#' and extends its rows with a novel variable
#' 
#' @param pt a  probability table for degrees
#' @param name a string
#' @param values a list of values
#' @return the table with additional rows minus the original rows
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variable.dpp_matching_probas.B <- function(pt, name, values) {

    if (class(pt) != "dpp_matching_probas")
        stop("pt should be a matching Probability Table.")

    pt$data <- expand_with_variable_row.data.frame(pt$data, name, values)

    pt$Bi[length(pt$Bi)+1] <- name

    pt
}

#' Expands the rows of matching probabilities with additional variables
#'
#' Takes the probability table for matching passed as a parameter, 
#' and extends it with various variables
#' 
#' @param pt a probability table for degrees
#' @param names a list of strings
#' @param dictionary a dictionary in which the domains of each variable of names is described
#' @return the table with additional rows minus the original rows
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
expand_with_variables.dpp_matching_probas.B <- function(pt, names, dictionary) {

    if (class(pt) != "dpp_matching_probas")
        stop("cpt should be a matching Probability Table.")

    pt2 <- pt
    for (n in names) {
        values <- dictionary$encoding[[n]] #levels
        if (is.null(values)) {
            stop("there is no variable ",n," in the dictionary")
        }
        cat("adding for ",n," values ",values,"\n")
        pt2 <- expand_with_variable.dpp_matching_probas.B(pt2, n, values)
    }
    pt2
}


#' Fix inputs parameters for DPP resolution
#' 
#' This function can fix what is not constitent in the naming 
#' of attributes, such as attributes not provided in pdi, pdj or pij
#' but provided for another item. 
#' 
#' This function is called by \code{\link{matching.prepare}} if the optional parameter \code{fix=TRUE} is used. 
#' 
#' @param sample.A the sample to use for population A created using \code{\link{create_sample}}
#' @param sample.B the sample to use for population B created using \code{\link{create_sample}}
#' @param pdi the distribution of degrees for population A created with \code{\link{create_degree_probabilities_table}}
#' @param pdj the distribution of degrees for population B created with \code{\link{create_degree_probabilities_table}}
#' @param pij the matching probabilities created with \code{\link{create_matching_probabilities_table}}
#' @return a named list containing the fixed variables
#' 
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.fix <- function(sample.A, sample.B, pdi, pdj, pij) {
    
    res_pdi <- pdi
    res_pdj <- pdj
    res_pij <- pij

    res_pdi <- remove_excess_probas_pd(res_pdi, sample.A$dictionary)
    res_pdj <- remove_excess_probas_pd(res_pdj, sample.B$dictionary)
    res_pij <- remove_excess_probas_pd(res_pij, sample.A$dictionary)
    res_pij <- remove_excess_probas_pd_lines(res_pij, sample.B$dictionary)
    # TODO pij

    mods_pdi <- res_pdi$attributes
    mods_pdj <- res_pdj$attributes
    mods_pij_A <- pij$Ai
    mods_pij_B <- pij$Bi

    # in case probability distributions do miss variables of pij, add them
    missing_pdi <- setdiff(mods_pij_A, mods_pdi)
    if (length(missing_pdi) > 0) {
        message(paste("some attributes are missing in pdi: ", paste(missing_pdi, collapse=","), "; we will add them", sep=""))
        res_pdi <- expand_with_variables.dpp_degree_cpt(res_pdi, missing_pdi, sample.A$dictionary)
    }

    missing_pdj <- setdiff(mods_pij_B,mods_pdj)
    if (length(missing_pdj) > 0) {
        message(paste("some attributes are missing in pdj: ", paste(missing_pdi, collapse=","), "; we will add them", sep=""))
        res_pdj <- expand_with_variables.dpp_degree_cpt(res_pdj, missing_pdj, sample.B$dictionary)
    }

    # in case one dimension of pij is missing variables of pij, add them
    missing_pij_A <- setdiff(mods_pdi,mods_pij_A)
    if (length(missing_pij_A) > 0) {
        message(paste("some attributes are missing in pij: ", paste(missing_pij_A, collapse=","), "; we will add them", sep=""))
        res_pij <- expand_with_variables.dpp_matching_probas.A(res_pij, missing_pij_A, sample.A$dictionary)
    }

    missing_pij_B <- setdiff(mods_pdj,mods_pij_B)
    if (length(missing_pij_B) > 0) {
        message(paste("some attributes are missing in pij: ", paste(missing_pij_B, collapse=","), "; we will add them", sep=""))
        res_pij <- expand_with_variables.dpp_matching_probas.B(res_pij, missing_pij_B, sample.B$dictionary)
    }

    # in case the order of variables in pdi and pdj is not the same as in pij, switch it 
    # TODO

    #print(names(list_variabes_modalities(colnames(res_pij$data)[1])))

    res_pdi <- reorder_column_modalities(res_pdi, names(list_variabes_modalities(colnames(res_pij$data)[1])))
    res_pdi$data <- reorder_columns(res_pdi$data, colnames(res_pij$data))

    res_pdj <- reorder_column_modalities(res_pdj, names(list_variabes_modalities(rownames(res_pij$data)[1])))
    res_pdj$data <- reorder_columns(res_pdj$data, rownames(res_pij$data))

    list(sample.A=sample.A, sample.B=sample.B, pdi=res_pdi, pdj=res_pdj, pij=res_pij)
}


#' Removes useless modalities from a probability distribution
#' 
#' Often, the user provides probabilities to deal with all the possible modalities
#' in its dataset; but sometimes the sample just does not contains some of the 
#' modalities which might exist. This function takes as a reference the dictionary 
#' assuming it contains only what is required. In removes all the columns related
#' to modalities not really found in the probability table.
#' 
#' If data is removed, a warning will be emitted.
#' 
#' @param pd a probability distribution
#' @param dictionary a dictionary
#' @return the same probability distribution with maybe less elements
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
remove_excess_probas_pd <- function(pd, dictionary) {

    for (n in colnames(pd$data)) {
        # print(n)
        k2v <- list_variabes_modalities(n)
        for (k in names(k2v)) {
            v <- k2v[k]
            # print("k,v")
            # print(k)
            # print(v)
            # print(dictionary$encoding[[k]])
            if (!(v %in% dictionary$encoding[[k]])) {
                warning(paste("modality ",k,"=",v," is present in parameters but not in the sample, removing it", sep=""))
                pd$data[n] <- NULL
            }
        }
    }
    pd
    
}

#' Removes useless modalities from a probability distribution
#' 
#' Often, the user provides probabilities to deal with all the possible modalities
#' in its dataset; but sometimes the sample just does not contains some of the 
#' modalities which might exist. This function takes as a reference the dictionary 
#' assuming it contains only what is required. In removes all the rows related
#' to modalities not really found in the probability table.
#' 
#' If data is removed, a warning will be emitted.
#' 
#' @param pd a probability distribution
#' @param dictionary a dictionary
#' @return the same probability distribution with maybe less elements
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
remove_excess_probas_pd_lines <- function(pd, dictionary) {
    
    lines_to_remove <- list()

    for (n in rownames(pd$data)) {
        # print(n)
        k2v <- list_variabes_modalities(n)
        for (k in names(k2v)) {
            v <- k2v[k]
            # print("k,v")
            # print(k)
            # print(v)
            # print(dictionary$encoding[[k]])
            if (!(v %in% dictionary$encoding[[k]])) {
                warning(paste("modality ",k,"=",v," is present in parameters but not in the sample, removing it", sep=""))
                #pd$data[n] <- NULL
                lines_to_remove[length(lines_to_remove)+1] <- n
            }
        }
    } 

    res <- pd
    res$data <- res$data[!(rownames(pd$data) %in% lines_to_remove),]
    
    res
}

#' Reorders the columns of a dataset so they are in the expected order
#' 
#' The user, or the automatic fixing process, might lead to situations where 
#' a pdi has for variable list(a=1, a=2) whilst pij has list(a=2,a=1). 
#' As all the solving process relies on the hypothesis the order is OK, 
#' checking of such a dataset would fail.
#' 
#' This function helps fixing these situations by just changing the order of
#' the dataset's columns
#' 
#' @param df a data frame
#' @param expected_order a list of strings corresponding to the same col names with a different order
#' @return the data frame in a different order
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
reorder_columns <- function (df, expected_order) {

    current_order <- colnames(df)
    find_idx <- function(n) {
        which(current_order==n)[1]
    }
    # print("expected")
    # print(expected_order)
    # print("current")
    # print(current_order)
    order <- sapply(expected_order, find_idx, simplify=T)
    order <- order[which(!is.na(order))]

    # print(order)

    # print("so before:")
    # print(colnames(df))
    # print("then")
    df[,order]
    # print(colnames(res))
    # res 
}

#' Reorders the modalities inside columns of a dataset
#' 
#' The user, or the automatic fixing process, might lead to situations where 
#' a pdi has for variable "b=1,a=2" whilst the order is "a=2,b=1" in other datasets.
#' As all the solving process relies on the hypothesis the order is OK, 
#' checking of such a dataset would fail.
#' 
#' This function helps fixing these situations by just changing the order of
#' modalities inside each column name
#' 
#' @param dataset a dataset containing a dataframe
#' @param expected_order a list of strings corresponding to the modalities names
#' @return the same dataset with different names for the same data
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
reorder_column_modalities <- function(dataset, expected_order) {
    
    reorder_values <- function(n) {

        k2v <- list_variabes_modalities(n)
        
        ppp <- function(k) {
            paste(k,"=",k2v[k],sep="")
        }

        r <- sapply(expected_order, ppp, simplify=T)
        paste(r, collapse=",")
    }

    colnames(dataset$data) <- sapply(colnames(dataset$data), reorder_values, simplify=T)
    dataset$attributes <- expected_order

    dataset

}

#' Prepares a case for Direct Probabilistic Peering
#'
#' Prepares a case based on two samples and two sets of probabilities.
#' It takes two samples, probabilities distributions for degrees 
#' (count of links per entity) and matching probabilities. 
#' It ensures all the parameters are consistent. In case \code{fix=TRUE}, a prior temptative 
#' to fix minor inconsistencies will be done using \code{\link{matching.fix}}.
#'
#' @param sample.A the sample to use for population A created using \code{\link{create_sample}}
#' @param sample.B the sample to use for population B created using \code{\link{create_sample}}
#' @param pdi the distribution of degrees for population A created with \code{\link{create_degree_probabilities_table}}
#' @param pdj the distribution of degrees for population B created with \code{\link{create_degree_probabilities_table}}
#' @param pij the matching probabilities created with \code{\link{create_matching_probabilities_table}}
#' @param fix if \code{TRUE}, the function will fix mismatches (notably missing variables in pdi, pdj or pij) by itself
#' @return a case ready to be solved with \code{\link{matching.solve}}
#'
#' @seealso \code{\link{matching.solve}} to use the result of this function
#' @examples
#' data(cas1)
#' prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#'
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.prepare <- function(sample.A, sample.B, pdi, pdj, pij, fix=TRUE) {

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
    fixed <- matching.fix(sample.A, sample.B, pdi, pdj, pij)
    sample.A_fixed <- fixed$sample.A
    sample.B_fixed <- fixed$sample.B
    pdi_fixed <- fixed$pdi
    pdj_fixed <- fixed$pdj
    pij_fixed <- fixed$pij

    # TODO ensure that every dimension value in dictionaries has a counterpart in our parameters!

    # ... the list of attributes tackled on sample.A, pdi_fixed and colname(pij_fixed) should be consistent 
    if (all.equal(pdi_fixed$attributes, pij_fixed$Ai) != TRUE) {
        stop("the attributes and values defined in pdi_fixed (",paste(pdi_fixed$attributes,collapse=","),") ",
            "should be the same as the columns of pij_fixed (", paste(pij_fixed$Ai,collapse=","), ").")
    }
    if (!all(pij_fixed$Ai %in% colnames(sample.A$sample))) {
         stop("the Ai attributes defined pdi_fixed and pij_fixed ",
            "(",paste(pij_fixed$Ai,collapse=","),")",
            " should be present in the sample of A (",
            paste(colnames(sample.A$sample),collapse=","),").")
    }
    if (!all.equal(pdj_fixed$attributes, pij_fixed$Bi)) {
        stop("the attributes and values defined in pdj_fixed (",paste(pdj_fixed$attributes,collapse=","),") ",
            "should be the same as the rows of pij_fixed (", paste(pij_fixed$Bi,collapse=","), ").")
    }
    if (!all(pij_fixed$Bi %in% colnames(sample.B_fixed$sample))) {
        stop("the Bj attributes defined pdj_fixed and pij_fixed ",
            "(",paste(pij_fixed$Bi,collapse=","),")",
            " should be present in the sample of B (",
            paste(colnames(sample.B_fixed$sample),collapse=","),").")
    }

    # ... for sure the list of indices should be the very same 
    if (!identical(colnames(pdi_fixed$data),colnames(pij_fixed$data))) {
        stop(paste("the variables and modalities for pdi_fixed (",paste(colnames(pdi_fixed$data),collapse=" "),") and pij_fixed (",paste(colnames(pij_fixed$data),collapse=" "),") should be exactly the same",sep=""))
    }
    if (!identical(colnames(pdj_fixed$data),rownames(pij_fixed$data))) {
        stop(paste("the variables and modalities for pdj_fixed (",paste(colnames(pdj_fixed$data),collapse=" "),") and pij_fixed (",paste(rownames(pij_fixed$data),collapse=" "),") should be exactly the same",sep=""))
    }

    # ... the list of values for Ai should be covered everywhere. 
    list_variabes_modalities(colnames(pdi_fixed$data))

    inputs <- list(
                pdi_fixed=pdi_fixed, 
                pdj_fixed=pdj_fixed,
                pij_fixed=pij_fixed,
                sample.A=list(dictionary=sample.A_fixed$dictionary),
                sample.B=list(dictionary=sample.B_fixed$dictionary)
                )

    # store as inputs
    
    inputs$pij_fixed <- pij_fixed
    
    # analyze sample frequencies
    # ... compute frequencies fi    
    inputs$ci <- measure_contigencies(sample.A_fixed, colnames(pdi_fixed$data))
    inputs$cj <- measure_contigencies(sample.B_fixed, colnames(pdj_fixed$data))
    
    # analyze degrees
    # ... compute average degrees
    inputs$di <- compute_average_degree(pdi_fixed)
    inputs$dj <- compute_average_degree(pdj_fixed)

    inputs$max.di <- assess.max.pd(inputs$pdi_fixed)
    inputs$max.dj <- assess.max.pd(inputs$pdj_fixed)
    inputs$min.di <- assess.min.pd(inputs$pdi_fixed)
    inputs$min.dj <- assess.min.pd(inputs$pdj_fixed)

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

