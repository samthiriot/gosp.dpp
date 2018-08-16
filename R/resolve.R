

#' Rounds an object by ensuring its totals remains equal before and after rounding
#'
#' Freely inspired by \url{https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum}
#' 
#' @param x the object to round
#' @param ... ignored
#' @return the rounded version of the object
#'
#' @export
round_sum <- function (x, ...) {
   UseMethod("round_sum", x)
 }

round_sum.numeric <- function(x) {
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

round_sum.matrix <- function(m) {

    # convert the matrix to a vector (by columns)
    x <- unlist(m)

    # round it 
    r <- round_sum.numeric(x)

    # convert it to a matrix
    mat <- matrix(r, ncol=ncol(m), nrow=nrow(m))

    # convert it back to a matrix
    res <- as.data.frame(mat, row.names=row.names(m), col.names=colnames(m))
    colnames(res) <- colnames(m)
    
    res
}
round_sum.data.frame <- round_sum.matrix

# for a distribution of degree, 
# sums the total degree, on the basis of 
# 0  
sum_degrees <- function(pdx) {

    # ... sum the weight we are playing with 
    power <- rep(0, ncol(pdx))
    for (i in 1:nrow(pdx)) {
        power <- power + (i-1)*pdx[i,]
    }

    power
}

#' Update a given distribution column of degrees probabilities 
#' so its matching a target average degree for this column.
#' 
#' It works by finding a pivot and then increasing values on one 
#' side and removing probabilities on the other side.
#' 
#' @param pdx the distribution of degrees 
#' @param t.target the vector of the expected average degrees
update_degree_distribution.col <- function(pdx, t.target) {

    np.orig <- pdx * 0:(length(pdx)-1)
    t.orig <- sum(np.orig)
    # cat("t.orig", t.orig, "\n")
    
    # quick exit if there is nothing to do
    if (abs(t.target - t.orig) <= 1e-10) {
        return(pdx)
    }

    p.potential <- rep(0,times=length(pdx))
    p.potential.sum.neg <- 0

    # what are the free margins?
    for (i in 1:length(pdx)) {
        if (pdx[i] == 0) {
            # nothing
        } else if (
                    ( (t.orig < t.target ) & (i-1 < t.target) )
                    |
                    ( (t.orig > t.target ) & (i-1 > t.target) )
                    ) {
            p.potential[i] <- -pdx[i]
            p.potential.sum.neg <- p.potential.sum.neg - pdx[i]
        } else {
            p.potential[i] <- 1 - pdx[i]
        }
    }

    # 
    # print("potential")
    # print(p.potential)

    # what can we gain ?
    np.potential <- pmin(p.potential,-p.potential.sum.neg) * 0:(length(pdx)-1)
    #print(np.potential)

    np.min <- min(np.potential)
    np.max <- max(np.potential)
    np.potential.positive.max <- max(np.potential) # max(p.potential[which(np.potential == np.max)])
    # cat("np.potential.positive.max",np.potential.positive.max,"\n") 
    
    np.potential.negative.sum <- sum(np.potential[which(np.potential < 0)])
    # cat("np.potential.negative.sum",np.potential.negative.sum,"\n") 
    
    np.potential.cumulated <- np.potential.negative.sum + np.potential.positive.max

    # create factors, 0 or 1
    factors <- pmin( (np.potential == np.max) + (p.potential < 0), 1)
    # print("factors")
    # print(factors)

    rat <- (t.target - t.orig) / np.potential.cumulated
    # cat("rat",rat,"\n")


    p.add <- factors * rat * pmin(p.potential, -p.potential.sum.neg)

    # cat("p.add",p.add,"=",sum(p.add),"\n")
    if (abs(sum(p.add)) >= 1e-10 ) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }

    # print("p.modified ------------------------------")
    p.modified <- p.add + pdx
    # cat("p.modified",p.modified,"=",sum(p.modified),"\n")
    if (length(which(p.modified < 0)) > 0) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }
    if (length(which(p.modified > 1)) > 0) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }

    
    if (abs(sum(p.modified)-1) > 1e-10) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }
     
    # now ensure this leads to the expected result 
    np.res.cumulated <- sum(p.modified * 0:(length(pdx)-1))
    # cat("np.res.cumulated",np.res.cumulated,"for target", t.target,"\n")

    p.modified
}

#' Update a given distribution of degrees probabilities so its matching target values.
#' 
#' This process is done column by column by calling \link{update_degree_distribution.col}.
#' '
#' @param pdx the distribution of degrees 
#' @param dx the vector of the expected average degrees
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
update_degree_distribution <- function(pdx, dx) {

    # ... sum the weight we are playing with 
    power <- sum_degrees(pdx)

    # print("working on:")
    # print(pdx)
    # print("original average degree: ")
    # print(power)
    # print("yet we are supposed to reach")
    # print(dx)
    
    factord <- dx/power

    pdx.reweighted <- pdx

    for (c in 1:ncol(pdx)) {
        pdx.reweighted[,c] <- update_degree_distribution.col(pdx[,c],dx[c])
    }

    colnames(pdx.reweighted) <- colnames(pdx)

    pdx.reweighted

}


#' Patches a degree distribution contingencies so its totals fit expectations
#' 
#' Adapts a distribution of degrees so that each column times n sums to expected nn
#' (basically solves rounding issues for this specific case)
#'
#' @param pdn the contigencies of degrees
#' @param nn the expected total degree (ni) that should be reached
#' @param cn the expected total count (ci) that should be reached
#' @param verbose print detailed information during the process
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
rectify.degree.counts <- function(pdn, nn, cn, verbose=FALSE) {

    if (verbose) {
        cat("should rectify\n") 
        print(pdn)
        cat("to reach\n")
        print(nn)
    }
    for (col in 1:ncol(pdn)) {

        total.current <- sum(pdn[,col]*0:(nrow(pdn)-1))
        total.expected <- nn[col]

        if (total.current > total.expected) {
            to.remove <- total.current - total.expected

            # no more warning: we anyway check at the end how well we fixed it
            # if (to.remove > 1) {
            #     warning("/!\\ the solution to fix rounding errors below 1 slot (here:",to.remove,") is not managed well\n.",sep="")
            # }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if ( (pdn[2,col] >= to.remove) && (pdn[1,col] > 0)) {
                pdn[2,col] <- pdn[2,col] - to.remove
                pdn[1,col] <- pdn[1,col] + to.remove
            } else {
                if (verbose)
                    warning("/!\\ found no good solution to fix this rounding error by removing ",to.remove,"\n") 
            }

        } else if (total.current < total.expected) {
            to.add <- total.expected - total.current 

            # no more warning: we anyway check at the end how well we fixed it
            # if (to.add > 1) {
            #     warning("/!\\ the solution to fix rounding errors above 1 slot (here:",to.add,") is not managed well\n", sep="")
            # }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if ( (pdn[1,col] >= to.add) && (pdn[2,col] > 0)) {
                pdn[1,col] <- pdn[1,col] - to.add
                pdn[2,col] <- pdn[2,col] + to.add
            } else {
                if (verbose)
                    warning("/!\\ found no good solution to fix this rounding error of an additional ",to.add,"\n") 
            }
        }

    }

    # adapt the cn s we did not checked until now
    to.add <- unname(cn) - unname(colSums(pdn))
    if (!all(to.add == 0) && (sum(to.add) == 0)) {
        
        # we can solve the problem by playing on the first raw (the one which induces no degree)
        pdn[1,] <- pdn[1,] + to.add

        # nota: this solution might lead to inconsistent result; 
        # typically if degree 0 is impossible, we are adding negative values here ! 
        # it will be checked later. 
    } 

    # ensure we did it well, that is we achieve to reach the expected result
    if (all.equal(
                unname(colSums(pdn*(0:(nrow(pdn)-1)))), 
                unname(nn)
                ) != TRUE) {

        stop("was unable to fix the rounding of degree distributions ",
            "so the total slots ",
            paste(unname(colSums(pdn*(0:(nrow(pdn)-1)))), collapse=","),
            " sums up to ",
            paste(nn, collapse=",")
            )
    }

    if ((all.equal(unname(cn), unname(colSums(pdn))) != TRUE) || any(pdn[1,] < 0)) {
        stop("was unable to fix the rounding of degree distributions ",
            "so the total count ", paste(unname(colSums(pdn)), collapse=","),
            " matches ",
            paste(unname(cn), collapse=",")
            )
    } 
    

    pdn
}

#' Normalise an object
#' 
#' takes a vector, matrix or list and updates its content
#' by dividing it by its sum.
#'
#' @param df any data to normalise
#' @return the same object normalised
#'
normalise <- function(df) {
    df / sum(df)
}


#' Replaces NaN by zeros
#' 
#' Replaces NaN by zeros. Used when we might divide 0 
#' by 0 and know the result should be 0.
#' 
#' @param x a vector, list or matrix
#' @param ... ignored
#' @return the same type without NaNs
#'
nan_to_zeros <- function (x, ...) {
   UseMethod("nan_to_zeros", x)
 }

nan_to_zeros.list <- function (x, ...) {
    x[which(is.nan(x))] <- 0
    x
 }
nan_to_zeros.numeric <- nan_to_zeros.list

nan_to_zeros.matrix <- function (x, ...) {
    for (i in 1:ncol(x)) {
        x[,i] <- nan_to_zeros(x[,i])
    }
    x
}
nan_to_zeros.data.frame <- nan_to_zeros.matrix


#' Replaces Infinite by zeros
#' 
#' Replaces Infinite by zeros. Used when we divide 
#' by a value and know the result should be 0 
#' when the diviser is 0.
#' 
#' @param vv a vector or list
#' @return the same vector without Infinite
#'
inf_to_zeros <- function(vv) {
    vv[which(is.infinite(vv))] <- 0
    vv
}


#' Completes a partial solution by inference  
#' 
#' Takes an incomplete solution and enriches by applying the equations.
#' For instance if an equation links \code{hat.di * hat.ci = hat.ni}, 
#' then this equation will be applied if two of these variables are available.   
#'
#'
#' @param sol the current solution (a named list)
#' @param case the case to solve
#' @param verbose if TRUE, will display detailed information on the console
#' @param indent the level of indentation for verbose messages
#' @return a solution with possibly more variables known
#'
#' @seealso this function relies on \code{\link{round_sum}} to round values 
#'          so their sums is preserved; 
#'          \code{\link{normalise}} to normalise variables;
#'          \code{\link{update_degree_distribution}} to compute the detail 
#'          of degree distributions from average degrees.
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
propagate.direct <- function(sol,case, verbose=FALSE, indent=1) {

    info.rule <- function(name, sol, verbose=FALSE, indent=3) {
        if (verbose) { 
            cat(paste(rep("\t",times=indent),collapse=""),name,"\n",sep="")
        }
        sol$inference <- list(sol$inference, list(name))
        sol
    } 

    if (verbose) {
        cat(paste(rep("\t",times=indent),collapse=""),"solving equations based on known variables: ",paste.known(sol),"\n",sep="")
    }

    if (is.null(sol$inference)) {
        sol$inference <- list()
    }

    while (TRUE) {

        changed <- FALSE


        # forward nA + fi -> ci 
        if ( 
            ("hat.nA" %in% names(sol)) && ("hat.fi" %in% names(sol)) && (!"hat.ci" %in% names(sol))
            ) {
            sol$hat.ci <- round_sum(sol$hat.nA * sol$hat.fi)
            sol$hat.fi <- sol$hat.ci / sol$hat.nA
            sol <- info.rule("hat.nA, hat.fi -> hat.ci", sol, verbose=verbose, indent=indent+1)
            #print("sol$hat.ci <- round(sol$hat.nA * sol$hat.fi)")
            changed <- TRUE
        }
        
        # forward nB + fj -> cj 
        if ( 
            ("hat.nB" %in% names(sol)) && ("hat.fj" %in% names(sol)) && (!"hat.cj" %in% names(sol))
            ) {
            sol$hat.cj <- round_sum(sol$hat.nB * sol$hat.fj)
            # rounding might introduce a bias; update this proba
            sol$hat.fj <- sol$hat.cj / sol$hat.nB
            sol <- info.rule("hat.nB, hat.fj -> hat.cj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }


        # forward ci + di -> ni 
        if ( 
            ("hat.ci" %in% names(sol)) && ("hat.di" %in% names(sol)) && (!"hat.ni" %in% names(sol))
            ) {
            sol$hat.ni <- round_sum(sol$hat.ci * sol$hat.di)
            # adapt rounding
            sol$hat.di <- nan_to_zeros(sol$hat.ni / sol$hat.ci)
            sol <- info.rule("hat.ci, hat.di -> hat.ni", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

         # forward cj + dj -> nj 
        if ( 
            ("hat.cj" %in% names(sol)) && ("hat.dj" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- round_sum(sol$hat.cj * sol$hat.dj)
            sol$hat.dj <- nan_to_zeros(sol$hat.nj / sol$hat.cj)
            sol <- info.rule("hat.cj, hat.dj -> hat.nj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        if ( 
            ("hat.ni" %in% names(sol)) && (!"hat.nL" %in% names(sol))
            ) {
            sol$hat.nL <- sum(sol$hat.ni)
            sol <- info.rule("hat.ni -> hat.nL", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        if ( 
            ("hat.nj" %in% names(sol)) && (!"hat.nL" %in% names(sol))
            ) {
            sol$hat.nL <- sum(sol$hat.nj)
            sol <- info.rule("hat.nj -> hat.nL", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # forward sum(ci) -> nA
        if ( 
            ("hat.ci" %in% names(sol)) && (!"hat.nA" %in% names(sol))
            ) {
            sol$hat.nA <- sum(sol$hat.ci)
            sol <- info.rule("hat.ci -> hat.nA", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.cj" %in% names(sol)) && (!"hat.nB" %in% names(sol))
            ) {
            sol$hat.nB <- sum(sol$hat.cj)
            sol <- info.rule("hat.cj -> hat.nB", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # based on ni and pij (either defined, or from inputs)
        if ( 
            ("hat.ni" %in% names(sol))
            && (!"hat.nij" %in% names(sol))
            ) {
            pij <- if ("hat.pij" %in% names(sol)) sol$hat.pij else case$inputs$pij$data
            sol$hat.nij <- nan_to_zeros( t(sol$hat.ni * t(pij) / colSums(pij)) )
            # round per column (to keep the sums consistent)
            for (i in 1:ncol(sol$hat.nij)) {
                sol$hat.nij[,i] <- round_sum(sol$hat.nij[,i])
            }
            # adapt th based on rounded values
            sol$hat.pij <- sol$hat.nij / sum(sol$hat.nij)
            sol <- info.rule("hat.ni -> hat.nij, hat.pij", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.nj" %in% names(sol))
            && (!"hat.nij" %in% names(sol))
            ) {
            pij <- if ("hat.pij" %in% names(sol)) sol$hat.pij else case$inputs$pij$data
            sol$hat.nij <- nan_to_zeros( sol$hat.nj * pij / rowSums(pij) )
            # round per line (to keep the sums consistent)
            for (i in 1:nrow(sol$hat.nij)) {
                sol$hat.nij[i,] <- round_sum(sol$hat.nij[i,])
            }
            sol$hat.pij <- sol$hat.nij / sum(sol$hat.nij)
            sol <- info.rule("hat.nj -> hat.nij, hat.pij", sol, verbose=verbose, indent=indent+1)

            changed <- TRUE
        }


        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.ni" %in% names(sol))
            ) {
            sol$hat.ni <- colSums(sol$hat.nij)
            sol <- info.rule("hat.nij -> hat.ni", sol, verbose=verbose, indent=indent+1)

            changed <- TRUE
        }
        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- rowSums(sol$hat.nij)
            sol <- info.rule("hat.nij -> hat.nj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
 
        # forward ni + di -> ci 
        if (
            ("hat.ni" %in% names(sol)) && ("hat.di" %in% names(sol)) && (!"hat.ci" %in% names(sol))
            ) {
            candidate_ci <- sol$hat.ni / sol$hat.di
            # if di = 0, then we cannot extrapolate ci based on it...
            # ... yet we might just assume it might be any other value
            # let's say we try to keep its relative frequency.             
            fi_for_missing <- if ("hat.fi" %in% names(sol)) sol$hat.fi else case$stats$fi
            indices_not_inf <- which((!is.infinite(candidate_ci)) && (!is.nan(candidate_ci)))
            average_ratio <- mean(candidate_ci[indices_not_inf] / fi_for_missing[indices_not_inf])
            indices_wrong <- which(is.infinite(candidate_ci) | is.nan(candidate_ci))
            candidate_ci[indices_wrong] <- fi_for_missing[indices_wrong] * average_ratio
            sol$hat.ci <- round_sum(candidate_ci)
            sol$hat.di[indices_not_inf] <- sol$hat.ni[indices_not_inf] / sol$hat.ci[indices_not_inf]
            sol <- info.rule("hat.ni, hat.di -> hat.ci", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE

        }
        if (
            ("hat.nj" %in% names(sol)) && ("hat.dj" %in% names(sol)) && (!"hat.cj" %in% names(sol))
            ) {
            candidate_cj <- sol$hat.nj / sol$hat.dj
            # if di = 0, then we cannot extrapolate ci based on it...
            # ... yet we might just assume it might be any other value
            # let's say we try to keep its relative frequency.             
            fj_for_missing <- if ("hat.fj" %in% names(sol)) sol$hat.fj else case$stats$fj
            indices_not_inf <- which((!is.infinite(candidate_cj)) && (!is.nan(candidate_cj)))
            average_ratio <- mean(candidate_cj[indices_not_inf] / fj_for_missing[indices_not_inf])
            indices_wrong <- which(is.infinite(candidate_cj) | is.nan(candidate_cj))
            candidate_cj[indices_wrong] <- fj_for_missing[indices_wrong] * average_ratio
            sol$hat.cj <- round_sum(candidate_cj)
            sol$hat.dj[indices_not_inf] <- sol$hat.nj[indices_not_inf] / sol$hat.cj[indices_not_inf]
            sol <- info.rule("hat.nj, hat.dj -> hat.cj", sol, verbose=verbose, indent=indent+1)            
            changed <- TRUE

        }
        # if we have card and nx, then we can assess degree
        if ( 
            ("hat.ni" %in% names(sol)) && ("hat.ci" %in% names(sol)) &&  (!"hat.di" %in% names(sol))
            ) {
            sol$hat.di <- pmax(
                            case$inputs$min.di, 
                            pmin(
                                    case$inputs$max.di,
                                    nan_to_zeros(sol$hat.ni / sol$hat.ci)
                                    )
                            )
            sol <- info.rule("hat.ni, hat.ci -> hat.di", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.nj" %in% names(sol)) && ("hat.cj" %in% names(sol)) &&  (!"hat.dj" %in% names(sol))
            ) {
            sol$hat.dj <- pmax(
                            case$inputs$min.dj, 
                            pmin(
                                    case$inputs$max.dj,
                                    nan_to_zeros(sol$hat.nj / sol$hat.cj)
                                    )
                            )
            sol <- info.rule("hat.nj, hat.cj -> hat.dj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # pij -> pi 
        if ( 
            ("hat.pij" %in% names(sol)) && (!"hat.pi" %in% names(sol))
            ) {
            sol$hat.pi <- colSums(sol$hat.pij)
            sol <- info.rule("hat.pij -> hat.pi", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.pij" %in% names(sol)) && (!"hat.pj" %in% names(sol))
            ) {
            sol$hat.pj <- rowSums(sol$hat.pij)
            sol <- info.rule("hat.pij -> hat.pj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # fi * pi -> di  
        if ( 
            ("hat.pi" %in% names(sol)) && ("hat.fi" %in% names(sol)) && (!"hat.di" %in% names(sol))
            ) {
            sol$hat.di <- nan_to_zeros(sol$hat.pi / sol$hat.fi)
            # fix potential NaNs due to di=0 by their native di counterparts (they do not matter anyway)
            sol <- info.rule("hat.pi, hat.fi -> hat.di", sol, verbose=verbose, indent=indent+1)
            
            # in factor we might change it of a factor, only the ratio matters
            factor <- 1
            if (any(sol$hat.di > case$inputs$max.di)) {
                # we are higher than what is allowed
                factor <- min(case$inputs$max.di / sol$hat.di, na.rm=TRUE)
            }  else if (any(sol$hat.di < case$inputs$min.di)) {
                factor <- min(case$inputs$min.di / sol$hat.di, na.rm=TRUE) # TODO verifier
            }
            sol$hat.di <- sol$hat.di * factor
            
            changed <- TRUE
        }

        if ( 
            ("hat.pj" %in% names(sol)) && ("hat.fj" %in% names(sol)) &&  (!"hat.dj" %in% names(sol))
            ) {
            sol$hat.dj <- nan_to_zeros(sol$hat.pj / sol$hat.fj)
            
            sol <- info.rule("hat.pj, hat.fj -> hat.dj", sol, verbose=verbose, indent=indent+1)

            # in factor we might change it of a factor, only the ratio matters
            factor <- 1
            if (any(sol$hat.dj > case$inputs$max.dj)) {
                # we are higher than what is allowed
                factor <- min(case$inputs$max.dj / sol$hat.dj, na.rm=TRUE)
            }  else if (any(sol$hat.di < case$inputs$min.di)) {
                factor <- min(case$inputs$min.dj / sol$hat.dj, na.rm=TRUE) # TODO verifier
            }
            sol$hat.dj <- sol$hat.dj * factor
            changed <- TRUE
        }

        # fj + dj -> pj
        if ( 
            ("hat.fi" %in% names(sol)) && ("hat.di" %in% names(sol)) &&  (!"hat.pi" %in% names(sol))
            ) {
            sol$hat.pi <- normalise(sol$hat.fi * sol$hat.di)
            sol <- info.rule("hat.fi, hat.di -> hat.pi", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.fj" %in% names(sol)) && ("hat.dj" %in% names(sol)) &&  (!"hat.pj" %in% names(sol))
            ) {
            sol$hat.pj <- normalise(sol$hat.fj * sol$hat.dj)
            sol <- info.rule("hat.fj, hat.dj -> hat.pj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # pj + pij -> pij 
        if ( 
            ("hat.pi" %in% names(sol)) && (!"hat.pij" %in% names(sol))
            ) {
            sol$hat.pij <- t(t(case$inputs$pij$data) / colSums(case$inputs$pij$data)) * sol$hat.pi
            sol <- info.rule("hat.pi -> hat.pij", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.pj" %in% names(sol)) && (!"hat.pij" %in% names(sol))
            ) {
            sol$hat.pij <- case$inputs$pij$data / rowSums(case$inputs$pij$data) * sol$hat.pj
            sol <- info.rule("hat.pj -> hat.pij", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # ci -> fi
        if ( 
            ("hat.ci" %in% names(sol)) && (!"hat.fi" %in% names(sol))
            ) {
            sol$hat.fi <- sol$hat.ci / sum(sol$hat.ci)
            sol <- info.rule("hat.ci -> hat.fi", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.cj" %in% names(sol)) && (!"hat.fj" %in% names(sol))
            ) {
            sol$hat.fj <- sol$hat.cj / sum(sol$hat.cj)
            sol <- info.rule("hat.cj -> hat.fj", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # pdi -> di
        if ( 
            ("hat.pdi" %in% names(sol)) && (!"hat.di" %in% names(sol))
            ) {
            sol <- info.rule("hat.pdi -> hat.di", sol, verbose=verbose, indent=indent+1)
            sol$hat.di <- colSums(sol$hat.pdi*(0:(nrow(sol$hat.pdi)-1)))
            changed <- TRUE
        }
        if ( 
            ("hat.pdj" %in% names(sol)) && (!"hat.dj" %in% names(sol))
            ) {
            sol <- info.rule("hat.pdj -> hat.dj", sol, verbose=verbose, indent=indent+1)
            sol$hat.dj <- colSums(sol$hat.pdj*(0:(nrow(sol$hat.pdj)-1)))
            changed <- TRUE
        }

        # ndi -> ci
        if ( 
            ("hat.ndi" %in% names(sol)) && (!"hat.ci" %in% names(sol))
            ) {
            sol <- info.rule("hat.ndi -> hat.ci", sol, verbose=verbose, indent=indent+1)
            sol$hat.ci <- colSums(sol$hat.ndi)
            changed <- TRUE
        }
        if ( 
            ("hat.ndj" %in% names(sol)) && (!"hat.cj" %in% names(sol))
            ) {
            sol <- info.rule("hat.ndj -> hat.cj", sol, verbose=verbose, indent=indent+1)
            sol$hat.cj <- colSums(sol$hat.ndj)
            changed <- TRUE
        }


        # pdi, ci -> ndi
        if ( 
            all(c("hat.pdi", "hat.ci", "hat.ni") %in% names(sol)) && (!"hat.ndi" %in% names(sol))
            ) {

            sol <- info.rule("hat.pdi, hat.ci -> hat.ndi", sol, verbose=verbose, indent=indent+1)
            # print(sol$hat.pdi)
            sol$hat.ndi <- round_sum(t(t(sol$hat.pdi) * sol$hat.ci))
            sol$hat.ndi <- rectify.degree.counts(sol$hat.ndi, sol$hat.ni, sol$hat.ci, verbose=FALSE)
            
            # update pdi after rounding (when divisible, else we keep the theorical version)
            indices <- which(sol$hat.ci!=0)
            sol$hat.pdi[indices] <- t(t(sol$hat.ndi[indices]) / sol$hat.ci[indices])

            # update our past knowledge according to rounding
            if ("hat.di" %in% names(sol)) {
                sol$hat.di <- colSums(sol$hat.pdi * (0:(nrow(sol$hat.pdi)-1)))
            }

            changed <- TRUE
        }
        if ( 
            all(c("hat.pdj", "hat.cj", "hat.nj") %in% names(sol)) && (!"hat.ndj" %in% names(sol))
            ) {
            
            sol <- info.rule("hat.pdj, hat.cj -> hat.ndj", sol, verbose=verbose, indent=indent+1)
            sol$hat.ndj <- round_sum(t(t(sol$hat.pdj) * sol$hat.cj))
            sol$hat.ndj <- rectify.degree.counts(sol$hat.ndj, sol$hat.nj, sol$hat.cj, verbose=FALSE)   

            # update pdi after rounding (when divisible, else we keep the theorical version)
            indices <- which(sol$hat.cj!=0)
            sol$hat.pdj[indices] <- as.data.frame(
                                        t(t(sol$hat.ndj[indices]) / sol$hat.cj[indices]),
                                        optional=FALSE)

            # update our past knowledge according to rounding
            if ("hat.dj" %in% names(sol)) {
                sol$hat.dj <- colSums(sol$hat.pdj * (0:(nrow(sol$hat.pdj)-1)))
            }

            changed <- TRUE
        }
    
    
        # di -> pdi
        if ( 
            ("hat.di" %in% names(sol)) && (!"hat.pdi" %in% names(sol))
            ) {
            sol <- info.rule("hat.di -> hat.pdi", sol, verbose=verbose, indent=indent+1)
            sol$hat.pdi <- tryCatch({
                        update_degree_distribution(case$inputs$pdi$data, sol$hat.di)
                    }, error = function(e) {
                        if (verbose) 
                            cat(rep("\t",times=indent+2),"unable to update hat.pdi to fit hat.di=",paste(sol$hat.di,collapse=","),"\n",sep="")                   
                        
                        stop("The case is too constrained to be solved (hat.di are not compliant with the original pdi)")
                    })           
            changed <- TRUE
        }

        if ( 
            ("hat.dj" %in% names(sol)) && (!"hat.pdj" %in% names(sol))
            ) {
            sol <- info.rule("hat.dj -> hat.pdj", sol, verbose=verbose, indent=indent+1)
            sol$hat.pdj <- tryCatch({
                        update_degree_distribution(case$inputs$pdj$data, sol$hat.dj)
                    }, error = function(e) {
                        if (verbose) 
                            cat(rep("\t",times=indent+2),"unable to update hat.pdj to fit hat.dj=",paste(sol$hat.dj,collapse=","),"\n",sep="")                   
                        stop("The case is too constrained to be solved (probabilities hat.dj are not compliant with the original pdj)")
                    })
            changed <- TRUE
        }

        # stop looping when we changed nothing 
        # (so we are sure we cannot apply any novel equation)
        if (!changed) {
            break
        }

    }

    sol

}

#' Asserts two values are equal.
#' 
#' In the context of checking the consistency of a solution,
#' ensures the content are similar.
#' 
#' @param v1 a scalar, vector, matrix or data frame 
#' @param v2 a scalar, vector, matrix or data frame 
#' @param msg the message to display in case of failure
#' @param indent the indentation (count of tabs) for verbose display
#' @param verbose if TRUE, assertion errors will be displayed in the console
#' @param tolerance the numeric tolerance to define two numeric are equal
#' @return 1 in case of error else 0
#' 
assert.equal <- function(v1,v2,msg, verbose=FALSE, indent=3, tolerance=1.5e-8) {

    if (!identical(v1,v2) && !(all.equal(v1,v2,tolerance=tolerance)==TRUE)) {  # && ((length(v1) > 1) & sum((v1-v2)^2) > 0)

        if (verbose) {
            
            if (is.numeric(v1) && is.numeric(v2)) {
                cat(rep("\t",times=indent),
                        "assert error on '",msg,"': ",
                        paste(v1,collapse=",")," != ",paste(v2,collapse=","),
                        "\n",sep="")            
            } else {
                cat(rep("\t",times=indent),
                        "assert error on '",msg,"' ",
                        "\n",sep="")
                cat(rep("\t",times=indent),all.equal(v1,v2),"\n",sep="")
            }
            # print(v1)
            # print(v2)
            
            #print( all.equal(v1,v2) )
        }
        
        return(1)
    }
    0
}

#' Checks the consistency of a solution
#' 
#' Takes a temptative solution and the initial case to resolve. 
#' Tests wether the equations all lead to consistent results of 
#' do contradict each other. Returns an error with \code{stop()}
#' in case of error if \code{fail=TRUE}; else just returns the 
#' count of problems. 
#'
#' @param sol the current solution (a named list)
#' @param case the case to be solved
#' @param fail if TRUE, an error will call stop()
#' @param verbose if TRUE, will ddetect.problemsisplay detailed information on the console
#' @param indent the indentation (count of tabs) for verbose display
#' @return the count of problems 
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
detect.problems <- function(sol, case, fail=TRUE, verbose=FALSE, indent=1) {

    if (verbose) {
        cat(rep("\t",times=indent),"checking the consistency of the current solution with ",paste.known(sol),"\n",sep="")
    }    
    problems <- 0

    # nA = sum(ni)
    if (!is.null(sol$hat.nA) && !is.null(sol$hat.ci)) {
        problems <- problems + assert.equal(
                                    sol$hat.nA, sum(sol$hat.ci), 
                                    "hat.nA = sum(hat.ci)", 
                                    verbose=verbose, 
                                    indent=indent+1)
    }
    if (!is.null(sol$hat.nB) && !is.null(sol$hat.cj)) {
        problems <- problems + assert.equal(
                                    sol$hat.nB, sum(sol$hat.cj), 
                                    "hat.nB = sum(hat.cj)", 
                                    verbose=verbose, 
                                    indent=indent+1)
    }

    # ni / nA = fi 
    if (!is.null(sol$hat.nA) && !is.null(sol$hat.ci) && !is.null(sol$hat.fi)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.ci/sol$hat.nA), 
                                        as.vector(sol$hat.fi), 
                                        "hat.ci/hat.nA = hat.fi", 
                                        verbose=verbose, indent=indent+1, tolerance=1)
    }
    if (!is.null(sol$hat.nB) && !is.null(sol$hat.cj) && !is.null(sol$hat.fj)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.cj/sol$hat.nB), 
                                        as.vector(sol$hat.fj), 
                                        "hat.cj/hat.nB = hat.fj", 
                                        verbose=verbose, indent=indent+1, tolerance=1)
    }  

    # sum(fi) = 1
    if (!is.null(sol$hat.fi)) {
        problems <- problems + assert.equal(
                                        1, sum(sol$hat.fi), 
                                        "sum(hat.fi)=1", 
                                        verbose=verbose, indent=indent+1)
    }
    if (!is.null(sol$hat.fj)) {
        problems <- problems + assert.equal(
                                        1, sum(sol$hat.fj), 
                                        "sum(hat.fj)=1", 
                                        verbose=verbose, indent=indent+1)
    }

    # ni = ci*di
    if (!is.null(sol$hat.ci) && !is.null(sol$hat.di) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                    as.vector(sol$hat.ni), 
                                    as.vector(sol$hat.ci*sol$hat.di), 
                                    "hat.ni = hat.ci*hat.di", 
                                    verbose=verbose, indent=indent+1,
                                    tolerance=1)
    }
    if (!is.null(sol$hat.cj) && !is.null(sol$hat.dj) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                    as.vector(sol$hat.nj), 
                                    as.vector(sol$hat.cj*sol$hat.dj), 
                                    "hat.nj = hat.cj*hat.dj", 
                                    verbose=verbose, indent=indent+1,
                                    tolerance=1)
    }

    # nL = sum(ni)
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                        sol$hat.nL, sum(sol$hat.ni), 
                                        "hat.nL = sum(hat.ni)", 
                                        verbose=verbose, indent=indent+1)
    }
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                        sol$hat.nL, sum(sol$hat.nj), 
                                        "hat.nL = sum(hat.nj)", 
                                        verbose=verbose, indent=indent+1)
    }

    # nL = sum(nij)
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.nij)) {
        problems <- problems + assert.equal(
                                        sol$hat.nL, sum(sol$hat.nij), 
                                        "hat.nL = sum(hat.nij)", 
                                        verbose=verbose, indent=indent+1)
    }

    # sum(hat.pij) = 1
    if (!is.null(sol$hat.pij)) {
        problems <- problems + assert.equal(
                                        1, sum(sol$hat.pij), 
                                        "sum(hat.pij) = 1", 
                                        verbose=verbose, indent=indent+1)
    }

    # hat.pij = hat.nij/sum(hat.nij)
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.pij)) {
        problems <- problems + assert.equal(
                                        sol$hat.nij/sum(sol$hat.nij), sol$hat.pij, 
                                        "hat.pij = hat.nij/sum(hat.nij)", 
                                        verbose=verbose, indent=indent+1)
    }
    
    # colSum(hat.nij) = hat.ni 
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                        as.vector(colSums(sol$hat.nij)), as.vector(sol$hat.ni), 
                                        "hat.ni = colSums(hat.nij)", 
                                        verbose=verbose, indent=indent+1)
    }
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                        as.vector(rowSums(sol$hat.nij)), as.vector(sol$hat.nj), 
                                        "hat.nj = rowSums(hat.nij)", 
                                        verbose=verbose, indent=indent+1)
    }

    # vsum(pdi) = 1
    if (!is.null(sol$hat.pdi)) {
        for (i in 1:ncol(sol$hat.pdi)) {
            problems <- problems + assert.equal(
                                        1, 
                                        sum(sol$hat.pdi[,i]), 
                                        paste("sum(hat.pdi[",i,"])=1",sep=""), 
                                        verbose=verbose, indent=indent+1)       
        }
    } 
    if (!is.null(sol$hat.pdj)) {
        for (i in 1:ncol(sol$hat.pdj)) {
            problems <- problems + assert.equal(
                                        1, 
                                        sum(sol$hat.pdj[,i]), 
                                        paste("sum(hat.pdj[",i,"])=1",sep=""), 
                                        verbose=verbose, indent=indent+1)       
        }
    } 

    # sum(hat.ndi) = hat.ci 
    if (!is.null(sol$hat.ndi) && !is.null(sol$hat.ci)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.ci), 
                                        as.vector(colSums(sol$hat.ndi)), 
                                        "hat.ci = sum(hat.ndi)", 
                                        verbose=verbose, indent=indent+1)  
    }    
    if (!is.null(sol$hat.ndj) && !is.null(sol$hat.cj)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.cj), 
                                        as.vector(colSums(sol$hat.ndj)), 
                                        "hat.cj = sum(hat.ndj)", 
                                        verbose=verbose, indent=indent+1)  
    }    

    # sum(hat.ndi * n) = hat.ni
    if (!is.null(sol$hat.ndi) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                as.vector(colSums(sol$hat.ndi * 0:(nrow(sol$hat.ndi)-1))), 
                                as.vector(sol$hat.ni), 
                                "hat.ni = sum( n * ndi[n])", 
                                verbose=verbose, indent=indent+1)  
    }
    if (!is.null(sol$hat.ndj) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                as.vector(colSums(sol$hat.ndj * 0:(nrow(sol$hat.ndj)-1))), 
                                as.vector(sol$hat.nj), 
                                "hat.nj = sum( n * ndj[n])", 
                                verbose=verbose, indent=indent+1)  
    }

    if (!is.null(sol$hat.di)) {
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.di < case$inputs$min.di), 
                                "hat.di >= min.di", 
                                verbose=verbose, indent=indent+1)     
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.di > case$inputs$max.di), 
                                "hat.di <= max.di", 
                                verbose=verbose, indent=indent+1)     
    }
    if (!is.null(sol$hat.dj)) {
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.dj < case$inputs$min.dj), 
                                "hat.dj >= min.dj", 
                                verbose=verbose, indent=indent+1)     
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.dj > case$inputs$max.dj), 
                                "hat.dj <= max.dj", 
                                verbose=verbose, indent=indent+1)     
    }

    if (verbose) {
        if (problems > 0) {
            cat(rep("\t",times=indent),"=> the solution is NOT consistent: ",problems," problems detected\n", sep="")
        } else {
            cat(rep("\t",times=indent),"=> the solution is consistent\n", sep="")
        }
    }
    if (fail && (problems > 0)) {
        stop("The case is too constrained to be solved (led to incoherent state)")
    }

    problems

}

#' Detects the missing chains of variables 
#' 
#' For a given temptative solution supposed to be consistent, 
#' checks the successive variables which depend on each other and 
#' are not yet solved. These chains of variables probably require 
#' using hypothesis for actual resolution.
#'
#' @param sol the current solution (a named list)
#' @param verbose if TRUE, will display detailed information on the console
#' @return a list of vectors (the chains) of strings 
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
detect.missing.chains <- function(sol, verbose=FALSE) {

    if (verbose) {
        cat("we know: ",paste.known(sol),"\n",sep="")
    }
    
    missing.chains <- list()
    current.chain <- NULL

    # "hat.ni", "hat.nj", 
    chain.n <- c("hat.nA", "hat.ci", "hat.fi", "hat.di", "hat.pij", "hat.dj", "hat.fj", "hat.cj", "hat.nB")

    for (i in 1:length(chain.n)) {
        var.name <- chain.n[[i]]
        
        #print(var.name)

        if (!(var.name %in% names(sol))) {
            #cat("missing",var.name,"\n")
            if (is.null(current.chain)) {
                # new chain 
                current.chain <- c(var.name)
            } else {
                # continue with the same chain 
                current.chain <- c(current.chain, var.name)
            }
        } else {
            # not missing
            if (!is.null(current.chain)) {
                missing.chains[[length(missing.chains)+1]] <- current.chain
                current.chain <- NULL
            }
        }
    }

    if (!is.null(current.chain)) {
        missing.chains[[length(missing.chains)+1]] <- current.chain
        current.chain <- NULL
    }

    if (verbose) {
        if (length(missing.chains) > 0) {
            cat("\tthere are ",length(missing.chains)," missing chains:\n")
            for (c in missing.chains) {
                cat("\t\t", paste(c,collapse=","))
            }
        } else {
            cat("\tthere is no missing chain in this solution.")
        }
    }
    missing.chains
}


#' Concatenates the known values
#' 
#' Only keeps in a solution (a list) the variables names which correspondond  to 
#' the chain of variables we need to solve
#'
#' @param sol the current solution (a list)
#' @return a string
#' 
paste.known <- function(sol) {
    paste(
        intersect(
            names(sol),
            c("hat.nA", "hat.ci", "hat.fi", "hat.di", "hat.pij", "hat.dj", "hat.fj", "hat.cj", "hat.nB")
            ),
        collapse=","
        )
}

#' Solves a chain of missing variables
#' 
#' For a given incomplete solution, a set of dependant variables to solve,
#'  a generation case, and targets defined by the user, tries to identify 
#' the best solution (if any).   
#'
#' It works by (i) trying to define hypothesis such as "let's try to respect 
#' the ideal value for this variable". It then (ii) relies on inference to 
#' determine the consequences of this hypothesis, and checks (iii) whether
#' this solution is consistent. 
#' 
#' At the end of this process, we might have no solution, one unique solution or 
#' several ones. If several solutions are available, the best one is taken, with best
#' being defined as the solution which minimizes the cumulated errors weighted by the 
#' user weights.   
#'
#' @param sol the current solution (a named list)
#' @param chain the chain to be solved (a list of strings containing variable names)
#' @param case the case to solve
#' @param nA the target population size for A
#' @param nB the target population size for B
#' @param nu.A control for nA: 0 means "respect nA", non-null "adapt it to solve the case"
#' @param phi.A control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param delta.A control for degree A: 0 means "respect the input parameters pdi", non-null "adapt them to solve the case"
#' @param gamma control for pij: 0 means "respect the matching probabilities pij", non-null "adapt them to solve the case"
#' @param delta.B control for degree B: 0 means "respect the input parameters pdj", non-null "adapt them to solve the case"
#' @param phi.B control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param nu.B control for nB: 0 means "respect nB", non-null "adapt it to solve the case"
#' @param verbose if TRUE, will display detailed information on the console
#' @return a list of vectors (the chains) of strings 
#'
#' @seealso \code{\link{propagate.direct}} for the inference of the consequences of the hypothesis,
#'          \code{\link{detect.problems}} to ensure potential solutions are consistent
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @importFrom utils combn
#' 
resolve.missing.chain <- function(sol, chain, case, 
                                    nA, nB, 
                                    nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma, 
                                    verbose=FALSE) {

    if (verbose)
        cat("\tstarting the investigation of the missing chain: ",paste(chain,collapse=","),"\n",sep="")

    solutions <- list()

    testable_variables <- intersect(
                            chain,
                            c("hat.di","hat.dj","hat.fi","hat.fj","hat.pij","hat.nA","hat.nB")
                            )
                          
    hypothesis <- unlist(lapply(1:length(testable_variables), function(x) combn(testable_variables,x,simplify=FALSE)), recursive=FALSE)

    if (verbose)
        cat("\t\twill test ",length(hypothesis)," hypothesis\n",sep="")

    # iterate all the missing elements of the chain that we might put hypothesis on
    for (tested_variables in hypothesis) {

        var_names <- lapply(tested_variables, function (x) substr(x, 5, nchar(x)))
        hypothesis_name <- paste(
                    paste(
                        tested_variables,
                        var_names,
                        sep="="
                        ),
                    collapse=","
                )
        if (verbose) {
            cat("\t\ttrying with hypothesis: ", hypothesis_name, "\n", sep="")
        }
        
        sol.hyp <- sol

        explored <- list()

        explored$investigated.var.name <- hypothesis_name
        explored$hypothesis <- setdiff(chain,tested_variables)

        for (var.hat.name in tested_variables) {

            if (var.hat.name == "hat.di") { sol.hyp[[var.hat.name]] <- case$inputs$di * case$masks$mask.di.all }
            else if (var.hat.name == "hat.dj") { sol.hyp[[var.hat.name]] <- case$inputs$dj * case$masks$mask.dj.all }
            else if (var.hat.name == "hat.fi") { sol.hyp[[var.hat.name]] <- normalise(case$stats$fi * case$masks$mask.fi.all) }
            else if (var.hat.name == "hat.fj") { sol.hyp[[var.hat.name]] <- normalise(case$stats$fj * case$masks$mask.fj.all)}
            else if (var.hat.name == "hat.pij") { sol.hyp[[var.hat.name]] <- normalise(case$inputs$pij$data * case$masks$mask.pij.all) }
            else if (var.hat.name == "hat.nA") { sol.hyp[[var.hat.name]] <- nA }
            else if (var.hat.name == "hat.nB") { sol.hyp[[var.hat.name]] <- nB }

            else { 
                    # we should never reach this step. 
                    # we screwed up on the list of variables at the head of the for loop
                    stop("/!\\ cannot create such an hypothesis\n")
            }

        }


        if (is.null(sol.hyp$inference)) {
            sol.hyp$inference <- list()
        }
        sol.hyp$inference <- list(sol.hyp$inference, list(paste("hypothesis on ",hypothesis_name,sep="")))
        
        # propagate on this basis
        sol.hyp <- tryCatch({
                        propagate.direct(sol.hyp, case, verbose=verbose, indent=3)
                    }, error = function(e) {
                        if (verbose) {
                            cat("\t\t\tfound an invalid solution.\n")
                        }
                        NULL
                    })     
        
        # so know we have an estimation
        if (!is.null(sol.hyp)) {
            if (!all(chain %in% names(sol.hyp))) {
                # we cannot infer anything from this hypothesis; let's ignore it
                if (verbose) {
                    cat("\t\t\tthese hypothesis do not enable the resolution of all the variables.\n")
                }
            } else {
                
                # do we achieve to solve the problem on this basis ?
                nb.problems <- detect.problems(sol.hyp, case, fail=FALSE, verbose=verbose, indent=3)

                if (nb.problems > 0) {
                    #cat("this variable does not provides a valid solution to our problem","\n")
                    if (verbose) {
                        cat("\t\t\twe found an invalid solution","\n")
                    }   
                } else {
                    if (verbose) {
                        cat("\t\t\twe found a valid solution which provides: ", paste.known(sol.hyp),"\n")
                    }
                    sol.hyp <- quantify.errors(sol.hyp, case, nA, nB)
                    if (verbose) {

                        cat(paste("here are the NRMSE errors of this solution: ",
                            "nA=",sol.hyp$nrmse.nA,
                            ",fi=",sol.hyp$nrmse.fi, 
                            ",di=",sol.hyp$nrmse.di,
                            ",pij=",sol.hyp$nrmse.pij, 
                            ",dj=",sol.hyp$nrmse.dj, 
                            ",fj=",sol.hyp$nrmse.fj,
                            ",nB=",sol.hyp$nrmse.nB,
                            "\n",
                            sep=""
                            ))
                    }
                    explored$sol <- sol.hyp
                    solutions[[length(solutions)+1]] <- explored
                }
                
            }
        }


    }

    # we did everything we could.
    # now "solutions" contains the valid solutions which were found 
    res <- NULL 
    if (length(solutions) == 0) {
        if (verbose) 
            cat("\t\t\tfound no solution for this chain :-(\n")
    } else if (length(solutions) == 1) {
        # return the unique solution
        res <- solutions[[1]]$sol
        if (verbose)
            cat("\t\t\tfound one valid solution\n")
    } else {

        # create a vector with errors
        val.or.0 <- function(x) {
            if (is.null(x) || is.nan(x)) {
                return(0)
            } else {
                return(x)
            }
        }
        cumulated.errors <- rep(0, times=length(solutions))
        weights <- c(nu.A, phi.A, delta.A, gamma, delta.B, phi.B, nu.B)
        weights.names <- c("nu.A", "phi.A", "delta.A", "gamma", "delta.B", "phi.B", "nu.B")
        indices_weights_not_null <- which(weights != 0)


        if (verbose)
            cat("\t\tfound ",length(solutions)," solutions, we have to select the best according to weights",weights.names[indices_weights_not_null],"\n")
            weights <- c(nu.A, phi.A, delta.A, gamma, delta.B, phi.B, nu.B)

        # print("weights")
        # print(weights)
        for (i in 1:length(solutions)) {

            s <- solutions[[i]]
            # print(names(s$sol))

            errors <- c(
                        val.or.0(s$sol$nrmse.nA), 
                        val.or.0(s$sol$nrmse.fi), 
                        val.or.0(s$sol$nrmse.di), 
                        val.or.0(s$sol$nrmse.pij), 
                        val.or.0(s$sol$nrmse.dj), 
                        val.or.0(s$sol$nrmse.fj), 
                        val.or.0(s$sol$nrmse.nB)
                        )


            # print("errors for this sol")
            # print(errors[indices_weights_not_null])
            # print("weighted")
            # print(errors[indices_weights_not_null] / weights[indices_weights_not_null])
            cumulated.errors[i] <- sum(errors[indices_weights_not_null] / weights[indices_weights_not_null])
            
            if (verbose)
                cat(paste("\t\t\tsolution (",i,") => ", errors[indices_weights_not_null], "\t weighted: ",cumulated.errors,"\n",sep=""))
        }

        best.solutions <- which(cumulated.errors == min(cumulated.errors)) 
        best.solution <- NULL
        if (length(best.solutions) > 1) {
            best.solution <- sample(best.solutions, 1)
            # TODO warning ? systematic message ?
            if (verbose) {
                cat("\t\t\tthe best solutions are solutions ",paste(best.solutions,collapse=",")," based on hypothesis:\n",sep="")
                for (idx in best.solutions)
                    cat("\t\t\t\t* solution ",idx," solved with:\t",solutions[[idx]]$investigated.var.name,"\n",sep="")
                #cat("\t\t\tthere are multiple best solutions (that is: ",length(best.solutions),") ; just selecting one randomly\n",sep="")
                cat("\t\t=> selected ",best.solution," which is one of the ",length(best.solutions)," solutions having lowest weighted error\n",sep="")
            }
        } else {
            best.solution <- best.solutions[[1]]
            if (verbose)
                cat("\t\t=> will keep solution ",best.solution," which is the one with the lowest weighted error\n",sep="")
        }
        if (verbose) {
            cat("\t\t\tthis solution is based on hypothesis: ",solutions[[best.solution]]$investigated.var.name,"\n",sep="")
        }
        res <- solutions[[best.solution]]$sol
        #print(names(solutions[[best.solution]]$sol))
    }

    res
}

#' Resolves a partial solution by inference and hypothesis
#' 
#' Resolves a given case starting with the base solution
#' accordng to user parameters. 
#' It first uses \code{\link{propagate.direct}} to infer 
#' results from available elements; then it detects the 
#' chains of variable which are not constrained enough to be solved this way 
#' using \code{\link{detect.missing.chains}}. Then it solves each chain 
#' iteratively until the entire problem is solved. 
#' 
#' @param sol the current solution (a named list)
#' @param case the case to solve
#' @param nA the target population size for A
#' @param nB the target population size for B
#' @param nu.A control for nA: 0 means "respect nA", non-null "adapt it to solve the case"
#' @param phi.A control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param delta.A control for degree A: 0 means "respect the input parameters pdi", non-null "adapt them to solve the case"
#' @param gamma control for pij: 0 means "respect the matching probabilities pij", non-null "adapt them to solve the case"
#' @param delta.B control for degree B: 0 means "respect the input parameters pdj", non-null "adapt them to solve the case"
#' @param phi.B control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param nu.B control for nB: 0 means "respect nB", non-null "adapt it to solve the case"
#' @param verbose if TRUE, will display detailed information on the console
#' 
#' @return a list of vectors (the chains) of strings 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
resolve <- function(sol, case, 
                        nA, nB, 
                        nu.A, phi.A, delta.A, 
                        gamma,
                        nu.B, phi.B, delta.B,  
                        verbose=FALSE) {

    # direct propagation of what we know
    sol.tmp <- propagate.direct(sol, case, verbose=verbose)

    missing.chains <- detect.missing.chains(sol.tmp)
    while (length(missing.chains) > 0) {

        chain <- missing.chains[[1]]
        missing.chains[[1]] <- NULL

        found <- resolve.missing.chain(sol.tmp, chain, case, 
                    nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma, 
                    verbose=verbose)

        if (!is.null(found)) {
            if (verbose) {
                cat("\tthis missing chain was solved\n", sep="")
            }
            sol.tmp <- found 
            missing.chains <- detect.missing.chains(sol.tmp)
        }

    }


    # print("propagated")
    # print(sol)

    return(sol.tmp)
}


#' Quantifies the errors of a solution
#' 
#' Computes the errors in a (possibly partial) solution 
#' for the variables which are defined and adds them 
#' to these solutions. For instance if \code{sol$hat.fi}
#' is available, \code{sol$mse.fi} will be computed. 
#' 
#' @param sol the solution to evaluate
#' @param case the reference case
#' @param nA the expected size of population A
#' @param nB the expected size of population B
#' @return the solution with additional variables for errors 
#'  
quantify.errors <- function(sol, case, nA, nB) {


    # measure errors
    # MSE fi
    if (!is.null(sol$hat.fi)) {
        sol$mse.fi <- mean( ( sol$hat.fi - case$stats$fi)^2  )
        sol$rmse.fi <- sqrt(sol$mse.fi)
        sol$nrmse.fi <- sol$rmse.fi
    }
    if (!is.null(sol$hat.fj)) {
        sol$mse.fj <- mean( ( sol$hat.fj - case$stats$fj)^2  )
        sol$rmse.fj <- sqrt(sol$mse.fj)
        sol$nrmse.fj <- sol$rmse.fj
    }

    # MSE di
    if (!is.null(sol$hat.di)) {
        sol$mse.di <- mean( ( sol$hat.di - case$inputs$di)^2  )
        sol$rmse.di <- sqrt(sol$mse.di)
        sol$nrmse.di <- sol$rmse.di
    }
    if (!is.null(sol$hat.dj)) {
        sol$mse.dj <- mean( ( sol$hat.dj - case$inputs$dj)^2  )
        sol$rmse.dj <- sqrt(sol$mse.dj)
        sol$nrmse.dj <- sol$rmse.dj
    }

    # MSE pdi 
    if (!is.null(sol$hat.pdi)) {
        sol$mse.pdi <- mean( ( sol$hat.pdi - case$inputs$pdi$data)^2  )
        sol$rmse.pdi <- sqrt(sol$mse.pdi)
        sol$nrmse.pdi <- sol$rmse.pdi
    }
    if (!is.null(sol$hat.pdj)) {
        sol$mse.pdj <- mean( ( sol$hat.pdj - case$inputs$pdj$data)^2  )
        sol$rmse.pdj <- sqrt(sol$mse.pdj)
        sol$nrmse.pdj <- sol$rmse.pdj
    }

    # MSE pij
    if (!is.null(sol$hat.pij)) {
        sol$mse.pij <- mean( ( sol$hat.pij - case$inputs$pij$data )^2  )
        sol$rmse.pij <- sqrt(sol$mse.pij)
        sol$nrmse.pij <- sol$rmse.pij
    }

    # error n
    if (!is.null(sol$hat.nA)) {
        sol$error.nA <- abs( sol$hat.nA - nA )
        sol$rmse.nA <- sol$error.nA
        sol$nrmse.nA <- sol$rmse.nA / max(sol$hat.nA, nA)
    }
    if (!is.null(sol$hat.nB)) {
        sol$error.nB <- abs( sol$hat.nB - nB )
        sol$rmse.nB <- sol$error.nB
        sol$nrmse.nB <- sol$rmse.nB / max(sol$hat.nB, nB)
    }

    sol
}

#' Ensures the content of a solution is in the expect format
#' 
#' Ensures a solution is formatted well: notably guarantees
#' the lists or data frames which have to be names are named 
#' as expected.
#' 
#' @param sol the solution (might be partial)
#' @param case the case
#' @return the same solution fixed formats
#' 
ensure.form <- function(sol, case) {

    if (!is.null(sol$hat.ni)) {
        names(sol$hat.ni) <- colnames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.ci)) {
        names(sol$hat.ci) <- colnames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.fi)) {
        names(sol$hat.fi) <- colnames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.di)) {
        names(sol$hat.di) <- colnames(case$inputs$pij$data)
    }

    if (!is.null(sol$hat.nj)) {
        names(sol$hat.nj) <- rownames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.cj)) {
        names(sol$hat.cj) <- rownames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.fj)) {
        names(sol$hat.fj) <- rownames(case$inputs$pij$data)
    }
    if (!is.null(sol$hat.dj)) {
        names(sol$hat.dj) <- rownames(case$inputs$pij$data)
    }

    sol
}

#' Solves the equations by arbitrating.
#' 
#' Solves the underlying equations based on data, the inputs parameters, 
#' and the control parameters. It distributes the errors in the places 
#' accepted by the user. 
#'
#' This function is deterministic.
#'
#' @param case a case prepared with \code{\link{matching.prepare}}
#' @param nA the count of entities expected for population A
#' @param nB the count of entities expected for population B
#' @param nu.A control for nA: 0 means "respect nA", non-null "adapt it to solve the case"
#' @param phi.A control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param delta.A control for degree A: 0 means "respect the input parameters pdi", non-null "adapt them to solve the case"
#' @param gamma control for pij: 0 means "respect the matching probabilities pij", non-null "adapt them to solve the case"
#' @param delta.B control for degree B: 0 means "respect the input parameters pdj", non-null "adapt them to solve the case"
#' @param phi.B control for frequencies: 0 means "respect the original frequencies as detected in the sample", non-null "adapt it to solve the case"
#' @param nu.B control for nB: 0 means "respect nB", non-null "adapt it to solve the case"
#' @param verbose if TRUE, the resolution will emit messages for debug
#' 
#' @return a case ready for generation
#'
#' @export
#'
#' @seealso \code{\link{matching.prepare}} to prepare a case for this function, \code{\link{matching.generate}} to use the result for actual generation 
#' 
#' @examples
#'
#' # load sample data
#' data(cas1)
#' # prepare the case  
#' case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#' # resolve tbe case
#' disc <- matching.solve(case.prepared, 
#'                     nA=50000,nB=40000, 
#'                     nu.A=0, phi.A=0, delta.A=1, 
#'                     gamma=1, 
#'                     delta.B=0, phi.B=0, nu.B=0)
#' # print the resolution information
#' print(disc)
#' # access the solved frequencies, distribution of degrees and matching probabilities
#' print(disc$gen$hat.fi)
#' print(disc$gen$hat.pdi)
#' print(disc$gen$hat.pij)
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.solve <- function(case, 
                                nA, nB, 
                                nu.A, phi.A, delta.A, 
                                gamma, 
                                delta.B, phi.B, nu.B, 
                                verbose = FALSE) {

    if (class(case) != "dpp_prepared") 
        stop("case should be the result of a preparation of data by matching.prepare")

    mixx <- list()

    # start with masks
    # TODO reuse
    case$masks <- list()
    case$masks$mask.fi <- as.integer(case$stats$fi > 0)
    case$masks$mask.fj <- as.integer(case$stats$fj > 0)
    case$masks$mask.di <- as.integer(case$inputs$di > 0)
    case$masks$mask.dj <- as.integer(case$inputs$dj > 0)
    case$masks$mask.pij <- (case$inputs$pij$data > 0)*1

    case$masks$mask.pij.all <-  t(t(case$masks$mask.pij) * case$masks$mask.fi * case$masks$mask.di) * case$masks$mask.fj * case$masks$mask.dj
    case$masks$mask.fi.all <- as.integer(colSums(case$masks$mask.pij.all) > 0)
    case$masks$mask.fj.all <- as.integer(rowSums(case$masks$mask.pij.all) > 0)
    case$masks$mask.di.all <- case$masks$mask.di * case$masks$mask.fi.all
    case$masks$mask.dj.all <- case$masks$mask.dj * case$masks$mask.fj.all

    #print(case$masks)

    if (verbose) 
        cat("\nstarting the resolution of the case\n")
    

    sol <- list()
    if (nu.A == 0)  {   sol$hat.nA <- nA }
    if (phi.A == 0) {   sol$hat.fi <- normalise(case$stats$fi * case$masks$mask.fi.all) }
    if (delta.A == 0) { sol$hat.di <- case$inputs$di * case$masks$mask.di.all }

    if (nu.B == 0) {    sol$hat.nB <- nB }
    if (phi.B == 0) {   sol$hat.fj <- normalise(case$stats$fj * case$masks$mask.fj.all) }
    if (delta.B == 0) { sol$hat.dj <- case$inputs$dj * case$masks$mask.dj.all }

    if (gamma == 0) {   sol$hat.pij <- normalise(case$inputs$pij$data * case$masks$mask.pij.all) }

    if (verbose) {
        cat("\taccording to user weights, we already know: ",paste(names(sol),collapse=","),"\n",sep="")
    }

    # detect issues right now; maybe the problem is overconstrained or badly constrainted
    detect.problems(sol, case, verbose=verbose)

    sol <- resolve(sol, case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma, verbose=verbose)

    # ensure all the variables found a solution during the process
     if (!all(c(
            "hat.pdi","hat.pdj",
            "hat.ci","hat.cj",
            "hat.nij",
            "hat.nA","hat.nB") %in% names(sol))) {
        stop("The case is too constrained to be solved: some of the variables were not solved (",
            paste(setdiff(c(
                "hat.pdi","hat.pdj",
                "hat.ci","hat.cj",
                "hat.nij",
                "hat.nA","hat.nB"),
                names(sol)),
                collapse=","),
            ")")
    } 

    detect.problems(sol, case, verbose=verbose)

    if (verbose) {
        cat("\ncase solved.\n")
    }
    # measure errors
    sol <- quantify.errors(sol, case, nA, nB)

    # ensure the form is ok
    sol <- ensure.form(sol, case)

    # prepare the result
    res <- case
    res$inputs$nA <- nA
    res$inputs$nB <- nB
    res$inputs$nu.A <- nu.A 
    res$inputs$nu.B <- nu.B 
    res$inputs$phi.A <- phi.A 
    res$inputs$phi.B <- phi.B
    res$inputs$delta.A <- delta.A 
    res$inputs$delta.B <- delta.B
    res$inputs$gamma <- gamma 

    res$gen <- sol


    class(res) <- "dpp_resolved"
       
    res
}

#' Display the result of a resolution
#' 
#' @param x the case to print
#' @param ... ignored
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_resolved <- function(x,...) {

    cat("case prepared and ready for generation\n")
    cat("control parameters: nu.A=",x$inputs$nu.A,
            ", phi.A=",x$inputs$phi.A,", delta.A=",x$inputs$delta.A,
            ", gamma=",x$inputs$gamma,
            ", delta.B=",x$inputs$delta.B,", phi.B=",x$inputs$phi.B,", nu.B=",x$inputs$nu.B,
            "\n\n",
            sep="")

    disp <- function(name,x) {
        cat("$",name,":\n",sep="")
        print(x$gen[[name]])
        cat("\n")
    }
    disperr <- function(name,x) {
        subname <- substr(name,5,nchar(name))
        msename <- paste("rmse.",subname,sep="")
        cat("$",name," [RMSE ",x$gen[[msename]],"]:\n",sep="")
        print(x$gen[[name]])
        cat("\n")
    }

    cat("$hat.nA: ",x$gen$hat.nA," [err ",x$gen$error.nA,"]\n\n",sep="")

    disperr("hat.fi",x)
    disp("hat.ci",x)
    disperr("hat.di",x)
    disperr("hat.pdi",x)
    disp("hat.ni",x)

    cat("$hat.nL: ",x$gen$hat.nL,"\n",sep="")
    disp("hat.nij",x)
    disperr("hat.pij",x)

    disp("hat.nj",x)
    disperr("hat.dj",x)
    disperr("hat.pdj",x)
    disp("hat.cj",x)
    disperr("hat.fj",x)

    cat("$hat.nB: ",x$gen$hat.nB," [err ",x$gen$error.nB,"]\n",sep="")
      
    cat("The problem was resolved with the following process:\n")
    cat("\t",paste(x$gen$inference,collapse="\n\t"),sep="")
    cat("\n")
} 

