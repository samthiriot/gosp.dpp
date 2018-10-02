

#' Rounds an object by ensuring its totals remains equal before and after rounding
#'
#' Freely inspired by \url{https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum}
#' 
#' @param x the object to round
#' @param ... ignored
#' @return the rounded version of the object
#'
#' @export
#'
round_sum <- function (x, ...) {
    UseMethod("round_sum", x)
}

round_sum.numeric <- function(x, ...) {
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

round_sum.matrix <- function(x, ...) {

    # convert the matrix to a vector (by columns)
    y <- unlist(x)

    # round it 
    r <- round_sum.numeric(y)

    # convert it to a matrix
    mat <- matrix(r, ncol=ncol(x), nrow=nrow(x))

    # convert it back to a matrix
    res <- as.data.frame(mat, row.names=row.names(x), col.names=colnames(x))
    colnames(res) <- colnames(x)
    
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
#' @param verbose more text if TRUE (FALSE by default)
#' @param precision for comparisons
#'
#' @author samuel.thiriot@res-ear.ch
#'
#' @keywords internal
#'
update_degree_distribution.col <- function(pdx, t.target, verbose=FALSE, precision=1e-10) {

    np.orig <- pdx * 0:(length(pdx)-1)
    t.orig <- sum(np.orig)
    # cat("t.orig", t.orig, "\n")
    
    # quick exit if there is nothing to do
    if (abs(t.target - t.orig) <= precision) {
        return(pdx)
    }

    # quick case: if the target is 0, then the solution is (1,0,.....,0) 
    if (t.target <= precision) {
        pdx[0] <- 1
        for (i in 2:length(pdx)) {
            pdx[i] <- 0
        }
        return(pdx)
    }


    if (verbose) {
        cat("should update a column of degree distribution pdx=", 
            paste(pdx, collapse=","), "such as it sums up to ", t.target, 
            "instead of the current sum ", t.orig, "\n")
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
    # print("np.potential.cumulated")
    # print(np.potential.cumulated)

    # create factors, 0 or 1 
    factors <- pmin( (np.potential == np.max) + (p.potential < 0), 1)
    # print("factors")
    # print(factors)

    rat <- (t.target - t.orig) / np.potential.cumulated
    # cat("rat",rat,"\n")


    p.add <- factors * rat * pmin(p.potential, -p.potential.sum.neg)

    # cat("p.add", p.add,"summing up to ", abs(sum(p.add)),"\n")
    # cat("p.add",p.add,"=",sum(p.add),"\n")
    if (abs(sum(p.add)) >= precision) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }

    # print("p.modified ------------------------------")
    p.modified <- p.add + pdx
    # cat("p.modified",p.modified,"=",sum(p.modified),"\n")
    if (length(which(p.modified < 0)) > 0) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far - would have probabilities below 0)")
    }
    if (length(which(p.modified > 1)) > 0) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far - would have probabilities greater than 1)")
    }

    # cat("p.modified", p.modified, abs(sum(p.modified)-1), "\n")
    if (abs(sum(p.modified)-1) > precision) {
        stop("The case is too constrained to be solved (cannot adapt the degree probabilities that far)")
    }
     
    # now ensure this leads to the expected result 
    np.res.cumulated <- sum(p.modified * 0:(length(pdx)-1))
    # cat("np.res.cumulated",np.res.cumulated,"for target", t.target,"\n")

    p.modified
}

#' Update a given distribution of degrees probabilities so its matching target degrees
#' 
#' This process is done column by column by calling \link{update_degree_distribution.col}.
#' '
#' @param pdx the distribution of degrees 
#' @param dx the vector of the expected average degrees
#' @param precision the precision for comparison with 0; change only if suggested and relevant
#' @param verbose more messages if TRUE
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
#' @keywords internal
#'
update_degree_distribution <- function(pdx, dx, precision=.Machine$double.eps, verbose=F) {

    if (class(pdx) != "data.frame") {
        stop("pdx should be a data.frame")
    }
    if (class(dx) != "numeric") {
        stop("dx should be numeric")
    }
    if (ncol(pdx) != length(dx)) {
        stop("length of dx (",length(dx),") should match the columns count of pdx (",ncol(pdx),")")
    }

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
        pdx.reweighted[,c] <- update_degree_distribution.col(pdx[,c],dx[c], verbose=verbose, precision=precision)
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
#' @keywords internal
#'
rectify.degree.counts <- function(pdn, nn, cn, verbose=FALSE) {

    if (verbose) {
        cat("should rectify\n") 
        print(pdn)
        cat("to reach\n")
        print(nn)
    }
    for (col in 1:ncol(pdn)) {

        # if (verbose) cat("col ", col, "\n")

        total.current <- sum(pdn[,col]*0:(nrow(pdn)-1))
        total.expected <- nn[col]

        if (total.current > total.expected) {
            to.remove <- total.current - total.expected

            if (verbose) cat("col ", col, ": should remove ", to.remove, "\n")
            # no more warning: we anyway check at the end how well we fixed it
            # if (to.remove > 1) {
            #     warning("/!\\ the solution to fix rounding errors below 1 slot (here:",to.remove,") is not managed well\n.",sep="")
            # }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if ( (pdn[2,col] >= to.remove) ) { # && (pdn[1,col] > 0)
                pdn[2,col] <- pdn[2,col] - to.remove
                pdn[1,col] <- pdn[1,col] + to.remove
            } else if ( ( nrow(pdn) >= 3) && (pdn[3,col] >= to.remove/2) ) { # && (pdn[1,col] > 0)
                # cat("removing 2 and adding 1\n")
                # to remove -1, we also might remove 2 and add 1 !
                # print(pdn[,col])
                pdn[3,col] <- pdn[3,col] - to.remove
                pdn[2,col] <- pdn[2,col] + to.remove
                # print(pdn[,col])
            } else {
                if (verbose)
                    warning("/!\\ found no good solution to fix this rounding error by removing ",to.remove,"\n") 
            }

        } else if (total.current < total.expected) {
            to.add <- total.expected - total.current 

            if (verbose) cat("col ", col, ": should add ", to.add, "\n")

            # no more warning: we anyway check at the end how well we fixed it
            # if (to.add > 1) {
            #     warning("/!\\ the solution to fix rounding errors above 1 slot (here:",to.add,") is not managed well\n", sep="")
            # }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if (pdn[1,col] >= to.add) { # (pdn[2,col] > 0
                pdn[1,col] <- pdn[1,col] - to.add
                pdn[2,col] <- pdn[2,col] + to.add
            } else if ( (nrow(pdn) >= 3) && (pdn[3,col] >= to.add) ) { # (pdn[2,col] > 0
                pdn[2,col] <- pdn[2,col] - to.add
                pdn[3,col] <- pdn[3,col] + to.add
            } else {
                if (verbose) {
                    warning("/!\\ found no good solution to fix this rounding error of an additional ",to.add,"\n") 
                    print(pdn[,col])
                }
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

        # print(pdn)
        # print(nn)
        # print(cn)
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
#' @param d any data to normalise
#' @param ... additional parameters are likely to be ignored
#' @return the same object normalised
#'
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'
normalise <- function(d, ...) {
   UseMethod("normalise", d)
}
normalise.data.frame <- function(d, ...) {
    d / sum(d)
}
normalise.numeric <- function(d, ...) {
    d / sum(d) 
}
normalise.list <- function(d, ...) {
    total <- sum.list(d)
    lapply(d, "/", total) 
}
normalise.dpp_degree_cpt <- function(d, ...) {
    for (i in 1:ncol(d$data)) {
        d$data[i] <- normalise(d$data[i])  
    }
    d
}

#' Sums applied on a list
#' 
#' sums the content of the list (required for recent R versions)
#'
#' @param l the list to sum
#' @return the same object summed with operator "+"
#'
#' @keywords internal
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'
sum.list <- function(l) {
    Reduce("+", l) 
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
#' @keywords internal
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
#' @keywords internal
#'
inf_to_zeros <- function(vv) {
    vv[which(is.infinite(vv))] <- 0
    vv
}

#' A simple Iterated Proportional Fitting implementation 
#' 
#' taken from https://stats.stackexchange.com/questions/59115/iterative-proportional-fitting-in-r
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
#' @seealso The entry point for IPF 2D: \code{\link{ipf.2d}}
#' 
#' @keywords internal
#'
ipf.2d.stackoverflow <- function(Margins_, seedAry, maxiter=100, closure=0.001, verbose=FALSE) {
    #Check to see if the sum of each margin is equal
    MarginSums. <- unlist(lapply(Margins_, sum))
    if(any(MarginSums. != MarginSums.[1])) warning("IPF: sum of each margin not equal")

    #Replace margin values of zero with 0.001
    Margins_ <- lapply(Margins_, function(x) {
        # TODO ???
        if(any(x == 0)) warning("IPF: zeros in marginsMtx replaced with 0.001") 
        x[x == 0] <- 0.001
        x
    })

    #Check to see if number of dimensions in seed array equals the number of
    #margins specified in the marginsMtx
    numMargins <- length(dim(seedAry))
    if(length(Margins_) != numMargins) {
        stop("IPF: number of margins in marginsMtx not equal to number of margins in seedAry")
    }

    #Set initial values
    resultAry <- seedAry
    iter <- 0
    marginChecks <- rep(1, numMargins)
    margins <- seq(1, numMargins)

    #Iteratively proportion margins until closure or iteration criteria are met
    while((any(marginChecks > closure)) & (iter < maxiter)) {
        for(margin in margins) {
            marginTotal <- apply(resultAry, margin, sum)
            marginCoeff <- Margins_[[margin]]/marginTotal
            marginCoeff[is.infinite(marginCoeff)] <- 0
            resultAry <- sweep(resultAry, margin, marginCoeff, "*")
            marginChecks[margin] <- sum(abs(1 - marginCoeff))
        }    
        iter <- iter + 1
    }

    #If IPF stopped due to number of iterations then output info
    if(verbose && (iter == maxiter)) cat("IPF stopped due to number of iterations\n")

    #Return balanced array
    resultAry
}


#' Applies Iterative Proportional Fitting on a 2d dataframe
#' 
#' Takes as inputs a dataframe representing a 2d matrix, 
#' and target marings for cols and rows. 
#' Returns a table reweighted so its sums fit margins.
#' 
#' Either falls back to a local implementation \code{\link{ipf.2d.stackoverflow}}, or to an
#' higher quality one from the mipfp package if it is installed. 
#'
#' @param df the dataframe to reweight
#' @param margins_cols a vector containing the target marginals for cols
#' @param margins_rows a vector containing the target marginals for rows
#' @param verbose displays more messages if TRUE (FALSE by default)
#' @param max.iterations the maximum iterations before stopping (1000 by default)
#' @param precision the target precision, will stop once reached (defaults to 1e-6)
#' @return the dataframe reweighted
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @export
#'
ipf.2d <- function(df, margins_cols, margins_rows, max.iterations=1000, precision=1e-10, verbose=FALSE) {

    if (class(df) != "data.frame") 
        stop("df should be a data frame")
    if (class(margins_cols) != "numeric") 
        stop("margins_cols should be a numeric vector")
    if (class(margins_rows) != "numeric") 
        stop("margins_rows should be a numeric vector")

    if (length(margins_rows) != nrow(df)) 
        stop(paste("the length of margins_rows (", length(margins_rows) , ") should match the count of rows (", nrow(df), ") from df", sep=""))
    if (length(margins_cols) != ncol(df)) 
        stop(paste("the length of margins_cols (", length(margins_cols) , ") should match the count of columns (", ncol(df), ") from df", sep=""))

    if (FALSE && verbose) {
        cat("should reweight the matrix\n")
        print(df)
        cat("in order to reach column marginals: ", paste(margins_cols, collapse=","), "\n")
        cat("in order to reach row marginals: ", paste(margins_rows, collapse=","), "\n")
    }

    if (requireNamespace("mipfp", quietly = TRUE)) {
        # the mipfp package is available, let's use it, its a reference implementation
        if (verbose) cat("using for IPF the implementation from package mipfp\n")

        res <- mipfp::Ipfp(
            seed=as.matrix(df), 
            target.list=list(2,1), 
            target.data=list(margins_cols, margins_rows), 
            print = verbose, 
            iter = max.iterations, 
            tol = 1e-10, tol.margins = precision , na.target = FALSE
            )

        if (!res$conv) {
            # detect failure
            print(res)
            stop("IPF did not converged properly. That's a bit unexpected. Stopping.")
        }
        # extract the result of interest
        as.data.frame(res$p.hat)

    } else {
        # fallback to a local, less powerful implementation
        res <- ipf.2d.stackoverflow(
            Margins_=list(margins_rows, margins_cols), 
            seedAry=as.matrix(df), 
            maxiter=max.iterations, 
            closure=precision,
            verbose=verbose
            )
        normalise(as.data.frame(res))
    }
    
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
#' @param precision.pd the precision to be used for probabilistic distributions of degrees comparisons
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
#' @keywords internal
#'
propagate.direct <- function(sol, case, precision.pd=.Machine$double.eps, verbose=FALSE, indent=1) {

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
            # print(sol$hat.ci)
            # print(sol$hat.di)
            # print(sol$hat.ci * sol$hat.di)
            sol$hat.ni <- round_sum(sol$hat.ci * sol$hat.di)
            # adapt rounding
            sol$hat.di <- nan_to_zeros(sol$hat.ni / sol$hat.ci)
            # adapt because of rounding
            sol$hat.pi <- normalise(sol$hat.ni)
            sol <- info.rule("hat.ci, hat.di -> hat.ni [hat.pi]", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        # forward cj + dj -> nj 
        if ( 
            ("hat.cj" %in% names(sol)) && ("hat.dj" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- round_sum(sol$hat.cj * sol$hat.dj)
            sol$hat.dj <- nan_to_zeros(sol$hat.nj / sol$hat.cj)
            # adapt because of rounding
            sol$hat.pj <- normalise(sol$hat.nj)
            sol <- info.rule("hat.cj, hat.dj -> hat.nj [hat.pj]", sol, verbose=verbose, indent=indent+1)
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
            && ("hat.pij" %in% names(sol))
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
            sol <- info.rule("hat.ni, hat.pij -> hat.nij [hat.pij]", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }
        if ( 
            ("hat.nj" %in% names(sol))
            && ("hat.pij" %in% names(sol))
            && (!"hat.nij" %in% names(sol))
            ) {
            pij <- if ("hat.pij" %in% names(sol)) sol$hat.pij else case$inputs$pij$data
            sol$hat.nij <- nan_to_zeros( sol$hat.nj * pij / rowSums(pij) )
            # round per line (to keep the sums consistent)
            for (i in 1:nrow(sol$hat.nij)) {
                sol$hat.nij[i,] <- round_sum(sol$hat.nij[i,])
            }
            sol$hat.pij <- sol$hat.nij / sum(sol$hat.nij)
            sol <- info.rule("hat.nj, hat.pij -> hat.nij [hat.pij]", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }

        if ( 
            ("hat.nL" %in% names(sol))
            && ("hat.pij" %in% names(sol))
            && (!"hat.nij" %in% names(sol))
            ) {

            sol$hat.nij <- round_sum(sol$hat.pij * sol$hat.nL)
            
            # round.
            # if the totals of ci are known, then round per column 
            # (same for cj and rows)
            # if ("hat.ci" %in% names(sol)) {
            #     for (i in 1:ncol(sol$hat.nij)) {
            #         sol$hat.nij[,i] <- round_sum(sol$hat.nij[,i])
            #     }
            # } else if ("hat.cj" %in% names(sol)) {
            #     for (i in 1:nrow(sol$hat.nij)) {
            #         sol$hat.nij[i,] <- round_sum(sol$hat.nij[i,])
            #     }
            # } 
            # TODO what if both ci and cj ?
            # recompute pij, in case rounding would make things inconsistent
            sol$hat.pij <- normalise(sol$hat.nij)
            sol <- info.rule("hat.nL, hat.pij -> hat.nij [hat.pij]", sol, verbose=verbose, indent=indent+1)
            changed <- TRUE
        }


        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.ni" %in% names(sol))
            ) {
            sol$hat.ni <- colSums(sol$hat.nij)
            if ("hat.pi" %in% names(sol)) {
                sol <- info.rule("hat.nij -> hat.ni [hat.pi]", sol, verbose=verbose, indent=indent+1)
                sol$hat.pi <- normalise(sol$hat.ni)
            } else {
                sol <- info.rule("hat.nij -> hat.ni", sol, verbose=verbose, indent=indent+1)
            }
            changed <- TRUE
        }
        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- rowSums(sol$hat.nij)
            if ("hat.pj" %in% names(sol)) {
                sol <- info.rule("hat.nij -> hat.nj [hat.pj]", sol, verbose=verbose, indent=indent+1)
                sol$hat.pj <- normalise(sol$hat.nj)
            } else {
                sol <- info.rule("hat.nij -> hat.nj", sol, verbose=verbose, indent=indent+1)
            }
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

        # TODO IPF

        # pi, pj -> pij
        # (IPF)
        if ( 
            ("hat.pi" %in% names(sol)) && ("hat.pj" %in% names(sol)) && (!"hat.pij" %in% names(sol))
            ) {

            sol <- info.rule("hat.pi, hat.pj -> hat.pij (IPF)", sol, verbose=verbose, indent=indent+1)
    
            reweighted <- ipf.2d(  
                    df=case$inputs$pij$data, 
                    margins_cols=sol$hat.pi, margins_rows=sol$hat.pj, 
                    max.iterations=1000, precision=1e-10, verbose=F) # verbose

            sol$hat.pij <- reweighted
            changed <- TRUE
        }

        # fi * pi -> di  
        if ( 
            ("hat.pi" %in% names(sol)) && ("hat.fi" %in% names(sol)) && (!"hat.di" %in% names(sol))
            ) {
            # print("current hat.pi")
            # print(sol$hat.pi)
            # print("current hat.fi")
            # print(sol$hat.fi)
            sol$hat.di <- nan_to_zeros(sol$hat.pi / sol$hat.fi)
            # fix potential NaNs due to di=0 by their native di counterparts (they do not matter anyway)
            # print("current hat.di")
            # print(sol$hat.di)
            

            # the degrees higher than the possible one are set to the maximum
            ids_too_high <- which(sol$hat.di > case$inputs$max.di)
            sol$hat.di[ids_too_high] <- case$inputs$max.di[ids_too_high]

            ids_too_low <- which(sol$hat.di < case$inputs$min.di)
            sol$hat.di[ids_too_low] <- case$inputs$min.di[ids_too_low]

            # now, maybe we cannot reach the expected values...
            if (length(ids_too_high) + length(ids_too_low) > 0) {
                if (case$inputs$phi.A > 0) { # case$inputs$gamma
                    # the user prefers to report the errors on fi !
                    sol <- info.rule("hat.pi, hat.fi -> hat.di [hat.fi]", sol, verbose=verbose, indent=indent+1)
                    sol$hat.fi <- normalise(nan_to_zeros(sol$hat.pi / sol$hat.di))
                    #print("hat.fi is now")
                    #print(sol$hat.fi)
                } 
                # else if (case$inputs$phi.A < case$inputs$gamma) {
                #     # the user prefers to report the errors on pi 
                #     sol <- info.rule("hat.pi, hat.fi -> hat.di [+ hat.pi + !hat.pij]", sol, verbose=verbose, indent=indent+1)
                #     sol$hat.pi <- normalise(nan_to_zeros(sol$hat.pi * sol$hat.di))
                #     sol$hat.pij <- NULL
                #     print("hat.pi")
                #     print(sol$hat.pi)
                # } 
                else {
                    stop("not enough freedom on pdi to respect both fi and pi. Try relaxing phi.A, or add more potential degrees to pdi")
                }
                
            } else {
                sol <- info.rule("hat.pi, hat.fi -> hat.di", sol, verbose=verbose, indent=indent+1)
            }

            # print("updated hat.di")
            # print(sol$hat.di)

            changed <- TRUE
        }

        if ( 
            ("hat.pj" %in% names(sol)) && ("hat.fj" %in% names(sol)) &&  (!"hat.dj" %in% names(sol))
            ) {
            sol$hat.dj <- nan_to_zeros(sol$hat.pj / sol$hat.fj)
            
            # the degrees higher than the possible one are set to the maximum
            ids_too_high <- which(sol$hat.dj > case$inputs$max.dj)
            sol$hat.dj[ids_too_high] <- case$inputs$max.dj[ids_too_high]

            ids_too_low <- which(sol$hat.dj < case$inputs$min.dj)
            sol$hat.dj[ids_too_low] <- case$inputs$min.dj[ids_too_low]

            # now, maybe we cannot reach the expected values...
            if (length(ids_too_high) + length(ids_too_low) > 0) {
                if (case$inputs$phi.B > 0) { # case$inputs$gamma
                    # the user prefers to report the errors on fi !
                    sol <- info.rule("hat.pj, hat.fj -> hat.dj [hat.fj]", sol, verbose=verbose, indent=indent+1)
                    sol$hat.fj <- normalise(nan_to_zeros(sol$hat.pj / sol$hat.dj))
                } else {
                    stop("not enough freedom on pdj to respect both fj and pj. Try relaxing phi.B, or add more potential degrees to pdj")
                }
                
            } else {
                sol <- info.rule("hat.pj, hat.fj -> hat.dj", sol, verbose=verbose, indent=indent+1)
            }

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

        # pj + dj -> fj
        # TODO this should exist, but it leads to diffult problems. If a degree is close to 0, then dividing by it leads to weird frequencies.
        # if ( 
        #     ("hat.pi" %in% names(sol)) && ("hat.di" %in% names(sol)) &&  (!"hat.fi" %in% names(sol))
        #     ) {
        #     print(sol$hat.pi)
        #     print(sol$hat.di)
        #     print(sol$hat.pi / sol$hat.di)
        #     sol$hat.fi <- normalise(sol$hat.pi / sol$hat.di)
        #     sol <- info.rule("hat.pi, hat.di -> hat.fi", sol, verbose=verbose, indent=indent+1)
        #     print(sol$hat.fi)
        #     print(sum(sol$hat.fi))
        #     changed <- TRUE
        # }
        # if ( 
        #     ("hat.pj" %in% names(sol)) && ("hat.dj" %in% names(sol)) &&  (!"hat.fj" %in% names(sol))
        #     ) {
        #     sol$hat.fj <- normalise(sol$hat.pj / sol$hat.dj)
        #     sol <- info.rule("hat.pj, hat.dj -> hat.fj", sol, verbose=verbose, indent=indent+1)
        #     print(sol$hat.fj)
        #     print(sum(sol$hat.fj))
        #     changed <- TRUE
        # }

        # pj + pij -> pij 
        # if ( 
        #     ("hat.pi" %in% names(sol)) && (!"hat.pij" %in% names(sol))
        #     ) {
        #     sol$hat.pij <- t(t(case$inputs$pij$data) / colSums(case$inputs$pij$data)) * sol$hat.pi
        #     sol <- info.rule("hat.pi -> hat.pij", sol, verbose=verbose, indent=indent+1)
        #     changed <- TRUE
        # }
        # if ( 
        #     ("hat.pj" %in% names(sol)) && (!"hat.pij" %in% names(sol))
        #     ) {
        #     sol$hat.pij <- case$inputs$pij$data / rowSums(case$inputs$pij$data) * sol$hat.pj
        #     sol <- info.rule("hat.pj -> hat.pij", sol, verbose=verbose, indent=indent+1)
        #     changed <- TRUE
        # }

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

            sol <- info.rule("hat.pdi, hat.ci, hat.ni -> hat.ndi", sol, verbose=verbose, indent=indent+1)
            sol$hat.ndi <- round_sum(t(t(sol$hat.pdi) * sol$hat.ci))
            sol$hat.ndi <- rectify.degree.counts(sol$hat.ndi, sol$hat.ni, sol$hat.ci, verbose=F)
            
            # update pdi after rounding (when divisible, else we keep the theorical version)
            indices <- which(sol$hat.ci!=0)
            sol$hat.pdi[indices] <- t(t(sol$hat.ndi[indices]) / sol$hat.ci[indices])

            # update our past knowledge according to rounding
            if ("hat.di" %in% names(sol)) {
                sol$hat.di <- colSums(sol$hat.pdi * (0:(nrow(sol$hat.pdi)-1)))
                #print("now hat.di")
                #print(sol$hat.di)
            }

            changed <- TRUE
        }
        if ( 
            all(c("hat.pdj", "hat.cj", "hat.nj") %in% names(sol)) && (!"hat.ndj" %in% names(sol))
            ) {
            
            sol <- info.rule("hat.pdj, hat.cj, hat.nj -> hat.ndj", sol, verbose=verbose, indent=indent+1)
            sol$hat.ndj <- round_sum(t(t(sol$hat.pdj) * sol$hat.cj))
            sol$hat.ndj <- rectify.degree.counts(sol$hat.ndj, sol$hat.nj, sol$hat.cj, verbose=F)   

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
                        update_degree_distribution(case$inputs$pdi$data, sol$hat.di, precision=precision.pd, verbose=F)
                    }, error = function(e) {
                        if (verbose) 
                            cat(rep("\t",times=indent+2),
                                "unable to update hat.pdi to fit hat.di=",paste(sol$hat.di,collapse=","),": ",e$message,"\n",sep="")                   
                        
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
                            cat(rep("\t",times=indent+2),"unable to update hat.pdj to fit hat.dj=",paste(sol$hat.dj,collapse=","),": ",e$message,,"\n",sep="")                   
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
#' @keywords internal
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
#' @param tolerance.pd the tolerance for probability distributions
#' @param tolerance.degree.max the tolerance for min/max average degrees
#' @param tolerance.pij the tolerance for pij probabilities
#' @return the count of problems 
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @keywords internal
#'
detect.problems <- function(sol, case, fail=TRUE, verbose=FALSE, indent=1, tolerance.pd=1.5e-8, tolerance.degree.max=1.5e-8, tolerance.pij=1e-4) {

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
    if (!is.null(sol$hat.nL) && !is.null(sol[["hat.ni", exact=TRUE]])) {
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
    if (!is.null(sol$hat.nij) && !is.null(sol[["hat.ni", exact=TRUE]])) {
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

    # colSum(hat.pij) = hat.pi 
    if (!is.null(sol$hat.pij) && !is.null(sol[["hat.pi", exact=TRUE]])) {
        problems <- problems + assert.equal(
                                        as.vector(colSums(sol$hat.pij)), as.vector(sol$hat.pi), 
                                        "hat.pi = colSums(hat.pij)", 
                                        verbose=verbose, indent=indent+1,
                                        tolerance=tolerance.pij)
    }
    if (!is.null(sol$hat.pij) && !is.null(sol$hat.pj)) {
        problems <- problems + assert.equal(
                                        as.vector(rowSums(sol$hat.pij)), as.vector(sol$hat.pj), 
                                        "hat.pj = rowSums(hat.pij)", 
                                        verbose=verbose, indent=indent+1,
                                        tolerance=tolerance.pij)
    }

    # vsum(pdi) = 1
    if (!is.null(sol$hat.pdi)) {
        for (i in 1:ncol(sol$hat.pdi)) {
            problems <- problems + assert.equal(
                                        1, 
                                        sum(sol$hat.pdi[,i]), 
                                        paste("sum(hat.pdi[",i,"])=1 (you can tune tolerance.pd)",sep=""), 
                                        verbose=verbose, indent=indent+1,
                                        tolerance=tolerance.pd)       
        }
    } 
    if (!is.null(sol$hat.pdj)) {
        for (i in 1:ncol(sol$hat.pdj)) {
            problems <- problems + assert.equal(
                                        1, 
                                        sum(sol$hat.pdj[,i]), 
                                        paste("sum(hat.pdj[",i,"])=1 (you can tune tolerance.pd)",sep=""), 
                                        verbose=verbose, indent=indent+1,
                                        tolerance=tolerance.pd)       
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
    if (!is.null(sol$hat.ndi) && !is.null(sol[["hat.ni", exact=TRUE]])) {
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
                                any(case$inputs$min.di - sol$hat.di > tolerance.degree.max), 
                                "hat.di >= min.di (try increasing tolerance.degree.max)", 
                                verbose=verbose, indent=indent+1)     
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.di - case$inputs$max.di > tolerance.degree.max), 
                                "hat.di <= max.di (try increasing tolerance.degree.max)", 
                                verbose=verbose, indent=indent+1)     
    }
    if (!is.null(sol$hat.dj)) {
        problems <- problems + assert.equal(
                                FALSE, 
                                any(case$inputs$min.dj - sol$hat.dj > tolerance.degree.max), 
                                "hat.dj >= min.dj (try increasing tolerance.degree.max)", 
                                verbose=verbose, indent=indent+1)     
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.dj - case$inputs$max.dj > tolerance.degree.max), 
                                "hat.dj <= max.dj (try increasing tolerance.degree.max)", 
                                verbose=verbose, indent=indent+1)     
    }

    if (verbose) {
        if (problems > 0) {
            #print(sol)

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
#' @keywords internal
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
#' @keywords internal
#'
paste.known <- function(sol) {
    paste(
        intersect(
            names(sol),
            c("hat.nA", "hat.ci", "hat.fi", "hat.di", "hat.pij", "hat.dj", "hat.fj", "hat.cj", "hat.nB")
            ),
        collapse=", "
        )
}

#' generates the hypothesis to explore a given solution
#'
#' Explores all the combinations of variables which are missing 
#' and might be set to user inputs.
#'
#' @param sol a current solution
#' @return a list of lists of strings
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'
#' @importFrom utils combn
#'
#' @keywords internal
#'
generate.hypothesis <- function(sol) {

    required_variables <- c("hat.nA","hat.fi",
                            "hat.di", # "hat.pdi",
                            "hat.pij",
                            "hat.dj", #"hat.pdj",
                            "hat.fj","hat.nB")

    missing_variables <- base::setdiff(required_variables, names(sol))

    # generate the combinations
    hypothesis <- unlist(lapply(1:length(required_variables), function(x) combn(required_variables,x,simplify=FALSE)), recursive=FALSE)

    hypothesis

}

#' A string representing an hypothesis
#'
#' For a given hypothesis l("hat.di","hat.fi"), 
#' returns a string like "hat.di=di, hat.fi=fi"
#'
#' @param hypothesis a list of strings representing variables, as generated by \code{\link{generate.hypothesis}}
#' @return a string
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'
#' @keywords internal
#'
name.hypothesis <- function(hypothesis) {
    var_names <- lapply(hypothesis, function (x) substr(x, 5, nchar(x)))
    paste(
                paste(
                    hypothesis,
                    var_names,
                    sep="="
                    ),
                collapse=","
            )

}

#' Asserts a given hypothesis on a given solution for a given case
#' 
#' For an hypothesis such as l("hat.fi","hat.pdi"), sets the corresponding
#' variables in the solution such as sol$hat.fi=case$fi and sol$hat.pdi=case$pdi.
#'
#' @param hypothesis a list of strings representing variables, as generated by \code{\link{generate.hypothesis}}
#' @param sol a current solution
#' @param case the case containing the initial values
#' @param verbose FALSE by default
#' @return the solution with more variables
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'
#' @importFrom utils combn
#'
#' @keywords internal
#'
assert.hypothesis <- function(hypothesis, sol, case, nA, nB, verbose=FALSE) {
    # find the initial variable names (without the prefix "hat.")
    
    sol.hyp <- sol
    for (var.hat.name in hypothesis) {

        if (var.hat.name == "hat.di") { sol.hyp[[var.hat.name]] <- case$inputs$di * case$masks$mask.di.all }
        else if (var.hat.name == "hat.dj") { sol.hyp[[var.hat.name]] <- case$inputs$dj * case$masks$mask.dj.all }
        else if (var.hat.name == "hat.pdi") { sol.hyp[[var.hat.name]] <- case$inputs$pdi$data * case$masks$mask.di.all }
        else if (var.hat.name == "hat.pdj") { sol.hyp[[var.hat.name]] <- case$inputs$pdj$data * case$masks$mask.dj.all }
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

    sol.hyp
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
#' @param tolerance.pd the tolerance for probability distributions
#' @param tolerance.degree.max the tolerance for min/max average degrees
#' @param tolerance.pij the tolerance for pij probabilities
#' @return a list of vectors (the chains) of strings 
#'
#' @seealso \code{\link{propagate.direct}} for the inference of the consequences of the hypothesis,
#'          \code{\link{detect.problems}} to ensure potential solutions are consistent
#'  
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @importFrom utils combn
#' 
#' @keywords internal
#'
resolve.missing.chain <- function(sol, chain, case, 
                                    nA, nB, 
                                    nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma, 
                                    verbose=FALSE,
                                    tolerance.pd=1.5e-8, tolerance.degree.max=1.5e-8, tolerance.pij=1e-6) {

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
                            cat(paste("\t\t\tfound an invalid solution (inference failed: ", e, ")\n", sep=""))
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
                nb.problems <- detect.problems(
                                        sol.hyp, case, 
                                        fail=FALSE, verbose=verbose, indent=3, 
                                        tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

                if (nb.problems > 0) {
                    #cat("this variable does not provides a valid solution to our problem","\n")
                    if (verbose) {
                        cat(paste("\t\t\twe found an invalid solution (",nb.problems," inconsistencies found)\n", sep=""))
                    }   
                } else {
                    if (verbose) {
                        cat("\t\t\twe found a valid solution which provides: ", paste.known(sol.hyp),"\n")
                    }
                    sol.hyp <- quantify.errors(sol.hyp, case, nA, nB)
                    if (verbose) {

                        cat(paste("\t\t\there are the NRMSE errors of this solution: ",
                            "nA=",formatC(sol.hyp$nrmse.nA, format="G"), ", fi=",formatC(sol.hyp$nrmse.fi, format="G"), 
                            ", di=",formatC(sol.hyp$nrmse.di, format="G"), ", pij=",formatC(sol.hyp$nrmse.pij, format="G"), 
                            ", dj=",formatC(sol.hyp$nrmse.dj, format="G"), ", fj=",formatC(sol.hyp$nrmse.fj, format="G"),
                            ", nB=",formatC(sol.hyp$nrmse.nB, format="G"), "\n",
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
        cumulated.errors.constrainsts <- rep(0, times=length(solutions))

        weights <- c(nu.A, phi.A, delta.A, gamma, delta.B, phi.B, nu.B)
        weights.names <- c("nu.A", "phi.A", "delta.A", "gamma", "delta.B", "phi.B", "nu.B")
        indices_weights_not_null <- which(weights > 0)
        indices_weights_null <- which(weights == 0)
        rnmse.names <- c("RNMSE.nA", "RNMSE.fi", "RNMSE.pdi", "RNMSE.pij", "RNMSE.pdj", "RNMSE.fj", "RNMSE.nB")

        if (verbose)
            cat("\t\tfound ",length(solutions)," solutions, we have to select the best according to weights", 
                    paste(weights.names, weights, collapse=", ", sep="="),"\n")

        for (i in 1:length(solutions)) {

            s <- solutions[[i]]

            errors <- c(
                        val.or.0(s$sol$nrmse.nA), 
                        val.or.0(s$sol$nrmse.fi), 
                        val.or.0(s$sol$nrmse.pdi), 
                        val.or.0(s$sol$nrmse.pij), 
                        val.or.0(s$sol$nrmse.pdj), 
                        val.or.0(s$sol$nrmse.fj), 
                        val.or.0(s$sol$nrmse.nB)
                        )

            weighted <- errors[indices_weights_not_null] / weights[indices_weights_not_null]
            cumulated.errors[i] <- sum(weighted)
            cumulated.errors.constrainsts[i] <- sum(errors[indices_weights_null])

            if (verbose)
                cat(paste(
                    "\t\t\tsolution (",i,") (", formatC(cumulated.errors[i],format="G"), ") => ", s$investigated.var.name, "\n",
                    "\t\t\t\t\t\t", paste(rnmse.names[indices_weights_not_null], formatC(errors[indices_weights_not_null],format="G"),collapse=", ",sep="="), "\n",
                    "\t\t\t\t weighted: \t",paste(rnmse.names[indices_weights_not_null], formatC(weighted,format="G"), collapse=", ",sep="="),"\n",
                    "\t\t\t\t constrainsts: \t",paste(rnmse.names[indices_weights_null], formatC(errors[indices_weights_null],format="G"), collapse=", ",sep="="),"\n",
                    sep=""))
        }

        # the best solution is either the solution which minimizes the errors on constrainsts (if there is any difference), 
        # or the solution minimizing cumulated errors on relaxed parameters.
        weight.by.constrainsts <- length(unique(cumulated.errors.constrainsts))>1
        best.solutions <-   if (weight.by.constrainsts) 
                                which(cumulated.errors.constrainsts == min(cumulated.errors.constrainsts)) 
                            else 
                                which(cumulated.errors == min(cumulated.errors)) 
                                
        best.solution <- NULL

        if (weight.by.constrainsts && (length(best.solutions) > 1)) {
            # the solutions are equivalent by constraints; why one is the best according to weights ?
            best.solutions <- which(
                    (cumulated.errors.constrainsts == min(cumulated.errors.constrainsts)) &
                    (cumulated.errors == min(cumulated.errors))) 
                
        }   
        if (length(best.solutions) > 1) {
            best.solution <- sample(best.solutions, 1)
            # TODO warning ? systematic message ?
            if (verbose) {
                cat("\t\t\tthe best solutions are solutions ",paste(best.solutions,collapse=",")," based on hypothesis:\n",sep="")
                for (idx in best.solutions)
                    cat("\t\t\t\t* solution ",idx," solved with:\t",solutions[[idx]]$investigated.var.name,"\n",sep="")
                #cat("\t\t\tthere are multiple best solutions (that is: ",length(best.solutions),") ; just selecting one randomly\n",sep="")
                cat("\t\t=> selected ",best.solution, " which is one of the ",length(best.solutions)," solutions having ", sep="")
                if (weight.by.constrainsts)
                    cat("lowest error for constrainsts\n",sep="")
                else 
                    cat("lowest weighted error\n",sep="")
            }
        } else {
            best.solution <- best.solutions[[1]]
            if (verbose) {
                cat("\t\t=> will keep solution ",best.solution," which is the one with the ",sep="")
                if (weight.by.constrainsts)
                    cat("lowest error for constrainsts\n",sep="")
                else 
                    cat("lowest weighted error\n",sep="")
            }

        }
        if (verbose) {
            cat("\t\t\tthis solution is based on hypothesis: ",solutions[[best.solution]]$investigated.var.name,"\n",sep="")
        }
        res <- solutions[[best.solution]]$sol
        #print(names(solutions[[best.solution]]$sol))
    }

    res
}

#' Selects one best solution among a list of possible solutions
#'
#' Relying on the weigthed error of each solution, selects one of the best ones.
#' 
#' @param solutions a list of solutions such as generated inside \code{\link{resolve}}
#' @inheritParams resolve
#' @return a solution solved
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
#' @keywords internal
#'
best.solution <- function(solutions,
                        nu.A, phi.A, delta.A, 
                        gamma,
                        nu.B, phi.B, delta.B,  
                        verbose=FALSE) {

    # no solution ? => no solution !
    if (length(solutions) == 0) {
        if (verbose) 
            cat("\t\t\tfound no solution for this chain :-(\n")
        stop("no solution found; the case seems overconstrained. You might try to relax relaxation parameters, or relax.zeros your input distributions, or relax the tolerance margins. Keep trying!")
    } 

    # one solution
    if (length(solutions) == 1) {
        # return the unique solution
        if (verbose)
            cat("\t\t\tfound one valid solution\n")
        return(solutions[[1]]$sol)
    }

    # several solutions
    weights <- c(nu.A, phi.A, delta.A, gamma, delta.B, phi.B, nu.B)
    weights.names <- c("nu.A", "phi.A", "delta.A", "gamma", "delta.B", "phi.B", "nu.B")
    indices_weights_not_null <- which(weights > 0)
    indices_weights_null <- which(weights == 0)

    rnmse.names <- c("RNMSE.nA", "RNMSE.fi", "RNMSE.pdi", "RNMSE.pij", "RNMSE.pdj", "RNMSE.fj", "RNMSE.nB")

    if (verbose)
        cat("\t\tfound ",length(solutions)," solutions, we have to select the best according to weights", 
                paste(weights.names, weights, collapse=", ", sep="="),"\n")

    # mesure the weighted error for each solution
    cumulated.errors <- rep(0, times=length(solutions))
    for (i in 1:length(solutions)) {

        s <- solutions[[i]]

        rnmse <- c(s$sol$nrmse.nA, s$sol$nrmse.fi, s$sol$nrmse.pdi, s$sol$nrmse.pij, s$sol$nrmse.pdj, s$sol$nrmse.fj, s$sol$nrmse.nB)
        # The weighted sum is made of the errors of the relaxed parameters divided by relaxation (the more relaxed, the less important the error is).
        # We add the errors on parameters supposed to be relaxed times 10000, that is they count a lot !
        cumulated.errors[i] <- sum(rnmse[indices_weights_not_null] / weights[indices_weights_not_null]) + 
                    sum(rnmse[indices_weights_null] * 10000)

        

        if (verbose)
            cat(paste(
                "\t\t\tsolution (",i,") (", formatC(cumulated.errors[i],format="G"), ") => ", s$hypothesis, "\n",
                "\t\t\t\t errors: \t ", paste(rnmse.names, formatC(rnmse,format="G"),collapse=", ",sep="="), "\n",
                #"\t\t\t\t weighted: \t",paste(rnmse.names[indices_weights_not_null], formatC(weighted,format="G"), collapse=", ",sep="="),"\n",
                #rnmse[indices_weights_null]
                sep=""))
    }

    # select the solution(s) having the lowest error
    best.solutions <- which(cumulated.errors == min(cumulated.errors)) 
                            
    best <- NULL

    if (length(best.solutions) == 1) {
        if (verbose) 
            cat("\t\t=> will keep solution ",best," which is the one with the lowest weighted error\n",sep="")
        best <- best.solutions[[1]]
    } else {
        best <- sample(best.solutions, 1)
        # TODO warning ? systematic message ?
        if (verbose) {
            cat("\t\t\tthe best solutions are solutions ",paste(best.solutions,collapse=",")," based on hypothesis:\n",sep="")
            for (idx in best.solutions)
                cat("\t\t\t\t* solution ",idx," solved with:\t",solutions[[idx]]$hypothesis,"\n",sep="")
            #cat("\t\t\tthere are multiple best solutions (that is: ",length(best.solutions),") ; just selecting one randomly\n",sep="")
            cat("\t\t=> randomly selected ",best, " which is one of the ",length(best.solutions)," solutions having lowest weighted error\n", sep="")
        }
    }

    if (verbose) {
        cat("\t\t\tthis solution is based on hypothesis: ",solutions[[best]]$hypothesis,"\n",sep="")
    }
    
    solutions[[best]]$sol

}

#' Resolves a partial solution by inference and hypothesis
#' 
#' Resolves a given case starting with the base solution
#' according to user parameters. 
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
#' @param tolerance.pd the tolerance for probability distributions
#' @param tolerance.degree.max the tolerance for min/max average degrees
#' @param tolerance.pij the tolerance for pij probabilities
#' 
#' @return a list of vectors (the chains) of strings 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
#' @keywords internal
#'
resolve <- function(sol, case, 
                        nA, nB, 
                        nu.A, phi.A, delta.A, 
                        gamma,
                        nu.B, phi.B, delta.B,  
                        verbose=FALSE,
                        tolerance.pd=1.5e-8, tolerance.degree.max=1.5e-8, tolerance.pij=1e-6) {

    # direct propagation of what we know
    sol.infered <- propagate.direct(sol, case, verbose=verbose)
    detect.problems(sol.infered, case, 
                    verbose=verbose, 
                    tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

    solutions <- list()
    hypothesis <- generate.hypothesis(sol.infered)
    blacklist <- list() # nothing blacklisted so far
    for (h in hypothesis) {

        # skip blacklisted hypothesis which are known not to work
        skipped <- F
        for (b in blacklist) {
            if (all(b %in% h)) {
                skipped <- T
                if (verbose) {
                    cat("\t\tskipping ", name.hypothesis(h), " because ", name.hypothesis(b), "failed in the past\n", sep="")
                }
                break
            }
        }
        if (skipped & verbose) {
            #cat("\t\tskipping ", name.hypothesis(h), "\n", sep="")
            next
        }

        # the current explored solution, which will host the error rates as well
        explored <- list()
        explored$hypothesis <- name.hypothesis(h)

        # add the hypothesis to the solution
        if (verbose) {
            cat("\t\ttrying with hypothesis: ", name.hypothesis(h), "\n", sep="")
        }
        sol.hyp <- assert.hypothesis(h, sol.infered, case, nA, nB, verbose=verbose)

        # inference on this basis
        sol.hyp <- tryCatch({
                        propagate.direct(sol.hyp, case, verbose=verbose, indent=3)
                    }, error = function(e) {
                        if (verbose) {
                            cat(paste("\t\t\tfound an invalid solution (inference failed: ", e, ")\n", sep=""))
                        }
                        NULL
                    })
        
        if (is.null(sol.hyp)) {
            # we have no solution
            # do not try more complex combinations
            blacklist[[length(blacklist)+1]] <- h
            # continue
            next
        }

        if (!all(c("hat.nA","hat.ci","hat.ndi","hat.ni","hat.nij","hat.ni","hat.ndj","hat.cj","hat.nB") %in% names(sol.hyp))) {
            # we cannot infer everything from this hypothesis; let's ignore it
            if (verbose)
                cat("\t\t\tthese hypothesis do not enable the resolution of all the required variables.\n")
            # do not explore this solution more
            next
        }
        
         # do we achieve to solve the problem on this basis ?
        nb.problems <- detect.problems(
                                sol.hyp, case, 
                                fail=FALSE, verbose=verbose, indent=3, 
                                tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

        if (nb.problems > 0) {
            if (verbose) 
                cat(paste("\t\t\twe found an invalid solution (",nb.problems," inconsistencies found)\n", sep=""))
            # do not try more complex combinations
            blacklist[[length(blacklist)+1]] <- h
            # stop there
            next
        } 

        if (verbose)
            cat("\t\t\twe found a valid solution which provides: ", paste.known(sol.hyp),"\n")
        sol.hyp <- quantify.errors(sol.hyp, case, nA, nB)
        if (verbose) {
            cat(paste("\t\t\there are the NRMSE errors of this solution: ",
                "nA=",formatC(sol.hyp$nrmse.nA, format="G"), ", fi=",formatC(sol.hyp$nrmse.fi, format="G"), 
                ", di=",formatC(sol.hyp$nrmse.di, format="G"), ", pij=",formatC(sol.hyp$nrmse.pij, format="G"), 
                ", dj=",formatC(sol.hyp$nrmse.dj, format="G"), ", fj=",formatC(sol.hyp$nrmse.fj, format="G"),
                ", nB=",formatC(sol.hyp$nrmse.nB, format="G"), "\n",
                sep=""
                ))
        }
        explored$sol <- sol.hyp
        solutions[[length(solutions)+1]] <- explored

    }

    best.solution(solutions, nu.A, phi.A, delta.A, 
                        gamma,
                        nu.B, phi.B, delta.B,  
                        verbose=verbose)


    # missing.chains.orig <- detect.missing.chains(sol.tmp)

    # permutations <- function(l) {
    #     if (length(l) < 2) return (list(l))
    #     else if (length(l) == 2) return (list(list(l[0],l[1]), list(l[1],l[0]) )
    # }
    # for (missing.chains in length(missing.chains.orig)) {
    #     # will try by starting with the ith element
    #     # TODO TODO TODO 
    #     missing.chains <- missing.chains.orig
    #     while (length(missing.chains) > 0) {

    #         chain <- missing.chains[[1]]
    #         missing.chains[[1]] <- NULL

    #         found <- resolve.missing.chain(sol.tmp, chain, case, 
    #                     nA=nA, nB=nB, nu.A=nu.A, phi.A=phi.A, delta.A=delta.A, nu.B=nu.B, phi.B=phi.B, delta.B=delta.B, gamma=gamma, 
    #                     verbose=verbose,
    #                     tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

    #         if (!is.null(found)) {
    #             if (verbose) {
    #                 cat("\tthis missing chain was solved\n", sep="")
    #             }
    #             sol.tmp <- found 
    #             missing.chains <- detect.missing.chains(sol.tmp)
    #         }

    #     }

    #     if (!is.null(sol.tmp)) return(sol.tmp)
    # }


    # # print("propagated")
    # # print(sol)

    # return(sol.tmp)
}

#' Quantifies the goodness of fit between two dataframes
#' 
#' Computes the squared of each, then the difference between them,
#' then elevated to power 2, then summed and mult 4.
#' 
#' @param df1 data ref
#' @param df2 data gen 
#' @return a scalar 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @export
#'
gof.freeman_tukey <- function(df1, df2) {

    4*sum((sqrt(df1) - sqrt(df2))^2)

}

#' Quantifies the Fisher goodness of fit between two dataframes
#' 
#' Computes fisher test using \code{\link{fisher.test}};
#' in case of error, returns NA
#' 
#' @param x data ref
#' @param y data gen 
#' @return a scalar or NA if not possible
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @importFrom stats fisher.test
#' 
#' @export
#'
fisher.pvalue.or.NA <- function(x, y) {
    tryCatch({
        fisher.test(x, y, simulate.p.value=T)$p.value
    }, error = function(e) {
        NA
    })
}

#' Quantifies the chi2 goodness of fit between two dataframes
#' 
#' Computes fisher test using \code{\link{chisq.test}};
#' in case of error, returns NA
#' 
#' @param x data ref
#' @param y data gen 
#' @return a scalar or NA if not possible
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#'  
#' @importFrom stats chisq.test
#' 
#' @export
#'
chisq.test.or.NA  <- function(x, y) {
    tryCatch({
        chisq.test(x, B=y)$p.value
    }, error = function(e) {
        NA
    })
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
#' @param compute.fisher if TRUE, the compute the Fisher p value (slow)
#' @param compute.chi2 if TRUE, the compute the Chi2 
#' @return the solution with additional variables for errors 
#' 
#' @keywords internal
#'
quantify.errors <- function(sol, case, nA, nB, compute.fisher=F, compute.chi2=F) {

    # measure errors
    # MSE fi
    if (!is.null(sol$hat.fi)) {
        sol$mse.fi <- mean( ( sol$hat.fi - case$stats$fi)^2  )
        sol$rmse.fi <- sqrt(sol$mse.fi)
        sol$nrmse.fi <- sol$rmse.fi
        sol$ft.fi <- gof.freeman_tukey(sol$hat.fi, case$stats$fi)
    }
    if (!is.null(sol$hat.fj)) {
        sol$mse.fj <- mean( ( sol$hat.fj - case$stats$fj)^2  )
        sol$rmse.fj <- sqrt(sol$mse.fj)
        sol$nrmse.fj <- sol$rmse.fj
        sol$ft.fj <- gof.freeman_tukey(sol$hat.fj, case$stats$fj)
    }
    # Pearson ci, cj
    if (!is.null(sol$hat.ci)) {
        sol$xi2.p.fi <- chisq.test.or.NA(sol$hat.ci, case$stats$fi)
        if (compute.fisher) sol$fisher.p.fi <- fisher.pvalue.or.NA(sol$hat.ci, case$stats$fi)
    }
    if (!is.null(sol$hat.cj)) {
        sol$xi2.p.fj <- chisq.test.or.NA(sol$hat.cj, case$stats$fj)
        if (compute.fisher) sol$fisher.p.fj <- fisher.pvalue.or.NA(sol$hat.cj, case$stats$fj)
    }

    # MSE di
    if (!is.null(sol$hat.di)) {
        sol$mse.di <- mean( ( sol$hat.di - case$inputs$di)^2  )
        sol$rmse.di <- sqrt(sol$mse.di)
        sol$nrmse.di <- sol$rmse.di / max(case$inputs$di)
        # non sense? sol$ft.di <- gof.freeman_tukey(sol$hat.di, case$inputs$di)
    }
    if (!is.null(sol$hat.dj)) {
        sol$mse.dj <- mean( ( sol$hat.dj - case$inputs$dj)^2  )
        sol$rmse.dj <- sqrt(sol$mse.dj) / max(case$inputs$dj)
        sol$nrmse.dj <- sol$rmse.dj
        # non sense? sol$ft.dj <- gof.freeman_tukey(sol$hat.dj, case$inputs$dj)
    }

    # MSE pdi 
    if (!is.null(sol$hat.pdi)) {
        sol$mse.pdi <- mean( ( sol$hat.pdi - case$inputs$pdi$data)^2  )
        sol$rmse.pdi <- sqrt(sol$mse.pdi)
        sol$nrmse.pdi <- sol$rmse.pdi # / ncol(case$inputs$pdi$data)
        sol$ft.pdi <- gof.freeman_tukey(sol$hat.pdi, case$inputs$pdi$data)
    }
    if (!is.null(sol$hat.pdj)) {
        sol$mse.pdj <- mean( ( sol$hat.pdj - case$inputs$pdj$data)^2  )
        sol$rmse.pdj <- sqrt(sol$mse.pdj)
        sol$nrmse.pdj <- sol$rmse.pdj # / ncol(case$inputs$pdj$data)
        sol$ft.pdj <- gof.freeman_tukey(sol$hat.pdj, case$inputs$pdj$data)
    }
    # Pearson ndi, ndj
    if (!is.null(sol$hat.ndi)) {
        sol$xi2.p.pdi <- chisq.test.or.NA(sol$hat.ndi, case$inputs$pdi_fixed$data)
        if (compute.fisher) sol$fisher.p.pdi <- fisher.pvalue.or.NA(sol$hat.ndi, case$inputs$pdi_fixed$data)
    }
    if (!is.null(sol$hat.ndj)) {
        sol$xi2.p.pdj <- chisq.test.or.NA(sol$hat.ndj, case$inputs$pdi_fixed$data)
        if (compute.fisher) sol$fisher.p.pdj <- fisher.pvalue.or.NA(sol$hat.ndj, case$inputs$pdj_fixed$data)
    }

    # MSE pij
    if (!is.null(sol$hat.pij)) {
        # if (verbose) {
        #     print("computing the NRMSE of pij, for hat.pij=")
        #     print(sol$hat.pij)
        #     print("original being")
        #     print(case$inputs$pij$data)
        #     print("difference is ")
        #     print(sol$hat.pij - case$inputs$pij$data)
        #     print("square is")
        #     print(( sol$hat.pij - case$inputs$pij$data )^2)
        #     print("so mean is ")
        #     print(mean.data.frame( ( sol$hat.pij - case$inputs$pij$data )^2  ))
        #     print("and thus")
        #     print(sqrt(mean.data.frame( ( sol$hat.pij - case$inputs$pij$data )^2  )))
        # }
        sol$mse.pij <- mean.data.frame( ( sol$hat.pij - case$inputs$pij$data )^2  )
        sol$rmse.pij <- sqrt(sol$mse.pij)
        sol$nrmse.pij <- sol$rmse.pij
        sol$ft.pij <- gof.freeman_tukey(sol$hat.pij, case$inputs$pij$data)
    }
    # Pearson ndi, ndj
    if (!is.null(sol$hat.nij)) {
        sol$xi2.p.pij <- chisq.test.or.NA(sol$hat.nij, case$inputs$pij$data)
        if (compute.fisher) sol$fisher.p.pij <- fisher.pvalue.or.NA(sol$hat.nij, case$inputs$pij_fixed$data)
    }
    
    # error n
    if (!is.null(sol$hat.nA)) {
        sol$error.nA <- abs( sol$hat.nA - nA )
        sol$rmse.nA <- sol$error.nA
        sol$nrmse.nA <- abs(sol$hat.nA - nA)/nA
    }
    if (!is.null(sol$hat.nB)) {
        sol$error.nB <- abs( sol$hat.nB - nB )
        sol$rmse.nB <- sol$error.nB
        sol$nrmse.nB <- abs(sol$hat.nB - nB)/nB
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
#' @keywords internal
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

#' Computes the mean of all the values inside a dataframe 
#' 
#' Sometimes the standard mean() function works on matrices,
#' but it does not on dataframes. In case the parameter is a dataframe,
#' the mean will be computed as the sum of all the values divided by the 
#' count of elements. Else another standard mean will be called.
#' 
#' @param x the dataframe 
#' @param ... additional parameters are quietly ignored
#' @return a double value 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
mean.data.frame <- function(x, ...) {

    if (class(x) != "data.frame") 
        mean(x)
    else 
        sum(x) / (nrow(x) * ncol(x))
}

#' Solves the equations by arbitrating.
#' 
#' Solves the underlying equations based on data, the inputs parameters, 
#' and the control parameters. It tries to distribute the errors in the places 
#' accepted by the user. 
#'
#' This function is deterministic as long as only one solution is found. 
#' In case several solutions with equivalent quality are found, one of them 
#' is randomly chosen.
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
#' @param tolerance.pd the tolerance for probability distributions
#' @param tolerance.degree.max the tolerance for min/max average degrees
#' @param tolerance.pij the tolerance for pij probabilities
#' 
#' @return a case ready for generation
#'
#' @export
#'
#' @seealso \code{\link{matching.prepare}} to prepare a case for this function, 
#'          \code{\link{matching.generate}} to use the result for actual generation, 
#'          \code{\link{plot.dpp_resolved}} to plot the quality of the solution
#' 
#' @examples
#'
#' # load sample data
#' data(dwellings_households)
#' # prepare the case  
#' case.prepared <- matching.prepare(
#'                      dwellings_households$sample.A, dwellings_households$sample.B, 
#'                      dwellings_households$pdi, dwellings_households$pdj, 
#'                      dwellings_households$pij)
#' # resolve tbe case
#' solved <- matching.solve(case.prepared, 
#'                      nA=50000,nB=40000, 
#'                      nu.A=0, phi.A=0, delta.A=1, 
#'                      gamma=1, 
#'                      delta.B=0, phi.B=0, nu.B=0)
#' # print the resolution information
#' print(solved)
#' # access the solved frequencies, distribution of degrees and matching probabilities
#' print(solved$gen$hat.fi)
#' print(solved$gen$hat.pdi)
#' print(solved$gen$hat.pij)
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.solve <- function(case, 
                                nA, nB, 
                                nu.A, phi.A, delta.A, 
                                gamma, 
                                delta.B, phi.B, nu.B, 
                                verbose = FALSE,
                                tolerance.pd = 1.5e-6, tolerance.degree.max=1.5e-8, tolerance.pij=1e-6) {

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

    case$inputs$nA <- nA
    case$inputs$nB <- nB
    case$inputs$nu.A <- nu.A 
    case$inputs$nu.B <- nu.B 
    case$inputs$phi.A <- phi.A 
    case$inputs$phi.B <- phi.B
    case$inputs$delta.A <- delta.A 
    case$inputs$delta.B <- delta.B
    case$inputs$gamma <- gamma 
    case$inputs$tolerance.pd <- tolerance.pd
    case$inputs$tolerance.degree.max <- tolerance.degree.max

    #print(case$masks)

    if (verbose) 
        cat("\nstarting the resolution of the case\n")
    

    sol <- list()
    if (nu.A == 0)      { sol$hat.nA <- nA }
    if (phi.A == 0)     { sol$hat.fi <- normalise(case$stats$fi * case$masks$mask.fi.all) }
    if (delta.A == 0)   { 
        sol$hat.di <- case$inputs$di * case$masks$mask.di.all 
        sol$hat.pdi <- case$inputs$pdi$data #* case$masks$mask.di.all 
    }
        
    if (nu.B == 0)      { sol$hat.nB <- nB }
    if (phi.B == 0)     { sol$hat.fj <- normalise(case$stats$fj * case$masks$mask.fj.all) }
    if (delta.B == 0)   { 
        sol$hat.dj <- case$inputs$dj * case$masks$mask.dj.all 
        sol$hat.pdj <- case$inputs$pdj$data
    }

    if (gamma == 0)     { sol$hat.pij <- normalise(case$inputs$pij$data * case$masks$mask.pij.all) }

    if (verbose) {
        cat("\taccording to user weights, we already know: ",paste(names(sol),collapse=","),"\n",sep="")

    }

    # detect issues right now; maybe the problem is overconstrained or badly constrainted
    detect.problems(sol, case, 
                    verbose=verbose, 
                    tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

    sol <- resolve(sol, case, 
                nA=nA, nB=nB, 
                nu.A=nu.A, phi.A=phi.A, delta.A=delta.A, nu.B=nu.B, phi.B=phi.B, delta.B=delta.B, gamma=gamma, 
                verbose=verbose, 
                tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

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

    detect.problems(sol, case, 
                    verbose=verbose, 
                    tolerance.pd=tolerance.pd, tolerance.degree.max=tolerance.degree.max, tolerance.pij=tolerance.pij)

    if (verbose) {
        cat("\ncase solved.\n")
    }
    # measure errors
    sol <- quantify.errors(sol, case, nA, nB)

    # ensure the form is ok
    sol <- ensure.form(sol, case)

    # prepare the result
    res <- case
    
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

