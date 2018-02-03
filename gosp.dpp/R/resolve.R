

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
#' This process is done column by column by calling @link{update_degree_distribution.col}.
#' '
#' @param pdx the distribution of degrees 
#' @param dx the vector of the expected average degrees
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


#' Adapts a distribution of degrees so that each column times n sums to expected nn
#' (basically solves rounding issues for this specific case)
#'
#' @param pdn the contigencies of degrees
#' @param nn the expected totals that should be reached 
rectify.degree.counts <- function(pdn, nn) {

    for (col in 1:ncol(pdn)) {
        total.current <- sum(pdn[,col]*0:(nrow(pdn)-1))
        total.expected <- nn[col]

        #print(pdn)
        #print(nn)

        if (total.current > total.expected) {
            to.remove <- total.current - total.expected

            if (to.remove > 1) {
                warning("/!\\ the case for fixing rounding errors below 1 slot (here:",to.remove,") is not managed well\n.",sep="")
            }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if ( (pdn[2,col] >= to.remove) && (pdn[1,col] > 0)) {
                pdn[2,col] <- pdn[2,col] - to.remove
                pdn[1,col] <- pdn[1,col] + to.remove
            } else {
                warning("/!\\ found no good solution to fix this rounding error of an additional R\n") 
            }

        } else if (total.current < total.expected) {
            to.add <- total.expected - total.current 

            if (to.add > 1) {
                warning("/!\\ the case for fixing rounding errors above 1 slot (here:",to.add,") is not managed well\n", sep="")
            }

            # we are creating more slots than expected
            # ... we have to remove it somewhere 
            if ( (pdn[1,col] >= to.add) && (pdn[2,col] > 0)) {
                pdn[1,col] <- pdn[1,col] - to.add
                pdn[2,col] <- pdn[2,col] + to.add
            } else {
                warning("/!\\ found no good solution to fix this rounding error of an additional R\n") 
            }
        }

    }

    pdn
}

normalise <- function(df) {
    df / sum(df)
}


propagate.direct <- function(sol,case) {


    info.rule <- function(name, sol) {
        cat("\t",name,"\n",sep="")
        sol$inference[[length(sol$inference)+1]] <- name
        sol
    } 

    if (is.null(sol$inference)) {
        sol$inference <- list()
    }
    # TODO mask
    # TODO normalize

    # TODO should detect when info are not coherent (case if too many zeros !)
    while (TRUE) {

        changed <- FALSE


        # forward nA + fi -> ci 
        if ( 
            ("hat.nA" %in% names(sol)) && ("hat.fi" %in% names(sol)) && (!"hat.ci" %in% names(sol))
            ) {
            sol$hat.ci <- round_sum(sol$hat.nA * sol$hat.fi)
            sol$hat.fi <- sol$hat.ci / sol$hat.nA
            sol <- info.rule("hat.nA, hat.fi -> hat.ci", sol)
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
            sol <- info.rule("hat.nB, hat.fj -> hat.cj", sol)
            #print("sol$hat.cj <- round(sol$hat.nB * sol$hat.fj)")
            changed <- TRUE
        }


        # forward ci + di -> ni 
        if ( 
            ("hat.ci" %in% names(sol)) && ("hat.di" %in% names(sol)) && (!"hat.ni" %in% names(sol))
            ) {
            sol$hat.ni <- round_sum(sol$hat.ci * sol$hat.di)
            # adapt rounding
            sol$hat.di <- sol$hat.ni / sol$hat.ci 
            #print("sol$hat.ni <- round(sol$hat.ci * sol$hat.di)")
            sol <- info.rule("hat.ci, hat.di -> hat.ni", sol)
            #print(sol$hat.ni)
            changed <- TRUE
        }

         # forward cj + dj -> nj 
        if ( 
            ("hat.cj" %in% names(sol)) && ("hat.dj" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- round_sum(sol$hat.cj * sol$hat.dj)
            sol$hat.dj <- sol$hat.nj / sol$hat.cj
            sol <- info.rule("hat.cj, hat.dj -> hat.nj", sol)
            changed <- TRUE
        }

        if ( 
            ("hat.ni" %in% names(sol)) && (!"hat.nL" %in% names(sol))
            ) {
            sol$hat.nL <- sum(sol$hat.ni)
            sol <- info.rule("hat.ni -> hat.nL", sol)
            changed <- TRUE
        }

        if ( 
            ("hat.nj" %in% names(sol)) && (!"hat.nL" %in% names(sol))
            ) {
            sol$hat.nL <- sum(sol$hat.nj)
            sol <- info.rule("hat.nj -> hat.nL", sol)
            changed <- TRUE
        }

        # forward sum(ci) -> nA
        if ( 
            ("hat.ci" %in% names(sol)) && (!"hat.nA" %in% names(sol))
            ) {
            sol$hat.nA <- sum(sol$hat.ci)
            sol <- info.rule("hat.ci -> hat.nA", sol)
            changed <- TRUE
        }
        if ( 
            ("hat.cj" %in% names(sol)) && (!"hat.nB" %in% names(sol))
            ) {
            sol$hat.nB <- sum(sol$hat.cj)
            sol <- info.rule("hat.cj -> hat.nB", sol)
            changed <- TRUE
        }

        # TODO probably not true !!!
        if ( 
            ("hat.ni" %in% names(sol)) #&& ("hat.pij" %in% names(sol)) 
            && (!"hat.nij" %in% names(sol))
            ) {
            sol$hat.nij <- t(sol$hat.ni * t(case$inputs$pij$data) / colSums(case$inputs$pij$data))
            # round per column (to keep the sums consistent)
            for (i in 1:ncol(sol$hat.nij)) {
                sol$hat.nij[,i] <- round_sum(sol$hat.nij[,i])
            }
            # adapt th based on rounded values
            sol$hat.pij <- sol$hat.nij / sum(sol$hat.nij)
            # ... and propagate 

            sol <- info.rule("hat.ni -> hat.nij, hat.pij", sol)

            #print("sol$hat.nij <- round( sol$hat.ni * case$inputs$pij$data / colSums(case$inputs$pij$data) )")
            #            print(sol$hat.nij)

            # print("computed hat.nij from ni and pij")
            # print(sol$hat.nij/sum(sol$hat.nij))
            changed <- TRUE
        }
        if ( 
            ("hat.nj" %in% names(sol)) #&& ("hat.pij" %in% names(sol)) 
            && (!"hat.nij" %in% names(sol))
            ) {
            sol$hat.nij <- sol$hat.nj * case$inputs$pij$data / rowSums(case$inputs$pij$data)
            # round per line (to keep the sums consistent)
            for (i in 1:nrow(sol$hat.nij)) {
                sol$hat.nij[i,] <- round_sum(sol$hat.nij[i,])
            }
            sol$hat.pij <- sol$hat.nij / sum(sol$hat.nij)
            sol <- info.rule("hat.nj -> hat.nij, hat.pij", sol)
            #print("sol$hat.nij <- round(sol$hat.nj * case$inputs$pij$data / rowSums(case$inputs$pij$data))")
            #print("computed hat.nij from nj and pij")
            #print(sol$hat.nij/sum(sol$hat.nij))
            changed <- TRUE
        }


        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.ni" %in% names(sol))
            ) {
            sol$hat.ni <- colSums(sol$hat.nij)
            #print("sol$hat.ni <- colSums(sol$hat.nij)")
            sol <- info.rule("hat.nij -> hat.ni", sol)

            changed <- TRUE
        }
        if ( 
            ("hat.nij" %in% names(sol)) && (!"hat.nj" %in% names(sol))
            ) {
            sol$hat.nj <- rowSums(sol$hat.nij)
            sol <- info.rule("hat.nij -> hat.nj", sol)
            changed <- TRUE
        }
 
        # forward ni + di -> ci 
        if (
            ("hat.ni" %in% names(sol)) && ("hat.di" %in% names(sol)) && (!"hat.ci" %in% names(sol))
            ) {
            sol$hat.ci <- round_sum(sol$hat.ni / sol$hat.di )
            sol$hat.di <- sol$hat.ni / sol$hat.ci
            sol <- info.rule("hat.ni, hat.di -> hat.ci", sol)
            changed <- TRUE

        }
        if (
            ("hat.nj" %in% names(sol)) && ("hat.dj" %in% names(sol)) && (!"hat.cj" %in% names(sol))
            ) {
            sol$hat.cj <- round_sum(sol$hat.nj / sol$hat.dj )
            sol$hat.dj <- sol$hat.nj / sol$hat.cj
            sol <- info.rule("hat.nj, hat.dj -> hat.cj", sol)            
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
                                    sol$hat.ni / sol$hat.ci
                                    )
                            )
            sol <- info.rule("hat.ni, hat.ci -> hat.di", sol)
            changed <- TRUE
        }
        if ( 
            ("hat.nj" %in% names(sol)) && ("hat.cj" %in% names(sol)) &&  (!"hat.dj" %in% names(sol))
            ) {
            sol$hat.dj <- pmax(
                            case$inputs$min.dj, 
                            pmin(
                                    case$inputs$max.dj,
                                    sol$hat.nj / sol$hat.cj
                                    )
                            )
            sol <- info.rule("hat.nj, hat.cj -> hat.dj", sol)
            changed <- TRUE
        }

        # pij -> pi 
        if ( 
            ("hat.pij" %in% names(sol)) && (!"hat.pi" %in% names(sol))
            ) {
            sol$hat.pi <- colSums(sol$hat.pij)
            sol <- info.rule("hat.pij -> hat.pi", sol)
            #print("sol$hat.pi <- colSums(sol$hat.pij)")
            changed <- TRUE
        }
        if ( 
            ("hat.pij" %in% names(sol)) && (!"hat.pj" %in% names(sol))
            ) {
            sol$hat.pj <- rowSums(sol$hat.pij)
            sol <- info.rule("hat.pij -> hat.pj", sol)
            #print("sol$hat.pj <- rowSums(sol$hat.pij)")
            changed <- TRUE
        }

        # fi * pi -> di  
        # WTF ???? is that even valid ???
        if ( 
            ("hat.pi" %in% names(sol)) && ("hat.fi" %in% names(sol)) && (!"hat.di" %in% names(sol))
            #&& (!"hat.ci" %in% names(sol)
            ) {
            # in factor we might change it of a factor, only the ratio matters
            sol$hat.di <- sol$hat.pi / sol$hat.fi
            
            sol <- info.rule("hat.pi, hat.fi -> hat.di", sol)
            
            #print("sol$hat.di <- x * sol$hat.pi / sol$hat.fi")
            # 
            factor <- 1
            if (any(sol$hat.di > case$inputs$max.di)) {
                # we are higher than what is allowed
                factor <- min(case$inputs$max.di / sol$hat.di)
            }  else if (any(sol$hat.di < case$inputs$min.di)) {
                factor <- min(case$inputs$min.di / sol$hat.di) # TODO verifier
            }
            sol$hat.di <- sol$hat.di * factor

            # print("et donc")
            #  print(sol$hat.di)
            #  print(normalise(sol$hat.di * sol$hat.fi))
            #  print(sol$hat.pi)
            changed <- TRUE
        }

        if ( 
            ("hat.pj" %in% names(sol)) && ("hat.fj" %in% names(sol)) &&  (!"hat.dj" %in% names(sol))
            ) {
            # in factor we might change it of a factor, only the ratio matters
            sol$hat.dj <- sol$hat.pj / sol$hat.fj
            
            sol <- info.rule("hat.pj, hat.fj -> hat.dj", sol)
            
            #print("sol$hat.di <- x * sol$hat.pi / sol$hat.fi")
            # 
            factor <- 1
            if (any(sol$hat.dj > case$inputs$max.dj)) {
                # we are higher than what is allowed
                factor <- min(case$inputs$max.dj / sol$hat.dj)
            }  else if (any(sol$hat.di < case$inputs$min.di)) {
                factor <- min(case$inputs$min.dj / sol$hat.dj) # TODO verifier
            }
            sol$hat.dj <- sol$hat.dj * factor

            # print("et donc")
            #  print(sol$hat.di)
            #  print(normalise(sol$hat.di * sol$hat.fi))
            #  print(sol$hat.pi)
            changed <- TRUE
        }
        # TODO j !!!!

        # fj + dj -> pj
        if ( 
            ("hat.fi" %in% names(sol)) && ("hat.di" %in% names(sol)) &&  (!"hat.pi" %in% names(sol))
            ) {
            sol$hat.pi <- normalise(sol$hat.fi * sol$hat.di)
            sol <- info.rule("hat.fi, hat.di -> hat.pi", sol)
            #print("sol$hat.pi <- normalise(sol$hat.fi * sol$hat.di)")
            changed <- TRUE
        }
        if ( 
            ("hat.fj" %in% names(sol)) && ("hat.dj" %in% names(sol)) &&  (!"hat.pj" %in% names(sol))
            ) {
            sol$hat.pj <- normalise(sol$hat.fj * sol$hat.dj)
            #print("sol$hat.pj <- normalise(sol$hat.fj * sol$hat.dj)")
            sol <- info.rule("hat.fj, hat.dj -> hat.pj", sol)
            changed <- TRUE
        }

        # pj + pij -> pij 
        if ( 
            ("hat.pi" %in% names(sol)) && (!"hat.pij" %in% names(sol))
            ) {
            sol$hat.pij <- t(t(case$inputs$pij$data) / colSums(case$inputs$pij$data)) * sol$hat.pi
            # print("sol$hat.pij <- case$inputs$pij$data / rowSums(case$inputs$pij$data) * sol$hat.pj")
            # print(case$inputs$pij$data)
            # print(sol$hat.pj)
            sol <- info.rule("hat.pi -> hat.pij", sol)
            changed <- TRUE
        }
        if ( 
            ("hat.pj" %in% names(sol)) && (!"hat.pij" %in% names(sol))
            ) {
            sol$hat.pij <- case$inputs$pij$data / rowSums(case$inputs$pij$data) * sol$hat.pj
            # print("sol$hat.pij <- case$inputs$pij$data / rowSums(case$inputs$pij$data) * sol$hat.pj")
            # print(case$inputs$pij$data)
            # print(sol$hat.pj)
            sol <- info.rule("hat.pj -> hat.pij", sol)
            changed <- TRUE
        }

        # ci -> fi
        if ( 
            ("hat.ci" %in% names(sol)) && (!"hat.fi" %in% names(sol))
            ) {
            sol$hat.fi <- sol$hat.ci / sum(sol$hat.ci)
            #print("sol$hat.fi <- sol$hat.ci / sum(sol$hat.ci)")
            sol <- info.rule("hat.ci -> hat.fi", sol)
            changed <- TRUE
        }
        if ( 
            ("hat.cj" %in% names(sol)) && (!"hat.fj" %in% names(sol))
            ) {
            sol$hat.fj <- sol$hat.cj / sum(sol$hat.cj)
            #print("sol$hat.fj <- sol$hat.cj / sum(sol$hat.cj)")
            sol <- info.rule("hat.cj -> hat.fj", sol)
            changed <- TRUE
        }


        # pdi, ci -> ndi
        if ( 
            ("hat.pdi" %in% names(sol)) && ("hat.ci" %in% names(sol)) && (!"hat.ndi" %in% names(sol))
            ) {
            
            sol <- info.rule("hat.pdi, hat.ci -> hat.ndi", sol)
            sol$hat.ndi <- round_sum(t(t(sol$hat.pdi) * sol$hat.ci))
            sol$hat.ndi <- rectify.degree.counts(sol$hat.ndi, sol$hat.ni)   
            changed <- TRUE
        }
        if ( 
            ("hat.pdj" %in% names(sol)) && ("hat.cj" %in% names(sol)) && (!"hat.ndj" %in% names(sol))
            ) {
            
            sol <- info.rule("hat.pdj, hat.cj -> hat.ndj", sol)
            sol$hat.ndj <- round_sum(t(t(sol$hat.pdj) * sol$hat.cj))
            sol$hat.ndj <- rectify.degree.counts(sol$hat.ndj, sol$hat.nj)   
            changed <- TRUE
        }
    
    
        # di -> pdi
        if ( 
            ("hat.di" %in% names(sol)) && (!"hat.pdi" %in% names(sol))
            ) {
            sol <- info.rule("hat.di -> hat.pdi", sol)
            sol$hat.pdi <- tryCatch({
                        update_degree_distribution(case$inputs$pdi$data, sol$hat.di)
                    }, error = function(e) {
                        stop("The case is too constrained to be solved (hat.di are not compliant with the original pdi)")
                    })            #print("sol$hat.fi <- sol$hat.ci / sum(sol$hat.ci)")
            changed <- TRUE
        }

        if ( 
            ("hat.dj" %in% names(sol)) && (!"hat.pdj" %in% names(sol))
            ) {
            sol <- info.rule("hat.dj -> hat.pdj", sol)
            sol$hat.pdj <- tryCatch({
                        update_degree_distribution(case$inputs$pdj$data, sol$hat.dj)
                    }, error = function(e) {
                        stop("The case is too constrained to be solved (probabilities hat.dj are not compliant with the original pdj)")
                    })            #print("sol$hat.fi <- sol$hat.ci / sum(sol$hat.ci)")
            changed <- TRUE
        }

    
        if (!changed) {
            break
        }

    }

    #print("propagated")
    #print(sol)

    sol

}

assert.equal <- function(v1,v2,msg) {

    if (!identical(v1,v2) && !(all.equal(v1,v2)==TRUE)) {  # && ((length(v1) > 1) & sum((v1-v2)^2) > 0)
        cat("\nASSERT ERROR:",msg,"\n")
        print(v1)
        print(v2)
        
        print( all.equal(v1,v2) )
        #cat(msg,"(",str(v1),",",str(v2),")\n")
        return(1)
    }
    0
}

detect.problems <- function(sol, case, fail=TRUE) {
    
    problems <- 0

    # nA = sum(ni)
    if (!is.null(sol$hat.nA) && !is.null(sol$hat.ci)) {
        problems <- problems + assert.equal(sol$hat.nA, sum(sol$hat.ci), "hat.nA = sum(hat.ci)")
    }
    if (!is.null(sol$hat.nB) && !is.null(sol$hat.cj)) {
        problems <- problems + assert.equal(sol$hat.nB, sum(sol$hat.cj), "hat.nB = sum(hat.cj)")
    }

    # ni / nA = fi 
    if (!is.null(sol$hat.nA) && !is.null(sol$hat.ci) && !is.null(sol$hat.fi)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.ci/sol$hat.nA), 
                                        as.vector(sol$hat.fi), 
                                        "hat.ci/hat.nA = hat.fi")
    }
    if (!is.null(sol$hat.nB) && !is.null(sol$hat.cj) && !is.null(sol$hat.fj)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.cj/sol$hat.nB), 
                                        as.vector(sol$hat.fj), 
                                        "hat.cj/hat.nB = hat.fj")
    }  

    # sum(fi) = 1
    if (!is.null(sol$hat.fi)) {
        problems <- problems + assert.equal(1, sum(sol$hat.fi), "sum(hat.fi)=1")
    }
    if (!is.null(sol$hat.fj)) {
        problems <- problems + assert.equal(1, sum(sol$hat.fj), "sum(hat.fj)=1")
    }

    # ni = ci*di
    if (!is.null(sol$hat.ci) && !is.null(sol$hat.di) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                    as.vector(sol$hat.ci*sol$hat.di), 
                                    as.vector(sol$hat.ni), 
                                    "hat.ni = hat.ci*hat.di")
    }
    if (!is.null(sol$hat.cj) && !is.null(sol$hat.dj) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                    as.vector(sol$hat.cj*sol$hat.dj), 
                                    as.vector(sol$hat.nj), 
                                    "hat.nj = hat.cj*hat.dj")
    }

    # nL = sum(ni)
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(sol$hat.nL, sum(sol$hat.ni), "hat.nL = sum(hat.ni)")
    }
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(sol$hat.nL, sum(sol$hat.nj), "hat.nL = sum(hat.nj)")
    }

    # nL = sum(nij)
    if (!is.null(sol$hat.nL) && !is.null(sol$hat.nij)) {
        problems <- problems + assert.equal(sol$hat.nL, sum(sol$hat.nij), "hat.nL = sum(hat.nij)")
    }

    # sum(hat.pij) = 1
    if (!is.null(sol$hat.pij)) {
        problems <- problems + assert.equal(1, sum(sol$hat.pij), "sum(hat.pij) = 1")
    }

    # hat.pij = hat.nij/sum(hat.nij)
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.pij)) {
        problems <- problems + assert.equal(sol$hat.nij/sum(sol$hat.nij), sol$hat.pij, "hat.pij = hat.nij/sum(hat.nij)")
    }
    
    # colSum(hat.nij) = hat.ni 
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(as.vector(colSums(sol$hat.nij)), as.vector(sol$hat.ni), "hat.ni = colSums(hat.nij)")
    }
    if (!is.null(sol$hat.nij) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(as.vector(rowSums(sol$hat.nij)), as.vector(sol$hat.nj), "hat.nj = rowSums(hat.nij)")
    }

    # vsum(pdi) = 1
    if (!is.null(sol$hat.pdi)) {
        for (i in 1:ncol(sol$hat.pdi)) {
            problems <- problems + assert.equal(1, sum(sol$hat.pdi[,i]), paste("sum(hat.pdi[",i,"])=1",sep=""))       
        }
    } 
    if (!is.null(sol$hat.pdj)) {
        for (i in 1:ncol(sol$hat.pdj)) {
            problems <- problems + assert.equal(1, sum(sol$hat.pdj[,i]), paste("sum(hat.pdj[",i,"])=1",sep=""))       
        }
    } 

    # sum(hat.ndi) = hat.ci 
    if (!is.null(sol$hat.ndi) && !is.null(sol$hat.ci)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.ci), 
                                        as.vector(colSums(sol$hat.ndi)), 
                                        "hat.ci = sum(hat.ndi)")   
    }    
    if (!is.null(sol$hat.ndj) && !is.null(sol$hat.cj)) {
        problems <- problems + assert.equal(
                                        as.vector(sol$hat.cj), 
                                        as.vector(colSums(sol$hat.ndj)), 
                                        "hat.cj = sum(hat.ndj)")   
    }    

    # sum(hat.ndi * n) = hat.ni
    if (!is.null(sol$hat.ndi) && !is.null(sol$hat.ni)) {
        problems <- problems + assert.equal(
                                as.vector(colSums(sol$hat.ndi * 0:(nrow(sol$hat.ndi)-1))), 
                                as.vector(sol$hat.ni), 
                                "hat.ni = sum( n * ndi[n])")   
    }
    if (!is.null(sol$hat.ndj) && !is.null(sol$hat.nj)) {
        problems <- problems + assert.equal(
                                as.vector(colSums(sol$hat.ndj * 0:(nrow(sol$hat.ndj)-1))), 
                                as.vector(sol$hat.nj), 
                                "hat.nj = sum( n * ndj[n])")   
    }

    # TODO ? min max di dj
    if (!is.null(sol$hat.di)) {
        #print("TESTING DI")
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.di < case$inputs$min.di), 
                                "hat.di >= min.di")      
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.di > case$inputs$max.di), 
                                "hat.di <= max.di")      
    }
    if (!is.null(sol$hat.dj)) {
        #print("TESTING DJ")
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.dj < case$inputs$min.dj), 
                                "hat.dj >= min.dj")      
        problems <- problems + assert.equal(
                                FALSE, 
                                any(sol$hat.dj > case$inputs$max.dj), 
                                "hat.dj <= max.dj")      
    }

    if (fail && (problems > 0)) {
        stop("The case is too constrained to be solved (led to incoherent state)")
    }
    problems

}

detect.missing.chains <- function(sol) {

    cat("we know: ",names(sol),"\n")

    missing.chains <- list()
    current.chain <- NULL

    # "hat.ni", "hat.nj", 
    chain.n <- c("hat.nA", "hat.ci", "hat.fi", "hat.di", "hat.pij", "hat.dj", "hat.fj", "hat.cj", "hat.nB")

    for (i in 1:length(chain.n)) {
        var.name <- chain.n[[i]]
        
        #print(var.name)

        if (!(var.name %in% names(sol))) {
            cat("missing",var.name,"\n")
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

    print("missing chains")
    print(missing.chains)
}

resolve.missing.chain <- function(sol, chain, case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma) {

    cat("starting the investigation of missing chain: ",paste(chain,collapse=","),"\n",sep="")

    solutions <- list()

    for (var.hat.name in chain) {

        sol.hyp <- sol

        explored <- list()

        cat("let's investigate", var.hat.name, "\n")

        explored$investigated.var.name <- var.hat.name 
        explored$hypothesis <- setdiff(chain,c(var.hat.name))

        # let's add the hypothesis we know all the other ones
        # for (other.var.hat.name in setdiff(chain,c(var.hat.name))) {

        #     var.name <- substr(other.var.hat.name, 5, nchar(other.var.hat.name))

        #     cat("start with the hypothesis", other.var.hat.name, "=", var.name, "\n")

        #     # TODO
        #     if (other.var.hat.name == "hat.di") { sol.hyp[[other.var.hat.name]] <- case$inputs$di }
        #     else if (other.var.hat.name == "hat.dj") { sol.hyp[[other.var.hat.name]] <- case$inputs$dj }
        #     else if (other.var.hat.name == "hat.fi") { sol.hyp[[other.var.hat.name]] <- case$stats$fi }
        #     else if (other.var.hat.name == "hat.fj") { sol.hyp[[other.var.hat.name]] <- case$stats$fj }
        #     else if (other.var.hat.name == "hat.pij") { sol.hyp[[other.var.hat.name]] <- case$inputs$pij$data }
        #     else if (other.var.hat.name == "hat.nA") { sol.hyp[[other.var.hat.name]] <- nA }
        #     else if (other.var.hat.name == "hat.nB") { sol.hyp[[other.var.hat.name]] <- nB }

        #     else { cat("/!\\ cannot create such an hypothesis\n") }

        # }

        if (var.hat.name == "hat.di") { sol.hyp[[var.hat.name]] <- case$inputs$di }
        else if (var.hat.name == "hat.dj") { sol.hyp[[var.hat.name]] <- case$inputs$dj }
        else if (var.hat.name == "hat.fi") { sol.hyp[[var.hat.name]] <- case$stats$fi }
        else if (var.hat.name == "hat.fj") { sol.hyp[[var.hat.name]] <- case$stats$fj }
        else if (var.hat.name == "hat.pij") { sol.hyp[[var.hat.name]] <- case$inputs$pij$data }
        else if (var.hat.name == "hat.nA") { sol.hyp[[var.hat.name]] <- nA }
        else if (var.hat.name == "hat.nB") { sol.hyp[[var.hat.name]] <- nB }

        else { cat("/!\\ cannot create such an hypothesis\n") }

        if (is.null(sol.hyp$inference)) {
            sol.hyp$inference <- list()
        }
        sol.hyp$inference[[length(sol.hyp$inference)+1]] <- paste("hypothesis on ",var.hat.name,"=",substr(var.hat.name,5,nchar(var.hat.name)),sep="")

        # propagate on this basis
        sol.hyp <- tryCatch({
                        propagate.direct(sol.hyp, case)
                    }, error = function(e) {
                        cat("/!\\ error in this exploration, this solution is not valid\n\n")
                    })     
        
        # so know we have an estimation
        #candidate.sol <- sol.hyp[[var.hat.name]]
        #explored$value <- candidate.sol

        if (!all(chain %in% names(sol.hyp))) {
            # we cannot infer anything from this hypothesis; let's ignore it
            cat("we can infer nothing from our hypothesis\n")
            
        } else {
            cat("we found a solution based on an hypothesis on ",var.hat.name,"\n")
        
            #print(sol.hyp[[var.hat.name]])


            # let's come back to the original problem
            #sol.hyp <- sol
            # add our solution
            #sol.hyp[[var.hat.name]] <- candidate.sol
            # propagate on this basis
            #sol.hyp <- propagate.direct(sol.hyp, case)

            # do we achieve to solve the problem on this basis ?
            nb.problems <- detect.problems(sol.hyp, case, fail=FALSE)

            if (nb.problems > 0) {
                cat("this variable does not provides a valid solution to our problem","\n")
                
            } else {
                cat("we have a valid solution for ",var.hat.name,"\n")
                #print(sol.hyp)
                sol.hyp <- quantify.errors(sol.hyp, case, nA, nB)
                explored$sol <- sol.hyp
                solutions[[length(solutions)+1]] <- explored
            }
            
        }


    }


    cat("found ",length(solutions)," solutions\n")
    #print(solutions)

    res <- NULL 
    if (length(solutions) == 0) {
        cat("found no solution for this chain :-(\n\n")
    } else if (length(solutions) == 1) {
        # return the unique solution
        res <- solutions[[1]]$sol
    } else {
        # TODO pls 
        #print("TODO select the best? a mix ?")
        
        # create a vector with errors

        val.or.0 <- function(x) {
            if (is.null(x)) {
                return(0)
            } else {
                return(x)
            }
        }
        cumulated.errors <- rep(0, times=length(solutions))
        weights <- c(nu.A, phi.A, delta.A, gamma, delta.B, phi.B, nu.B)
        print("weights")
        print(weights)
        for (i in 1:length(solutions)) {

            s <- solutions[[i]]
            print(names(s$sol))

            errors <- c(
                        val.or.0(s$sol$error.nA), 
                        val.or.0(s$sol$mse.fi), 
                        val.or.0(s$sol$mse.di), 
                        val.or.0(s$sol$mse.pij), 
                        val.or.0(s$sol$mse.dj), 
                        val.or.0(s$sol$mse.fj), 
                        val.or.0(s$sol$error.nB)
                        )

            print("errors for this sol")
            print(errors)
            print("weighted")
            print(errors * weights)
            cumulated.errors[i] <- sum(errors * weights)
            
        }
        cat("\t=>", cumulated.errors,"\n")
        best.solution <- which.min(cumulated.errors)
        cat("\t=> will keep solution ",best.solution," which is the one with the lowest weighted error\n")
        res <- solutions[[best.solution]]$sol
        #print(names(solutions[[best.solution]]$sol))
    }

    res
}


resolve <- function(sol, case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma) {

    cat("we know: ",names(sol),"\n")

    # direct propagation of what we know
    sol.tmp <- propagate.direct(sol,case)

    # detect what to solve 

    #complete <- c("hat.nA", "hat.fi", "hat.ni", "hat.di", "hat.pij", "hat.nij", "hat.dj", "hat.fj", "hat", "mse.di", "mse.pij","mse.dj","mse.fj","err.nB")

    #chain.p <- c("hat.fi", "hat.di", "hat.pij", "hat.dj", "hat.fj")

    
    missing.chains <- detect.missing.chains(sol.tmp)
    while (length(missing.chains) > 0) {

    #for (chain in missing.chains) {
        chain <- missing.chains[[1]]
        missing.chains[[1]] <- NULL

        found <- resolve.missing.chain(sol.tmp, chain, case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma)

        if (!is.null(found)) {
            cat("found a solution which we will use. It provides ",paste(names(found),collapse=","),"\n", sep="")
            sol.tmp <- found 
            missing.chains <- detect.missing.chains(sol.tmp)
        }

        cat("still ",length(missing.chains)," missing\n")
        #print(missing.chains)
    }

    # print("propagated")
    # print(sol)

    return(sol.tmp)
}


quantify.errors <- function(sol, case, nA, nB) {

    # measure errors
    # MSE fi
    if (!is.null(sol$hat.fi)) {
       sol$mse.fi <- mean( ( sol$hat.fi - case$stats$fi)^2  )
    }
    if (!is.null(sol$hat.fj)) {
        sol$mse.fj <- mean( ( sol$hat.fj - case$stats$fj)^2  )
    }

    # MSE di
    if (!is.null(sol$hat.di)) {
        sol$mse.di <- mean( ( sol$hat.di - case$inputs$di)^2  )
    }
    if (!is.null(sol$hat.dj)) {
        sol$mse.dj <- mean( ( sol$hat.dj - case$inputs$dj)^2  )
    }

    # MSE pdi 
    if (!is.null(sol$hat.pdi)) {
        sol$mse.pdi <- mean( ( sol$hat.pdi - case$inputs$pdi$data)^2  )
    }
    if (!is.null(sol$hat.pdj)) {
        sol$mse.pdj <- mean( ( sol$hat.pdj - case$inputs$pdj$data)^2  )
    }

    # MSE pij
    if (!is.null(sol$hat.pij)) {
        sol$mse.pij <- mean( ( sol$hat.pij - case$inputs$pij$data )^2  )
    }

    # error n
    if (!is.null(sol$hat.nA)) {
        sol$error.nA <- abs( sol$hat.nA - nA )
    }
    if (!is.null(sol$hat.nB)) {
        sol$error.nB <- abs( sol$hat.nB - nB )
    }

    sol
}

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

matching.arbitrate <- function(case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma) {

    
    mixx <- list()

    # start with masks
    # TODO reuse
    mask.fi <- as.integer(case$stats$fi > 0)
    mask.fj <- as.integer(case$stats$fj > 0)
    mask.di <- as.integer(case$inputs$di > 0)
    mask.dj <- as.integer(case$inputs$di > 0)
    mask.pij <- (case$inputs$pij$data > 0)*1

    mask.pij.all <-  t(t(mask.pij) * mask.fi * mask.di) * mask.fj * mask.dj
    mask.fi.all <- as.integer(colSums(mask.pij.all) > 0)
    mask.fj.all <- as.integer(rowSums(mask.pij.all) > 0)



    cat("\n\n=============== QUE VIVA LA RESOLUTION !\n")
    

    sol <- list()
    if (nu.A == 0)  { sol$hat.nA <- nA }
    if (phi.A == 0) { sol$hat.fi <- case$stats$fi }
    if (delta.A == 0) { sol$hat.di <- case$inputs$di }

    if (nu.B == 0) { sol$hat.nB <- nB }
    if (phi.B == 0) { sol$hat.fj <- case$stats$fj }
    if (delta.B == 0) { sol$hat.dj <- case$inputs$dj }

    if (gamma == 0) { sol$hat.pij <- case$inputs$pij$data }

    #parameters <- list(nA=nA, nB=nB, nu.A=nu.A, phi.A=phi.A, delta.A=delta.A, nu.B=nu.B, phi.B=phi.B, delta.B=delta.B, gamma=gamma)

    # detect issues right now; maybe the problem is overconstrained or badly constrainted
    detect.problems(sol, case)

    sol <- resolve(sol, case, nA, nB, nu.A, phi.A, delta.A, nu.B, phi.B, delta.B, gamma)

    
    detect.problems(sol, case)

    # measure errors
    sol <- quantify.errors(sol, case, nA, nB)

    # ensure the form is ok
    sol <- ensure.form(sol, case)

    print("end of resolution")
    print(sol)

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
        msename <- paste("mse.",subname,sep="")
        cat("$",name," [MSE ",x$gen[[msename]],"]:\n",sep="")
        print(x$gen[[name]])
        cat("\n")
    }

    cat("$hat.nA: ",x$gen$hat.nA," [err ",x$gen$error.nA,"]\n",sep="")

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
# TODO we need to first round the factors for weights !

