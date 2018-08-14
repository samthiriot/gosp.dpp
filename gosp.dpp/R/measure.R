
#' Display the result of measures on a generated population
#'
#' @param x the measure to display
#' @param ... ignored 
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_measure <- function(x,...) {

    cat("measures on the generated population\n\n")

    cat("$nL=",x$nL,"\n",sep="")

    cat("\n$fi [MSE ",x$fi$mse.f.orig,"]\n",sep="")
    print(x$fi$hat.f)
    
    cat("\n$pdi [MSE ",x$pdi$mse.orig,"]\n",sep="")
    print(x$pdi$hat.pd)

    cat("\n$pij [MSE ",x$pij$mse.orig,"]\n",sep="")
    print(x$pij$hat.pij)

    cat("\n$pdj [MSE ",x$pdj$mse.orig,"]\n",sep="")
    print(x$pdj$hat.pd)
    
    cat("\n$fj [MSE ",x$fj$mse.f.orig,"]\n",sep="")
    print(x$fj$hat.f)
    
}


#' Display a generated population
#'
#' @param x the generation result to display
#' @param ... ignored 
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_result <- function(x,...) {
    cat("generation result\n")

    cat("$pop\n")
    print(x$pop)

    cat("\n$measure\n")
    print(x$measure)
}

#' Measures the matching probabilities from a generated population
#'
#' Measures the matching probabilities from a generation result. 
#' 
#' @param pop the generated population
#' @param sample.A the original sample for population A
#' @param sample.B the original sample for population B
#' @param pij the original matching probabilities
#' @param mix.pij the solved matching probabilities
#' @param A2l2B the joined population A / links / population B
#' @param verbose if TRUE, detailed messages are printed
#'
#' @return a dataframe with the measured probabilities
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
measure.pij <- function(pop, sample.A, sample.B, pij, mix.pij, A2l2B, verbose=FALSE) {


	hat.nij <- pij$data

	for (cA in colnames(pij$data)) {

		k2v.A <- extract_attributes_values(cA)

		for (cB in rownames(pij$data)) {

			k2v.B <- extract_attributes_values(cB)

			if (verbose)
				cat("\tcount how many links connect ", cA, " and ", cB, "\n", sep="")

			criteriaAraw <- 1:nrow(A2l2B)
			for (k in names(k2v.A)) {
				if (!(k %in% colnames(A2l2B))) {
					k <- paste(k,"_A",sep="")
				}
	        	criteriaAraw <- intersect(criteriaAraw, which(A2l2B[k] == k2v.A[[k]]))
	        }

			criteriaBraw <- criteriaAraw
			for (k in names(k2v.B)) {
				if (!(k %in% colnames(A2l2B))) {
					k <- paste(k,"_B",sep="")
				}
	        	criteriaBraw <- intersect(criteriaBraw, which(A2l2B[k] == k2v.B[[k]]))
	        }

			count <- nrow(A2l2B[criteriaBraw,])
			if (verbose)
				cat("\t=>", count, "\n")
			
			hat.nij[cB,cA] <- count

		}

	}

    hat.pij <- hat.nij / nrow(pop$links)
	
	list( 
			hat.nij=hat.nij, 
			hat.pij=hat.pij, 
			mse.orig=mean( (pij$data-hat.pij)^2 ) , 
			mse.target=mean( (mix.pij-hat.pij)^2 ) 
			)
}

# TODO cleanup
#' Measures the the contigencies and frequencies from a generated population
#'
#' Measures the contigencies and frequencies from a generated population. 
#' 
#' @param target.ni the original ni 
#' @param orig.fi frequencies 
#' @param colname the name of the column on which fi depends on 
#' @param dico the original dictionary 
#' @param pop the generated population 
#'
#' @return a list containing the measures
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
measure.ci_fi <- function(target.ni, orig.fi, colname, dico, pop) {

    hat.ni <- c()
    for (codeA in dico$encoding[[colname]]) {
        hat.ni <- c(hat.ni, nrow(pop[which(pop[colname]==codeA),]))
    }

    hat.fi <- hat.ni / sum(hat.ni)

    list(
        hat.n=hat.ni,
        mse.n.target=mean( (target.ni - hat.ni)^2 ),
        hat.f=hat.fi,
        mse.f.orig=mean( (orig.fi - hat.fi)^2 )
        )

}

#' Measure the probability distribution
#' 
#' Measure the distributions of probabilities from the 
#' generated population.
#' 
#' @param dico the dictionary
#' @param pdn.orig the original expected distributions
#' @param sample the sample to measure
#' @param verbose detailed messages printed if TRUE
#' 
#' @return a list containing the measures
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
measure.pdn <- function(dico, pdn.orig, sample, verbose=FALSE) {


	hat.pdn <- pdn.orig$data
	for (codeA in dico$encoding[[pdn.orig$attributes]]) {
		#codeA <- dico$encoding[[colname]][[cA]]

		total <- 0
		for (x in 1:nrow(pdn.orig$data)) {
			degree <- x-1

			if (verbose)
				cat("\tmeasuring how many entities have",pdn.orig$attributes,"=",codeA,"and current degree=",degree,"")
			count <- nrow(sample[which( (sample[pdn.orig$attributes]==codeA) & (sample$current.degree==degree) ),])
			
			if (verbose)
				cat("\t=>", count, "\n")
			hat.pdn[x,codeA] <- count 
			total <- total + count 
		}
		hat.pdn[,codeA] <- hat.pdn[,codeA] / total 
	}

	list(
		hat.pd=hat.pdn,
		mse.orig=mean( (pdn.orig$data-hat.pdn)^2 ),
		count.over.degree=nrow(sample[which( sample$current.degree>sample$target.degree),]),
		count.under.degree=nrow(sample[which( sample$current.degree<sample$target.degree),])
		)
}

#' Merge entities and attributes
#' 
#' Merges a generated population: returns the join of population A,
#' links and population B. 
#' 
#' The underlying operation is done by \code{\link{merge}}.
#'
#' @param pop the generated population
#' @return a dataframe from the join 
#' 
#' @export 
#' 
merge_links <- function(pop) {

	# merge the datasets	
	A2l <- merge(pop$A, pop$links, by="id.A", suffixes=c("_A", "_l"))
	A2l2B <- merge(A2l, pop$B, by="id.B", suffixes=c("_A", "_B"))

	A2l2B
}

#' Measures the characteristics of a generated population
#'
#' Measures the characteristics of a generated population:
#' the contigencies and frequencies with \code{\link{measure.ci_fi}},
#' the degree distributions with \code{\link{measure.pdn}},
#' the matching probabilities with \code{\link{measure.pij}}.
#' 
#' @param case the original case
#' @param pop the generated population
#' @param sample.A the original sample A
#' @param sample.B the original sample B
#' @param pij the target matching probabilities
#' @param verbose detailed messages are printed if TRUE
#' 
#' @return a list containg the measures 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
measure.population <- function(case, pop, sample.A, sample.B, pij, verbose=FALSE) {

	# merge the datasets	
	A2l2B <- merge_links(pop)

	if (verbose)
		cat("total links measured: ", nrow(pop$links),"\n")

	res <- list(
				nL=nrow(A2l2B),
				
                fi=measure.ci_fi(case$gen$hat.ni, case$stats$fi, case$inputs$pij$Ai, case$inputs$sample.A$dictionary, pop$A),
                fj=measure.ci_fi(case$gen$hat.nj, case$stats$fj, case$inputs$pij$Bi, case$inputs$sample.B$dictionary, pop$B),
                
				pdi=measure.pdn(case$inputs$sample.A$dictionary, case$inputs$pdi, pop$A),
				pdj=measure.pdn(case$inputs$sample.B$dictionary, case$inputs$pdj, pop$B),

				pij=measure.pij(pop, sample.A, sample.B, pij, case$gen$hat.pij, A2l2B, verbose)
				)

	res

}
