
#' @importFrom utils head
#' @importFrom utils tail
NULL

#' Takes an original sample, and returns \code{count.required} entities 
#' selected according to the weight column \code{colname.weight} of this sample (if any).
#'
#' If the population is too small, it will be made bigger by replication of entities. 
#' If the population is too big, only some entities will be selected depending to weights. 
#' If the population is exactly of the right size, then it will be entirely kept as 
#' in the original sample.
#'
#' The actual sampling is entirely delegated to the \link{dplyr}{sample_n} function.
#' 
#'
#' @param sample the sample to use for loading 
#' @param count.required the count of entities expected as a result 
#' @param colname.weight a string describing the column to be used as a weight in the sample
#'
#' @importFrom dplyr sample_n
#'
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.generate.resize_population <- function(sample, count.required, colname.weight) {
	# weights.factor <- count.required / weights.available
	#print(weights.factor)

	#print(count.required)
	
	selected <- sample_n(sample, count.required, weight=sample[,colname.weight], replace=TRUE)
	#print(nrow(selected))

	# return everything but the weights
	selected[,names(selected) != colname.weight]
}


matching.generate.copy_population <- function(n, sample, verbose=FALSE) {

	if (verbose)
    	cat("copying population based on attributes '",paste(names(n),collapse=","),"' with ",
    		paste(n,collapse=","),"\n",sep="")
    
    target <- sample$sample[FALSE,] 	# empty dataframe with the same properties
    colnames(target) <- colnames(sample$sample)

	for (name in names(n)) {

        count.required <- n[name]

      	if (verbose)
    		cat("\tcopying ", count.required, " having ",name,"\n",sep="")
        
        k2v <- extract_attributes_values(name)
        #print(k2v)

        # accumulate the conditions on each key / value
        selected_ids <- 1:nrow(sample$sample) 
        for (k in names(k2v)) {
        	selected_ids <- intersect(selected_ids, which(sample$sample[k] == k2v[[k]]))
        }

        # identify the subpart of the population having these values
        available <- sample$sample[selected_ids,]
        
		count.available <- nrow(available)
		weights.available <- sum(available[,sample$dictionary$colname.weight])

		if (verbose)
			cat("\tshould copy ", count.required, " individuals for ", name,
					" and found ",count.available, " individuals with weights summing to ", weights.available,"\n",sep="")

		toAdd <- matching.generate.resize_population(available, count.required, sample$dictionary$colname.weight)

		# extend the original target population
		target <- rbind(target, toAdd)
		
	}

    if (sum(n) != nrow(target)) {   
        stop(paste("we should have copied ",sum(n)," rows but we only copied",nrow(target),"instead\n",sep=""))
    }
	target
}


#' Adds a target degree column to a population 
#' based on the distribution of contigency for degrees passed as a parameter.
#'
#' @param samp the original sample (used for its dictionary) 
#' @param pop the population to which the column will be added to
#' @param n the expected contigencies (used as a control)
#' @param ndx the contigencies for each degree
#' @param verbose if TRUE, detailed messages are printed
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.generate.add_degree <- function(samp, pop, n, ndx, verbose) {


	# add the column with NA
	pop$target.degree <- rep.int(NA, nrow(pop))

	# for each different value  
	#  colnames(pdi$data)
	for (name in names(ndx)) {

		k2v <- extract_attributes_values(name)
		
		totalForDegree <- 0
		lastDegreeNonNull <- 0

		for (degree in 0:(nrow(ndx)-1)) {

			count.required <- ndx[degree+1,name]
			totalForDegree <- totalForDegree + count.required*degree

			# accumulate selection criteria: write set of attributes, and also no defined degree
			criteriaRaw <- which( is.na(pop["target.degree"]) )
			for (k in names(k2v)) {
	        	criteriaRaw <- intersect(criteriaRaw, which(pop[k] == k2v[[k]]))
	        }

			if (verbose)
				cat("\tset target degree ", degree, " for ", count.required, " over ", length(criteriaRaw), " having ", name,"\n", sep="")
			
			if (count.required == 0) {
				next 
			}
			
			lastDegreeNonNull <- degree

			criteria <- NULL 
			if (count.required < length(criteriaRaw)) {
				criteria <- sample(criteriaRaw, count.required)
			} else {
				criteria <- criteriaRaw
			}

			pop[criteria,"target.degree"] <- degree
		}


		if (totalForDegree < n[name]) {
	
			warning(paste("oops, not created enough slots here:",totalForDegree," for ",n[name]," expected; ",
					"defining the last ",length(criteriaRaw), "to degree ",lastDegreeNonNull,"\n",sep=""))

		}
		
	}

	still.missing <- pop[ which( is.na(pop["target.degree"]) ),]
	if (nrow(still.missing)>0) {
		warning("some entities were not given any degree, they will be fixed at zero target degree")
		print(head(still.missing))
		pop[ which( is.na(pop["target.degree"]) ),"target.degree"] <- 0
	}
	

	# so we might assess how much of each we need to precisely meet the expectation 
	pop
}

#' Generates a synthetic population based on the resolution of the problem,
#' applied on two samples \code{sample.A} and \code{sample.B}. 
#' It first resizes the two populations so they reach the expected sizes, 
#' the adds degrees to them, then creates links based on degrees 
#' and characteristics, before measuring the result.  
#'
#' @param case the case resolved by the \code{\link{matching.solve}} function. 
#' @param sample.A the sample to use as population A 
#' @param sample.B the sample to use as population B
#' @param verbose displays detailed messages if TRUE
#' @return the generated population
#' 
#' @export 
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
matching.generate <- function(case, sample.A, sample.B, verbose=FALSE) {

	if (class(case) != "dpp_resolved") 
		stop("case should be the result of a solving by matching.solve")

	if (verbose)
    	cat("starting generation\n")

	# copy the individuals required for A
	targetA <- sample.A$sample[FALSE,] 	# empty with the same properties

	nni <- case$gen$hat.ci 
	nnj <- case$gen$hat.cj 

	# copy individuals
	targetA <- matching.generate.copy_population(n=nni, sample=sample.A, verbose=verbose)
	targetB <- matching.generate.copy_population(n=nnj, sample=sample.B, verbose=verbose)

    # cat("creating ids...\n")
	# create unique ids (unique at the level of both populations)
	targetA$id.A <- seq.int(1, nrow(targetA))
	targetB$id.B <- seq.int(nrow(targetA)+1, nrow(targetA)+nrow(targetB))
	
	# create a matching table with links (we init it with nothing so we avoid costly resizing later)
	links <- data.frame(
			id.A=c(),
			id.B=c()
			)

    if (verbose)
    	cat("\nadding target degrees for population A...\n")
	# add columns for target degree 
	targetA <- matching.generate.add_degree(sample.A, targetA, case$gen$hat.ni, case$gen$hat.ndi, verbose=verbose)
	
	if (verbose)
    	cat("\nadding target degrees for population B...\n")
	targetB <- matching.generate.add_degree(sample.B, targetB, case$gen$hat.nj, case$gen$hat.ndj, verbose=verbose)
	
    # add columns for obtained degrees (0 for now)
	targetA$current.degree <- rep.int(0, nrow(targetA))
	targetB$current.degree <- rep.int(0, nrow(targetB))

    if (verbose)
    	cat("\nmatching A and B to create total",case$gen$hat.nL,"links...\n")

	# match them 
	for (cA in colnames(case$gen$hat.nij)) {

		k2v.A <- extract_attributes_values(cA)
		
		for (cB in rownames(case$gen$hat.nij)) {

			k2v.B <- extract_attributes_values(cB)

			count.required <- case$gen$hat.nij[cB,cA] 

			if (verbose)
				cat("\tshould create", count.required, "links for A:\t", cA, "\tand B:\t", cB, "\n")

			pass.remaining <- 1

			links.toadd.A <- c()
			links.toadd.B <- c()
			

			while ( (pass.remaining > 0) & (count.required > 0) ) {

				pass.remaining <- pass.remaining - 1

				criteriaAraw <- which( targetA$current.degree < targetA$target.degree )
				for (k in names(k2v.A)) {
		        	criteriaAraw <- intersect(criteriaAraw, which(targetA[k] == k2v.A[[k]]))
		        }

				criteriaBraw <- which( targetB$current.degree < targetB$target.degree )
				for (k in names(k2v.B)) {
		        	criteriaBraw <- intersect(criteriaBraw, which(targetB[k] == k2v.B[[k]]))
		        }
				
				available.degree.A <- max(targetA[criteriaAraw,"target.degree"] - targetA[criteriaAraw,"current.degree"])
				available.degree.B <- max(targetB[criteriaBraw,"target.degree"] - targetB[criteriaBraw,"current.degree"])

				# first pass for largely connectable ones
				if ( (available.degree.A > 1) | (available.degree.B > 1) ) {
					pass.remaining <- 1
				}

				count.found <- min(length(criteriaAraw), length(criteriaBraw), count.required)

				criteriaA <- sample(criteriaAraw, count.found)
				criteriaB <- sample(criteriaBraw, count.found)

				# add the links
				links.toadd.A <- c(links.toadd.A, targetA[criteriaA,"id.A"])
				links.toadd.B <- c(links.toadd.B, targetB[criteriaB,"id.B"])

				# update the degree
				targetA[criteriaA,"current.degree"] <- targetA[criteriaA,"current.degree"] + 1
				targetB[criteriaB,"current.degree"] <- targetB[criteriaB,"current.degree"] + 1

				count.required <- count.required - count.found
			}


			if ( count.required > 0)  {
				warning(paste("\t/!\\ failed to create ", count.required," link(s) ",
						"for ",case$inputs$pij$Ai, "=", cA, 
						"\tand B:\t",  case$inputs$pij$Bi, "=", cB,"\n",sep=""))
			}

			links.toadd <- data.frame(
						id.A=links.toadd.A,
						id.B=links.toadd.B
						)
			links <- rbind(links, links.toadd)
			
		}

	}

	# add a column for current degree
	res <- case
	res$pop <- list(A=targetA, B=targetB, links=links)
    class(res$pop) <- "dpp_population" 
	res$measure <- measure.population(case, res$pop, sample.A, sample.B, case$inputs$pij)
    class(res$measure) <- "dpp_measure"

    class(res) <- "dpp_result"

	res
}

#' Display a generated population
#' 
#' @param x the population to print
#' @param ... ignored
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_population <- function(x,...) {
    cat("Synthetic population with ",nrow(x$A)," parents, ",nrow(x$B)," children, ",nrow(x$links)," links\n",sep="")
    
    cat("\n$A (first lines)\n")
    print(head(x$A))
    
    cat("\n$B (first lines)\n")
    print(head(x$B))

    cat("\n$links (first lines)\n")
    print(head(x$links,n=15))
}
