
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
#' The actual sampling is entirely delegated to the \link{dplyr::sample_n} function.
#' 
#' TODO ensure it is able to load exactly the same entities  
#'
#' @param sample the sample to use for loading 
#' @param count.required the count of entities expected as a result 
#' @param colname.weight a string describing the column to be used as a weight in the sample
#'
#' @importFrom dplyr sample_n
#'
matching.generate.resize_population <- function(sample, count.required, colname.weight) {
	# weights.factor <- count.required / weights.available
	#print(weights.factor)
	
	selected <- sample_n(sample, count.required, weight=sample[,colname.weight], replace=TRUE)
	#print(nrow(selected))

	selected
}

matching.generate.copy_population <- function(n, select.colname, sample) {

    cat("copying population based on attribute '",select.colname,"' with ",n,"\n")

    #cols_without_weight <- -which(names(sample$sample) == sample$dictionnary$colname.weight)

    target <- sample$sample[FALSE,] 	# empty with the same properties


	for (name in names(n)) {

        count.required <- n[name]

        #cat("* copying", count.required, " having ",select.colname,"=",name,"\n")

        code <- sample$dictionnary$encoding[[select.colname]][[name]] # .table
        
        #cat("* copying", count.required, " having ",select.colname,"=",code,"\n")
        
        if (is.null(code)) {
            stop(paste("was unable to find the encoded value for ",select.colname,"=",name,sep=""))
        }
        
        # TODO do not copy weight
        #cols <- (colnames(sample$sample) %in% c(sample$dictionnary$colname.weight))
		#print(cols_without_weight)
        available <- sample$sample[which(sample$sample[select.colname] == code),]

		count.available <- nrow(available)
		weights.available <- sum(available[,sample$dictionnary$colname.weight])

		cat("should copy ", count.required, " individuals for ", select.colname, "=", code, "(for",name,") ",
					"and found ",count.available, " individuals with weights summing to ", weights.available,"\n")

		toAdd <- matching.generate.resize_population(available, count.required, sample$dictionnary$colname.weight)

		# extend the original target population
		target <- rbind(target, toAdd)
		
	}

    if (sum(n) != nrow(target)) {   
        cat("hum, we should have copied ",sum(n),"rows but we only copied",nrow(target),"instead\n")
        stop("wrong size")
    }
	target
}

#' Adds a target degree column to a population 
#' based on the distribution of contigency for degrees passed as a parameter.
#'
#' @param samp the original sample (used for its dictionnary) 
#' @param pop the population to which the column will be added to
#' @param n the expected contigencies (used as a control)
#' @param ndx the contigencies for each degree
#' @param colname the name of the column on which the degrees depend on
matching.generate.add_degree <- function(samp, pop, n, ndx, colname) {

	# we have to reweight the distribution of degrees
	# TODO do that earlier at discretisation tile

    cat("adding degress based on ndx", "\n")
    print(ndx)

	# add the column with NA
	pop$target.degree <- rep.int(NA, nrow(pop))

	# for each different value  
	#  colnames(pdi$data)
	for (name in names(samp$dictionnary$encoding[[colname]])) {

		code <- samp$dictionnary$encoding[[colname]][[name]]
		
		totalForDegree <- 0
		lastDegreeNonNull <- 0

		for (degree in 0:(nrow(ndx)-1)) {

			count.required <- ndx[degree+1,make.names(name)]
			totalForDegree <- totalForDegree + count.required*degree

			criteriaRaw <- which( (pop[colname] == code) & is.na(pop["target.degree"]) )
			
			cat("set target degree", degree, "for", count.required, "over", length(criteriaRaw), "having", colname, "=", code, "(",name,")\n")
			
			if (count.required == 0) {
				next 
			}
			
			lastDegreeNonNull <- degree

			#print(criteriaRaw)
			criteria <- NULL 
			if (count.required < length(criteriaRaw)) {
				criteria <- sample(criteriaRaw, count.required)
			} else {
				criteria <- criteriaRaw
			}

			#print("criteria")
			#print(criteria)
			
			pop[criteria,"target.degree"] <- degree
		}


		if (totalForDegree < n[name]) {
	
			cat("oops, not created enough slots here:",totalForDegree,"for",n[name],"expected\n")

			criteriaRaw <- which( (pop[colname] == code) & is.na(pop["target.degree"]) )
			
			cat("defining the last ",length(criteriaRaw), "to degree ",lastDegreeNonNull,"\n")
			
			pop[criteriaRaw,"target.degree"] <- lastDegreeNonNull
	
		}
		
	}

	still.missing <- pop[ which( is.na(pop["target.degree"]) ),]
	if (nrow(still.missing)>0) {
		print("oooops some values are still missing, they will be fixed at zero target degree")
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
#' @param case the case resolved by the \link{resolve} function. 
#' @param sample.A the sample to use as population A 
#' @param sample.B the sample to use as population B
#' @return the generated population
#' 
#' @export 
matching.generate <- function(case, sample.A, sample.B) {

    cat("starting generation\n")

	# copy the individuals required for A
	targetA <- sample.A$sample[FALSE,] 	# empty with the same properties

	nni <- case$gen$hat.ci 
	nnj <- case$gen$hat.cj 

	# copy individuals
	targetA <- matching.generate.copy_population(n=nni, select.colname=case$inputs$pij$Ai, sample=sample.A)
	targetB <- matching.generate.copy_population(n=nnj, select.colname=case$inputs$pij$Bi, sample=sample.B)

    cat("creating ids...\n")
	# create unique ids (unique at the level of both populations)
	targetA$id.A <- seq.int(1, nrow(targetA))
	targetB$id.B <- seq.int(nrow(targetA)+1, nrow(targetA)+nrow(targetB))
	
	# create a matching table with links (we init it with nothing so we avoid costly resizing later)
	links <- data.frame(
			id.A=c(),
			id.B=c()
			)

	# remove the weight columns (which have no meaning here anymore)
    # TODO !!!

    cat("\nadding target degrees...\n")
	# add columns for target degree 
	targetA <- matching.generate.add_degree(sample.A, targetA, case$gen$hat.ni, case$gen$hat.ndi, case$inputs$pdi$attributes)
	targetB <- matching.generate.add_degree(sample.B, targetB, case$gen$hat.nj, case$gen$hat.ndj, case$inputs$pdj$attributes)
	
    # add columns for obtained degrees (0 for now)
	targetA$current.degree <- rep.int(0, nrow(targetA))
	targetB$current.degree <- rep.int(0, nrow(targetB))

    cat("\nmatching A and B to create total",case$gen$hat.nL,"links...\n")

	# match them 
	for (cA in colnames(case$gen$hat.nij)) {

		codeA <- sample.A$dictionnary$encoding[[case$inputs$pij$Ai]][[cA]]

		for (cB in rownames(case$gen$hat.nij)) {

			codeB <- sample.B$dictionnary$encoding[[case$inputs$pij$Bi]][[cB]]

			count.required <- case$gen$hat.nij[cB,cA] 

			cat("should create", count.required, "links for A:\t", case$inputs$pij$Ai, "=", cA, "\tand B:\t",  case$inputs$pij$Bi, "=", cB, "\n")

			pass.remaining <- 1

			links.toadd.A <- c()
			links.toadd.B <- c()
			

			while ( (pass.remaining > 0) & (count.required > 0) ) {

				pass.remaining <- pass.remaining - 1

				#print(cat("should find ", count.required, "for ", case$inputs$pij$Ai, "=", cA, "and", cB))

				criteriaAraw <- which( (targetA[case$inputs$pij$Ai] == codeA) & (targetA$current.degree < targetA$target.degree) )
				criteriaBraw <- which( (targetB[case$inputs$pij$Bi] == codeB) & (targetB$current.degree < targetB$target.degree) )
				
				available.degree.A <- max(targetA[criteriaAraw,"target.degree"] - targetA[criteriaAraw,"current.degree"])
				available.degree.B <- max(targetB[criteriaBraw,"target.degree"] - targetB[criteriaBraw,"current.degree"])

				#print(cat("\tavailable degree A", available.degree.A, "available.degree.B", available.degree.B))

				# first pass for largely connectable ones
				if (available.degree.A > 1) {
					#print("\tdoing a pass for high degree availability (A) first")
					pass.remaining <- 1
					
					#criteriaAraw <- which( (targetA[case$inputs$pij$Ai] == codeA) & (targetA$current.degree+1 < targetA$target.degree) )
					
				}
				if (available.degree.B > 1) {
					#print("\tdoing a pass for high degree availability (B) first")
					pass.remaining <- 1
		
					#criteriaBraw <- which( (targetB[case$inputs$pij$Bi] == codeB) & (targetB$current.degree+1 < targetB$target.degree) )		
				}

				count.found <- min(length(criteriaAraw), length(criteriaBraw), count.required)

				# TODO keep it ???
				#print(targetA[criteriaAraw,"target.degree"]-targetA[criteriaAraw,"current.degree"])
				cat("\tfound", length(criteriaAraw), "in A and", length(criteriaBraw), "in B")
				cat("\t=> creating", count.found, "links\n")
				
				criteriaA <- sample(criteriaAraw, count.found) #, prob=targetA[criteriaAraw,"target.degree"]-targetA[criteriaAraw,"current.degree"]
				criteriaB <- sample(criteriaBraw, count.found)

				# add the links
				links.toadd.A <- c(links.toadd.A, targetA[criteriaA,"id.A"])
				links.toadd.B <- c(links.toadd.B, targetB[criteriaB,"id.B"])

				# update the degree

				#print(head(targetA[criteriaA,"current.degree"], n=50))
				targetA[criteriaA,"current.degree"] <- targetA[criteriaA,"current.degree"] + 1
				#print(head(targetA[criteriaA,"current.degree"], n=50))

				targetB[criteriaB,"current.degree"] <- targetB[criteriaB,"current.degree"] + 1

				#print(targetA[criteriaA,"current.degree"])

				count.required <- count.required - count.found
			}


			if ( count.required > 0)  {

				cat("\t/!\\ failed to create", count.required,"link (rounding effect ?)\n")

				if (count.required > 1) {
					# TODO for pop B ?
					criteriaAllAvailable <- which(targetA$current.degree < targetA$target.degree)
					actualCount <- sum(targetA[criteriaAllAvailable,"target.degree"]-targetA[criteriaAllAvailable,"current.degree"])

					actualCountUs <- sum(targetA[criteriaAraw,"target.degree"]-targetA[criteriaAraw,"current.degree"])


					print(cat("not enough individuals found in the population A ", actualCountUs," to fullfill the demand ", count.required,
								". We created ", nni[cA], " entities for", sum(case$gen$hat.nij[,cA]), " required for matching" ))
					print(cat("we were supposed to match based on "))
					print(case$gen$hat.nij)
					print(cat("with degrees being "))
					print(nni)

					print(cat("this total of slots are available:",actualCount))

					print(cat("but only this total of slots are available for us:",actualCountUs))

					print("population A is currently")
					print(head(targetA[criteriaAraw,], n=100))
				}
			}

			links.toadd <- data.frame(
						id.A=links.toadd.A,
						id.B=links.toadd.B
						)
			links <- rbind(links, links.toadd)
			
		}

	}

	# compute degree (?)

	# print("should add links A")
	#print(links)


	# add a column for current degree

	res <- case
	res$pop <- list(A=targetA, B=targetB, links=links)
    class(res$pop) <- "dpp_population" 
	res$measure <- measure.population(case, res$pop, sample.A, sample.B, case$inputs$pij)
    class(res$measure) <- "dpp_measure"

    class(res) <- "dpp_result"

	res
}


print.dpp_population <- function(x,...) {
    cat("result population with ",nrow(x$A)," parents, ",nrow(x$B)," children, ",nrow(x$links)," links\n",sep="")
    
    cat("\n$A (first lines)\n")
    print(head(x$A))
    
    cat("\n$B (first lines)\n")
    print(head(x$B))

    cat("\n$links (first lines)\n")
    print(head(x$links,n=15))
}
