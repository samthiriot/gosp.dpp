
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

print.dpp_result <- function(x,...) {
    cat("generation result\n")

    cat("$pop\n")
    print(x$pop)

    cat("\n$measure\n")
    print(x$measure)
}

measure.pij <- function(pop, sample.A, sample.B, pij, mix.pij, A2l2B) {


	hat.nij <- pij$data

	for (cA in colnames(pij$data)) {

		codeA <- sample.A$dictionnary$encoding[[pij$Ai]][[cA]]

		for (cB in rownames(pij$data)) {

			codeB <- sample.B$dictionnary$encoding[[pij$Bi]][[cB]]

			cat("\tcount how many links connect", pij$Ai, "=", codeA, "(", cA, ") and", pij$Bi, "=", codeB, "(", cB, ")")

			colname.A <- pij$Ai
			if (!(colname.A %in% colnames(A2l2B))) {
				colname.A <- paste(colname.A,"_A",sep="")
			}

			colname.B <- pij$Bi
			if (!(colname.B %in% colnames(A2l2B))) {
				colname.B <- paste(colname.B,"_B",sep="")
			}

			criteria <- which( (A2l2B[colname.A] == codeA) & (A2l2B[colname.B] == codeB ) )
			count <- nrow(A2l2B[criteria,])
			cat("\t=>", count, "\n")
			
			hat.nij[cB,cA] <- count

		}

	}

    hat.pij <- hat.nij / nrow(pop$links)
	
    #print(pij$data-hat.pij)

	list( 
			hat.nij=hat.nij, 
			hat.pij=hat.pij, 
			mse.orig=mean( (pij$data-hat.pij)^2 ) , 
			mse.target=mean( (mix.pij-hat.pij)^2 ) 
			)
}

measure.ci_fi <- function(target.ni, orig.fi, colname, dico, pop) {

    hat.ni <- c()
    for (codeA in dico$encoding[[colname]]) {
        #codeA <- dico$encoding[[colname]][[cA]]
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

# measures the distribution of degrees from the population
measure.pdn <- function(dico, pdn.orig, sample) {


	hat.pdn <- pdn.orig$data
	for (codeA in dico$encoding[[pdn.orig$attributes]]) {
		#codeA <- dico$encoding[[colname]][[cA]]

		total <- 0
		for (x in 1:nrow(pdn.orig$data)) {
			degree <- x-1
			cat("\tmeasuring how many entities have",pdn.orig$attributes,"=",codeA,"and current degree=",degree,"")
			count <- nrow(sample[which( (sample[pdn.orig$attributes]==codeA) & (sample$current.degree==degree) ),])
			cat("\t=>", count, "\n")
			hat.pdn[x,codeA] <- count 
			total <- total + count 
		}
		# print(hat.pdn[,codeA])
		hat.pdn[,codeA] <- hat.pdn[,codeA] / total 
		# print(hat.pdn[,codeA])
		
	}

	list(
		hat.pd=hat.pdn,
		mse.orig=mean( (pdn.orig$data-hat.pdn)^2 ),
		count.over.degree=nrow(sample[which( sample$current.degree>sample$target.degree),]),
		count.under.degree=nrow(sample[which( sample$current.degree<sample$target.degree),])
		)
}
measure.population <- function(case, pop, sample.A, sample.B, pij) {

	# merge the datasets	
	A2l <- merge(pop$A, pop$links, by="id.A", suffixes=c("_A", "_l"))
	A2l2B <- merge(A2l, pop$B, by="id.B", suffixes=c("_A", "_B"))

	cat("total links measured: ", nrow(pop$links),"\n")

	res <- list(
				nL=nrow(A2l2B),
				
				#ni=measure.ni_pi(case$gen$ni, case$stats$pi, case$gen$hat.pi, case$inputs$pij$Ai, case$inputs$sample.A$dictionnary, A2l2B),
				#nj=measure.ni_pi(case$gen$nj, case$stats$pj, case$gen$hat.pj, case$inputs$pij$Bi, case$inputs$sample.B$dictionnary, A2l2B),
				
                fi=measure.ci_fi(case$gen$hat.ni, case$stats$fi, case$inputs$pij$Ai, case$inputs$sample.A$dictionnary, pop$A),
                fj=measure.ci_fi(case$gen$hat.nj, case$stats$fj, case$inputs$pij$Bi, case$inputs$sample.B$dictionnary, pop$B),
                
				# TODO add mix
				pdi=measure.pdn(case$inputs$sample.A$dictionnary, case$inputs$pdi, pop$A),
				pdj=measure.pdn(case$inputs$sample.B$dictionnary, case$inputs$pdj, pop$B),

				pij=measure.pij(pop, sample.A, sample.B, pij, case$gen$hat.pij, A2l2B)
				)

	#print(res)

	res

}
