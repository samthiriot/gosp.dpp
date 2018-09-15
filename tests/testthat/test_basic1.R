
# tests based on case1 


library(gosp.dpp)

context("tests on generation with one attribute only (case 1)")

test_that("case 1 properly loaded", {

	data(cas1)

	expect_is(cas1$sample.A, "dpp_sample")
	expect_is(cas1$sample.B, "dpp_sample")

	expect_is(cas1$pdi, "dpp_degree_cpt")
	expect_is(cas1$pdj, "dpp_degree_cpt")

	expect_is(cas1$pij, "dpp_matching_probas")

})


test_that("case 1 preparation", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	expect_is(prepared, "dpp_prepared")
	
	# TODO more tests
  
})

test_that("resolution with nA, phi.A, delta.B, phi.B and nB", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
							nA=50000,nB=40000, 
							nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0,
							verbose=FALSE
							)

	expect_is(solved, "dpp_resolved")

	expect_equal(50000, solved$gen$hat.nA)
	expect_equal(40000, solved$gen$hat.nB)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)

	expect_equal(prepared$inputs$pdj$data, solved$gen$hat.pdj, tolerance=1e-5)
  
  	# TODO finish this block and reuse it on every test
  	# the actual generation is long and uses memory and CPU; so we avoid it on CRAN
  	skip_on_cran()
  	case <- matching.generate(solved, cas1$sample.A, cas1$sample.B, verbose=FALSE)

  	expect_equal(nrow(case$pop$A), 50000)
  	expect_equal(nrow(case$pop$B), 40000)

  	# ensure the weight column in not present
  	expect_false(cas1$sample.A$dictionary$colname.weight %in% colnames(case$pop$A))
  	expect_false(cas1$sample.B$dictionary$colname.weight %in% colnames(case$pop$B))

  	# TODO more tests

})



# SMALL CHAINS
# small chains which can be solved step by step with one unique hypothesis

test_that("constraints: nA, phi.A, phi.B", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
							nA=50000,nB=40000, 
							nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=1,
							verbose=FALSE) 

	expect_is(solved, "dpp_resolved")

	expect_equal(50000, solved$gen$hat.nA)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	
  
})


test_that("constraints: nA, phi.A, delta.B, phi.B, nu.B", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
								nA=50000, nB=40000, 
								nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0)
	expect_is(solved, "dpp_resolved")

	expect_equal(50000, solved$gen$hat.nA)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	
  
})


test_that("constraints: phi.A,phi.B, nu.B", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
								nA=50000, nB=40000, 
								nu.A=1, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=0
								)

	expect_is(solved, "dpp_resolved")

	expect_equal(40000, solved$gen$hat.nB)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	
})

test_that("constraints: phi.A, phi.B, nu.B", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
								nA=50000, nB=40000, 
								nu.A=1, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=0
								)

	expect_is(solved, "dpp_resolved")

	expect_equal(40000, solved$gen$hat.nB)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	
})

test_that("constraints: phi.A, delta.A (free on matching and B)", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
								nA= 50000, nB=40000, 
								nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=1, phi.B=1, nu.B=1,
								verbose=FALSE
								)

	expect_is(solved, "dpp_resolved")

	expect_equal(50000, solved$gen$hat.nA)
	
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$pdi, unname(solved$gen$hat.pdi$data), tolerance=1e-5)
	
})

test_that("constraints: gamma (free on A and B)", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	# solved <- matching.solve(prepared, nA=50000, nB=40000, nu.A=1, phi.A=1, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1, verbose=T)
	
	solved <- matching.solve(prepared, 
		nA=50000, nB=40000, 
		nu.A=1, phi.A=1, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1,
		verbose=F)
	
	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	
	expect_equal(as.matrix(cas1$pij$data), solved$gen$hat.pij, tolerance=0.1)
	
})



# LONG CHAINS
# test the resolution of long chains for which several hypothesis have to be piled to be solved.

context("tests on case 1 with the exploration of several hypothesis")


test_that("constraints: A free (case 1) with equal weights", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
							nA=50000,nB=40000, 
							nu.A=1, phi.A=1, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0, 
							verbose=FALSE)

	#print(solved)

	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	expect_false(is.null(solved$gen$hat.dj))
	expect_false(is.null(solved$gen$hat.ci))
	expect_false(is.null(solved$gen$hat.cj))
	expect_false(is.null(solved$gen$hat.pij))
	expect_false(is.null(solved$gen$hat.nij))
	expect_false(is.null(solved$gen$hat.ndi))
	expect_false(is.null(solved$gen$hat.ndj))

	# based on how the algo is defined and how the resolution is weighted,
	# we expect the selected solution to be one with hat.fi=fi and hat.di=di
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(prepared$stats$pdi, unname(solved$gen$hat.pdi$data), tolerance=1e-5)
	
	# we know these elements are wrong:
	# ... the error has to be transfered into hat.pij
	expect_false(all(TRUE==all.equal(cas1$pij$data, solved$gen$hat.pij, tolerance=1e-5)))

	# this has to be because the weight = 0
	expect_equal(unname(prepared$inputs$dj), unname(solved$gen$hat.dj), tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	expect_equal(40000, solved$gen$hat.nB, tolerance=1)

})


test_that("constraints: A free (case 1) weighting nu.A", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
							nA=50000,nB=40000, 
							nu.A=1, phi.A=10, delta.A=10, gamma=10, delta.B=0, phi.B=0, nu.B=0, 
							verbose=FALSE)

	#print(solved)

	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	expect_false(is.null(solved$gen$hat.dj))
	expect_false(is.null(solved$gen$hat.ci))
	expect_false(is.null(solved$gen$hat.cj))
	expect_false(is.null(solved$gen$hat.pij))
	expect_false(is.null(solved$gen$hat.nij))
	expect_false(is.null(solved$gen$hat.ndi))
	expect_false(is.null(solved$gen$hat.ndj))

	# based on how the algo is defined and how the resolution is weighted,
	# we expect the selected solution to be one with hat.nA=nA, hat.fi=fi 
	expect_equal(prepared$stats$fi, solved$gen$hat.fi, tolerance=1e-5)
	expect_equal(50000, solved$gen$hat.nA, tolerance=1)
	# also pdi is respected, as a side effect (fi is respected, and pij is free, so pdi can be preserved)
	expect_equal(prepared$stats$pdi, unname(solved$gen$hat.pdi$data), tolerance=1e-5)

	# this cannot be true, the error has to be reported into di and pij
	expect_false(all(TRUE == all.equal(cas1$pij$data, solved$gen$hat.pij, tolerance=1e-5)))

	# this has to be because the weight = 0
	expect_equal(unname(prepared$inputs$dj), unname(solved$gen$hat.dj), tolerance=1e-5)
	expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-5)
	expect_equal(40000, solved$gen$hat.nB, tolerance=1)

})


test_that("constraints: nothing (totally free - long chain)", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	solved <- matching.solve(prepared, 
							nA=50000,nB=40000, 
							nu.A=1, phi.A=1, delta.A=1, gamma=1, delta.B=1, phi.B=1, nu.B=1,
							verbose=FALSE)

	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	expect_false(is.null(solved$gen$hat.dj))
	expect_false(is.null(solved$gen$hat.ci))
	expect_false(is.null(solved$gen$hat.cj))
	expect_false(is.null(solved$gen$hat.pij))
	expect_false(is.null(solved$gen$hat.nij))
	expect_false(is.null(solved$gen$hat.ndi))
	expect_false(is.null(solved$gen$hat.ndj))

})

# SMALL COUNTS
context("tests on case 1 with small sizes")

{
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
	for (factor in c(5000,
		#2000,
		1000,500,333,
		300,
		200,101,100
		,50
		)) {

		nA <- 5*factor
		nB <- 4*factor

		test_that(paste("resolution with small values for nA=",nA," and nB=",nB,"(test factor ",factor,")",sep=""), {
				
			#cat("test with size nA=",nA," and nB=",nB,"\n",sep="")

			solved <- matching.solve(prepared, 
										nA=nA,nB=nB, 
										nu.A=0, phi.A=1, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0,
										verbose=FALSE)

			expect_is(solved, "dpp_resolved")

			# we are quiet tolerant here: 1% 
			expect_equal(nA, solved$gen$hat.nA, tolerance=0.01*nA)
			expect_equal(nB, solved$gen$hat.nB, tolerance=0.01*nB)

			# very tolerant...
			expect_equal(prepared$stats$fj, solved$gen$hat.fj, tolerance=1e-2)
			expect_equal(prepared$stats$pdj, unname(solved$gen$hat.pdj$data), tolerance=1e-5)

		
		})
	}
}


context("tests on case 1 with zero cells")

# ZERO CELLS
# ensure expected failures do fail

test_that("constraints: pdi with zero (p(di=0)=1.0)", {
	
	data(cas1)

	cas1.zero.di <- cas1 
	cas1.zero.di$pdi <- create_degree_probabilities_table(
								probabilities=data.frame(
				                    'surface=1'=c(0.2, 0.8, 0, 0, 0),
				                    'surface=2'=c(1.0, 0.0, 0.0, 0, 0),
				                    'surface=3'=c(0.05, 0.8, 0.1, 0.05, 0),
									check.names=FALSE
				                    )
								)

	prepared <- matching.prepare(cas1.zero.di$sample.A, cas1.zero.di$sample.B, cas1.zero.di$pdi, cas1.zero.di$pdj, cas1.zero.di$pij)

	solved <- matching.solve(prepared, 
		nA=50000, nB=40000, 
		nu.A=1, phi.A=1, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1,
		verbose=FALSE
		)
	
	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	
})


test_that("constraints: pdj with zero (p(dj=0)=1.0)", {
	
	data(cas1)

	cas1.zero.dj <- cas1 
	cas1.zero.dj$pdj <- create_degree_probabilities_table(
			                probabilities=data.frame(
			                    'size=1'=c(0, 1),
			                    'size=2'=c(0, 1),
			                    'size=3'=c(1, 0),
			                    'size=4'=c(0, 1),
								check.names=FALSE
								)
			                )

	prepared <- matching.prepare(cas1.zero.dj$sample.A, cas1.zero.dj$sample.B, cas1.zero.dj$pdi, cas1.zero.dj$pdj, cas1.zero.dj$pij)

	solved <- matching.solve(prepared, 
		nA=50000, nB=40000, 
		nu.A=0, phi.A=0, delta.A=0, gamma=0, delta.B=0, phi.B=1, nu.B=1,
		verbose=FALSE
		)
	
	expect_is(solved, "dpp_resolved")

	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	
})

test_that("constraints: pij with zero", {
	
	data(cas1)

	# in the example case, replace one of the matching probabilities by a zero
	cas1.zero.dj <- cas1
	cas1.zero.dj$pij <- create_matching_probabilities_table(
                data.frame(
                    'surface=1'=c(0.2, 0.1, 0.05, 0.025),
                    'surface=2'=c(0.0375, 0.125, 0.0, 0.05),
                    'surface=3'=c(0.0125, 0.025, 0.2, 0.175),
                    row.names=c("size=1", "size=2", "size=3", "size=4"),
                    check.names=FALSE
                    )
                )

	prepared <- matching.prepare(cas1.zero.dj$sample.A, cas1.zero.dj$sample.B, cas1.zero.dj$pdi, cas1.zero.dj$pdj, cas1.zero.dj$pij)

	solved <- matching.solve(prepared, 
		nA=50000, nB=40000, 
		nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=0, phi.B=1, nu.B=1,
		verbose=FALSE
		)

	#print(solved)
	
	expect_is(solved, "dpp_resolved")

	expect_equal(0, solved$gen$hat.pij["size=3","surface=2"])
	expect_false(is.null(solved$gen$hat.nA))
	expect_false(is.null(solved$gen$hat.nB))
	expect_false(is.null(solved$gen$hat.di))
	
})


context("tests on case 1 with expected failures")

# FAILURES
# ensure expected failures do fail

test_that("constraints: phi.A, gamma (too constrained)", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

   	expect_error(do.call(
   		matching.solve,
   		list(prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=0, phi.B=0, nu.B=0)),
   		"The case is too constrained to be solved.*")
	
})

test_that("constraints: nu.A, delta.A, gamma, delta.B, phi.B, nu.B (too constrained)", {
	
	data(cas1)

	prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

   	expect_error(do.call(
   		matching.solve,
   		list(prepared, nA=50000,nB=40000, nu.A=0, phi.A=1, delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=0)),
   		"The case is too constrained to be solved.*")
	
})

