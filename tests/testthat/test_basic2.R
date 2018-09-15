
# tests based on case1 


library(gosp.dpp)

data(dwellings_households)

context("tests on generation with several attributes (case dwellings households)")

test_that("case 1 properly loaded", {

	expect_is(dwellings_households$sample.A, "dpp_sample")
	expect_is(dwellings_households$sample.B, "dpp_sample")

	expect_is(dwellings_households$pdi, "dpp_degree_cpt")
	expect_is(dwellings_households$pdj, "dpp_degree_cpt")

	expect_is(dwellings_households$pij, "dpp_matching_probas")

})


test_that("case preparation", {
	
	case.prepared <- matching.prepare(
						dwellings_households$sample.A, 
						dwellings_households$sample.B, 
						dwellings_households$pdi, 
						dwellings_households$pdj, 
						dwellings_households$pij
						)
	expect_is(case.prepared, "dpp_prepared")
	
	expect_equal(length(case.prepared$stats$fi), ncol(dwellings_households$pdi$data))
	expect_equal(length(case.prepared$stats$fj), ncol(dwellings_households$pdj$data))

	# TODO more tests
  
})

test_that("resolution with nA, phi.A, delta.B, phi.B and nB", {
	
	case.prepared <- matching.prepare(
						dwellings_households$sample.A, 
						dwellings_households$sample.B, 
						dwellings_households$pdi, 
						dwellings_households$pdj, 
						dwellings_households$pij)
	factor <- 1

	disc <- matching.solve(case.prepared, 
						nA=50000*factor,nB=40000*factor, 
						nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0,
						verbose=FALSE)

	expect_is(disc, "dpp_resolved")

	expect_equal(disc$gen$hat.nA, 50000*factor)
	expect_equal(disc$gen$hat.nB, 40000*factor)
	
	expect_equal(case.prepared$stats$fi, disc$gen$hat.fi, tolerance=3e-5)
	expect_equal(case.prepared$stats$fj, disc$gen$hat.fj, tolerance=3e-5)

	expect_equal(case.prepared$inputs$pdj$data, disc$gen$hat.pdj, tolerance=1e-5)
  
  	# TODO finish this block and reuse it on every test
  	# the actual generation is long and uses memory and CPU; so we avoid it on CRAN
  	skip_on_cran()
  	case <- matching.generate(disc, dwellings_households$sample.A, dwellings_households$sample.B, verbose=FALSE)

  	expect_equal(nrow(case$pop$A), 50000*factor)
  	expect_equal(nrow(case$pop$B), 40000*factor)

  	# ensure the weight column in not present
  	expect_false(dwellings_households$sample.A$dictionary$colname.weight %in% colnames(case$pop$A))
  	expect_false(dwellings_households$sample.B$dictionary$colname.weight %in% colnames(case$pop$B))

  	# TODO more tests

})
