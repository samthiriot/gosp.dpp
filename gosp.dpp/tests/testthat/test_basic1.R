
# tests based on case1 


library(gosp.dpp)
data(cas1)

context("tests on case 1")

test_that("case 1 properly loaded", {

	expect_is(cas1$sample.A, "dpp_sample")
	expect_is(cas1$sample.B, "dpp_sample")

	expect_is(cas1$pdi, "dpp_degree_cpt")
	expect_is(cas1$pdj, "dpp_degree_cpt")

	expect_is(cas1$pij, "dpp_matching_probas")

})


test_that("case 1 preparation", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	expect_is(case.prepared, "dpp_prepared")
	
  
})

test_that("resolution with nA, phi.A, delta.B, phi.B and nB", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0)

	expect_is(disc, "dpp_resolved")

	expect_equal(50000, disc$gen$hat.nA)
	expect_equal(40000, disc$gen$hat.nB)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)

	expect_equal(case.prepared$inputs$pdj$data, disc$gen$hat.pdj, tolerance=1e-5)
  
})



# SMALL CHAINS
# small chains which can be solved step by step with one unique hypothesis

test_that("constraints: nA, phi.A, phi.B", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=1)

	expect_is(disc, "dpp_resolved")

	expect_equal(50000, disc$gen$hat.nA)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)
	
  
})


test_that("constraints: nA, phi.A, delta.B, phi.B, nu.B", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0)

	expect_is(disc, "dpp_resolved")

	expect_equal(50000, disc$gen$hat.nA)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)
	
  
})


test_that("constraints: phi.A,phi.B, nu.B", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=1, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=0)

	expect_is(disc, "dpp_resolved")

	expect_equal(40000, disc$gen$hat.nB)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)
	
})

test_that("constraints: phi.A, phi.B, nu.B", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=1, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=0)

	expect_is(disc, "dpp_resolved")

	expect_equal(40000, disc$gen$hat.nB)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)
	
})

test_that("constraints: phi.A, delta.A (free on matching and B)", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000, nB=40000, nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=1, phi.B=1, nu.B=1)

	expect_is(disc, "dpp_resolved")

	expect_equal(50000, disc$gen$hat.nA)
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$pdi, unname(disc$gen$hat.pdi$data), tolerance=1e-5)
	
})

test_that("constraints: phi.A, gamma (free on A and B)", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1)
	
	expect_is(disc, "dpp_resolved")

	expect_false(is.null(disc$gen$hat.nA))
	expect_false(is.null(disc$gen$hat.nB))
	expect_false(is.null(disc$gen$hat.di))
	
	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(as.matrix(cas1$pij$data), disc$gen$hat.pij, tolerance=1e-5)
	
})


# LONG CHAINS
# test the resolution of long chains for which several hypothesis have to be piled to be solved.


test_that("constraints: nothing (totally free - long chain)", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	disc <- matching.arbitrate(case.prepared, nA=50000,nB=40000, nu.A=1, phi.A=1, delta.A=1, gamma=1, delta.B=1, phi.B=1, nu.B=1)

	expect_is(disc, "dpp_resolved")

	# TODO expect nA and nB inbetween

	expect_equal(case.prepared$stats$fi, unname(disc$gen$hat.fi), tolerance=1e-5)
	expect_equal(case.prepared$stats$fj, unname(disc$gen$hat.fj), tolerance=1e-5)
	
  	expect_equal(case.prepared$inputs$di, disc$gen$hat.di, tolerance=1e-5)
	expect_equal(case.prepared$inputs$dj, unname(disc$gen$hat.dj), tolerance=1e-5)
	
})

# FAILURES
# ensure expected failures do fail

test_that("constraints: phi.A, gamma (too constrained)", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

   	expect_error(do.call(
   		matching.arbitrate,
   		list(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=0, phi.B=0, nu.B=0)),
   		"The case is too constrained to be solved.*")
	
})

test_that("constraints: nu.A, delta.A, gamma, delta.B, phi.B, nu.B (too constrained)", {
	
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

   	expect_error(do.call(
   		matching.arbitrate,
   		list(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=1, delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=0)),
   		"The case is too constrained to be solved.*")
	
})

