
# tests on 

library(gosp.dpp)


# TODO make igraph optional

library("igraph")

# import::from(igraph, is.directed)
# import::from(igraph, as.igraph)


# if ("igraph" %in% installed.packages()) {

# } else {
# 	skip("igraph not installed")
# }


data(cas1)


context("tests on igraph conversion")

test_that("export a population as an igraph (without attributes)", {

	# prepare the case
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	# resolve the case
	disc <- matching.arbitrate(case.prepared, nA=5000,nB=4000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0)

	# generate
	generated <- matching.generate(disc, cas1$sample.A, cas1$sample.B)

	# convert to an igraph graph 
	g <- as.igraph(generated$pop)

	# is it the right type ?
	expect_is(g, "igraph")

	expect_true(is.directed(g))

	# is it the right size ? 
	# ... in the count of links (edges)
	expect_equal(gsize(g), nrow(generated$pop$links))
	# ... in the count of entities (vertices)
	expect_equal(vcount(g), nrow(generated$pop$A) + nrow(generated$pop$B))
	
	# does it has attributes ?
	# ... by default it should not have attributes
	expect_equal(0, length(vertex_attr_names(g)))

})


test_that("export a population as an igraph (with attributes)", {

	# prepare the case
	case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)

	# resolve the case
	disc <- matching.arbitrate(case.prepared, nA=5000,nB=4000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=0, phi.B=0, nu.B=0)

	# generate
	generated <- matching.generate(disc, cas1$sample.A, cas1$sample.B)

	# convert to an igraph graph 
	g <- as.igraph(generated$pop, with.attributes=TRUE)

	# is it the right type ?
	expect_is(g, "igraph")

	expect_true(is.directed(g))

	# is it the right size ? 
	# ... in the count of links (edges)
	expect_equal(gsize(g), nrow(generated$pop$links))
	# ... in the count of entities (vertices)
	expect_equal(vcount(g), nrow(generated$pop$A) + nrow(generated$pop$B))
	
	# does it has attributes ?
	# ... it should have the right count of attributes
	expect_equal(length(union(colnames(generated$pop$A),colnames(generated$pop$B))), length(vertex_attr_names(g)))

	# let's check a few attributes of a few individuals
	# ... of A (test 50)
	for (idx in sample(1:nrow(generated$pop$A),20)) { 
		id <- generated$pop$A[idx,"id.A"]
		for (name in colnames(generated$pop$A)) {

			expect_equal(
				generated$pop$A[idx,name],
				vertex_attr(g, name, id)
				)

		}
	}
	# ... of B (test 50)
	for (idx in sample(1:nrow(generated$pop$B),20)) { 
		id <- generated$pop$B[idx,"id.B"]
		for (name in colnames(generated$pop$B)) {

			expect_equal(
				generated$pop$B[idx,name],
				vertex_attr(g, name, id)
				)

		}
	}

	
})

