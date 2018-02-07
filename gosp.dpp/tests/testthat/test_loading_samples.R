
# tests on 

library(gosp.dpp)


context("tests on loading samples")


test_that("samples: detection of one attribute for a pdi", {

	pdi <- create_degree_probabilities_table(
                probabilities=data.frame(
                    'surface=1'=c(0.2, 0.8, 0, 0, 0),
                    'surface=2'=c(0.15, 0.8, 0.05, 0, 0),
                    'surface=3'=c(0.05, 0.8, 0.1, 0.05, 0),
                    check.names=FALSE
                    )
                )
	
	# standard checks
	expect_is(pdi, "dpp_degree_cpt")

	# we should now have a weight column 
	expect_equivalent(pdi$attributes, c("surface"))

})


test_that("samples: detection of several attributes for a pdi", {

	pdi <- create_degree_probabilities_table(
                probabilities=data.frame(
                    'surface=1,other=1'=c(0.2, 0.8, 0, 0, 0),
                    'surface=2,other=1'=c(0.15, 0.8, 0.05, 0, 0),
                    'surface=3,other=1'=c(0.05, 0.8, 0.1, 0.05, 0),
                    check.names=FALSE
                    )
                )
	
	# standard checks
	expect_is(pdi, "dpp_degree_cpt")

	# we should now have a weight column 
	expect_equivalent(pdi$attributes, c("surface","other"))

})


test_that("samples: error when malformed attributes for a pdi", {

   	expect_error(do.call(
   		create_degree_probabilities_table,
   		list(data.frame(
                    'surface=1,other1'=c(0.2, 0.8, 0, 0, 0),
                    'surface=2,other=1'=c(0.15, 0.8, 0.05, 0, 0),
                    'surface=3,other=1'=c(0.05, 0.8, 0.1, 0.05, 0),
                    check.names=FALSE
                    ))),
   		"invalid name.*")

   	expect_error(do.call(
   		create_degree_probabilities_table,
   		list(data.frame(
                    'surface=1'=c(0.2, 0.8, 0, 0, 0),
                    'surface=2,other=1'=c(0.15, 0.8, 0.05, 0, 0),
                    'surface=3,other=1'=c(0.05, 0.8, 0.1, 0.05, 0),
                    check.names=FALSE
                    ))),
   		"all the column names should contain.*")

})


test_that("samples: load sample with weights and dictionary", {

	# create a useless dataframe without weight 
	df <- data.frame(
				useless=sample(1:2, size=100, replace=TRUE),
				myweights=runif(100)
				)
	s <- create_sample(
			data=df,
			encoding=list("useless"=list("male"=1,"female"=2)),
			weight.colname="myweights"
			)
	
	# standard checks
	expect_is(s, "dpp_sample")

	# we should now have a weight column 
	expect_equal(s$dictionary$colname.weight, "myweights")

})


test_that("samples: automatic creation of weights", {

	# create a useless dataframe without weight 
	df <- data.frame(useless=sample(1:2, size=100, replace=TRUE))
	s <- create_sample(data=df,encoding=list("useless"=list("male"=1,"female"=2)))
	
	# standard checks
	expect_is(s, "dpp_sample")

	# we should now have a weight column 
	expect_false(is.null(s$dictionary$colname.weight))
	# whose name is a string 
	expect_is(s$dictionary$colname.weight, "character")

	# this column should exist in the sample
	expect_false(is.null(s$sample[,s$dictionary$colname.weight]))

	# this column should only contain "1"
	expect_true(all(s$sample[,s$dictionary$colname.weight]==1))

})

test_that("samples: automatic creation of encoding (one int attribute)", {

	# create a useless dataframe without weight 
	df <- data.frame(
				useless=sample(1:2, size=100, replace=TRUE),
				somedouble=runif(100)
				)
	s <- create_sample(data=df)
	
	# standard checks
	expect_is(s, "dpp_sample")

	# we should now have a weight column 
	expect_false(is.null(s$dictionary))
	expect_false(is.null(s$dictionary$encoding))
	expect_false(is.null(s$dictionary$decoding))

	# expect only one column to be created (just for "useless")
	expect_equal(length(names(s$dictionary$encoding)),1)

	# expect those columns to exist:
	expect_false(is.null(s$dictionary$encoding$useless))
	expect_equal(s$dictionary$encoding$useless[["1"]],1)
	expect_equal(s$dictionary$encoding$useless[["2"]],2)
	expect_equal(length(names(s$dictionary$encoding$useless)),2)

	# expect the weight column not to have been created 
	expect_true(is.null(s$dictionary$encoding[[s$dictionary$colname.weight]]))

})


test_that("samples: fail if the weight column is missing", {

	# create a useless dataframe without weight 
	df <- data.frame(useless=sample(1:2, size=100, replace=TRUE))
	
   	expect_error(do.call(
   		create_sample,
   		list(
   				data=df,
   				encoding=list("useless"=list("male"=1,"female"=2)), 
   				# this is the faulty parameter
   				weight.colname="notthere")),
   		"There is no column weight.colname.*"
		)

})


test_that("samples: automatic creation of missing parts in dictionary", {

	# create a useless dataframe without weight 
	df <- data.frame(	gender=sample(1:2, size=100, replace=TRUE),
						useless_more=sample(1:5, size=100, replace=TRUE)
						)
	s <- create_sample(data=df,encoding=list("gender"=list("male"=1,"female"=2)))

	expect_false(is.null(s$dictionary$encoding$gender))
	expect_false(is.null(s$dictionary$encoding$useless_more))
	expect_equal(length(s$dictionary$encoding$useless_more), 5)

})


test_that("samples: automatic creation of a missing dictionary", {

	# create a useless dataframe without weight 
	df <- data.frame(	gender=sample(1:2, size=100, replace=TRUE),
						useless_more=sample(1:5, size=100, replace=TRUE)
						)
	s <- create_sample(data=df)

	expect_false(is.null(s$dictionary$encoding$gender))
	expect_false(is.null(s$dictionary$encoding$useless_more))
	expect_equal(length(s$dictionary$encoding$useless_more), 5)

})
