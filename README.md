[![Travis-CI Build Status](https://travis-ci.org/samthiriot/gosp.dpp.svg?branch=master)](https://travis-ci.org/samthiriot/gosp.dpp)

# gosp.dpp 

Generation of Synthetic Populations: Direct Probabilistic Pairing


# user install

From R, you can install it in 2 steps only:

install the devtools package 

    install.packages("devtools")
	
then use it to install the package from github

	install_github("samthiriot/gosp.dpp")


# developer install

clone the repository

    install.packages("devtools")
	library(devtools)
	devtools::install()
	devtools::load_all()

enjoy!


# releasing

generate the data

	library(devtools)
	devtools::install()
	source("data-raw/cas1.R")
	source("data-raw/case2.R")

run local tests

	library(devtools)
	devtools::test()

check the package locally

	library(devtools)
	devtools::check(manual=TRUE)

check on various platforms

before release, we test the package on Windows, MacOSX and Linux

	library(rhub)
	check()

update comments for CRAN: if relevant, update the comments in cran-comments.md


## TODO

actual upload:
