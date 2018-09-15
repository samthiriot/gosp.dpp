[![Travis-CI Build Status](https://travis-ci.org/samthiriot/gosp.dpp.svg?branch=master)](https://travis-ci.org/samthiriot/gosp.dpp)

## gosp.dpp 

Generation of Synthetic Populations: Direct Probabilistic Pairing


## user install

From R, you can install it in 2 steps only:

install the devtools package 

    install.packages("devtools")
	
then use it to install the package from github

	library(devtools)
	install_github("samthiriot/gosp.dpp")

you would also better install the optional dependancies:

    install.packages(c("ggplot2", "gridExtra", "igraph", "mipfp"))


## first steps

To now how to start, you might have a look to the vignettes:
* simple example of generation of a synthetic populations made of dwellings and households: http://htmlpreview.github.io/?https://raw.githubusercontent.com/samthiriot/gosp.dpp/master/inst/doc/compose_dwellings_and_households.html


## developer install

clone the repository

from inside the clone

    install.packages(c("devtools", "rhub", "knitr"))
	library(devtools)
	devtools::install()
	devtools::load_all()

enjoy!


## releasing

generate the data

	library(devtools)
	devtools::install()
	source("data-raw/dwellings_households.R")

run local tests

	library(devtools)
	devtools::test()

check the package locally

	library(devtools)
	devtools::check(manual=TRUE, vignette=TRUE)

check on various platforms

before release, we test the package on Windows, MacOSX and Linux

	library(rhub)
	check()

update comments for CRAN: if relevant, update the comments in cran-comments.md

