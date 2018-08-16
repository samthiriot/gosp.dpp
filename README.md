[![Travis-CI Build Status](https://travis-ci.org/samthiriot/direct-probabilistic-pairing.svg?branch=master)](https://travis-ci.org/samthiriot/direct-probabilistic-pairing)

# gosp.dpp 

Generation of Synthetic Populations: Direct Probabilistic Pairing


# test the git version 

If you want to test or use the development version, you can clone this repository, and:

install the devtools package 

    install.packages("devtools")
	
the first time, from the gosp.dpp directory

	library(devtools)
	devtools::install()
	devtools::load_all()
	source("data-raw/cas1.R")
	source("data-raw/case2.R")

then from the gosp.dpp directory:

    library(devtools)
    devtools::check()


# more info

see http://r-pkgs.had.co.nz/check.html

# releasing:

## generate the data

	library(devtools)
	devtools::install()
	source("data-raw/cas1.R")
	source("data-raw/case2.R")


use devtools checking: http://r-pkgs.had.co.nz/check.html

see https://www.r-project.org/nosvn/pandoc/devtools.html