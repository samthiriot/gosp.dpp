# direct-probabilistic-pairing
experiments on generation of structured synthetic populations

# test the git version 

If you want to test or use the development version, you can clone this repository, and:

install the devtools package 

    install.packages("devtools")
	install.packages("dplyr")
	install.packages("igraph")
	
the first time, from the gosp.dpp directory

	library(devtools)
	devtools::install()
	devtools::load_all()
	source("data-raw/cas1.R")
	source("data-raw/case2.R")

No !!!
	devtools::use_data_raw("gosp.dpp")

then from the root of the clone:

    library(devtools)
    devtools::check("gosp.dpp")


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