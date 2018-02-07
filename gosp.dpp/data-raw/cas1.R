#
# generates the cas1 data


library(gosp.dpp)


cas1 = list()

# parents
cas1$sample.A <- gosp.dpp::create_sample(
                data=data.frame(
                        read.csv(
                            "data-raw/logements.csv", 
                            sep=";", 
                            dec=",")
                        ),
                encoding = list(
                        'surface'=list(
                            'surface=small'=1,
                            'surface=medium'=2,
                            'surface=large'=3
                            )
                       ),
                weight.colname="weight"
                )



# children
cas1$sample.B <- create_sample(
                data=data.frame(
                        read.csv(
                            "data-raw/foyers.csv", 
                            sep=";", 
                            dec=",")
                        ),
                encoding = list(
                        'size'=list(
                            'size=1'=1,
                            'size=2'=2,
                            'size=3'=3,
                            'size=4'=4
                            )
                       ),
                weight.colname="weight"
                )


# pdi
cas1$pdi <- create_degree_probabilities_table(
                probabilities=data.frame(
                    'surface=1'=c(0.2, 0.8, 0, 0, 0),
                    'surface=2'=c(0.15, 0.8, 0.05, 0, 0),
                    'surface=3'=c(0.05, 0.8, 0.1, 0.05, 0),
                    check.names=FALSE
                    ),
                attributes.names=c("surface")
                )


# pdj
cas1$pdj <- create_degree_probabilities_table(
                probabilities=data.frame(
                    'size=1'=c(0, 1),
                    'size=2'=c(0, 1),
                    'size=3'=c(0, 1),
                    'size=4'=c(0, 1),
                    check.names=FALSE
                    ),
                attributes.names=c("size")
                )

cas1$pij <- create_matching_probabilities_table(
                data=data.frame(
                    'surface=1'=c(0.2, 0.1, 0.05, 0.025),
                    'surface=2'=c(0.0375, 0.125, 0.1, 0.05),
                    'surface=3'=c(0.0125, 0.025, 0.1, 0.175),
                    row.names=c("size=1", "size=2", "size=3", "size=4"),
                    check.names=FALSE
                    ),
                Ai=c("surface"),
                Bi=c("size")
                )

devtools::use_data(cas1, overwrite=TRUE)