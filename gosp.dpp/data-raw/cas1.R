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
                            'small'=1,
                            'medium'=2,
                            'large'=3
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
                            '1 person'=1,
                            '2 persons'=2,
                            '3 persons'=3,
                            '4 and more'=4
                            )
                       ),
                weight.colname="weight"
                )


# pdi
cas1$pdi <- create_degree_probabilities_table(
                probabilities=data.frame(
                    'small'=c(0.2, 0.8, 0, 0, 0),
                    'medium'=c(0.15, 0.8, 0.05, 0, 0),
                    'large'=c(0.05, 0.8, 0.1, 0.05, 0)
                    ),
                attributes.names=c("surface")
                )


# pdj
cas1$pdj <- create_degree_probabilities_table(
                probabilities=data.frame(
                    '1 person'=c(0, 1),
                    '2 persons'=c(0, 1),
                    '3 persons'=c(0, 1),
                    '4 and more'=c(0, 1)
                    ),
                attributes.names=c("size")
                )

cas1$pij <- create_matching_probabilities_table(
                data=data.frame(
                    'small'=c(0.2, 0.1, 0.05, 0.025),
                    'medium'=c(0.0375, 0.125, 0.1, 0.05),
                    'large'=c(0.0125, 0.025, 0.1, 0.175),
                    row.names=c("1 person", "2 persons", "3 persons", "4 and more")
                    ),
                Ai=c("surface"),
                Bi=c("size")
                )

devtools::use_data(cas1, overwrite=TRUE)