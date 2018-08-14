
#' Plots the differences between two populations
#' 
#' 
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'

#' Identifies the columns present it two datasets 
#'
#' @param t1 a dataframe
#' @param t2 a dataframe
#' @return a list of strings
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'  
compute_common_columns <- function(t1, t2) {
    common_cols <- intersect(names(t1), names(t2)) 
    common_cols[common_cols != "id"]
}

plot.difference.populations <- function(p1, p2, entitytype="entity?", par=T) {

    cols_to_process <- compute_common_columns(p1, p2)
    if (par) { 
        par(mfrow=c(2,length(cols_to_process)))
    }

    for (n in cols_to_process){
        print(n)

        titleA <- if (n == cols_to_process[1]) paste("reference", entitytype) else NULL
        titleB <- if (n == cols_to_process[1]) paste("generated", entitytype) else NULL
        
        hist(p1[[n]], col="gray", breaks=unique(p1[[n]]), include.lowest=T, main=titleA, xlab=n, freq=F)
        hist(p2[[n]], col="lightblue", breaks=unique(p2[[n]]), include.lowest=T, main=titleB, xlab=n, freq=F)

    }

}

plot.differences <- function(sp, sampleA, sampleB, nameA="A", nameB="B") {

    

    pAref <- sampleA # sp$inputs$sample.A$sample
    pAgen <- sp$pop$A 
    
    pBref <- sampleB # sp$inputs$sample.B$sample
    pBgen <- sp$pop$B 
    
    cols_A <- compute_common_columns(pAref, pAgen)
    cols_B <- compute_common_columns(pBref, pBgen)

    par(mfrow=c(length(cols_A)+length(cols_B),2))

    plot.difference.populations(pAref, pAgen, entitytype=nameA, par=F)
    plot.difference.populations(pBref, pBgen, entitytype=nameB, par=F)

}

#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 geom_histogram
#' 
plot.unconstrained <- function(sp, sampleA, sampleB, nameA="A", nameB="B") {

    # TODO read https://stackoverflow.com/questions/38507729/recast-in-r-with-sumproduct

    pAref <- sampleA # sp$inputs$sample.A$sample
    pAgen <- sp$pop$A 
    
    pBref <- sampleB # sp$inputs$sample.B$sample
    pBgen <- sp$pop$B 
    
    cols_A <- compute_common_columns(pAref, pAgen)
    cols_B <- compute_common_columns(pBref, pBgen)


    plots_list <- list()

    cols_to_process <- cols_A
    for (n in cols_to_process){
        print(n)

        pAref[[n]] <- factor(pAref[[n]])
        p <- ggplot(pAref, aes(pAref[[n]], weight=pAref$weight)) + geom_histogram() + xlab(n)
        plots_list <- append(plots_list, list(p))
    }

    print("list:")
    print(plots_list)
    do.call("grid.arrange", c(plots_list, ncol=1, nrow=length(plots_list)))


}

#' Plots a visual synthese of the errors induced by the generation process. 
#'
#' Plots all the values controlled by the algorithm, and compares their expected value (passed as data or parameter)
#' and their actual value (as measured in the resulting synthetic population).
#' The resulting graphs contain: 
#' a graph showing the relaxation parameters passed to the solving function; 
#' a graph showing the Normalized Root Mean Square Error (NRMSE) for the each control variable 
#' a graph showing the population sizes asked for and generated for populations A and B
#' a graph showing the difference between the input and observed peering probabilities pij
#' two graphs showing the initial and observed frequencies of control variables in both populations A and B
#' two graphs showing the initial and observed average degrees in both populations A and B
#' two graphs showing the difference between the expected and measured distribution of probability of degrees for both A and B
#' 
#' @param x an object returned by the matching.generate method
#' @param sampleA the original sample for population A
#' @param sampleB the original sample for population B
#' @param nameA a meaningfull label for the entity type of population A, such as "dwellings" (default to "A")
#' @param nameB a meaningfull label for the entity type of population B, such as "households" (default to "B")
#' @param colorRef the color to be used to plot values passed as parameters (defaults to "darkgray")
#' @param colorSynthetic the color to be used to plot values measured in the synthetic population (defaults to "blue")
#' @param ... other parameters will be ignored quietly
#' @return nothing
#' 
#' @examples 
#' data(cas1)
#' case.prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#' disc <- matching.arbitrate(case.prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                            delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' sp <- matching.generate(case=disc, sample.A=cas1$sample.A, sample.B=cas1$sample.B, verbose=TRUE)
#' plot(sp, cas1$sample.A$sample, cas1$sample.B$sample, "dwellings", "households")
#'
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab ggtitle scale_fill_gradient2 geom_tile scale_fill_manual
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#'
plot.dpp_result <- function(x, sampleA, sampleB, nameA="A", nameB="B", colorRef="darkgray", colorSynthetic="blue", ...) {

    sp <- x

    # TODO ensure this object is of the right type
    if (class(sp) != "dpp_result") stop("the data to analyze x should be the result of a matching.generate call")
    if (class(sampleA) != "data.frame") stop("sampleA should be a data frame")
    if (class(sampleB) != "data.frame") stop("sampleB should be a data frame")
    
    # plot relaxation parameters
    all_relaxation <- data.frame(
        parameter=c("nA","fi","pdi/di","pij","pdj/dj","fj","nB"),
        relaxation=c(sp$inputs$nu.A, sp$inputs$phi.A, sp$inputs$delta.A, sp$inputs$gamma, sp$inputs$delta.B, sp$inputs$phi.B, sp$inputs$nu.B)
        )
    all_relaxation$parameter <- factor(all_relaxation$parameter, levels=c("nA","fi","pdi/di","pij","pdj/dj","fj","nB"))
    plot_all_relaxation <- ggplot(all_relaxation, aes(factor(parameter), factor(relaxation))) + geom_bar(stat="identity", fill=colorRef)

    # plot errors on each value
    all_errors <- data.frame(
        error=c("nA","fi","pdi","di","pij","dj","pdj","fj","nB"),
        NRMSE=c(sp$gen$nrmse.nA, sp$gen$nrmse.fi, sp$gen$nrmse.pdi, sp$gen$nrmse.di, sp$gen$nrmse.pij, sp$gen$nrmse.dj, sp$gen$nrmse.pdj, sp$gen$nrmse.fj, sp$gen$nrmse.nB)
        )
    all_errors$error <- factor(all_errors$error, levels=c("nA","fi","pdi","di","pij","dj","pdj","fj","nB"))
    plot_all_errors <- ggplot(all_errors, aes(factor(error), factor(NRMSE))) + 
                        geom_bar(stat="identity", fill=colorSynthetic)

    # this scale will be used for various heatmaps    
    scale_gray_blue <- scale_fill_manual(values=c("darkgray","blue"))

    #population_sizes <- data.frame(type=c(paste(nameA, "(nA)"), paste(nameB,"(nB)")), parameter=c(sp$inputs$nA,sp$inputs$nB), synthetic=c(sp$gen$hat.nA,sp$gen$hat.nB))
    population_sizes <- data.frame(
        type=c(paste(nameA, "(nA)"), paste(nameB,"(nB)")), 
        count=c(sp$inputs$nA,sp$inputs$nB,sp$gen$hat.nA,sp$gen$hat.nB), 
        state=c("theoretical","theoretical","synthetic","synthetic"))
    population_sizes$type <- factor(population_sizes$type, c(paste(nameA, "(nA)"), paste(nameB,"(nB)")))
    population_sizes$state <- factor(population_sizes$state, levels=c("theoretical","synthetic"))

    # ggplot(population_sizes, aes(x=type)) + geom_bar(data=parameter, stat="identity", fill=colorRef) 
    # geom_bar(stat="identity", aes(y=parameter, fill=colorRef), position="dodge") +
    plot_population_sizes <- ggplot(population_sizes, aes(x=factor(type), y=factor(count), fill=factor(state))) + 
                        geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue


    # compare expected and actual degree
    degrees_A <- data.frame(
        attributes=c(names(sp$inputs$di),names(sp$gen$hat.di)),
        average.degree=c(sp$inputs$di,sp$gen$hat.di),
        state=c(rep("theoretical",length(sp$inputs$di)), rep("synthetic",length(sp$gen$hat.di)))
        )
    degrees_A$state <- factor(degrees_A$state, levels=c("theoretical","synthetic"))
    plot_degrees_A <- ggplot(degrees_A, aes(x=factor(attributes), y=factor(average.degree), fill=factor(state))) + 
                        geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ylab("average degree") + 
                        ggtitle(paste("average degree",nameA))

    degrees_B <- data.frame(
        attributes=c(names(sp$inputs$dj),names(sp$gen$hat.dj)),
        average.degree=c(sp$inputs$dj,sp$gen$hat.dj),
        state=c(rep("theoretical",length(sp$inputs$dj)), rep("synthetic",length(sp$gen$hat.dj)))
        )
    degrees_B$state <- factor(degrees_B$state, levels=c("theoretical","synthetic"))
    plot_degrees_B <- ggplot(degrees_B, aes(x=factor(attributes), y=factor(average.degree), fill=factor(state))) + 
                        geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ylab("average degree") + 
                        ggtitle(paste("average degree",nameB))


    # compare expected and actual frequencies
    frequencies_A <- data.frame(
        attributes=c(names(sp$stats$fi),names(sp$gen$hat.fi)),
        frequency=c(sp$stats$fi,sp$gen$hat.fi),
        state=c(rep("theoretical",length(sp$stats$fi)), rep("synthetic",length(sp$gen$hat.fi)))
        ) 
    frequencies_A$state <- factor(frequencies_A$state, levels=c("theoretical","synthetic"))
    plot_frequencies_A <- ggplot(frequencies_A, aes(x=factor(attributes), y=factor(frequency), fill=factor(state))) + 
                        geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ylab("freq") + 
                        ggtitle(paste("frequencies",nameA))

    frequencies_B <- data.frame(
        attributes=c(names(sp$stats$fj),names(sp$gen$hat.fj)),
        frequency=c(sp$stats$fj,sp$gen$hat.fj),
        state=c(rep("theoretical",length(sp$stats$fj)), rep("synthetic",length(sp$gen$hat.fj)))
        ) 
    frequencies_B$state <- factor(frequencies_B$state, levels=c("theoretical","synthetic"))
    plot_frequencies_B <- ggplot(frequencies_B, aes(x=factor(attributes), y=factor(frequency), fill=factor(state))) + 
                        geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ylab("freq") + 
                        ggtitle(paste("frequencies",nameB))


    heat_map_gradient <- scale_fill_gradient2(limits=c(-1,1)) # , trans="log"

    # plot pdi
    diff_pdi <- sp$inputs$pdi$data - sp$gen$hat.pdi
    diff_pdi$degree <- factor(seq(0,nrow(diff_pdi)-1))
    data_hm_pdi <- melt(diff_pdi)
    plot_pdi <- ggplot(data_hm_pdi, aes(factor(variable), factor(degree))) + 
                        geom_tile(aes(fill=value)) + 
                        heat_map_gradient + 
                        ggtitle(paste("difference degree for",nameA))

    diff_pdj <- sp$inputs$pdj$data - sp$gen$hat.pdj
    diff_pdj$degree <- factor(seq(0,nrow(diff_pdj)-1))
    data_hm_pdj <- melt(diff_pdj)
    plot_pdj <- ggplot(data_hm_pdj, aes(factor(variable), factor(degree))) + 
                        geom_tile(aes(fill=value)) + 
                        heat_map_gradient + 
                        ggtitle(paste("difference degree for",nameB))
    
    # plot pij
    diff_pij <- sp$inputs$pij$data - sp$gen$hat.pij
    diff_pij$attributesB <- row.names(diff_pij)
    data_hm_pij <- melt(diff_pij)
    plot_pij <- ggplot(data_hm_pij, aes(factor(variable), factor(attributesB))) + 
                        geom_tile(aes(fill=value)) + 
                        heat_map_gradient + 
                        ggtitle(paste("difference peering"))

    # actual plot over a grid
    grid.arrange(
        plot_all_relaxation, plot_all_errors, 
        plot_population_sizes, plot_pij,
        plot_frequencies_A, plot_frequencies_B,
        plot_degrees_A, plot_degrees_B,
        plot_pdi, plot_pdj,
        ncol=2, nrow=5)

}