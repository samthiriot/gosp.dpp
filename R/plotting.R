
#' Identifies the columns present it two datasets 
#'
#' @param t1 a dataframe
#' @param t2 a dataframe
#' @return a list of strings
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'  
#' @keywords internal
#'
compute_common_columns <- function(t1, t2) {
    common_cols <- intersect(names(t1), names(t2)) 
    common_cols[common_cols != "id"]
}

#' Stops if ggplot2 is not installed
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
ensure_presence_ggplot <- function() {

    if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Packages \"ggplot2\" and \"gridExtra\" a required for plotting features. You might install them with install.packages(c(\"ggplot2\",\"gridExtra\"))", call. = FALSE)
    }

}


# TODO
# Plots the differences between two populations
# 
# 
#
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
#' @importFrom graphics hist par
#' @importFrom stats frequency
#' 
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

#' @importFrom graphics par
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

# @importFrom gridExtra grid.arrange
# @importFrom ggplot2 geom_histogram ggplot xlab aes
#
plot.unconstrained <- function(sp, sampleA, sampleB, nameA="A", nameB="B") {

    ensure_presence_ggplot()

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
        p <- ggplot2::ggplot(pAref, ggplot2::aes(pAref[[n]], weight=pAref$weight)) + 
                    ggplot2::geom_histogram() + 
                    ggplot2::xlab(n)
        plots_list <- append(plots_list, list(p))
    }

    do.call("gridExtra::grid.arrange", c(plots_list, ncol=1, nrow=length(plots_list)))

}

#' Adds line breaks to forge axes labels
#' 
#' In case labels refer to several attributes, such as "a=1,b=2", replaces the commas by
#' carriage returns
#'
#' @param labels a vector of strings 
#' @return a vector of strings
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
#' @keywords internal
#'
add_linebreaks_attributes <- function(labels) {
    gsub(",", "\n", labels)
}


#' Plots the relaxation parameters
#' 
#' Plots as a bar chart all the relaxation parameters as they were passed by the user
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_relaxation <- function(sp, colorRef="darkgray") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    # plot relaxation parameters
    all_relaxation <- data.frame(
        parameter=c("nA","fi","pdi/di","pij","pdj/dj","fj","nB"),
        relaxation=c(sp$inputs$nu.A, sp$inputs$phi.A, sp$inputs$delta.A, sp$inputs$gamma, sp$inputs$delta.B, sp$inputs$phi.B, sp$inputs$nu.B)
        )
    all_relaxation$parameter <- factor(all_relaxation$parameter, levels=c("nA","fi","pdi/di","pij","pdj/dj","fj","nB"))
    plot_all_relaxation <- ggplot2::ggplot(all_relaxation, ggplot2::aes(parameter, relaxation)) + 
                                ggplot2::geom_bar(stat="identity", fill=colorRef)

    plot_all_relaxation
}



#' Plots the error measures
#' 
#' Plots as a bar chart all the NRMSE (Normalized Rootsquared Mean Squarred Error) error measures 
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_errors <- function(sp, colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")

    ensure_presence_ggplot()

    all_errors <- data.frame(
        error=c("nA","fi","pdi","di","pij","dj","pdj","fj","nB"),
        NRMSE=c(sp$gen$nrmse.nA, sp$gen$nrmse.fi, sp$gen$nrmse.pdi, sp$gen$nrmse.di, sp$gen$nrmse.pij, sp$gen$nrmse.dj, sp$gen$nrmse.pdj, sp$gen$nrmse.fj, sp$gen$nrmse.nB)
        )
    all_errors$error <- factor(all_errors$error, levels=c("nA","fi","pdi","di","pij","dj","pdj","fj","nB"))
    plot_all_errors <- ggplot2::ggplot(all_errors, ggplot2::aes(error, NRMSE)) + 
                            ggplot2::geom_bar(stat="identity", fill=colorSynthetic)

    plot_all_errors
}

#' Plots the expected and generated population sizes
#' 
#' Plots as a bar chart the sizes expected and generated of the populations
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile scale_fill_manual
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_population_sizes <- function(sp, nameA="A", nameB="B", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")

    ensure_presence_ggplot()
    
    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    population_sizes <- data.frame(
        type=c(paste(nameA, "(nA)"), paste(nameB,"(nB)")), 
        count=c(sp$inputs$nA,sp$inputs$nB,sp$gen$hat.nA,sp$gen$hat.nB), 
        state=c("theoretical","theoretical","synthetic","synthetic"))
    population_sizes$type <- factor(population_sizes$type, c(paste(nameA, "(nA)"), paste(nameB,"(nB)")))
    population_sizes$state <- factor(population_sizes$state, levels=c("theoretical","synthetic"))

    # ggplot(population_sizes, ggplot2::aes(x=type)) + ggplot2::geom_bar(data=parameter, stat="identity", fill=colorRef) 
    # ggplot2::geom_bar(stat="identity", ggplot2::aes(y=parameter, fill=colorRef), position="dodge") +
    res_plot <- ggplot2::ggplot(population_sizes, ggplot2::aes(x=type, y=factor(count), fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue

    res_plot
}


#' Plots average degrees, expected and generated, as bar charts
#' 
#' Plots as a bar charts the avarage degree as they were expected and generated.
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile scale_fill_manual
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_average_degree_A <- function(sp, nameA="A", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    degrees_A <- data.frame(
        attributes=c(names(sp$inputs$di),names(sp$gen$hat.di)), # add_linebreaks_attributes(
        average.degree=c(sp$inputs$di,sp$gen$hat.di),
        state=c(rep("theoretical",length(sp$inputs$di)), rep("synthetic",length(sp$gen$hat.di)))
        )
    degrees_A$state <- factor(degrees_A$state, levels=c("theoretical","synthetic"))
    plot_degrees_A <- ggplot2::ggplot(degrees_A, ggplot2::aes(x=attributes, y=average.degree, fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ggplot2::ylab("average degree") + 
                        ggplot2::ggtitle(paste("average degree",nameA," (NRMSE ",sp$gen$nrmse.di,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_degrees_A
}

#' @rdname plot_average_degree_A
plot_average_degree_B <- function(sp, nameB="B", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    degrees_B <- data.frame(
        attributes=c(names(sp$inputs$dj),names(sp$gen$hat.dj)), # add_linebreaks_attributes(
        average.degree=c(sp$inputs$dj,sp$gen$hat.dj),
        state=c(rep("theoretical",length(sp$inputs$dj)), rep("synthetic",length(sp$gen$hat.dj)))
        )
    degrees_B$state <- factor(degrees_B$state, levels=c("theoretical","synthetic"))
    plot_degrees_B <- ggplot2::ggplot(degrees_B, ggplot2::aes(x=attributes, y=average.degree, fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ggplot2::ylab("average degree") + 
                        ggplot2::ggtitle(paste("average degree",nameB," (NRMSE ",sp$gen$nrmse.dj,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_degrees_B 
}


#' Plots expected and generated frequencies 
#' 
#' Plots as a bar charts the frequencies as they were expected and generated.
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile scale_fill_manual
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_frequencies_A <- function(sp, nameA="A", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    frequencies_A <- data.frame(
        attributes=c(names(sp$stats$fi),names(sp$gen$hat.fi)), #add_linebreaks_attributes(
        frequency=c(sp$stats$fi,sp$gen$hat.fi),
        state=c(rep("theoretical",length(sp$stats$fi)), rep("synthetic",length(sp$gen$hat.fi)))
        ) 
    frequencies_A$state <- factor(frequencies_A$state, levels=c("theoretical","synthetic"))
    res_plot <- ggplot2::ggplot(frequencies_A, ggplot2::aes(x=attributes, y=frequency, fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ggplot2::ylab("freq") + 
                        ggplot2::ggtitle(paste("frequencies",nameA," (NRMSE ",sp$gen$nrmse.fi,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

    res_plot
}
#' @rdname plot_frequencies_A
plot_frequencies_B <- function(sp, nameB="B", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    frequencies_B <- data.frame(
        attributes=c(names(sp$stats$fj),names(sp$gen$hat.fj)), # add_linebreaks_attributes
        frequency=c(sp$stats$fj,sp$gen$hat.fj),
        state=c(rep("theoretical",length(sp$stats$fj)), rep("synthetic",length(sp$gen$hat.fj)))
        ) 
    frequencies_B$state <- factor(frequencies_B$state, levels=c("theoretical","synthetic"))
    res_plot <- ggplot2::ggplot(frequencies_B, ggplot2::aes(x=attributes, y=frequency, fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ggplot2::ylab("freq") + 
                        ggplot2::ggtitle(paste("frequencies",nameB," (NRMSE ",sp$gen$nrmse.fj,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))


    res_plot
}


#' Plots expected and generated probability distribution of degrees 
#' 
#' Plots as a heatmap the difference between expected probabilities and 
#' generated ones. Locations in red have a lower probability than expected; in blue, it's the opposite.
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile scale_fill_manual scale_fill_gradient2
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_errors_pdi <- function(sp, nameA="A", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    heat_map_gradient <- ggplot2::scale_fill_gradient2(limits=c(-1,1)) # , trans="log"

    diff_pdi <- errors.pdi(sp)
    colnames(diff_pdi) <- colnames(diff_pdi) # add_linebreaks_attributes
    diff_pdi$degree <- factor(seq(0,nrow(diff_pdi)-1))
    data_hm_pdi <- melt(diff_pdi)
    plot_pdi <- ggplot2::ggplot(data_hm_pdi, ggplot2::aes(variable, degree)) + 
                        ggplot2::geom_tile(ggplot2::aes(fill=value)) + 
                        heat_map_gradient + 
                        ggplot2::ggtitle(paste("difference degree for",nameA," (NRMSE ",sp$gen$nrmse.pdi,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_pdi
}
#' @rdname plot_errors_pdi
plot_errors_pdj <- function(sp, nameB="B", colorRef="darkgray", colorSynthetic="blue") {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    heat_map_gradient <- ggplot2::scale_fill_gradient2(limits=c(-1,1)) # , trans="log"

    diff_pdj <- errors.pdj(sp)
    colnames(diff_pdj) <- colnames(diff_pdj) # add_linebreaks_attributes(
    diff_pdj$degree <- factor(seq(0,nrow(diff_pdj)-1))
    data_hm_pdj <- melt(diff_pdj)
    plot_pdj <- ggplot2::ggplot(data_hm_pdj, ggplot2::aes(variable, degree)) + 
                        ggplot2::geom_tile(ggplot2::aes(fill=value)) + 
                        heat_map_gradient + 
                        ggplot2::ggtitle(paste("difference degree for",nameB," (NRMSE ",sp$gen$nrmse.pdj,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))
    plot_pdj
}

#' Plots expected and generated pairing probability 
#' 
#' Plots as a heatmap the difference between expected pairing probabilities and 
#' generated ones. Locations in red have a lower probability than expected; in blue, it's the opposite.
#'
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}} or \code{\link{matching.solve}}.
#' @inheritParams plot.dpp_result
#' @return a ggplot ready to display
#' 
#' @seealso \code{\link{plot.dpp_result}} for the plot of all the results in one call
#' 
#' @export 
#' 
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle geom_tile scale_fill_manual scale_fill_gradient2
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
plot_errors_pij <- function(sp) {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    heat_map_gradient <- ggplot2::scale_fill_gradient2(limits=c(-1,1)) # , trans="log"

    diff_pij <- errors.pij(sp)
    colnames(diff_pij) <- colnames(diff_pij) # add_linebreaks_attributes(
    diff_pij$attributesB <- row.names(diff_pij)
    #diff_pij$variable <- add_linebreaks_attributes(diff_pij$variable)
    data_hm_pij <- melt(diff_pij)
    plot_pij <- ggplot2::ggplot(data_hm_pij, ggplot2::aes(variable, attributesB)) + 
                        ggplot2::geom_tile(ggplot2::aes(fill=value)) + 
                        heat_map_gradient + 
                        ggplot2::ggtitle(paste("difference pairing (NRMSE ",sp$gen$nrmse.pij,")")) +
                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_pij
}

#' Plots a visual synthese of the errors induced by the resolution process. 
#'
#' Plots all the values controlled by the algorithm, and compares their expected value (passed as data or parameter)
#' and their actual value (as measured in the resulting synthetic population).
#' The resulting graphs contain: 
#' \itemize{
#'  \item a graph showing the relaxation parameters passed to the solving function, as computed by \code{\link{plot_relaxation}}
#'  \item a graph showing the Normalized Root Mean Square Error (NRMSE) for the each control variable, as computed by \code{\link{plot_errors}}
#'  \item a graph showing the population sizes asked for and generated for populations A and B, as computed by \code{\link{plot_population_sizes}}
#'  \item a graph showing the difference between the input and observed pairing probabilities pij, as computed by \code{\link{plot_errors_pij}} 
#'  \item two graphs showing the initial and observed frequencies of control variables in both populations A and B, as computed by \code{\link{plot_frequencies_A}}
#'  \item two graphs showing the initial and observed average degrees in both populations A and B, as computed by \code{\link{plot_average_degree_A}}
#'  \item two graphs showing the difference between the expected and measured distribution of probability of degrees for both A and B, as computed by \code{\link{plot_errors_pdi}}
#' }
#' 
#' Note you can call the individual plotting functions for rendering for papers or zooming.
#'
#' @param x an object returned by the \code{\link{matching.solve}} or the \code{\link{matching.generate}} method
#' @param nameA a meaningfull label for the entity type of population A, such as "dwellings" (default to "A")
#' @param nameB a meaningfull label for the entity type of population B, such as "households" (default to "B")
#' @param colorRef the color to be used to plot values passed as parameters (defaults to "darkgray")
#' @param colorSynthetic the color to be used to plot values measured in the synthetic population (defaults to "blue")
#' @param ... other parameters will be ignored quietly
#' @return nothing
#' 
#' @examples 
#' data(cas1)
#' prepared <- matching.prepare(cas1$sample.A, cas1$sample.B, cas1$pdi, cas1$pdj, cas1$pij)
#' solved <- matching.solve(prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                            delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' plot(solved, "dwellings", "households")
#' # we might as well drive the same on a generated case:
#' # sp <- matching.generate(solved, sample.A=cas1$sample.A, sample.B=cas1$sample.B, verbose=TRUE)
#' # plot(sp, "dwellings", "households")
#' 
#' @export
#'
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#'
# @importFrom ggplot2 ggplot aes ggplot2::geom_bar xlab ylab ggtitle scale_fill_gradient2 geom_tile scale_fill_manual
# @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#'
plot.dpp_resolved <- function(x, nameA="A", nameB="B", colorRef="darkgray", colorSynthetic="blue", ...) {

    sp <- x

    # TODO ensure this object is of the right type
    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    ensure_presence_ggplot()

    # plot relaxation parameters
    plot_all_relaxation <- plot_relaxation(sp, colorRef)

    # plot errors on each value
    plot_all_errors <- plot_errors(sp, colorSynthetic)

    # this scale will be used for various heatmaps    
    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    plot_population_sizes <- plot_population_sizes(sp, nameA=nameA, nameB=nameB, colorRef=colorRef, colorSynthetic=colorSynthetic)

    plot_degrees_A <- plot_average_degree_A(sp, nameA=nameA, colorRef=colorRef, colorSynthetic=colorSynthetic)
    plot_degrees_B <- plot_average_degree_B(sp, nameB=nameB, colorRef=colorRef, colorSynthetic=colorSynthetic)

    # compare expected and actual frequencies
    plot_frequencies_A <- plot_frequencies_A(sp, nameA=nameA, colorRef=colorRef, colorSynthetic=colorSynthetic)
    plot_frequencies_B <- plot_frequencies_B(sp, nameB=nameB, colorRef=colorRef, colorSynthetic=colorSynthetic)

    plot_pdi <- plot_errors_pdi(sp, nameA=nameA, colorRef=colorRef, colorSynthetic=colorSynthetic)
    plot_pdj <- plot_errors_pdj(sp, nameB=nameB, colorRef=colorRef, colorSynthetic=colorSynthetic)

    plot_pij <- plot_errors_pij(sp)

    # actual plot over a grid
    gridExtra::grid.arrange(
        plot_all_relaxation, plot_all_errors, 
        plot_population_sizes, plot_pij,
        plot_frequencies_A, plot_frequencies_B,
        plot_degrees_A, plot_degrees_B,
        plot_pdi, plot_pdj,
        ncol=2, nrow=5)

}

#' @rdname plot.dpp_resolved
plot.dpp_result <- plot.dpp_resolved




# TODO example
#
#' Plots the difference between proportions in the original sample and synthethic population
#'
#' Measures the difference between the total weight of the various levels of a variable in the original sample,
#' and the actual count of each level in the synthetic population.
#'
#' @param sample a sample as produced by \code{\link{create_sample}} 
#' @param generated a population as produced by \code{\link{matching.generate}} 
#' @param var.name the name of the variable to measure and plot
#' @param colorRef the color to be used to plot values passed as parameters (defaults to "darkgray")
#' @param colorSynthetic the color to be used to plot values measured in the synthetic population (defaults to "blue")
#' @param as.factor if TRUE (default), plots the variable as a factor instead of numerical
#'
#' @return the plot
#'
#' @importFrom reshape2 dcast
#'
plot_variable <- function(sample, generated, var.name, colorRef="darkgray", colorSynthetic="blue", as.factor=TRUE) {

    ensure_presence_ggplot()

    if (class(sample) != "dpp_sample") 
        stop("the parameter 'sample' should be the result of a create_sample call")
    if (class(generated) != "data.frame")
        stop("the parameter 'sample' should be the result of a matching.generate call")
  

    # TODO check input types
    sample <- dcast(sample$sample, paste(var.name,"~.",sep=""), fun.aggregate=sum, value.var=sample$dictionary$colname.weight)
    synthetic <- dcast(generated, paste(var.name,"~.",sep=""), fun.aggregate=length, value.var="current.degree")

    prop_sample <- normalise(sample$.)
    prop_synthetic <- normalise(synthetic$.)
    df <- data.frame(
            modalities=sample[[var.name]], 
            proportion=c(prop_sample,prop_synthetic), 
            state=factor(c(rep("theoretical",nrow(sample)), rep("synthetic",nrow(synthetic))), levels=c("theoretical","synthetic"))
            )
    if (as.factor) {
        df$modalities <- as.factor(df$modalities)
    }
    rmse <- sqrt( mean( (prop_sample - prop_synthetic)^2 ) )
    
    scale_gray_blue <- ggplot2::scale_fill_manual(values=c(colorRef,colorSynthetic))

    ggplot2::ggplot(df, ggplot2::aes(x=modalities, y=proportion, fill=state)) + 
                        ggplot2::geom_bar(stat="identity", position = 'dodge2') + 
                        scale_gray_blue + 
                        ggplot2::ylab("freq") + 
                        ggplot2::ggtitle(paste("proportions of attribute ", var.name, " (RMSE:",rmse,")",sep=""))

}


# TODO doc
create_probability_view_init <- function(sp, maxcol=NULL) {

    if ( (class(sp) != "dpp_prepared") && (class(sp) != "dpp_resolved") ) 
        stop("the parameter 'sp' should be the result of a 'matching.prepare' or 'matching.solve' call")

    count_cols_left <- 4 
    count_cols_right <- if (is.null(maxcol)) ncol(sp$inputs$pij$data) else maxcol
 
    count_rows_left <- nrow(sp$inputs$pij$data)

    missingcol <- if (is.null(maxcol)) "" else " & ..."

    end_line <- "\\\\\n"

    # define the header of the table
    s <- paste(
        "\\begin{tabular}{r|cccc|", 
        paste(rep("c",count_cols_right),collapse=""), 
        "|c}\n",
        sep="")
    
    # Mod A
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\text{Mod}_i^A$} ",
        " & \\rot{",
        paste(names(sp$stats$fi[1:count_cols_right-1]), collapse="} & \\rot{"),
        "} ", missingcol,
        "& ",
        end_line,
        "\\hline\n",
        sep="")
    # fi
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$f_i$} ",
        " & ",
        paste(formatC(sp$stats$fi[1:count_cols_right-1],format="G"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$stats$fi), format="G"),
        end_line,
        sep="")
    # di
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\tilde{d}_i$} ",
        " & ",
        paste(formatC(sp$inputs$di[1:count_cols_right-1], format="G"), collapse=" & "),
        missingcol, " & ",
        end_line,
        sep="")

    # pi
    s <-  paste(
        s,
        "\\multicolumn{5}{r|}{$p_i$} ",
        " & ",
        paste(formatC(colSums(sp$inputs$pij$data[1:count_cols_right-1]), format="G"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$inputs$pij$data), format="G"),
        end_line,
        sep="")

    s <- paste(s, "\n", sep="")

    # define the header of the left part
    s <- paste(s, 
            "$\\text{Mod}_j^B$ & $f_j$ & $\\tilde{d}_j$ & $p_j$ & $j/i$ & ", 
            paste(seq(1, count_cols_right-1), collapse="&"), 
            missingcol, "&",
            end_line, 
            "\\hline\n",
            sep="")

    for (j in 1:count_rows_left) {
        s <- paste(s,
                paste(names(sp$stats$fj[j]), collapse=""), " & ",
                formatC(sp$stats$fj[j], format="G"), " & ",
                formatC(sp$inputs$dj[j], format="G"), " & ",
                formatC(sum(sp$inputs$pij$data[j,1:count_cols_right-1]), format="G"), " & ",
                j, " & ",
                paste(sapply(sp$inputs$pij$data[j,1:count_cols_right-1], formatC, format="G"), collapse=" & "), missingcol, " & ",
                end_line,
                sep="")
    }

    s <- paste(s, 
        "\\hline\n", 
        " & ", formatC(sum(sp$stats$fj), format="G"), 
        " & ", 
        " & ", formatC(sum(sp$inputs$pij$data), format="G"), 
        " & ",
        paste(rep("&",count_cols_right), collapse=""),
        " & ", formatC(sum(sp$inputs$pij$data), format="G"), 
        end_line,
        sep="")


    s <- paste(s, "\\end{tabular}\n", sep="")

    s

 }

create_probability_view_solved  <- function(sp, maxcol=NULL) {

    if (class(sp) != "dpp_resolved") 
        stop("the parameter 'sp' should be the result of a 'matching.solve' call")

    count_cols_left <- 4 
    count_cols_right <- if (is.null(maxcol)) ncol(sp$inputs$pij$data) else maxcol
 
    count_rows_left <- nrow(sp$inputs$pij$data)
    missingcol <- if (is.null(maxcol)) "" else " & ..."

    end_line <- "\\\\\n"

    # define the header of the table
    s <- paste(
        "\\begin{tabular}{r|cccc|", 
        paste(rep("c",count_cols_right),collapse=""), 
        "|c}\n",
        sep="")
    
    # Mod A
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\text{Mod}_i^A$} ",
        " & \\rot{",
        paste(names(sp$stats$fi[1:count_cols_right-1]), collapse="} & \\rot{"),
        "} ", missingcol, "& ",
        end_line,
        "\\hline\n",
        sep="")
    # fi
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{f_i}$} ",
        " & ",
        paste(formatC(sp$gen$hat.fi[1:count_cols_right-1],format="G"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$gen$hat.fi), format="G"),
        end_line,
        sep="")
    # di
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{\\tilde{d}}_i$} ",
        " & ",
        paste(formatC(sp$gen$hat.di[1:count_cols_right-1], format="G"), collapse=" & "),
        missingcol, " & ",
        end_line,
        sep="")

    # pi
    s <-  paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{p}_i$} ",
        " & ",
        paste(formatC(colSums(sp$gen$hat.pij[,1:count_cols_right-1]), format="G"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$gen$hat.pij), format="G"),
        end_line,
        sep="")

    s <- paste(s, "\n", sep="")

    # define the header of the left part
    s <- paste(s, 
            "$\\text{Mod}_j^B$ & $\\hat{f}_j$ & $\\hat{\\tilde{d}}_j$ & $\\hat{p}_j$ & $j/i$ & ", 
            paste(seq(1, count_cols_right-1), collapse="&"), 
            missingcol, "&",
            end_line, 
            "\\hline\n",
            sep="")

    for (j in 1:count_rows_left) {
        s <- paste(s,
                paste(names(sp$gen$hat.fj[j]), collapse=""), " & ",
                formatC(sp$gen$hat.fj[j], format="G"), " & ",
                formatC(sp$gen$hat.dj[j], format="G"), " & ",
                formatC(sum(sp$gen$hat.pij[j,]), format="G"), " & ",
                j, " & ",
                paste(sapply(sp$gen$hat.pij[j,1:count_cols_right-1], formatC, format="G"), collapse=" & "), missingcol, " & ",
                end_line,
                sep="")
    }

    s <- paste(s, 
        "\\hline\n", 
        " & ", formatC(sum(sp$gen$hat.fj), format="G"), 
        " & ", 
        " & ", formatC(sum(sp$gen$hat.pij), format="G"), 
        " & ",
        paste(rep("&",count_cols_right), collapse=""),
        " & ", formatC(sum(sp$gen$hat.pij), format="G"), 
        end_line,
        sep="")


    s <- paste(s, "\\end{tabular}\n", sep="")

    s

 }

create_discrete_view_solved  <- function(sp, maxcol=NULL) {

    if (class(sp) != "dpp_resolved") 
        stop("the parameter 'sp' should be the result of a 'matching.solve' call")

    count_cols_left <- 4 
    count_cols_right <- if (is.null(maxcol)) ncol(sp$inputs$pij$data) else maxcol
 
    count_rows_left <- nrow(sp$inputs$pij$data)
    missingcol <- if (is.null(maxcol)) "" else " & ..."

    end_line <- "\\\\\n"

    # define the header of the table
    s <- paste(
        "\\begin{tabular}{r|cccc|", 
        paste(rep("c",count_cols_right),collapse=""), 
        "|c}\n",
        sep="")
    
    # Mod A
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\text{Mod}_i^A$} ",
        " & \\rot{",
        paste(names(sp$stats$fi[1:count_cols_right-1]), collapse="} & \\rot{"),
        "} ", missingcol ," & ",
        end_line,
        "\\hline\n",
        sep="")
    # fi
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{c_i}$} ",
        " & ",
        paste(formatC(sp$gen$hat.ci[1:count_cols_right-1],format="d"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$gen$hat.ci), format="d"),
        end_line,
        sep="")
    # di
    s <- paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{\\tilde{d}}_i$} ",
        " & ",
        paste(formatC(sp$gen$hat.di[1:count_cols_right-1], format="G"), collapse=" & "),
        missingcol, " & ",
        end_line,
        sep="")

    # pi
    s <-  paste(
        s,
        "\\multicolumn{5}{r|}{$\\hat{n}_i$} ",
        " & ",
        paste(formatC(colSums(sp$gen$hat.nij[,1:count_cols_right-1]), format="d"), collapse=" & "),
        missingcol, " & ",
        formatC(sum(sp$gen$hat.nij), format="d"),
        end_line,
        sep="")

    s <- paste(s, "\n", sep="")

    # define the header of the left part
    s <- paste(s, 
            "$\\text{Mod}_j^B$ & $\\hat{c}_j$ & $\\hat{\\tilde{d}}_j$ & $\\hat{n}_j$ & $j/i$ & ", 
            paste(seq(1, count_cols_right-1), collapse="&"), 
            missingcol, "&",
            end_line, 
            "\\hline\n",
            sep="")

    for (j in 1:count_rows_left) {
        s <- paste(s,
                paste(names(sp$gen$hat.cj[j]), collapse=""), " & ",
                formatC(sp$gen$hat.cj[j], format="d"), " & ",
                formatC(sp$gen$hat.dj[j], format="G"), " & ",
                formatC(sum(sp$gen$hat.nij[j,]), format="d"), " & ",
                j, " & ",
                paste(sapply(sp$gen$hat.nij[j,1:count_cols_right-1], formatC, format="d"), collapse=" & "), missingcol, " & ",
                end_line,
                sep="")
    }

    s <- paste(s, 
        "\\hline\n", 
        " & ", formatC(sum(sp$gen$hat.cj), format="d"), 
        " & ", 
        " & ", formatC(sum(sp$gen$hat.nij), format="d"), 
        " & ",
        paste(rep("&",count_cols_right), collapse=""),
        " & ", formatC(sum(sp$gen$hat.nij), format="d"), 
        end_line,
        sep="")


    s <- paste(s, "\\end{tabular}\n", sep="")

    s

 }
