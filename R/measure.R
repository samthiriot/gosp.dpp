
#' Display the result of measures on a generated population
#'
#' @param x the measure to display
#' @param ... ignored 
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_measure <- function(x,...) {

    cat("measures on the generated population\n\n")

    cat("$nL=",x$nL,"\n",sep="")

    cat("\n$fi [MSE ",x$fi$mse.f.orig,"]\n",sep="")
    print(x$fi$hat.f)
    
    cat("\n$pdi [MSE ",x$pdi$mse.orig,"]\n",sep="")
    print(x$pdi$hat.pd)

    cat("\n$pij [MSE ",x$pij$mse.orig,"]\n",sep="")
    print(x$pij$hat.pij)

    cat("\n$pdj [MSE ",x$pdj$mse.orig,"]\n",sep="")
    print(x$pdj$hat.pd)
    
    cat("\n$fj [MSE ",x$fj$mse.f.orig,"]\n",sep="")
    print(x$fj$hat.f)
    
}


#' Display a generated population
#'
#' @param x the generation result to display
#' @param ... ignored 
#'
#' @export
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch>
#' 
print.dpp_result <- function(x,...) {
    cat("generation result\n")

    cat("$pop\n")
    print(x$pop)

    cat("\n$measure\n")
    print(x$measure)
}


#' Merge entities and attributes
#' 
#' Merges a generated population: returns the join of population A,
#' links and population B. 
#' 
#' The underlying operation is done by \code{\link{merge}}.
#'
#' @param pop the generated population
#' @return a dataframe from the join 
#' 
#' @export 
#' 
#' @keywords internal
#'
merge_links <- function(pop) {

	# merge the datasets	
	A2l <- merge(pop$A, pop$links, by="id.A", suffixes=c("_A", "_l"))
	A2l2B <- merge(A2l, pop$B, by="id.B", suffixes=c("_A", "_B"))

	A2l2B
}


#' Measure the errors in the distributions of degrees pdi
#' 
#' Returns the difference between the actual and input distributions of degrees. 
#' For instance a "minus one" means we reached one less than expected.
#' 
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}} or \code{\link{matching.solve}}.
#' @return a dataframe containing the difference between the expected and actual distribution of degrees
#' 
#' @seealso \code{\link{errors.pdj}} for the measure of distributions of degrees pdj
#' 
#' @export 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
errors.pdi <- function(sp) {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    sp$gen$hat.pdi - as.data.frame(sp$inputs$pdi_fixed)

}

#' Measure the errors in the distributions of degrees pdj
#' 
#' Returns the difference between the actual and input distributions of degrees. 
#' For instance a "minus one" means we reached one less than expected.
#' 
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}} or \code{\link{matching.solve}}.
#' @return a dataframe containing the difference between the expected and actual distribution of degrees
#' 
#' @seealso \code{\link{errors.pdi}} for the measure of distributions of degrees pdi
#' 
#' @export 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
errors.pdj <- function(sp) {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    sp$gen$hat.pdj - as.data.frame(sp$inputs$pdj_fixed)

}

#' Measure the errors in the pairing probabilities  pij
#' 
#' Returns the difference between the actual and input pairing probabilities. 
#' For instance a "minus one" means we reached one less than expected.
#' 
#' @param sp a synthetic population, as produced by \code{\link{matching.generate}} or \code{\link{matching.solve}}.
#' @return a dataframe containing the difference between the expected and actual pairing probabilities
#' 
#' @seealso \code{\link{errors.pdi}} for the measure of distributions of degrees pdi
#' 			\code{\link{errors.pdj}} for the measure of distributions of degrees pdj
#' 
#' @export 
#' 
#' @author Samuel Thiriot <samuel.thiriot@res-ear.ch> 
#' 
errors.pij <- function(sp) {

    if ((class(sp) != "dpp_result") && (class(sp) != "dpp_resolved")) 
        stop("the data to analyze x should be the result of a matching.solve or matching.generate call")
    
    sp$gen$hat.pij - as.data.frame(sp$inputs$pij_fixed)

}
