
#' Creates a LaTeX view with the probabilistic values after resolution
#'
#' Create a table view of the solved problem, with probabilistic values. 
#' Ready to be written into a tex file.
#' 
#' @param sp a solved case produced by \code{\link{matching.solve}} 
#' @param maxcol an optional integer containing the highest count of columns to display (useful for large tables)
#' @return a string containing a LaTeX table
#' 
#' @examples 
#' data(dwellings_households)
#' prepared <- matching.prepare(
#'                      dwellings_households$sample.A, dwellings_households$sample.B, 
#'                      dwellings_households$pdi, dwellings_households$pdj, 
#'                      dwellings_households$pij)
#' solved <- matching.solve(
#'                      prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                      delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' # create the string
#' pv <- as.latex.table.probabilistic.init(solved)
#' print(pv)
#' # write it into a file
#' # write(pv, file="case1_probabilistic_init.tex")
#' 
#' @export
#'
#' @keywords latex
#'
as.latex.table.probabilistic.init <- function(sp, maxcol=NULL) {

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
                formatC(sum(sp$inputs$pij$data[j,]), format="G"), " & ",
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

#' Creates a LaTeX view with the probabilistic values after resolution
#'
#' Create a table view of the solved problem, with probabilistic values. 
#' Ready to be written into a tex file.
#' 
#' @param sp a solved case produced by \code{\link{matching.solve}} 
#' @param maxcol an optional integer containing the highest count of columns to display (useful for large tables)
#' @return a string containing a LaTeX table
#' 
#' @examples 
#' data(dwellings_households)
#' prepared <- matching.prepare(
#'                      dwellings_households$sample.A, dwellings_households$sample.B, 
#'                      dwellings_households$pdi, dwellings_households$pdj, 
#'                      dwellings_households$pij)
#' solved <- matching.solve(
#'                      prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                      delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' # create the string
#' pv <- as.latex.table.probabilistic.solved(solved)
#' print(pv)
#' # write it into a file
#' # write(pv, file="case1_probabilistic_solved.tex")
#' 
#' @export
#'
#' @keywords latex
#'
as.latex.table.probabilistic.solved  <- function(sp, maxcol=NULL) {

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

#' Creates a LaTeX view with the discrete values
#'
#' Create a table view of the solved problem, with integer values. 
#' Ready to be written into a tex file.
#' 
#' @param sp a solved case produced by \code{\link{matching.solve}} 
#' @param maxcol an optional integer containing the highest count of columns to display (useful for large tables)
#' @return a string containing a LaTeX table
#' 
#' @examples 
#' data(dwellings_households)
#' prepared <- matching.prepare(
#'                      dwellings_households$sample.A, dwellings_households$sample.B, 
#'                      dwellings_households$pdi, dwellings_households$pdj, 
#'                      dwellings_households$pij)
#' solved <- matching.solve(
#'                      prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                      delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' # create the string
#' dv <- as.latex.table.discrete(solved)
#' print(dv)
#' # write it into a file
#' # write(dv, , file="case1_discrete.tex")
#' 
#' @export
#'
#' @keywords latex
#'
as.latex.table.discrete <- function(sp, maxcol=NULL) {

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


#' Writes relaxation parameters as LaTeX
#'
#' Create LaTeX code that depicts the relaxation parameters.
#' 
#' @param sp a solved case produced by \code{\link{matching.solve}} 
#' @return a string containing a LaTeX string.
#' 
#' @examples 
#' data(dwellings_households)
#' prepared <- matching.prepare(
#'                      dwellings_households$sample.A, dwellings_households$sample.B, 
#'                      dwellings_households$pdi, dwellings_households$pdj, 
#'                      dwellings_households$pij)
#' solved <- matching.solve(
#'                      prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, 
#'                      delta.A=0, gamma=0, delta.B=0, phi.B=0, nu.B=1, verbose=TRUE)
#' # create the string
#' rel <- as.latex.relaxation.constrainsts(solved)
#' print(rel)
#' # write it into a file
#' # write(rel, , file="relaxation_parameters.tex")
#' 
#' @export
#'
#' @keywords latex
#'
as.latex.relaxation.constrainsts <- function(sp) {

    weights <- c(sp$inputs$nu.A, sp$inputs$phi.A, sp$inputs$delta.A, sp$inputs$gamma, sp$inputs$delta.B, sp$inputs$phi.B, sp$inputs$nu.B)
    weights.names <- c("\\nu^A", "\\phi^A", "\\delta^A", "\\gamma", "\\delta^B", "\\phi^B", "\\nu^B")

    index.0 <- which(weights == 0)
    index.1 <- which(weights == 1)
    index.others <- which((weights != 0) & (weights != 1))
    
    elems <- list()
    if (length(index.0) > 0) 
        elems[[length(elems)+1]] <- paste("$", paste(weights.names[index.0], collapse="=\\allowbreak"), "=0$", sep="") 

    if (length(index.1) > 0) 
        elems[[length(elems)+1]] <- paste("$", paste(weights.names[index.1], collapse="=\\allowbreak"), "=1$", sep="")

    if (length(index.others) > 0) 
        elems[[length(elems)+1]] <- paste("$", paste(weights.names[index.others], weights[index.others], collapse="=\\allowbreak"), "=1$", sep="")
            

    paste("Relaxation parameters are ", paste(elems, collapse=", "), ". ", sep="")
    
}