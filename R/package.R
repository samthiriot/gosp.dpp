.onLoad <- function(libname = find.package("gosp.dpp"), pkgname = "gosp.dpp") {

	# avoid warnings about "no visible binding for local variable XXX" due to ggplot2 functions
	# see https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
	missingvariables.names <- c(	
			"classes", "average.degree", "state", 
			"NRMSE", "variable", "degree", "attributesB", 
			"value", "count", "modalities", 
			"proportion", "type")

	# if (getRversion() >= "3.1.0") {
	# 	# remove warnings for foreign checks in plotting.R (due to ggplot2 functions)
	# utils::suppressForeignCheck(missingvariables.names, add=TRUE)
	# } else 
	if (getRversion() >= "2.15.1") {
		# cannot remove the warnings in the most beautiful way, but we still can do it somehow
		utils::globalVariables(missingvariables.names)
	}

}