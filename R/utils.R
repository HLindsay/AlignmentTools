#'@title .checkForMcols
#'@description Internal function for checking required metadata columns
#'are present.
#'@author Helen Lindsay
#'@param alns GAlignments object containing metadata columns to check
#'@param mcols character(n) Name(s) of metadata columns to check for
#'@param func.nm character(1) The name of the function that called .checkForMcols
#'@param err.func function The error function to invoke if test fails,
#'e.g. warning or stop. (Default: warning)
#'@rdname checkForMcols
.checkForMcols <- function(alns, mcols, func.nm = NULL, err.func = warning){
    absent_mcols <- !mcols %in% names(mcols(alns))
    if (! any(absent_mcols)){ return(TRUE) }

    if (! is.null(func.nm)){
      warn_str <- paste("%s requires that alignments have",
                         "the following missing metadata columns:\n%s")
    } else {
      func.nm <- ""
      warn_str <- paste("%sAlignments must have the following",
                        "missing metadata columns:\n%s")
    }

    err.func(sprintf(warn_str, func.nm,
                     paste(mcols[absent_mcols], collapse = "\n")))
    return(FALSE)
}
