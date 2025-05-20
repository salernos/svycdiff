
#=== METHODS ===================================================================

#--- PRINT ---------------------------------------------------------------------

#' @export
print.svycdiff <- function(x, digits = 4L, ...) {

  cat(

    "\nOutcome Model:  \n",

    paste(deparse(x$fit_y$call), sep = "\n", collapse = "\n"), "\n",

    "\nTreatment Model:  \n",

    paste(deparse(x$fit_a$call), sep = "\n", collapse = "\n"), "\n",

    "\nSelection Model:  \n",

    paste(deparse(x$fit_s$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  cat("\nCDIFF:  \n")

  cdiff <- x$cdiff

  names(cdiff) <- c("CDIFF", "SE", "LCL", "UCL", "P-Value")

  print(round(cdiff, digits))
}

#--- SUMMARY -------------------------------------------------------------------

#' @export
summary.svycdiff <- function(object, digits = 4L, ...) {

  cat("\nCDIFF:  \n")

  cdiff <- object$cdiff

  names(cdiff) <- c("CDIFF", "SE", "LCL", "UCL", "P-Value")

  print(round(cdiff, digits))

  cat("\nOutcome Model:  \n")

  summary(object$fit_y)

  cat("\nTreatment Model:  \n")

  summary(object$fit_a)

  cat("\nSelection Model:  \n")

  summary(object$fit_s)
}
