#' SIP-package internal setup
#'
#' Declares global variables and imports common functions
#' used throughout the SIP package. This file is executed
#' when the package is loaded.
#'
#' @name SIP-package
#' @keywords internal
#' @importFrom stats na.omit setNames
NULL

utils::globalVariables(c(
  "VAR1", "VAR2", "aes", "correlation",
  "element_blank", "element_text", "na.omit", "setNames"
))

