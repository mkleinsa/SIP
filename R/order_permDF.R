#' This function reorders the permuted sample-by-variable data frame
#' to match the ID order of the primary sample-by-variable data frame.
#'
#' @param df,perm Data frames
#' @param id.var String
#' @returns Data frame

order_permDF <- function(df = df, perm = perm, id.var = id.var) {
  perm <- perm[, colnames(df)]
  perm <- perm[match(df[[id.var]], perm[[id.var]]), ]
  return(perm)
}
