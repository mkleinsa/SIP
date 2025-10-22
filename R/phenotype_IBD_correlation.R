#' @title
#' Phenotype-IBD correlation
#'
#' @description
#' This function first calculates the correlation between phenotypes
#' for sample pairs. Then it calculates the correlation between the
#' phenotype correlation and identity-by-descent for the sample pairs.
#' (N.B., Please omit missing values before running this function.)
#'
#' @details
#' The function requires the following inputs:
#'
#' 1) A sample-by-variable data frame with sample ID and phenotypes.
#' Columns should include an individual ID variable and phenotype names.
#'
#' 2) A string identifying the ID variable name (e.g., id.var="ID").
#'
#' 3) A vector of phenotype names.
#' (e.g., pheno.vars=c("PHENO1","PHENO2",...)).
#'
#' 4) A relatedness data frame containing identity-by-descent (IBD) for
#' pairs of individuals in the sample-by-variable data frame
#' Column names should include two ID variables and an IBD variable.
#'
#' 5) A vector of the 2 ID variable names in the relatedness data frame.
#' (e.g., rid.vars=c("ID1","ID2")).
#'
#' 6) A string identifying the IBD variable name (e.g., ibd.var="PropIBD").
#'
#' @param df,rel.df Data frame
#' @param id.var,ibd.var Strings
#' @param pheno.vars,rid.vars Character vectors
#' @returns Correlation
#' @export
#' @examples
#' phenotype_IBD_correlation(df = sipPair_exampleData,
#' rel.df = sipPair_relatednessData, id.var = "IID",
#' rid.vars = c("IID1","IID2"), ibd.var = "PropIBD",
#' pheno.vars = paste0("PHENO",1:300))

phenotype_IBD_correlation <- function(df = NULL, rel.df = NULL,
                                          id.var = NULL, rid.vars = NULL,
                                          ibd.var = NULL, pheno.vars = NULL) {

  tryCatch(

    expr = {

      if (is.null(df) | is.null(id.var) |
          is.null(pheno.vars) | is.null(rel.df) |
          is.null(rid.vars) | is.null(ibd.var)) {
        stop("Error: One or more required inputs have null values.
             Please check your variable assignments.")
      }

      # input checks #
      if (is.integer(ncol(df))==FALSE | is.integer(nrow(df))==FALSE) {
        stop("Error: df file is not properly formatted. Please input a
             sample-by-variable dataframe.")
      } else {
        df <- as.data.frame(df)
      }

      if (is.integer(ncol(rel.df))==FALSE | is.integer(nrow(rel.df))==FALSE) {
        stop("Error: Relatedness file is not properly formatted. Please input
             a sample-by-variable dataframe.")
      } else {
        rel.df <- as.data.frame(rel.df)
      }

      if (FALSE %in% (rid.vars %in% colnames(rel.df))) {
        stop("Error: ID variable pairs are not in relatedness dataframe column
             names. Please provide two ID variables as a vector of strings.")
      }

      if (length(rid.vars) != 2) {
        stop("Error: Improper number of relatedness dataframe ID variable names.
             Please provide two ID variables as a vector of strings.")
      }

      if (FALSE %in% (ibd.var %in% colnames(rel.df))) {
        stop("Error: Identity-by-descent variable is not in relatedness
             dataframe column names. Please provide one IBD variable name
             as a string.")
      }

      rel.df[[ibd.var]] <- as.numeric(rel.df[[ibd.var]])
      if (TRUE %in% is.na(rel.df[[ibd.var]])
          | TRUE %in% is.null(rel.df[[ibd.var]])
          | TRUE %in% (rel.df[[ibd.var]]) < 0
          | TRUE %in% (rel.df[[ibd.var]]) > 1) {
        stop("Error: One of more pairs of individuals do not have valid
             IBD values.")
      }

      if (FALSE %in% (id.var %in% colnames(df))) {
        stop("Error: ID variable is not in df column names. Please input
             one ID variable as a string.")
      }

      if (!(TRUE %in% is.na(pheno.vars))
          & (FALSE %in% (pheno.vars %in% colnames(df)))) {
        stop("Error: Phenotype names are not in df column names.
             Please input phenotype variable names as a vector of strings.")
      }

      # calculate correlation #
      pcor <- stats::setNames(as.data.frame(
        stats::cor(t(as.data.frame(df)[,pheno.vars]))), df[[id.var]])

      # get unique pairs #
      pcor[upper.tri(pcor)] <- NA

      # reformat #
      pcor[["VAR1"]] <- df[[id.var]]
      pcor <- setNames(data.table::melt(data.table::setDT(pcor),
                                        id.vars = "VAR1"),
                       c("VAR1","VAR2","pheno_correlation"))
      pcor <- na.omit(pcor)
      pcor <- pcor[pcor[["VAR1"]] != pcor[["VAR2"]],]

      # get correlation with IBD status #
      picor.x <- merge(pcor, rel.df,
                       by.x = c("VAR1","VAR2"),
                       by.y = c(rid.vars[1], rid.vars[2]))
      picor.y <- merge(pcor, rel.df,
                       by.x = c("VAR1","VAR2"),
                       by.y = c(rid.vars[2], rid.vars[1]))
      picor <- rbind(picor.x, picor.y)

      if (nrow(picor) == nrow(rel.df)) {
        return(stats::cor(as.data.frame(picor)[,c(ibd.var,
                                                  "pheno_correlation")])[2])
      } else {
        print("Error")
      }

    }, error = function(e){
      print(
        sprintf("An error occurred at %s : %s",
                Sys.time(),
                e)
      )
    }
  )
}



