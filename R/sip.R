#' @title
#' Single-iteration permutation for large-scale biobank data
#'
#' @description
#' This function performs the basic permutation for permuting phenotype
#' vectors in biobank data.
#'
#' @details
#' The function requires the following inputs:
#'
#' 1) A sample-by-variable dataframe with phenotypes and covariates.
#' Column names should include an ID variable, sex variable,
#' genotypic covariate names, phenotypic covariate names,
#' and phenotype names. (N.B. a secondary ID variable can be included
#' in the genotypic covariate names.)
#'
#' 2) A string identifying the ID variable name (e.g., id.var="ID").
#'
#' 3) A vector of genotypic covariates
#' (e.g., geno.vars=c("ID2","Batch","PC1","PC2",...)).
#'
#' 4) Optional: within.sex = FALSE. Default is within.sex = TRUE
#' and will permute males and females separately.
#'
#' 5) If within.sex = TRUE (the default), a string identifying the
#' sex variable name (e.g., sex.var="Inferred_Sex").
#'
#' 6) If within.sex = TRUE (the default), male and female values
#' in the sex vector (e.g., male.val=1, female.val=2).
#'
#' 7) Optional: a seed for sampling. If a seed is not provided, one will be
#' chosen randomly during the sampling process (e.g., seed=123).
#'
#' 8) N.B. Any column names not specified in (2)-(6) are assumed to be
#' phenotypes or phenotypic covariates.
#'
#' @param df Data frame
#' @param id.var,sex.var Strings
#' @param male.val,female.val Strings or integers
#' @param geno.vars Character vectors
#' @param seed Number
#' @param within.sex Boolean (TRUE/FALSE)
#' @returns Data frame
#' @export
#' @examples
#' sip(df = sip_exampleData, id.var = "IID", sex.var = "SEX", male.val = 1,
#' female.val = 2, geno.vars = c("FID","ANCESTRY","BATCH",paste0("PC",1:4)))

sip <- function(df = NULL, id.var = NULL, sex.var = NULL, male.val = NULL,
                female.val = NULL, geno.vars = NULL,
                within.sex = TRUE, seed = NULL) {

  tryCatch(

    expr = {

      # input checks #

      if (missing(df) || is.null(df) ||
          missing(id.var) || is.null(id.var) ||
          missing(geno.vars) || is.null(geno.vars)) {
        stop("`df`, `id.var`, and `geno.vars` must be provided.")
      }

      if (!is.data.frame(df)) {
        stop("Error: df file is not properly formatted.
             Please input a sample-by-variable dataframe.")
      } #else {
        #df <- as.data.frame(df)
      #}

      if (!(id.var %in% colnames(df))) {
        stop("Error: ID variable is not in df column names.
             Please input one ID variables as a string.")
      }

      if (!all(geno.vars %in% names(df))) {
        stop("Error: Genotypic variables are not in df column names.
             Please input genotypic variable names as a vector of strings.")
      }

      # get phenotype and phenotypic covariate names #
      pheno.vars <- colnames(df)[!(colnames(df) %in%
                                      c(id.var, sex.var, geno.vars))]

      # divide df by sex #
      if (within.sex) {

        if (missing(sex.var) || missing(male.val) || missing(female.val)) {
          stop("Error: One or more sex-related inputs have null values.
             Please check your variable assignments or choose
             'within.sex = FALSE' to ignore sex while permuting.")
        }

        if (!(sex.var %in% colnames(df))) {
          stop("Error: Sex variable is not in df column names.
             Please input a sex variable as a string or choose
             'within.sex = FALSE' to ignore sex while permuting.")
        }

        if (length(unique(df[[sex.var]])) > 2) {
          stop("Error: more than two sexes detected in the sex column.")
        }

        if (!all(unique(df[[sex.var]]) %in% c(male.val, female.val))) {
          stop("Error: sex values do not match the male
             or female values provides. Please input proper sex values
             (e.g., male.val=1, female.val=2).")
        }

        if (female.val %in% unique(df[[sex.var]])) {
          fdf <- get_sexDF(df = df, sex.var = sex.var, sex.val = female.val)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mdf <- get_sexDF(df = df, sex.var = sex.var, sex.val = male.val)
        }

        # separate df into fixed and permutable vectors #
        if (female.val %in% unique(df[[sex.var]])) {
          ffix <- get_fixed(df = fdf, id.var = id.var, sex.var = sex.var,
                            geno.vars = geno.vars)
          fpheno <- get_permute(df = fdf, id.var = id.var,
                                pheno.vars = pheno.vars)
          fpheno[[id.var]] <- NULL
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mfix <- get_fixed(df = mdf, id.var = id.var, sex.var = sex.var,
                            geno.vars = geno.vars)
          mpheno <- get_permute(df = mdf, id.var = id.var,
                                pheno.vars = pheno.vars)
          mpheno[[id.var]] <- NULL
        }

        # set seed #
        if (is.null(seed)) {
          seed <- sample(seq(999999), 1)
        }
        message(paste("Seed:",seed))

        # get permutation index #
        if (female.val %in% unique(df[[sex.var]])) {
          fidx <- get_permIdx(df = fdf, seed = seed)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          midx <- get_permIdx(df = mdf, seed = seed)
        }

        # reintegrate into permuted datasets #
        if (female.val %in% unique(df[[sex.var]])) {
          fperm <- get_permDF(fix.df = ffix, perm.df = fpheno, perm.idx = fidx)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mperm <- get_permDF(fix.df = mfix, perm.df = mpheno, perm.idx = midx)
        }

        # recombine permuted sex-specific datasets #
        if (female.val %in% unique(df[[sex.var]])
            & male.val %in% unique(df[[sex.var]])) {
          perm <- rbind(fperm, mperm)
        } else if (female.val %in% unique(df[[sex.var]])
                   & !(male.val %in% unique(df[[sex.var]]))) {
          perm <- fperm
        } else if (!(female.val %in% unique(df[[sex.var]]))
                   & male.val %in% unique(df[[sex.var]])) {
          perm <- mperm
        }

      } else {

        # separate df into fixed and permutable vectors #

        fix <- get_fixed(df = df, id.var = id.var, sex.var = sex.var,
                          geno.vars = geno.vars)
        pheno <- get_permute(df = df, id.var = id.var,
                              pheno.vars = pheno.vars)
        pheno[[id.var]] <- NULL

        # set seed #
        if (missing(seed)) {
          seed <- sample(seq(999999), 1)
        }
        message(paste("Seed:",seed))

        # get permutation index #
        idx <- get_permIdx(df = df, seed = seed)

        # reintegrate into permuted datasets #
        perm <- get_permDF(fix.df = fix, perm.df = pheno, perm.idx = idx)

      }

      # reorder permuted data #
      operm <- order_permDF(df = df, perm = perm, id.var = id.var)
      return(operm)

    }, error = function(e){
      message(
        sprintf("An error occurred at %s : %s",
                Sys.time(),
                e)
      )
    }
  )
}
