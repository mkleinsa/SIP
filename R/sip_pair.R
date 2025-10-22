#' @title
#' Single-iteration paired permutation for large-scale biobank data
#' with relatedness
#'
#' @description
#'# This function performs paired permutation for permuting phenotype
#' vectors among related individuals in biobank data.
#'
#' @details
#' The function requires the following inputs:
#'
#' 1) A sample-by-variable data frame with phenotypes and covariates.
#' Column names should include an ID variable, sex variable,
#' genotypic covariate names, phenotypic covariate names, and phenotype names.
#' (N.B. a secondary ID variable can be included
#' in the genotypic covariate names.)
#'
#' 2) A string identifying the ID variable name (e.g., id.var="ID").
#'
#' 3) A vector of genotypic covariates
#' (e.g., geno.vars=c("ID2","Batch","PC1","PC2",...)).
#'
#' 4) A relatedness data frame containing identity-by-descent (IBD) for
#' pairs of individuals in the sample-by-variable data frame.
#' Column names should include two ID variables and an IBD variable.
#'
#' 5) A vector of the 2 ID variable names in the relatedness data frame.
#' (e.g., rid.vars=c("ID1","ID2")).
#'
#' 6) A string identifying the IBD variable name (e.g., ibd.var="PropIBD").
#'
#' 7) Optional: within.sex = FALSE. Default is within.sex = TRUE
#' and will permute males and females separately.
#'
#' 8) If within.sex = TRUE (the default), a string identifying the
#' sex variable name (e.g., sex.var="Inferred_Sex").
#'
#' 9) If within.sex = TRUE (the default), male and female values
#' in the sex vector (e.g., male.val=1, female.val=2).
#'
#' 10) Optional: a seed for sampling. If a seed is not provided, one will be
#' chosen randomly during the sampling process (e.g., seed=123).
#'
#' 11) N.B. Any column names not specified in (2)-(9) are assumed to be
#' phenotypes or phenotypic covariates.
#'
#' @param df,rel.df Data frame
#' @param id.var,sex.var,ibd.var Strings
#' @param male.val,female.val Strings or integers
#' @param geno.vars,rid.vars Character vectors
#' @param seed Number
#' @param within.sex Boolean (TRUE/FALSE)
#' @returns Data frame
#' @export
#' @examples
#' sip_pair(df = sipPair_exampleData, id.var = "IID",
#' sex.var = "SEX", male.val = "M", female.val = "F",
#' geno.vars = c("FID","BATCH",paste0("PC",1:4)),
#' rel.df = sipPair_relatednessData, rid.vars=c("IID1","IID2"),
#' ibd.var="PropIBD")

sip_pair <- function(df = NULL, id.var = NULL,
                     sex.var = NULL, male.val = NULL,
                    female.val = NULL, geno.vars = NULL,
                    within.sex = TRUE, seed = NULL,
                    rel.df = NULL, rid.vars = NULL, ibd.var = NULL) {

  tryCatch(

    expr = {

      if (is.null(df) | is.null(id.var) |
          is.null(geno.vars) | is.null(rel.df) |
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

      if (!(TRUE %in% is.na(geno.vars))
          & (FALSE %in% (geno.vars %in% colnames(df)))) {
        stop("Error: Genotypic variables are not in df column names.
             Please input genotypic variable names as a vector of strings.")
      }

      # get phenotype and phenotypic covariate names #
      pheno.vars <- colnames(df)[!(colnames(df) %in%
                                      c(id.var, sex.var, geno.vars))]

      # divide df by sex #
      if (within.sex) {

        if (is.null(sex.var) | is.null(male.val) | is.null(female.val)) {
          stop("Error: One or more sex-related inputs have null values.
             Please check your variable assignments or choose
             'within.sex = FALSE' to ignore sex while permuting.")
        }

        if (FALSE %in% (sex.var %in% colnames(df))) {
          stop("Error: Sex variable is not in df column names.
             Please input a sex variable as a string or choose
             'within.sex = FALSE' to ignore sex while permuting.")
        }

        if (length(unique(df[[sex.var]])) > 2) {
          stop("Error: more than two sexes detected in the sex column.")
        }

        if (FALSE %in% ((unique(df[[sex.var]]))
                        %in% c(male.val, female.val))) {
          stop("Error: sex values do not match the male
             or female values provides. Please input proper sex values
             (e.g., male.val=1, female.val=2).")
        }

        if (female.val %in% unique(df[[sex.var]])) {
          fdf <- get_sexDF(df, sex.var, sex.val = female.val)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mdf <- get_sexDF(df, sex.var, sex.val = male.val)
        }

        # separate df into fixed and permutable vectors #
        if (female.val %in% unique(df[[sex.var]])) {
          ffix <- get_fixed(df = fdf, id.var = id.var, sex.var = sex.var,
                            geno.vars = geno.vars)
          fpheno <- get_permute(df = fdf, id.var = id.var,
                                pheno.vars = pheno.vars)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mfix <- get_fixed(df = mdf, id.var = id.var, sex.var = sex.var,
                            geno.vars = geno.vars)
          mpheno <- get_permute(df = mdf, id.var = id.var,
                                pheno.vars = pheno.vars)
        }

        # pair individuals by most closely related within sex #
        if (female.val %in% unique(df[[sex.var]])) {
          fpair <- get_relPairs(df = fdf, id.var = id.var, rel.df = rel.df,
                                rid.vars = rid.vars, ibd.var = ibd.var)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mpair <- get_relPairs(df = mdf, id.var = id.var, rel.df = rel.df,
                                rid.vars = rid.vars, ibd.var = ibd.var)
        }

        # get relatedness degree categories to permute within #
        if (female.val %in% unique(df[[sex.var]])) {
          fdeg <- get_relCat(pairs = fpair, rid.vars = rid.vars,
                             ibd.var = ibd.var)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mdeg <- get_relCat(pairs = mpair, rid.vars = rid.vars,
                             ibd.var = ibd.var)
        }

        # check for uneven number of individuals and create a dataset for #
        # permutation by single sample #
        if (female.val %in% unique(df[[sex.var]])) {
          fex <- get_singSamp(df = fdf, deglst = fdeg, id.var = id.var)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mex <- get_singSamp(df = mdf, deglst = mdeg, id.var = id.var)
        }

        # set seed #
        if (is.null(seed)) {
          seed <- sample(seq(999999), 1)
        }
        message(paste("Seed:",seed))
        set.seed(seed)

        # get paired permutations #
        if (female.val %in% unique(df[[sex.var]])) {
          fmtch1 <- get_pairPerm(deglst = fdeg, rid.vars = rid.vars,
                                 ibd.var = ibd.var, seed = seed)
        }

        if (male.val %in% unique(df[[sex.var]])) {
          mmtch1 <- get_pairPerm(deglst = mdeg, rid.vars = rid.vars,
                                 ibd.var = ibd.var, seed = seed)
        }

        # get single permutations #
        if (exists("fex")) {
          fex <- fex[order(fex[[id.var]]),]
          fmtch2 <- data.frame(fixid=fex[[id.var]],
                               permid=fex[[id.var]][get_permIdx(df = fex,
                                                                seed = seed)])
        }

        if (exists("mex")) {
          mex <- mex[order(mex[[id.var]]),]
          mmtch2 <- data.frame(fixid=mex[[id.var]],
                               permid=mex[[id.var]][get_permIdx(df = mex,
                                                                seed = seed)])
        }

        # reintegrate into permuted datasets #
        if (female.val %in% unique(df[[sex.var]])) {
          if (exists("fmtch2")) {
            fmtch <- rbind(fmtch1, fmtch2)
          } else {
            fmtch <- fmtch1
          }

          fperm <- get_pairPermDF(fix.df = ffix, perm.df = fpheno,
                                  perm.map = fmtch, id.var = id.var)
        }

        if (male.val %in% unique(df[[sex.var]])) {

          if (exists("mmtch2")) {
            mmtch <- rbind(mmtch1, mmtch2)
          } else {
            mmtch <- mmtch1
          }

          mperm <- get_pairPermDF(fix.df = mfix, perm.df = mpheno,
                                  perm.map = mmtch, id.var = id.var)
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

        # pair individuals by most closely related within sex #
        pair <- get_relPairs(df = df, id.var = id.var, rel.df = rel.df,
                              rid.vars = rid.vars, ibd.var = ibd.var)

        # get relatedness degree categories to permute within #
        deg <- get_relCat(pairs = pair, rid.vars = rid.vars,
                           ibd.var = ibd.var)

        # check for uneven number of individuals and create a dataset for #
        # permutation by single sample #
        ex <- get_singSamp(df = df, deglst = deg, id.var = id.var)

        # set seed #
        if (is.null(seed)) {
          seed <- sample(seq(999999), 1)
        }
        message(paste("Seed:",seed))
        set.seed(seed)

        # get paired permutations #
        mtch1 <- get_pairPerm(deglst = deg, rid.vars = rid.vars,
                               ibd.var = ibd.var, seed = seed)

        # get single permutations #
        if (exists("ex")) {
          ex <- ex[order(ex[[id.var]]),]
          mtch2 <- data.frame(fixid=ex[[id.var]],
                               permid=ex[[id.var]][get_permIdx(df = ex,
                                                                seed = seed)])
        }

        # reintegrate into permuted datasets #
        if (exists("mtch2")) {
          mtch <- rbind(mtch1, mtch2)
        } else {
          mtch <- mtch1
        }

        perm <- get_pairPermDF(fix.df = fix, perm.df = pheno,
                                perm.map = mtch, id.var = id.var)

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
