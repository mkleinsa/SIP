#' @title
#' Phenotype correlation plot
#'
#' @description
#' This function returns a plot showing correlation
#' between phenotypes. (N.B., please omit missing values
#' before running this function.)
#'
#' @details
#' The function requires the following inputs:
#'
#' 1) A sample-by-variable data frame with phenotypes.
#' Column names should include phenotype names.
#'
#' 2) A vector of phenotype names.
#' (e.g., pheno.vars=c("PHENO1","PHENO2",...)).
#'
#' @param df Data frame
#' @param pheno.vars Character vector
#' @returns Plot
#' @export
#' @examples
#' plot_phenotype_correlations(df = sip_exampleData,
#' pheno.vars = paste0("PHENO",1:500))

plot_phenotype_correlations <- function(df = NULL, pheno.vars = NULL) {

  tryCatch(

    expr = {

        if (is.null(df) | is.null(pheno.vars)) {
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

        if (!(TRUE %in% is.na(pheno.vars))
            & (FALSE %in% (pheno.vars %in% colnames(df)))) {
          stop("Error: Phenotype names are not in df column names.
             Please input phenotype variable names as a vector of strings.")
        }

      # calculate correlation #
      pcor <- stats::setNames(as.data.frame(
        stats::cor(as.data.frame(df)[,pheno.vars])),pheno.vars)

      # reformat #
      pcor[["VAR1"]] <- pheno.vars
      pcor <- stats::setNames(data.table::melt(data.table::setDT(pcor),
                                        id.vars = "VAR1"),
                       c("VAR1","VAR2","correlation"))

      # plot #
      pcor <- pcor[order(pcor[["correlation"]], decreasing = TRUE),]
      pcor[["VAR1"]] <- factor(pcor[["VAR1"]], levels = unique(pcor[["VAR1"]]))
      pcor[["VAR2"]] <- factor(pcor[["VAR2"]], levels = unique(pcor[["VAR1"]]))

      ggplot2::ggplot(data = pcor,
                      aes(x = VAR1, y = VAR2, fill = correlation)) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw(base_size = 12,base_family="Arial") +
        ggplot2::theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.title = element_text(hjust=0.5, size = 14),
              plot.subtitle = element_text(hjust=0.5),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              legend.title = element_text(size = 10),
              legend.position = "right") +
        ggplot2::scale_fill_viridis_c()

      }, error = function(e){
        message(
          sprintf("An error occurred at %s : %s",
                  Sys.time(),
                  e)
       )
    }
  )
}

