#' @title Look up genes in DEGs
#' @description Check if any gene in your list is a DEG for any trait
#'
#' @param DEG named list of data tables containing DEGs for a trait (columns: gene, ID, logFC, pval.adj)
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/DEG.csv"
#' @param logFC.theta threshold for abs(logFC), default=log(1.5, base=2)
#' @param alpha significance threshold for p-value, default=0.05
#' @return a data table of genes with 1 or 0 in each DEG column
#' @export
#'
#' @importFrom data.table fwrite
pairwise <- function(){

}

#' @title Look up genes in DEGs
#' @description Check if any gene in your list is a DEG for any trait
#'
#' @param DEG named list of data tables containing DEGs for a trait (columns: gene, ID, logFC, pval.adj)
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/DEG.csv"
#' @param logFC.theta threshold for abs(logFC), default=log(1.5, base=2)
#' @param alpha significance threshold for p-value, default=0.05
#' @return a data table of genes with 1 or 0 in each DEG column
#' @export
#'
#' @importFrom data.table fwrite
multi.trait <- function(){

}
