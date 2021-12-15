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
#' @examples
lookup.DEG <- function(DEG, genes, file.path = "lookups/DEG.csv", logFC.theta = log(1.5, base = 2), alpha = 0.05) {
  stopifnot(is.data.frame(genes), length(DEG) >= 1, is.numeric(logFC.theta), is.numeric(alpha))
  genes <- unique(genes, by = "ID")
  for (i in length(DEG)) {
    DEG[[i]][abs(logFC) > logFC.theta & pval.adj <= alpha]
    dt <- DEG[[i]][ID %in% genes$ID, .(gene, ID, direction = ifelse(logFC > 0, "+", "-"))]
    genes <- genes[, eval(names(DEG)[i]) := ifelse(ID %in% dt$ID, 1, 0)]
  }
  fwrite(genes, file.path)
}


#' @title Search OMIM
#' @description Get OMIM's entry for a gene
#'
#' @param gene gene symbol
#' @return data table of OMIM's entry for the gene
#' @export
#'
#' @importFrom httr GET content stop_for_status content_type
#' @examples
search.OMIM <- function(gene) {
  stopifnot(is.character(gene), length(gene)==1)
  gene <- tolower(gene)
  server <- "https://api.omim.org"
  r <- GET(paste(server, "/api/geneMap/search?search=", gene, "&include=geneMap&apiKey=cBg3cWK0TxibU0wvtkwtMA&format=json&start=0&limit=100", sep = ""), content_type("application/json"))
  stop_for_status(r)
  if (length(jsonlite::fromJSON(jsonlite::toJSON(content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList) != 0) {
    all.phen <- jsonlite::fromJSON(jsonlite::toJSON(content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList[[1]]$phenotypeMap$phenotype
    phen <- all.phen[sapply(all.phen, function(i) {
      !grepl("\\{|\\[|\\?", i)
    })]
    if (length(phen) == 0) {
      return(NA)
    } else {
      return(unlist(phen))
    }
  }
  else {
    return(NA)
  }
}

#' @title Look up genes in OMIM
#' @description if any gene in your list has an entry in OMIM showing relevant phenotypes for a trait of interest
#'
#' @param terms named list of terms related to each trait of interest
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/OMIM.csv"
#' @return a data table of genes with 1 or 0 for each trait column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples
lookup.OMIM <- function(terms.lst, genes, file.path = "lookups/OMIM.csv") {
  stopifnot(length(terms.lst)>=1, is.data.frame(genes), nrow(genes)>0)
  genes <- unique(genes, by = "ID")
  genes$OMIM <- sapply(genes$name, search.OMIM)

  for (trait in names(terms.lst)) {
    pattern <- paste(terms.lst[[trait]], collapse = "|")
    genes <- genes[, eval(trait) := as.integer(grepl(pattern, OMIM, ignore.case = TRUE)), by = name]
  }

  fwrite(genes[, OMIM := NULL], file.path)
}

#' @title Look up genes in high confidence genes data
#' @description Check if any gene in your list is already a high confidence likely effector gene for any of the traits of interest
#'
#' @param hc named list of high confidence genes for each trait of interest
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/HC.csv"
#' @return a data table of genes with 1 or 0 for each trait column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples
lookup.hc <- function(hc.genes, genes, file.path = "lookups/HC.csv") {
  stopifnot(is.data.frame(genes), length(hc.genes) >= 1, nrow(genes)>0)
  genes <- unique(genes, by = "ID")

  for (trait in names(hc.genes)) {
    genes <- genes[, eval(trait) := ifelse(name %in% hc.genes[[trait]], 1, 0), by = name]
  }

  fwrite(genes, file.path)
}
