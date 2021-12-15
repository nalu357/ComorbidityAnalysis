#' Define regions aroung distinct association signals from all traits of interest
#'
#' @param file.lst list of paths for excel file containing indpendent signals for the traits of interest (columns: SNPID, CHR, POS)
#' @param path filepath to save regions for trait
#' @param wd window around signals, defines regions' size
#' @return regions to perform colocalization
#' @export
get.regions <- function(filename.lst, path, wd = 1e+6) {
  dt <- data.table::data.table(trait = traits.lst, file = file.lst)
  signals <- data.table::data.table()
  for (file in file.lst) {
    signals <- rbind(signals, data.table::as.data.table(openxlsx::read.xlsx(file, startRow = 2)))
  }
  regions <- unique(na.omit(signals))[, ":="(start = POS - wd, end = POS + wd)]
  data.table::fwrite(regions, paste0(path, "/GWAS_regions.csv"))
  return(regions)
}

#' Define regions aroung distinct association signals from all traits of interest
#'
#' @param gwas.file excluding indels (columns: snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases)
#' @param wd window around signals, defines regions' size
#' @return regions to perform colocalization
#' @export
read_GWAS <- function(gwas.file, dict, snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases) {
  data.table::fread(gwas.file, select=c(snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases), verbose = FALSE)
  data.table::setcolorder(data, c(snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases))
  data.table::setnames(data, c("snp", "chr", "pos", "ea", "nea", "eaf", "maf", "beta", "se", "pval", "n", "ncases"))
  if (snp == ".") { # add rsID
    return()
  }
}

# Extract GWAS variants within 1Mb from merged independent signals
select_GWAS <- function(gwas, regions) {
  clumped.gwas <- data.table::data.table()
  for (i in 1:nrow(signals)) {
    clumped.gwas <- rbind(clumped.gwas, gwas[CHR %in% regions$CHR[i] & POS %between% .(regions$start[i], regions$end[i])])
    clumped.gwas <- unique(clumped.gwas)
  }
  return(clumped.gwas)
}

# Add ID column and get list of IDs with corresponding traits
add_ID <- function(dict, trait, listID) {
  dict[[trait]] <- dict[[trait]][, ID := paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]
  dict[[trait]][, add_key_value(listID, ID, trait), by = ID]
}

# Flip alleles per ID
flip.alleles <- function(dict_GWAS, IDlist, logFile) {
  for (id in keys(IDlist)) {
    if (length(IDlist[[id]]) > 1) {
      ref <- IDlist[[id]][1]
      # dict_ref <- ifelse(has.key(ref, dict_eQTL), dict_eQTL[[ref]], dict_GWAS[[ref]])
      refEA <- dict_GWAS[[ref]][ID == id, EA]
      lapply(IDlist[[id]][-1], function(i) {
        # dict <- ifelse(i %in% TISSUE$tissue, dict_eQTL[[i]], dict_GWAS[[i]])
        altEA <- dict_GWAS[[i]][ID == id, EA]
        if (altEA != refEA) {
          cat(paste0("FLIPPING beta of ID=", id, " and trait=", i, " with REF=", altEA, " using ref=", ref, " with REF=", refEA),
            file = logFile, append = TRUE, sep = "\n"
          )
          dict_GWAS[[i]][ID == id, beta := -beta]
        }
      })
    }
  }
}

preprocess.gwas <- function(signals.lst, gwas.file, traits) {
  print("Loading independent GWAS signals")
  regions <- get.regions(signals.lst)

  print("Loading GWAS data")
  read.gwas()

  print("Selecting GWAS signals around merged independent signals")
  select.gwas()

  print("Adding ID for GWAS datasets")
  add.ID()

  print("Flipping alleles when needed")
  flip.alleles()

  print("Outputing GWAS files for colocalization")
  for (trait in traits) {
    data.table::fwrite(GWAS[[paste(trait, "signals", sep = "_")]], paste0("/project_data/processed_data/GWAS", paste(traits, collapse = "_"), "/GWAS_", trait, "_precoloc_regions.csv"))
  }
}
