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

#' Perform pairwise colocalization analysis on multiple loci using coloc.abf
#'
#' @param regions regions to perform colocalization (columns: CHR, start end)
#' @param trait1 summary stats of trait1 as a data.table (essential columns: )
#' @param trait2 summary stats of trait2 as a data.table (essential columns: )
#' @param cc_ratio1 case vs. control ratio for trait1
#' @param cc_ratio2 case vs. control ratio for trait2
#' @param pp4.thres PP4 threshold for defining colocalized regions (default=0.8)
#' @return result of colocalization for each region along with 95% credible set for regions that reach PP4 threshold
#' @export
GWAS.coloc <- function(trait1, trait2, cc_ratio1, cc_ratio2, regions, pp4.thres=0.8, output.path=paste0("/project_data/out/GWAS_", trait1, "_", trai2, "/credible_set.csv")){
  data <- data.table::merge(trait1, trait2, by=c("ID, chr, pos"), allow.cartesian = T)
  data[rsID.x=="", rsID.x:=rsID.y]
  data[rsID.y=="", rsID.y:=rsID.x]

  print("Colocalization analysis using coloc.abf")
  res <- vector(mod ="list", length=nrow(GWAS_regions))
  names(res) <- seq(1:nrow(regions))
  credset <- vector(mod ="list", length=nrow(regions))
  names(credset) <- seq(1:nrow(regions))

  for (i in 1:nrow(regions)){
    data_tmp <- data[chr==regions[i]$CHR & pos %between% c(regions[i]$start, regions[i]$end)]

    res[[i]] <- coloc.abf(dataset1=list(beta=data_tmp$beta.x, varbeta=(data_tmp$se.x)^2, type="cc", s=cc_ratio1, N=nrow(trait1)),
                          dataset2=list(beta=data_tmp$beta.y, varbeta=(data_tmp$se.y)^2, type="cc", s=cc_ratio2, N=nrow(trait2)))

    # Get 95% credible set for regions where PP4 > pp4.thres
    if (res[[i]]$summary[6]>pp4.thres){
      o <- order(res[[i]]$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(res[[i]]$results$SNP.PP.H4[o])
      w <- which(cs > 0.95)[1]
      credset[[i]] <- data_tmp[as.numeric(gsub("(SNP.)","",res[[i]]$results[o,][1:w,]$snp)), ID]
    }
  }
  # check how to handle multiple returns
  return(res, credset)
}

#' Output a concise data table for the 95% credible set of colocalized regions
#'
#' @param res list of colocalization result for each region
#' @param credset list of 95% credible set for each region
#' @return data table for 95% credible set with PP4 for each variant
#' @export
output.credset <- function(res, credset) {
  # output credible set as data.table for GWAS eQTL colocalization
  credset.dt <- data.table()
  for(i in 1:length(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))])){
    reg <- names(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][i])
    credible.set <- Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))][[i]]
    o <- order(res[[reg]]$results$SNP.PP.H4, decreasing=TRUE)

    credset.dt <- rbind(credset.dt, data.table(region=rep(reg, region.n), SNP=credible.set,
                                               rsID=data[, rsID.x[match(credible.set, ID)]],
                                               CHR=data[, CHR[match(credible.set, ID)]],
                                               POS=data[, POS[match(credible.set, ID)]],
                                               PP4=rep(round(unname(res[[reg]]$summary[paste0("PP.H4.abf")]), digits=3), region.n),
                                               SNP.PP4=res[[reg]]$results$SNP.PP.H4[o] %>% .[1:region.n]))
  }
  fwrite(credset.dt, output.path)
}

# named list of gwas (name=trait)
# named list of mQTL for tissues of interest (name=tissue)
GWAS.mQTL.coloc <- function(gwas, mqtl, credset, regions, range=1e+6){
  res <- vector(mod ="list", length=length(unique(credset$region)))
  names(res) <- unique(credset$region)

  for (reg in unique(credset.dt$region)){
    # ----------- Select only variants included in the credible set -----------
    sel.snp <- credset[region==reg, max(PP4)]$rsID
    lapply(gwas, function(i){i <- i[ID %in% sel.snp & POS %between% c(sel.snp-range, sel.snp+range)]})
    lapply(mQTL, function(i){i <- i[ID %in% sel.snp & geneID %in% genes & POS %between% c(sel.snp-range, sel.snp+range)]})
    # TRAIT[, select_variants(trait, GWAS, credset.dt[region==reg], mQTL=FALSE), by=trait]
    # TISSUE[, select_variants(tissue, eQTL, credset.dt[region==reg], mQTL=TRUE), by=tissue]

    for (tissue in names(mQTL)){
      # mQTL[[tissue]] <- mQTL[[tissue]][ID %in% credset[region==reg, SNP]]
      for (gene in unique(mQTL[[tissue]]$geneID)){
        # ----------- Prepare beta and ses matrices -----------
        print("Preparing input matrices")
        ses <- gwas[[1]][, .(SNPID=ID, rsID=rsID, se=se)]
        data.table::setnames(ses, "se", paste0("GWAS_", names(gwas)[1]))
        betas <- gwas[[1]][, .(SNPID=ID, rsID=rsID, beta=beta)]
        data.table::setnames(betas, "beta", paste0("GWAS_", names(gwas)[1]))

        for (trait in names(gwas)[-1]) {
          ses <- data.table::merge(ses, gwas[[trait]][, .(SNPID=ID, se=se)], by="SNPID", allow.cartesian=T)
          data.table::setnames(ses, "se", paste0("GWAS_", trait))
          betas <- data.table::merge(betas, gwas[[trait]][, .(SNPID=ID, beta=beta)], by="SNPID", allow.cartesian=T)
          data.table::setnames(betas, "beta", paste0("GWAS_", trait))
        }

        ses <- data.table::merge(ses, mQTL[[tissue]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), se=se)], by="SNPID", allow.cartesian=T)
        ses <- data.table::dcast(ses, as.formula(paste(paste(colnames(ses)[c(-ncol(ses), -ncol(ses)+1)], collapse="+"), "~ geneID")), value.var = "se")
        betas <- data.table::merge(betas, mQTL[[tissue]][geneID==gene, .(SNPID=ID, geneID=paste0("eQTL_", tissue, "_", geneID), beta=beta)], by="SNPID", allow.cartesian=T)
        betas <- data.table::dcast(betas, as.formula(paste(paste(colnames(betas)[c(-ncol(betas), -ncol(betas)+1)], collapse="+"), "~ geneID")), value.var = "beta")

        betas <- data.table::as.data.table(betas)
        ses <- data.table::as.data.table(ses)

        # ----------- Run HyPrColoc -----------
        print("Colocalization analysis using HyPrColoc")
        id <- betas$SNPID
        rsid <- betas$rsID
        betas_mat <- as.matrix(betas[, c('SNPID', 'rsID'):=NULL])
        rownames(betas_mat) <- id
        ses_mat <- as.matrix(ses[, c('SNPID', 'rsID'):=NULL])
        rownames(ses_mat) <- id
        binary.traits = c(rep(1,nrow(TRAIT)), rep(0,ncol(betas_mat)-nrow(TRAIT)))
        traits <- colnames(betas_mat)

        betas_mat <- na.omit(betas_mat)
        ses_mat <- na.omit(ses_mat)

        if (nrow(betas_mat)>1){
          res[[reg]] <- hyprcoloc::hyprcoloc(betas_mat, ses_mat, trait.names=traits, snp.id=id, bb.alg=FALSE,
                                             binary.outcomes=binary.traits, prior.1=1e-10, prior.2=0.7, snpscores=T)
          print(res[[reg]])
          res.sen = hyprcoloc::sensitivity.plot(betas_mat, ses_mat, trait.names = traits, snp.id=id, bb.alg=FALSE,
                                                prior.2 = c(0.1, 0.3, 0.5, 0.7, 0.9), equal.thresholds = FALSE,
                                                binary.outcomes=binary.traits, prior.1=1e-10)

        }
      }
    }
  }
}

select_variants <- function(dt, credset, mQTL=FALSE) {
  credset.id <- credset[, SNP]
  if (length(credset.id)!=1) {
    dict[[paste0(i, "_coloc")]] <- dict[[i]][ID %in% credset.id]
  }
  else{
    if(mQTL==TRUE) {
      genes <- dict[[i]][ID==credset.id, geneID]
      # dict[[paste0(i, "_coloc")]] <- dict[[i]][geneID %in% genes & POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
      dict[[paste0(i, "_coloc")]] <- dict[[i]][geneID %in% genes & POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
    }
    else {
      dict[[paste0(i, "_coloc")]] <- dict[[i]][POS %between% c(credset[, POS]-1000000, credset[, POS]+1000000)]
    }
  }
}

###################################################################
# --------------- Get LD between SNPs in each region --------------
###################################################################
# USE SOME R PACKAGE OR DO IT IN PLINK
for (i in 1:length(Filter(Negate(is.null), credset)[!duplicated(Filter(Negate(is.null), credset))])){
  region.data <- data[CHR==GWAS_regions[region]$CHR & POS %between% c(GWAS_regions[region]$POS - range, GWAS_regions[region]$POS + range) & rsID.x!=".", rsID.x]
  chr <- data[rsID.x %in% region.data, unique(CHR)]
  fwrite(as.list(region.data), paste0("/project_data/processed_data/LDvariants/", TRAIT$trait[2], "/region", region, ".chr", chr), sep="\n")
}


###############################################################
# --------------- Plot regional association plot --------------
###############################################################
source("~/gwas_eqtl_colocalization/scripts/PlotFunctions.R")

ld.regions <- list("3965689"=4291928, "53501946"=53800954, "150521096"=150537635, "133414622"=133864599, "124468572"=123732769, # "124509177"=123450765,
                   "124509177"=123732769, "51180765"=50788778, "44938870"=45411941, "422144"=653575,"10808687"=9974824) # , "9974824"=10808687, "10808687"=9974824)

# read cred.set data table
# data = merge of trait1 and trait2
# function(credset, data, range=1e+6, ld.file, GRCh37_Genes)
# check if any logP==inf --> reduce it
ld.pos <- GWAS_regions[as.numeric(region.nr)]$POS
ld.file <- fread(paste0("/project_data/processed_data/LDvariants/chr", sel.chr , "_", ifelse(ld.pos %in% names(ld.regions),
                                                                                             as.numeric(unname(ld.regions[match(ld.pos, names(ld.regions))])),
                                                                                             ld.pos), ".ld"))[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]

for (reg in unique(credset$regions)){
  sel.snp <- credset[region==reg, max(PP4)]$rsID
  if (is.null(sel.snp)) print(paste0("ERROR: no rsID available for the lead SNP ", sel.snp))
  sel.chr <- credset[SNP==sel.snp, CHR]
  sel.pos <- credset[SNP==sel.snp, POS]
  sel.PP4 <- credset[SNP==sel.snp, PP4]
  credible.set <- as.vector(credset[region==reg, rsID])

  # Get data
  trait1.region <- trait1[CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range) & !is.null(rsID.x),
                          .(CHR=CHR, SNP=rsID.x, P=pval.y, BP=POS, logP=-log10(as.numeric(pval.y)))]
  trait2.region <- unique(trait2[CHR==sel.chr & between(POS, sel.pos-range, sel.pos+range),
                                 .(CHR=CHR, SNP=rsID.x, P=pval.x, BP=POS, logP=-log10(as.numeric(pval.x)))],
                          by="SNP")

  # Regional association plot
  data.lst <- list(trait2, trait1)
  names(data.lst) <-  c(trait2, trait1)
  locus.zoom(data = data.lst,
             offset_bp = range,
             genes.data = GRCh37_Genes,
             file.name = paste0(trait1, "_", trait2, "_T2D_", region.nr, "_beta.png"),
             secondary.snp = ifelse(credible.set<2, NA, credible.set),
             snp=sel.rsid,
             ignore.lead=TRUE,
             ld.file=ld.file,
             pp="PP4",
             pp.value=round(unname(res[[region.nr]]$summary["PP.H4.abf"]), digits=3),
             nplots=TRUE)

  # Plot PP4 of credible set
  locus.zoom(data = credset[region==reg, .(CHR=CHR, SNP=SNP, PP4=PP4, BP=POS, logP=PP4)],
             region = c(sel.chr, credset[region==reg, min(POS)]-100,credset[region==reg, max(POS)]+100),
             genes.data = GRCh37_Genes,
             file.name = paste0(trait1, "_", trait2, "_PP4_", region.nr, "_beta.png"),
             snp=sel.rsid,
             ignore.lead=TRUE,
             ld.file=ld.file,
             sig.type="PP4")
}

