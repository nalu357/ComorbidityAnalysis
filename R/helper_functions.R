sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
})

write_output <- function(var, dict, filename){
  # var can be one variable or a list
  fwrite(dict[[paste(var, collapse="_")]], filename)
}

add_key_value <- function(hash, key, new.value){
  if(has.key(key, hash)){
    hash[key] <- append(hash[[key]], new.value)
  }
  else{
    .set(hash, key, new.value)
  }
}

annotate_rsids <- function(trait, trait_dt, rsid){
  # read data
  # rsid <- fread(pste0("/project_data/processed_data/Variants2rsID/", trait,".id.rsid.bed"), verbose=FALSE)
  # trait_dt <- fread(paste0("/project_data/processed_data/GWASKneeOA/GWAS_", trait,"_precoloc_regions"), verbose=FALSE)

  # rename columns
  colnames(rsid) <- c("CHR:POS", "CHR", "POS", "rsID", "A1", "A2", "ID")

  # remove indels
  new_rsid <- rsid[nchar(A1)==1]

  # separate A2 column when there is a komma
  new_rsid <- separate(new_rsid, "A2", c("A2a", "A2b"), sep=",")
  new_rsid <- melt(new_rsid, measure.vars=c("A2a","A2b"), variable.name="variable", value.name="A2")

  # add IDs to trait_rsid
  new_rsid[, ID:=paste(`CHR:POS`, sort_alleles(A1, A2), sep="_")]

  # merge rsid into trait data.table
  new_trait_dt <- merge(trait_dt[, .SD, .SDcols = !"rsID"], new_rsid[, .(ID, POS=as.integer(POS), CHR=as.integer(CHR), rsID)], by=c("ID", "CHR", "POS"), all.x=TRUE)

  # save anotated dataset
  fwrite(new_trait_dt[, `CHR:POS`:=NULL], paste0("/project_data/processed_data/GWASKneeOA/GWAS_", trait,"_rsid_precoloc_regions.csv"), verbose=FALSE)
}


################################################
#--------------- Add rsID to T2D --------------
################################################
# T2D <- fread("/storage/hmgu/projects/OA_T2D/data_original/Mahajan.NatGenet2018b.T2D.European.txt.gz")
# T2D.UKBB <- fread("/storage/hmgu/projects/OA_T2D/data_original/Mahajan.NatGenet2018b.UKBB.HRC.T2D.European.txt.gz")
# T2D.rsid <- fread(paste0("/project_data/processed_data/Variants2rsID/T2D.rsid"))
# T2D.sumstats <- T2D.UKBB[, .(CHR=Chr, POS=Pos, rsID=rep(".", nrow(T2D.UKBB)), pval=as.numeric(Pvalue), beta=as.numeric(Beta),
#                      se=as.numeric(SE), MAF=fifelse(EAF < (1-EAF), EAF, (1-EAF)), EAF=EAF, N=N,
#                      Ncases=as.integer(74124), EA=EA, NEA=NEA, ID=paste(`#SNP`, sort_alleles(EA, NEA), sep="_"))]
# tst <- merge(T2D.sumstats, T2D.rsid, by="ID", all.x=TRUE) %>% .[, rsID:=V4] %>% .[, .(CHR,POS,rsID,pval,beta,se,MAF,EAF,N,Ncases=".",EA,NEA,ID)]
# fwrite(tst, paste0("/project_data/processed_data/T2D_UKBB_rsid.csv"))

###########################################################
#--------------- Add rsID to GO of UKBB only --------------
###########################################################
# OAphen <- c("AllOA", "KneeOA", "KneeHipOA", "HipOA", "TJR", "THR") # "TKR"
# for (phen in OAphen){
#   # phen <- "TKR"
#   OA <- fread(get.GO.filename(phen, subset="UKBB")) %>% .[, ID:=paste(CPTID, sort_alleles(EA, NEA), sep="_")]
#   rsid <- fread(paste0("/project_data/", phen, ".simple.bed")) %>% .[, ID:=paste(V1, sort_alleles(V5, V6), sep="_")]
#   result <- merge(OA, rsid, by="ID", all.x=TRUE) %>% .[, .(CHR=CHR,POS=POS,rsID=V4,pval=P,beta=BETA,se=SE,EAF=EAF,N=NCASES+NCONTROLS,Ncases=NCASES,Ncontrols=NCONTROLS,EA=EA,NEA=NEA,ID)]
#   fwrite(result, paste0("/project_data/processed_data/", phen, "_UKBB_rsid.csv"))
# }

######################################################################
#--------------- Add EAF eQTL summary stats --------------
######################################################################
# for (tissue in c("PancreaticIslets")) {  # "Synovium", "HighGradeCartilage", "LowGradeCartilage",
#   dt <- fread(paste0("/project_data/processed_data/eQTL_", tissue,".csv"))
#   if (tissue=="PancreaticIslets") {
#     raw <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/InsPIRE_PancreaticIslets_Gene_eQTLs_Nominal_Pvalues.txt.gz"))
#     col_names <- colnames(raw)[-1]
#     raw <- raw[, !"FreqALT"]
#     setnames(raw, col_names)
#     raw <- raw[, `:=` (ID=paste(paste(SNPchr, SNPposition, sep=":"), sort_alleles(REF, ALT), sep="_"), geneID=GeneID)]
#     result <- merge(dt, raw, by=c("ID", "geneID"), all.x=TRUE) %>% .[, .(CHR,POS,rsID.x,geneID,pval,beta,se,MAF,EAF=FreqALT,EA,NEA,ID)]
#   }
#   else{
#     raw <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/FunGen_eQTL_", tissue, ".FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz")) %>%
#       .[, `:=` (ID=paste(genotype_id, sort_alleles(REF, ALT), sep="_"), geneID=sub(".*_", "", phenotype_id))]
#     result <- merge(dt, raw, by=c("ID", "geneID"), all.x=TRUE) %>%
#       .[, EAF:=as.numeric(str_match(INFO, "AF=\\s*(.*?)\\s*;")[,2])] %>%
#       .[, .(CHR,POS,rsID,geneID,pval=pval.x,beta,se,MAF,EAF,N,EA,NEA,ID)]
#   }
#   fwrite(result, paste0("/project_data/processed_data/eQTL_", tissue,"_EAF.csv"))
#  }


######################################################################
#--------------- Add EAF eQTL summary stats --------------
######################################################################
# for (tissue in c("PancreaticIslets")) {  # "Synovium", "HighGradeCartilage", "LowGradeCartilage",
#   dt <- fread(paste0("/project_data/processed_data/eQTL_", tissue,".csv"))
#   if (tissue=="PancreaticIslets") {
#     raw <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/InsPIRE_PancreaticIslets_Gene_eQTLs_Nominal_Pvalues.txt.gz"))
#     col_names <- colnames(raw)[-1]
#     raw <- raw[, !"FreqALT"]
#     setnames(raw, col_names)
#     # raw <- raw[, `:=` (ID=paste(paste(SNPchr, SNPposition, sep=":"), sort_alleles(REF, ALT), sep="_"), geneID=sub("\\..*", "",GeneID))]
#     raw <- raw[, geneID:=sub("\\..*", "",GeneID)]
#     result <- merge(dt, raw, by=c("rsID", "geneID"), all.x=TRUE) %>% .[rsID!=".", .(CHR,POS,rsID,geneID,pval,beta,se,MAF,EAF=FreqALT,EA,NEA,ID)]
#   }
#   else{
#     raw <- fread(paste0("/storage/hmgu/projects/OA_T2D/data_original/FunGen_eQTL_", tissue, ".FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz")) %>%
#       .[, `:=` (ID=paste(genotype_id, sort_alleles(REF, ALT), sep="_"), geneID=sub(".*_", "", phenotype_id))]
#     result <- merge(dt, raw, by=c("ID", "geneID"), all.x=TRUE) %>%
#       .[, EAF:=as.numeric(str_match(INFO, "AF=\\s*(.*?)\\s*;")[,2])] %>%
#       .[, .(CHR,POS,rsID,geneID,pval=pval.x,beta,se,MAF,EAF,N,EA,NEA,ID)]
#   }
#   fwrite(result, paste0("/project_data/processed_data/eQTL_", tissue,"_EAF.csv"))
# }
