

downloadSupplementaryFiles <- function(gse) {
  if (!file.exists("/home/rstudio/projects/GSE152939/GSE152939_I138fs_raw_counts.txt.gz")) {
    supplementaryFiles <- GEOquery::getGEOSuppFiles(gse)
  } else {
    supplementaryFiles <- data.frame(row = c(1, 2, 3, 4))
    rownames(supplementaryFiles) <- c(
      "/home/rstudio/projects/GSE152939/GSE152939_I138fs_norm_counts_5CPM.txt.gz",
      "/home/rstudio/projects/GSE152939/GSE152939_I138fs_raw_counts.txt.gz",
      "/home/rstudio/projects/GSE152939/GSE152939_K88N_norm_counts_5CPM.txt.gz",
      "/home/rstudio/projects/GSE152939/GSE152939_K88N_raw_counts.txt.gz"
    )
  }

  return(supplementaryFiles)
}


readFile <- function(supplementaryFiles) {
  filenames <- rownames(supplementaryFiles)
  CRX_ExperimentRawCount <- read.delim(filenames[2], header=TRUE, check.names = FALSE)
  return (CRX_ExperimentRawCount)
}

filterOutLowCount <- function(CRX_ExperimentRawCount) {
  countPerMillion <- edgeR::cpm(CRX_ExperimentRawCount[3:28])
  rownames(countPerMillion) <- CRX_ExperimentRawCount[,1]
  keep = rowSums(countPerMillion >1) >= 3
  CRX_ExperimentRawCountFiltered <- CRX_ExperimentRawCount[keep,]
  return (CRX_ExperimentRawCountFiltered)
}

normalizeTheCounts <- function(CRX_ExperimentRawCountFiltered) {
  samples <- data.frame(lapply(colnames(CRX_ExperimentRawCountFiltered)[3:28], FUN=function(x){unlist(strsplit(x, split = "\\_"))[c(1, 2, 3)]}))
  colnames(samples) <- colnames(CRX_ExperimentRawCountFiltered)[3:28]
  rownames(samples) <- c("condition", "time","patient")
  samples <- data.frame(t(samples))

  filtered_data_matrix <- as.matrix(CRX_ExperimentRawCountFiltered[,3:28])
  rownames(filtered_data_matrix) <- CRX_ExperimentRawCountFiltered$ens.id
  d = edgeR::DGEList(counts=filtered_data_matrix, group=samples$condition)
  d = edgeR::calcNormFactors(d)
  normalized_counts <- edgeR::cpm(d)
  return (normalized_counts)
}

getEnsemblbiomart <- function() {
  ensemblMart <- biomaRt::useMart("ensembl")
  ensemblDataSet <-  biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensemblMart)

  return(ensemblDataSet)
}

mapTheData <- function(CRX_ExperimentRawCountFiltered, normalized_counts, ensemblDataSet) {
  #print(any(duplicated(CRX_ExperimentRawCountFiltered$ens.id)))

  conversion_stash <- "CRX_id_conversion.rds"
  if(file.exists(conversion_stash)){
    CRX_id_conversion <- readRDS(conversion_stash)
  } else {
    CRX_id_conversion <- biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                                        filters = c("ensembl_gene_id"),
                                        values = CRX_ExperimentRawCountFiltered$ens.id,
                                        mart = ensemblDataSet)
    saveRDS(CRX_id_conversion, conversion_stash)
  }

  normalized_counts_annot <- merge(CRX_id_conversion, normalized_counts,
                                   by.x = 1, by.y = 0, all.y=TRUE)

  #print(any(duplicated(normalized_counts_annot$ensembl_gene_id)))

  return(normalized_counts_annot)
}











