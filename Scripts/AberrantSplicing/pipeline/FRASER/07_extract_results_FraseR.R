#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "07_results.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'   - outputDir: '`sm cfg.getProcessedResultsDir() + "/aberrant_splicing/datasets/"`'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'   - assemblyVersion: '`sm cfg.genome.getBSGenomeVersion()`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - spliceTypeSetup: '`sm str(projectDir / ".drop" / "helpers" / "spliceTypeConfig.R")`'
#'   - calculatePhi: '`sm str(projectDir / ".drop" / "helpers" / "calculatePhi.R")`'
#'   - addSpliceType: '`sm str(projectDir / ".drop" / "helpers" / "spliceType_frameshift_annotation.R")`'
#'   - subtypes: '`sm str(projectDir / ".drop" / "helpers" / "subtypes_exonSkipping_inconclusive.R")`'
#'   - annotate_blacklist: '`sm str(projectDir / ".drop" / "helpers" / "annotate_blacklist.R")`'
#'   - blacklist_19: '`sm str(projectDir / ".drop" / "helpers" / "resource" / "hg19-blacklist.v2.bed.gz")`'
#'   - blacklist_38: '`sm str(projectDir / ".drop" / "helpers" / "resource" / "hg38-blacklist.v2.bed.gz")`'
#'   - fdsin: '`sm cfg.getProcessedDataDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                 "padjBetaBinomial_theta.h5"`'
#'   - txdb: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/txdb.db"`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/gene_name_mapping_{annotation}.tsv"`'
#'   - OMIM: 'Data/project_data/external/omim_genes.tsv'
#'  output:
#'   - resultTableJunc: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv"`'
#'   - resultTableGene: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv"`'
#'   - fds: '`sm cfg.getProcessedResultsDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS"`'
#'  type: script
#'---

#'
# snakemake <- readRDS('/data/nasif12/home_if12/yepez/workspace/cagi6_sk/.drop/tmp/AS/all--cagi6/07_results.Rds')
saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
source(snakemake@input$add_HPO_cols)
source(snakemake@input$calculatePhi)
source(snakemake@input$spliceTypeSetup, echo=FALSE)
source(snakemake@input$addSpliceType, echo=FALSE)
source(snakemake@input$subtypes, echo=FALSE)
source(snakemake@input$annotate_blacklist)
library(AnnotationDbi)

getResultByGene <- function(res){
  res[, genomicLocation := paste0(seqnames, ":", start, "-", end, ":", as.character(strand))]
  res[, c("start", "end", "width", "strand", "STRAND", "PAIRED_END") := NULL]
  
  res_gene <- res[, lapply(.SD, function(x){
    paste(x, collapse=",")
    }), by="hgncSymbol,sampleID,seqnames,numSamplesPerGene,numEventsPerGene"]
  return(res_gene)
}

opts_chunk$set(fig.width=12, fig.height=8)

annotation    <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
outputDir  <- snakemake@params$outputDir
assemblyVersion <- snakemake@params$assemblyVersion

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds_input <- loadFraserDataSet(dir=workingDir, name=dataset)
fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(dataset, annotation, sep = '--'), rewrite = TRUE)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)

# Read annotation and match the chr style
txdb <- loadDb(snakemake@input$txdb)
orgdb <- fread(snakemake@input$gene_name_mapping)

seqlevelsStyle(orgdb$seqnames) <- seqlevelsStyle(fds)
seqlevelsStyle(txdb) <- seqlevelsStyle(fds)

# Annotate the fds with gene names
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, feature = 'gene_name',
                              featureName = 'hgnc_symbol', keytype = 'gene_id')

# remove globin and NA genes
fds <- fds[grep('^HB', rowData(fds)$hgnc_symbol, invert = T),]
fds <- fds[!is.na(rowData(fds)$hgnc_symbol),]

# Compute phi and extract rho
fds <- calculatePhi(fds)
phi_melt <- melt(cbind(as.data.table(granges(fds))[,1:3], as.data.table(assays(fds)$phi)), 
                 id.vars = c('seqnames', 'start', 'end'), variable.name = 'sampleID', value.name = 'phi')
rho_psi <-  as.data.table(rowRanges(fds, type = 'j'))[, .(seqnames, start, end, rho_psi5, rho_psi3)]
rho_theta <- as.data.table(rowRanges(fds, type = 'ss'))[, .(seqnames, start, end, rho_theta)]


# Add the junction annotations to the fds
fds <- createFDSAnnotations(fds, txdb)
mcols(fds, type="ss")$annotatedJunction <- "none"

saveFraserDataSet(fds, dir=outputDir)

# Extract and save results per junction
res_junc_dt <- results(fds, padjCutoff=snakemake@params$padjCutoff, 
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff,
                    additionalColumns = 'annotatedJunction') %>% as.data.table()
print('Results per junction extracted')
fwrite(res_junc_dt, gsub('.tsv', '_all.tsv', snakemake@output$resultTableJunc), sep = '\t', quote = F)

# Remove phi
res_junc_dt <- merge(res_junc_dt, phi_melt, by = c('sampleID', 'seqnames', 'start', 'end'), all.x = T, sort = F)
res_junc_dt <- res_junc_dt[!(psiValue > .8 & (phi < .2 | is.na(phi)))]

# Remove rho
res_junc_dt <- merge(res_junc_dt, rho_theta, by = c('seqnames', 'start', 'end'), all.x = T, sort = F)
res_junc_dt <- merge(res_junc_dt, rho_psi, by = c('seqnames', 'start', 'end'), all.x = T, sort = F)
res_junc_dt[, rho_filtered := ((type == 'psi3' & rho_psi3 > .1) | (type == 'psi5' & rho_psi5 > .1) | (type == 'theta' & rho_theta > .1))]
res_junc_dt <- res_junc_dt[rho_filtered == F]

# Add blacklist regions
blacklist_input <- ifelse(assemblyVersion %in% c('hg19', 'hs37d5', '37'), 
                          snakemake@input$blacklist_19,
                          snakemake@input$blacklist_38)
blacklist_gr <- import(blacklist_input, format = "BED")
res_junc_dt <- addBlacklistLabels(res_junc_dt, blacklist_gr, "splicing")
res_junc_dt <- res_junc_dt[blacklist == F]

# Add OMIM genes
ot <- fread(snakemake@input$OMIM)
res_junc_dt <- merge(res_junc_dt, ot[, .(SYMBOL, OMIM, OMIM_Inheritance)], 
                     by.x = 'hgncSymbol', by.y = 'SYMBOL', all.x = T, sort = F)


# Calculate splice types and frameshift
res_junc_dt <- aberrantSpliceType(res_junc_dt, fds, txdb)

# Add the subtypes for exonSkipping and inconclusive
res_junc_dt <- checkExonSkipping(res_junc_dt, txdb)
res_junc_dt <- checkInconclusive(res_junc_dt, txdb)

# Add UTR labels
# res_junc_dt <- addUTRLabels(res_junc_dt, txdb)

# Add N of outliers
res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]

# Aggregate results by gene
if(nrow(res_junc_dt) > 0){
  res_genes_dt <- getResultByGene(res_junc_dt)
  
  # add HPO overlap information
  sa <- fread(snakemake@config$sampleAnnotation)
  if(!is.null(sa$HPO_TERMS)){
    if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
      res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
    }
  }
} else res_genes_dt <- data.table()



# Write results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)
