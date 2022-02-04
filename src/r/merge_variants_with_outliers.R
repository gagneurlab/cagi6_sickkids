library(data.table)
library(BiocParallel)
library(GenomeInfoDb)
library(magrittr)

if("snakemake" %in% ls()){
    threads    <- snakemake@threads
    file_exp   <- snakemake@input$EXP
    file_var   <- snakemake@input$VAR
    file_eve   <- snakemake@input$EVE
    out_all    <- snakemake@output$ALL
    out_browse <- snakemake@output$BROWSE
} else {
    threads  <- 40
    file_exp <- "Data/project_data/processed_results/combined_rna_results_cagi6.tsv"
    file_var <- list.files(path="Data/project_data/vep_annotation", pattern="rare_vars.tsv.gz", full.names = TRUE)
    file_eve <- "Data/project_data/EVE_dt_all.Rds" # contains EVE scores for CAGI6 variants
    out_all  <- "Data/project_data/processed_results/combined_rna_results_cagi6_with_variants_all.tsv.gz"
    out_browse <- "Data/project_data/processed_results/combined_rna_results_cagi6_with_variants_browseable.tsv.gz"
}


#' read sample annotation
dir <- 'Data/project_data'
sa_file <- file.path(dir, 'sample_annotation.tsv')
sa <- fread(sa_file)

names(file_var) <- gsub("(_genome)?.rare_vars.tsv.gz", "", basename(file_var))
file_var <- file_var[intersect(union(sa$DNA_ID, sa$RNA_ID), names(file_var))] # needs to be sorted
names(file_var) <- sa$RNA_ID
file_var[1:5]
sa[, .(DNA_ID, RNA_ID)][1:5]

register(MulticoreParam(threads, length(file_var), progressbar=TRUE))

length(file_var)
file_eve

extract_variants_per_sample <- function(i, files=file_var, symbols=NULL, eve){
    # read variants
    dt_var <- fread(files[i])
    dt_var[, c('Allele', 'EXON', 'INTRON', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
               'ALLELE_NUM', 'DISTANCE', 'STRAND', 'FLAGS', 'MINIMISED', 'SYMBOL_SOURCE', 'HGNC_ID', 'MANE', 'TSL', 'APPRIS', 'TREMBL', 'REFSEQ_MATCH', 'SOURCE',
               'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'miRNA', 'HGVS_OFFSET', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 
               'SOMATIC', 'PHENO', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'Condel', 'BLOSUM62', 'miRNA.1',
               'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 
               'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL',
               'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs') := NULL]
    dt_var <- dt_var[BIOTYPE == 'protein_coding']
    if(!is.null(symbols)) dt_var <- dt_var[SYMBOL %in% symbols]
    dt_var[, sampleID := i]
    
    seqlevelsStyle(dt_eve$CHR) <- seqlevelsStyle(dt_var$seqnames)
    # Add EVE score
    dt_var <- merge(dt_var, dt_eve[,.(seqnames=CHR, start=START, EVE)], 
            all.x=TRUE, by=c("seqnames", "start"))
    
    # fix prioritization
    dt_var[SELECTION_TYPE == "high_impact", SELECTION_TYPE:="medium_impact"]
    dt_var[EVE > 0.5, SELECTION_TYPE:="medium_impact"]
    dt_var[
        grepl("(likely_)?pathogenic(&|$)", CLIN_SIG) | 
        IMPACT %in% c("HIGH") |
        (!is.na(five_prime_UTR_variant_consequence) & five_prime_UTR_variant_consequence != "") |
        CADD_PHRED >= 20 | 
        SpliceAI_MAX >= 0.5 | 
        EVE > 0.64,
        SELECTION_TYPE:="high_impact"]
    
    # set correct ordering
    dt_var[,NumRareVarsPerSymbol:=length(unique(VCFRowID)), by="SYMBOL"]
    dt_var[,NumTranscriptsPerVar:=.N, by="VCFRowID,SYMBOL"]
    dt_var[,IMPACT:=factor(IMPACT, levels=c("HIGH", "MODERATE", "MODIFIER", "LOW"))]
    dt_var[,SELECTION_TYPE:=factor(SELECTION_TYPE, levels=c("high_impact", "medium_impact", "rare"))]
    dt_var[,CANONICAL:=factor(CANONICAL, levels=c("YES", ""))]
    
    # boil it down to 2 variants per gene/sample pair
    dt2merge <- dt_var[order(CANONICAL, SELECTION_TYPE, IMPACT, EVE)]
    dt2merge[,VarTransOrder:=1:.N,by="VCFRowID,SYMBOL"]
    dt2merge <- dt2merge[VarTransOrder == 1]
    dt2merge[,VarSymbolOrder:=1:.N,by="SYMBOL"]
    dt2merge <- dt2merge[VarSymbolOrder <= 2]
    
    # return it
    dt2merge
}


#' 
#' read RNA outlier results
#' 
dt_out <- fread(file.path(dir, 'processed_results', 'combined_rna_exp_spl_results.tsv'))

#' 
#' read EVE database
#'
dt_eve <- readRDS(file_eve)
dt_eve[, SAMPLE := NULL]
dt_eve <- unique(dt_eve)

#' extract all variants to be merged with outlier results
res_vars_dt <- bplapply(names(file_var), extract_variants_per_sample, 
                     files=file_var, symbols=NULL, 
                     eve=dt_eve) %>% rbindlist(fill = TRUE)
res_vars_dt[!grepl('chr', seqnames), seqnames := paste0('chr', seqnames)]
print('variants merged')


#' Add MAE
mae_file <- file.path(dir, 'processed_results', 'mae', 'mae', 'MAE_results_cagi6_simplified.tsv')
res_mae <- fread(mae_file)
seqlevelsStyle(res_mae$contig) <- seqlevelsStyle(res_vars_dt$seqnames)
res_vars_dt <- merge(res_vars_dt, res_mae[, .(sampleID, contig, position, refCount, altCount, altRatio, mae_padjust)], 
                     by.x = c('sampleID', 'seqnames', 'start'), by.y = c('sampleID', 'contig', 'position'), all.x = T, sort = F)

#' 
#' Add RNA expression and splicing outliers with variant data
#' 
dt_all <- merge(dt_out, res_vars_dt, 
                by.x=c("sampleID", "hgncSymbol"), by.y=c("sampleID", "SYMBOL"), all=TRUE)
dt_all <- dt_all[BIOTYPE != '']
dt_all[, outlier:=""]
dt_all[!is.na(foldChange) & deltaPsi == "" & mae_padjust == "", outlier:="e"]
dt_all[ is.na(foldChange) & deltaPsi != "" & mae_padjust == "", outlier:="s"]
dt_all[ is.na(foldChange) & deltaPsi == "" & mae_padjust != "", outlier:="m"]
dt_all[!is.na(foldChange) & deltaPsi != "" & mae_padjust == "", outlier:="es"]
dt_all[!is.na(foldChange) & deltaPsi == "" & mae_padjust != "", outlier:="em"]
dt_all[ is.na(foldChange) & deltaPsi != "" & mae_padjust != "", outlier:="sm"]
dt_all[!is.na(foldChange) & deltaPsi != "" & mae_padjust != "", outlier:="esm"]
dt_all[, outlier := factor(outlier, levels = c("esm", "es", "em", "sm", "e", "s", "m"))]
print('outliers added')

# Add HPO and OMIM
source('.drop/helpers/add_HPO_cols.R')
dt_all <- add_HPO_cols(dt_all)

ot <- fread('Data/project_data/external/omim_genes.tsv')
dt_all <- merge(dt_all, ot[, .(SYMBOL, OMIM, OMIM_Inheritance)], 
                by.x = 'hgncSymbol', by.y = 'SYMBOL', all.x = T, sort = F)

# Add semantic similarity
st <- lapply(list.files(paste0('Data/project_data/hpo/samples_gene/'), full.names = T), fread) %>% rbindlist()
setnames(st, old = c('geneID', 'SAMPLE_ID'), new = c('hgncSymbol', 'sampleID'))
dt_all <- merge(dt_all, st, by = c('hgncSymbol', 'sampleID'), all.x = T, sort = F)
print('semantic similarity added')

#' 
#' subset to a browseable set
#' 
browseable_dt <- dt_all[, .(sampleID, hgncSymbol, var=paste0(seqnames, ":", start, ":", REF, ">", ALT), 
    GT, CLIN_SIG, CADD_PHRED, REVEL, SpliceAI_MAX, EVE, SELECTION_TYPE, outlier, exp_effect=foldChange,
    splice_effect=deltaPsi, MAE=altRatio, Consequence, IMPACT, MAX_AF,
    PUBMED=gsub("[&].*", " &X", PUBMED),
    NumVars=NumRareVarsPerSymbol, Semantic=round(semantic_sim, 2),
    OMIM, OMIM_Inheritance, HPO_label_overlap)][order(outlier)]

browseable_dt


#' 
#' save table to tsv
#' 
fwrite(dt_all, file=out_all)
fwrite(browseable_dt, file=out_browse)

