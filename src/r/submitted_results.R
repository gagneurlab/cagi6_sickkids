library(data.table)
library(magrittr)

vars_cagi <- fread('Data/project_data/processed_results/combined_rna_results_cagi6_with_variants_all.tsv.gz')
model_res <- readRDS('Data/project_data/processed_results/model/model_prokisch_causal_vars.Rds')[, .(sampleID, hgncSymbol, model)]
boot_cast <- readRDS('Data/project_data/processed_results/model/bootstrap_prokisch_causal_vars.Rds')[, .(sampleID, hgncSymbol, se)]

all_res <- merge(vars_cagi, merge(model_res, boot_cast, 
        by = c('sampleID', 'hgncSymbol'), sort=FALSE),
        by = c('sampleID', 'hgncSymbol'), sort=FALSE)


all_res[, var := paste0(seqnames, ":", start, ":", REF, ":", ALT)]
submission_res <- all_res[, .(
        sampleID, 
        "HGNC Symbol"=hgncSymbol, 
        P=round(model, 6),
        SD=round(se, 8),
        Expression=ifelse(is.na(exp_padjust), "no", "yes"), 
        Splicing=ifelse(spl_padjust == "", "no", "yes"),
        ASE=ifelse(is.na(mae_padjust), "no", "yes"),
        Other="no",
        "RNA impact"=NA_character_,
        "Refseq transcript version"=NA_character_,
        "DNA variant"=var,
        "Additional comments"=NA_character_,
        foldChange,
        exp_padjust,
        aberrantSpliceType,
        genomicLocation,
        altRatio)]


# join variants into one line and take unique only considering MAE events
submission_res[, `DNA variant` := paste(paste(`DNA variant`, P, SD, sep = ':'), collapse=","), by=c("sampleID", "HGNC Symbol")]
dt_cagi <- submission_res[order(sampleID, -P, -ASE)][!duplicated(paste(sampleID, `HGNC Symbol`))]

# create nice RNA_impact column
dt_cagi[foldChange < 1, exp_impact := paste0('decreased expression, fold change: ', foldChange, ', FDR: ', round(exp_padjust, 4))]
dt_cagi[foldChange > 1, exp_impact := paste0('increased expression, fold change: ', foldChange, ', FDR: ', round(exp_padjust, 4))]

dt_cagi[, aberrantSpliceType := gsub('trunc, elong', 'exonTruncation&Elongation', aberrantSpliceType)]
dt_cagi[, aberrantSpliceType := gsub('NA', 'inconclusive', aberrantSpliceType)]
dt_cagi[Splicing == 'yes' & aberrantSpliceType == "", aberrantSpliceType:='inconclusive']
dt_cagi[, aberrantSpliceType := gsub('beyondGene', 'splicingBeyondGene', aberrantSpliceType)]
dt_cagi[, aberrantSpliceType := gsub('multigenic', 'multigenicSplicing', aberrantSpliceType)]
dt_cagi[, genomicLocation := gsub(':\\+|:-', '', genomicLocation)] # remove strand from coordinates
dt_cagi[Splicing == 'yes', spl_impact := paste0('splicing change: ', aberrantSpliceType, ', genomic location: ', genomicLocation)]

dt_cagi[!is.na(altRatio), mae_impact := paste0('Alt allele ratio: ', altRatio)]

dt_cagi[Expression == 'yes' | Splicing == 'yes' | ASE == 'yes', `RNA impact` := paste(exp_impact, spl_impact, mae_impact, sep = ' | ')]
dt_cagi[, `RNA impact` := gsub('NA [|] ', '', `RNA impact`)]
dt_cagi[, `RNA impact` := gsub(' [|] NA', '', `RNA impact`)]


# remove unwanted columns
dt_cagi[, c('exp_padjust', 'foldChange', 'aberrantSpliceType', 'genomicLocation', 
        'altRatio', "exp_impact", "spl_impact", "mae_impact") := NULL]


# Create and save first submission
#  --> Subset for 100 highest scores
setorder(dt_cagi, sampleID, -P)
dt_100 <- dt_cagi[, .SD[1:100,], by = sampleID]
fwrite(dt_100, 'Data/project_data/submission/droppers_model_1.tsv', sep = ';', quote = F)


# create second submission
manual_res <- fread('Data/project_data/manually_curated_cases.tsv')
dt_manual <- merge(dt_cagi, manual_res, by.x=c("sampleID", "HGNC Symbol"), by.y=c("sampleID", "hgncSymbol"), all=TRUE)
setnames(dt_manual, "P.x", "P")
setnames(dt_manual, "SD.x", "SD")
setnames(dt_manual, "Additional comments.x", "Additional comments")
dt_manual[!is.na(P.y), c("P", "SD", "Additional comments") := list(P.y, SD.y, `Additional comments.y`)]

dt_manual[,c("P.y", "SD.y", "Additional comments.y") := NULL]

setorder(dt_manual, sampleID, -P)
dt_manual <- dt_manual[, .SD[1:100,], by = sampleID]
fwrite(dt_manual, 'Data/project_data/submission/droppers_model_2.tsv', sep = ';', quote = F)
