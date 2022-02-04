library(data.table)
library(magrittr)
library(tidyr)

sa <- fread('Data/project_data/sample_annotation.tsv')

## Read aberrant expression results (CAGI alone and with GTEx) and combine them
res_cagi <- fread('Data/project_data/processed_results/aberrant_expression/cagi6/outrider/all/OUTRIDER_results.tsv')
res_gtex <- fread('Data/project_data/gtex_comb/processed_results/aberrant_expression/cagi6/outrider/gtex_comb/OUTRIDER_results.tsv')
res_gtex <- res_gtex[!grepl('^SRR', sampleID)] # remove GTEx samples
res_cagi[, aux := paste(sampleID, hgncSymbol, sep = '--')]
res_gtex[, aux := paste(sampleID, hgncSymbol, sep = '--')]
res_gtex <- res_gtex[! aux %in% res_cagi$aux] # do not include overlapping results twice

res_exp <- rbind(res_cagi, res_gtex, fill = T)
res_exp[, c('geneID', 'zScore', 'l2fc', 'aberrant', 'AberrantBySample', 'AberrantByGene', 'padj_rank') := NULL]
res_exp[, c('OMIM', 'OMIM_Inheritance', 'HPO_label_overlap', 'HPO_id_overlap', 'Ngenes_per_pheno', 'Nphenos_gene') := NULL]
setnames(res_exp, old = c('pValue', 'padjust', 'rawcounts', 'normcounts'),
         new = c('exp_pValue', 'exp_padjust', 'exp_rawcounts', 'exp_normcounts'))


## Read aberrant splicing results
res_spl <- fread('Data/project_data/processed_results/aberrant_splicing/results/cagi6/fraser/all/results.tsv')
res_spl <- res_spl[!is.na(hgncSymbol)]
res_spl[, c('zScore', 'phi', 'rho_theta', 'rho_psi5', 'rho_psi3', 'rho_filtered', 'blacklist', 
            'numSamplesPerJunc', 'numSamplesPerGene', 'numEventsPerGene', 'addHgncSymbols', 'seqnames') := NULL]
res_spl[, c('OMIM', 'OMIM_Inheritance', 'HPO_label_overlap', 'HPO_id_overlap', 'Ngenes_per_pheno', 'Nphenos_gene') := NULL]
setnames(res_spl, old = c('type', 'pValue', 'padjust', 'counts', 'totalCounts'),
         new = c('psiType', 'spl_pValue', 'spl_padjust', 'spl_counts', 'spl_totalCounts'))

## Combine results and save them
res_spl_out <- merge(res_exp, res_spl, by = c('sampleID', 'hgncSymbol'), all = T, sort = F)
fwrite(res_spl_out, 'Data/project_data/processed_results/combined_rna_exp_spl_results.tsv', 
       sep = '\t', quote = F)

# Read MAE results
# res_mae <- fread('Data/project_data/processed_results/mae/mae/MAE_results_cagi6.tsv')
# res_mae <- res_mae[contig != 'Y' & rare == T & gene_type == 'protein_coding']
# res_mae[, c('variantID', 'log2FC', 'AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'other_names', 'N_var', 'MAE') := NULL]
# res_mae[, c('OMIM', 'OMIM_Inheritance', 'gene_type') := NULL]
# setnames(res_mae, old = c('padj', 'gene_name'),
#          new = c('mae_padjust', 'hgncSymbol'))
# res_mae <- separate(res_mae, 'ID', into = c('DNA_ID', 'sampleID'), sep = '--')
# res_mae[, DNA_ID := NULL]
# fwrite(res_mae, 'Data/project_data/processed_results/mae/mae/MAE_results_cagi6_simplified.tsv')
res_mae <- fread('Data/project_data/processed_results/mae/mae/MAE_results_cagi6_simplified.tsv')
res_mae[, variant := paste0(contig, ":", position, refAllele, ">", altAllele)]
res_mae[, c("contig", "position", "refAllele", "altAllele") := NULL]
res_mae[, Ngene := .N, by = .(hgncSymbol, sampleID)]
collapse_res <- res_mae[Ngene > 1, lapply(.SD, function(x){
  paste(x, collapse=",")
  }), by = .(hgncSymbol, sampleID)]
res_mae <- rbind(res_mae[Ngene == 1], collapse_res)
res_mae[, Ngene := NULL]

res_all <- merge(res_spl_out, res_mae, by = c('sampleID', 'hgncSymbol'), all = T, sort = F)

# Add HPO
source('.drop/helpers/add_HPO_cols.R')
res_all <- add_HPO_cols(res_all)

# Add OMIM genes
ot <- fread('Data/project_data/external/omim_genes.tsv')
res_all <- merge(res_all, ot[, .(SYMBOL, OMIM, OMIM_Inheritance)], 
                 by.x = 'hgncSymbol', by.y = 'SYMBOL', all.x = T, sort = F)

# Add DNA_ID
res_all <- merge(res_all, sa[, .(RNA_ID, DNA_ID)], by.x = 'sampleID', by.y = 'RNA_ID', sort = F)

# Add semantic similarity
st <- lapply(list.files('Data/project_data/hpo/samples_gene', full.names = T), fread) %>% rbindlist()
setnames(st, old = c('geneID', 'SAMPLE_ID'), new = c('hgncSymbol', 'sampleID'))
res_all <- merge(res_all, st, by = c('hgncSymbol', 'sampleID'), all.x = T, sort = F)

# Save results
res_all[!is.na(exp_padjust) & !is.na(spl_padjust) & !is.na(mae_padjust)]
fwrite(res_all, 'Data/project_data/processed_results/combined_rna_results_cagi6.tsv', 
       sep = '\t', quote = F)

