# Read and export the expressed genes in the CAGI6 datasets

ods <- readRDS('Data/project_data/processed_results/aberrant_expression/cagi6/outrider/all/ods.Rds')
cagi_genes_dt <- fread('Data/project_data/processed_data/preprocess/cagi6/gene_name_mapping_cagi6.tsv')

fwrite(cagi_genes_dt[gene_id %in% rownames(ods), .(gene_name)], 
       'Data/project_data/processed_results/aberrant_expression/exp_genes_cagi6.tsv', sep = '\t', quote = F)

