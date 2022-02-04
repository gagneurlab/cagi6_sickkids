library(ontologyIndex)
library(ontologySimilarity)
library(data.table)
library(magrittr)
library(PCAN)

hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")

# Read HPO file from website and simplify it
hpo_gene <- fread("https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt", skip = 1)
colnames(hpo_gene) <- c("HPO_ID", "HPO_name", "entrezID", "geneID", "Additional_Info",  "source", "disease-ID")
hpo_gene <- hpo_gene[, .(geneID, HPO_ID)]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]
hpo_gene <- hpo_gene[HPO_ID %in% unique(hpo$id)]
hpByGene <- unstack(hpo_gene, HPO_ID~geneID)


dir <- 'Data/processed/hpo/hpo_sim/'
sa <- fread('Data/processed/sample_annotation.tsv')

lapply(1:nrow(sa), function(i){
  x <- sa[i, ]
  rna_id <- x$RNA_ID
  file <- paste0('Data/processed/hpo/samples_gene/', rna_id, '.tsv')
  
  if(!file.exists(file)){
    # Extract HPO terms
    patient_hpos <- strsplit(x$HPO_TERMS, split = ',') %>% unlist() %>% sort()
    patient_hpos <- paste0('HP:', gsub(".*:", '', patient_hpos)) # needed bc of some special characters
    hpo_files <- paste0(dir, patient_hpos, '.Rds')
    
    # Read hpGene and group by gene
    hpGeneResnik <- do.call(cbind, lapply(hpo_files, readRDS))
    hpMatByGene <- lapply(hpByGene, function(x){
                          hpGeneResnik[x, , drop=FALSE]
                        })
    
    # Compute the semantic similarity scores
    semantic_sim <- lapply(hpMatByGene,
                         hpSetCompSummary, method="bma", direction="symSim") %>% unlist()
    
    ss_dt <- as.data.table(semantic_sim, keep.rownames = 'geneID')
    ss_dt[, SAMPLE_ID := x$RNA_ID]
    ss_dt <- ss_dt[!is.na(semantic_sim)]
    setorder(ss_dt, -semantic_sim)
    fwrite(ss_dt, file, sep = '\t', quote = F)
  }
})


