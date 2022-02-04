library(ontologyIndex)
library(ontologySimilarity)
library(PCAN)
library(data.table)
library(magrittr)
threads <- 20

hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")

# Get ancestors and descendants terms
hp_ancestors <- hpo$ancestors
ic <- descendants_IC(hpo)


sa <- fread('Data/processed/sample_annotation.tsv')
patient_hpos <- lapply(sa$HPO_TERMS, strsplit, split = ',') %>% unlist() %>% unique() %>% sort()
patient_hpos <- gsub(".*:", '', patient_hpos) %>% unique() # needed bc of some special characters
patient_hpos <- paste0('HP:', patient_hpos)


dir <- 'Data/processed/hpo/hpo_sim/'
for(hpOfInterest in patient_hpos){
  file <- paste0(dir, hpOfInterest, '.Rds')
  if(!file.exists(file)){
    hpGeneResnik <- compareHPSets(
      hpSet1=names(ic), hpSet2=hpOfInterest,
      IC=ic,
      ancestors=hp_ancestors,
      method="Resnik",
      BPPARAM= MulticoreParam(threads)
    )
    saveRDS(hpGeneResnik, file)
  }
}

