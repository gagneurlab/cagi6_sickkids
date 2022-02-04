library(readxl)
library(data.table)
library(magrittr)


sa <- fread('Data/project_data/raw/CAGI6_Sickkids_Metadata - Sheet1.tsv')
sa[, Person := NULL]
# sa[, `...1` := NULL]
colnames(sa) <- gsub(' ', '_', colnames(sa))
# setnames(sa, 'Sample_ID', 'RNA_ID')
sa$Cohort %>% table()

sa[, RNA_BAM_FILE := paste0('Data/project_data/raw/bamFiles/', RNA_ID, '.mdup.bam')]
sa[, RNA_exists := file.exists(RNA_BAM_FILE)]
sa[RNA_exists == F]

# Add DNA info
sa[Cohort == 'Genome Clinic', DNA_VCF_FILE := paste0('Data/project_data/raw/vcfs_genomeClinic/vcfs_genomeClinic/',
                                                     RNA_ID, '_genome.vcf.gz')]
sa[Cohort == 'Complex Care', DNA_VCF_FILE := paste0('Data/project_data/raw/vcfs_complex_care/vcfs_complex_care/',
                                                     RNA_ID, '.vcf.gz')]
sa[, DNA_exists := file.exists(DNA_VCF_FILE)]
sa[DNA_exists == F, DNA_VCF_FILE := NA]

library(SeqArray)
sa$DNA_ID <- sapply(sa$DNA_VCF_FILE, function(vcf){
  if(!is.na(vcf)){
    return(seqVCF_SampID(vcf))
  } else return(NA)
})


sa[, DROP_GROUP := 'all']
sa[DNA_exists == T, DROP_GROUP := 'all,mae']

sa[, STRAND := 'reverse']
sa[, COUNT_MODE := "IntersectionStrict"]
sa[, PAIRED_END := TRUE]
sa[, COUNT_OVERLAPS := TRUE]

sa[, GENE_COUNTS_FILE := NA]
sa[, GENE_ANNOTATION := NA]

fwrite(sa[RNA_exists == T], 'Data/project_data/raw/sample_annotation.tsv', sep = '\t', quote = F)

# Create another sample annotation with the 50 blood samples with the highest coverage
sa_gtex <- fread('Data/gtex_analysis/sample_annotation.tsv')
sa_gtex <- sa_gtex[TISSUE_CLEAN == 'Whole_Blood', .(RNA_ID, sex, INDIVIDUAL_ID, DNA_ID, 
                                                    TISSUE_CLEAN, DROP_GROUP, COUNT_MODE, PAIRED_END, 
                                                    COUNT_OVERLAPS, STRAND, RNA_BAM_FILE, DNA_VCF_FILE,
                                                    GENE_COUNTS_FILE, GENE_ANNOTATION)]
cov_dt <- fread('Data/gtex_analysis/processed_data/aberrant_expression/v29/outrider/Whole_Blood/bam_coverage.tsv')
setorder(cov_dt, - record_count)
sa_gtex <- sa_gtex[RNA_ID %in% cov_dt[3:52,sampleID]]

sa_comb <- rbind(sa, sa_gtex, fill = T)
sa_comb[, DROP_GROUP := 'gtex_comb']
sa_comb[, STRAND := 'no']
fwrite(sa_comb, 'Data/project_data/gtex_comb/sample_annotation.tsv', sep = '\t', quote = F)
