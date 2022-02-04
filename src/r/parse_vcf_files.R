library(VariantAnnotation)
library(ensemblVEP)
library(data.table)
library(GenomicFiles)

# save.image(file = "./snakemake_debug.RData")

threads <- snakemake@threads
file <- snakemake@input$VCF
fileOut <- snakemake@output$TSV

register(MulticoreParam(threads, progressbar=TRUE))

cast2Numeric <- function(dt, cols){
    for(i in cols){
        dt[,c(i):=as.numeric(get(i))]
    }
    dt    
}

getSpliceAIMax <- function(dt, cols){
    saim <- dt[,get(cols[1])]
    for(i in 2:length(cols)){
        saim <- pmax(saim, dt[,get(cols[i])], na.rm=TRUE)
    }
    dt[,SpliceAI_MAX:=saim]
    dt
}

spliceAI_scores <- c(
    "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", 
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
numCols <- c("CADD_PHRED", "CADD_RAW", "gnomAD_AF", "MAX_AF", spliceAI_scores)


curSampleID <- gsub("(_genome)?.vep.vcf.gz", "", basename(file))
yieldSize <- 100000
vcfFile <- VcfFile(file, yieldSize=yieldSize)

# open(vcfFile)
# x <- yield(vcfFile)

#' 
#' Reading in the data
#' 
yield <- function(x, genome="hg19", ...){
    readVcf(x, genome=genome, ...)
}

#' 
#' Parsing a block of the input
#' 
map <- function(x, ...) {
    # transform CSQ into a VEP table of type data.table
    vep <- as.data.table(parseCSQToGRanges(x, VCFRowID=rownames(x)))
    
    # convert columns to correct data type 
    vep <- cast2Numeric(vep, numCols)
    vep <- getSpliceAIMax(vep, grep("_DS_", spliceAI_scores, value=TRUE))
    vep[,MAX_AF:=pmax(gnomAD_AF, MAX_AF, na.rm=TRUE)]
    
    # filter for 
    #   * rare variants
    #   * HGNC symbol (we have more on ENSEMBL IDs)
    #   * non intergenic
    res <- vep[
        CLIN_SIG == "pathogenic" | (
        (is.na(MAX_AF) | MAX_AF < 0.01) & 
        !is.na(SYMBOL) &
        Consequence != "intergenic")]
    res[,SELECTION_TYPE:="rare"]
    
    res[
        IMPACT %in% c("HIGH", "MODERATE") | 
        grepl("splice", Consequence) | 
        (!is.na(five_prime_UTR_variant_consequence) & five_prime_UTR_variant_consequence != "") |
        CADD_PHRED >= 10 | 
        SpliceAI_MAX >= 0.2 | 
        grepl("(likely_)?pathogenic(&|$)", CLIN_SIG),
    SELECTION_TYPE:="medium_impact"]
    
    res[
        grepl("(likely_)?pathogenic(&|$)", CLIN_SIG) | 
        IMPACT %in% c("HIGH") |
        (!is.na(five_prime_UTR_variant_consequence) & five_prime_UTR_variant_consequence != "") |
        CADD_PHRED >= 20 | 
        SpliceAI_MAX >= 0.5,
    SELECTION_TYPE:="high_impact"]
    
    # 
    # add variant level information
    # 
    x2res <- x[res[,VCFRowID],]
    res[,REF:=as.character(ref(x2res))]
    res[,ALT:=sapply(alt(x2res), function(x) paste(as.character(x), collapse=","))]
    res[,QUAL:=qual(x2res)]
    res[,FILTER:=filt(x2res)]
    res[,GT:=geno(x2res)$GT]
    res[,DP:=geno(x2res)$DP]
    res[,MQ:=info(x2res)$MQ]
    res[,VCFRowID:=names(x2res)]
    
    res
}

reduce <- rbind

# run the extraction on the VCF file
allres <- reduceByYield(vcfFile, YIELD=yield, MAP=map, REDUCE=reduce, parallel=TRUE)

# save the extraction to disc
fwrite(allres, file=fileOut)

