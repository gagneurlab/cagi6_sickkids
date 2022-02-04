##--------------------------------------------
## required packages for the aberrantSpliceTypes
message("Load aberrant splice type packages")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(rtracklayer) #to import that blacklist file
  library(GenomicRanges)
})


addUTRLabels <- function(junctions_dt, txdb){
  psi_positions <- which(junctions_dt$type != "theta")
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions])
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  ### UTR labels based on txdb file
  ### add 5' 3' UTR labels
  print("find UTR overlap")
  threes <- unique(from(findOverlaps(junctions_gr, threeUTRsByTranscript(txdb, use.names = TRUE))))
  fives <- unique(from(findOverlaps(junctions_gr, fiveUTRsByTranscript(txdb, use.names = TRUE))))
  junctions_dt[psi_positions, UTR := "no"]
  #print("UTRSSSSSSSSSS:")
  #print(threes)
  #print(fives)
  junctions_dt[psi_positions[threes], UTR := "3"]
  junctions_dt[psi_positions[fives], UTR := "5"]
  colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
  print("UTR labels done")
  #print(junctions_dt)
  return(junctions_dt)
}

### basic annotations (start, end, none, both) for full fds
createFDSAnnotations <- function(fds, txdb, BPPARAM=bpparam()){
  print("loading introns")
  #seqlevelsStyle(fds) <- seqlevelsStyle(txdb)[1]
  introns <- unique(unlist(intronsByTranscript(txdb)))
  # reduce the introns to only the actually expressed introns
  fds_known <- fds[unique(to(findOverlaps(introns, rowRanges(fds, type = "j"), type = "equal"))),]
  anno_introns <- as.data.table(rowRanges(fds_known, type="j"))[,.(seqnames, start, end, strand)]
  
  #calculate extra columns with mean/median intron expression count
  #add the new columns
  print("adding median count to introns")
  sampleCounts <- K(fds_known, type = "j")
  anno_introns[, meanCount := rowMeans(sampleCounts)]
  anno_introns[, medianCount := rowMedians(as.matrix(sampleCounts))]
  setorderv(anno_introns, "medianCount", order=-1) # order by medianCount (highest first)
  
  anno_introns_ranges <- makeGRangesFromDataFrame(anno_introns, keep.extra.columns = TRUE)
  
  ### get all fds junctions
  fds_junctions <- rowRanges(fds, type = "j")
  
  ### Do the annotation just for the most used intron (highest median expression)
  print("start calculating annotations")    
  overlaps <- findOverlaps(fds_junctions, anno_introns_ranges, select="first")
  annotations <- bplapply(1:length(fds_junctions), function(i){
    # only select first intron as already ordered by medianCount beforehand
    overlap <- overlaps[i]
    if(is.na(overlap)) return("none") #no overlap with any intron
    
    hit_equal <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap], type="equal"))
    if(length(hit_equal) > 0) return("both")
    
    hit_start <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap], type="start"))
    if(length(hit_start) > 0) return("start")
    hit_end   <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap], type="end"))
    if(length(hit_end) > 0) return("end")
    
    return("none") # overlaps but no start/end match
  }, BPPARAM=BPPARAM) # %>% unlist()
  annotations <- unlist(annotations)
  
  #table(annotations)
  rowRanges(fds)$annotatedJunction <- annotations
  print("annotations done")
  return(fds)
}