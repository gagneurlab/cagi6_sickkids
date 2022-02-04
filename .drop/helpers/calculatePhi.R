calculatePhi <- function(fds){
  # retrieve junction and splice site annotation with count information
  junction_dt <- data.table(
    as.data.table(rowRanges(fds, type="j"))[,.(seqnames, start, end, 
                                               strand, startID, endID)], 
    as.matrix(K(fds, type="psi5")))
  ss_dt <- data.table(
    as.data.table(rowRanges(fds, type="ss"))[,.(seqnames, start, end, 
                                                strand, spliceSiteID)], 
    as.matrix(N(fds, type="theta")))
  
  # melt to have one row per sample - junction combination
  junction_dt <- melt(junction_dt, variable.name="sampleID", value.name="k", 
                      id.vars=c("seqnames", "start", "end", "strand", 
                                "startID", "endID"))
  ss_dt <- melt(ss_dt, variable.name="sampleID", value.name="n", 
                id.vars=c("seqnames", "start", "end", "strand", 
                          "spliceSiteID"))
  
  # merge junction information with splice site annotation (theta)
  junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                       all.x=TRUE, by.x=c("startID", "sampleID"), 
                       by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
  setnames(junction_dt, "n", "n_donor")
  junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                       all.x=TRUE, by.x=c("endID", "sampleID"), 
                       by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
  setnames(junction_dt, "n", "n_acceptor")
  rm(ss_dt)
  gc()
  
  # deal with missing endIDs
  junction_dt_nacceptor <- data.table(
    as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
    as.matrix(N(fds, type="psi3")))
  junction_dt_nacceptor <- melt(junction_dt_nacceptor, 
                                variable.name="sampleID", value.name="n_psi3", 
                                id.vars=c("startID", "endID"))
  junction_dt[, n_psi3:=junction_dt_nacceptor[,n_psi3]]
  rm(junction_dt_nacceptor)
  gc()
  junction_dt_ndonor <- data.table(
    as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
    as.matrix(N(fds, type="psi5")))
  junction_dt_ndonor <- melt(junction_dt_ndonor, 
                             variable.name="sampleID", value.name="n_psi5", 
                             id.vars=c("startID", "endID"))
  junction_dt[, n_psi5:=junction_dt_ndonor[,n_psi5]]
  rm(junction_dt_ndonor)
  gc()
  
  # replace n (with non-split counts) by n_psi3/n_psi5 if NA
  junction_dt[is.na(n_acceptor), n_acceptor:=n_psi3]
  junction_dt[is.na(n_donor), n_donor:=n_psi5]
  
  # calculate phi
  junction_dt[, phi:= k / ((n_donor + n_acceptor)/2)]
  junction_dt[is.nan(phi), phi:=0] # n_donor = n_acceptor = 0
  
  # convert to matrix to store it as assay in the fds
  phi <- matrix(junction_dt[,phi], nrow=nrow(fds), ncol=ncol(fds), byrow=FALSE)
  rownames(phi) <- rownames(fds)
  colnames(phi) <- colnames(fds)
  assay(fds, "phi", type="psi5", withDimnames=FALSE) <- phi
  return(fds)
}
