# Contains the different functions needed to run the XGBoost model

library(BBmisc)
library(BiocParallel)

add_features_variants <- function(DT, expressed_genes){
  pt <- copy(DT)
  pt <- pt[!grepl('NA', var)] # remove combinations with no variants
  
  # remove common variants
  pt <- pt[MAX_AF < .01]
  pt[is.na(MAX_AF), MAX_AF:=0]
  pt[, var_freq := .N, by = c("var", "hgncSymbol")]
  pt <- pt[var_freq <= 10]
  
  # Add outlier T/F columns
  pt[, overexp := exp_effect > 1]
  pt[, underexp := exp_effect < 1]
  pt[is.na(exp_effect), c('overexp', 'underexp') := F]
  pt[, splice_effect := abs(as.numeric(gsub(",.*", "", splice_effect)))]
  pt[, splice := !is.na(splice_effect)]
  
  pt[, expressed_in_tissue := hgncSymbol %in% expressed_genes$gene_name]
  pt[, OMIM_gene := !is.na(OMIM)]
  pt[is.na(Semantic), Semantic := 0]
  pt[, SELECTION_TYPE := factor(SELECTION_TYPE, levels = c('high_impact', 'medium_impact', 'rare'))]
  setorder(pt, SELECTION_TYPE)
  
  return(pt)
}

# aggregate by sample/gene pair with 2 most severe values
getModelData <- function(PT){
  bplapply(unique(PT$sampleID), pt=PT,
           BPPARAM=MulticoreParam(10, 30, progressbar=TRUE), 
           FUN=function(i, pt){
             # If we have a homozygous stop variant the scores should go in both positions (1 and 2)
             # Hence we duplicate all homozygous variants
             het_string <- c("0|1", "1|0","1/0", "0/1", "0|2", "2|0", "0/2", "0|3", "3|0", "0/3")
             tmpdt <- rbind(pt[sampleID == i], pt[sampleID == i & !GT %in% het_string])
             tmpdt[,.(
               CAUSAL_GENE =unique(CAUSAL_GENE),
               underexp    =unique(underexp),
               overexp     =unique(overexp),
               splice      =unique(splice),
               MAE         =max(MAE, na.rm = T),
               Semantic    =unique(Semantic),
               OMIM        =unique(OMIM_gene),
               expressed   =unique(expressed_in_tissue),
               SpliceAI_1  =sort(c(SpliceAI_MAX, 0.1, 0.1), decreasing=TRUE)[1],
               SpliceAI_2  =sort(c(SpliceAI_MAX, 0.1, 0.1), decreasing=TRUE)[2],
               CADD_PHRED_1=sort(c(CADD_PHRED  , 9, 9),     decreasing=TRUE)[1],
               CADD_PHRED_2=sort(c(CADD_PHRED  , 9, 9),     decreasing=TRUE)[2],
               EVE_1       =sort(c(EVE         , 0.4, 0.4), decreasing=TRUE)[1],
               EVE_2       =sort(c(EVE         , 0.4, 0.4), decreasing=TRUE)[2],
               MAX_AF_1    =sort(c(MAX_AF      , 0.2, 0.2), decreasing=FALSE)[1],
               MAX_AF_2    =sort(c(MAX_AF      , 0.2, 0.2), decreasing=FALSE)[2]),
               by=c("sampleID", "hgncSymbol")]
           }) %>% rbindlist()
}


cv <- function(n, dt, params, BPPARAM=bpparam(), seed=.Random.seed, ...){
    cvs <- bplapply(1:n, cv_iter_xgboost, niter=n, dt=dt, seed=seed,
            params=params, BPPARAM=BPPARAM)
    all <- rbindlist(lapply(cvs, "[[", "predictions"))[order(-model)]
    all[,sampleRank:=rank(-model), by="sampleID"]
    all
    ans <- list(all=all, params=params, seed=seed, per_iteration=cvs)
    ans
}

    
cv_iter_xgboost <- function(iter, niter, dt, params, seed=42, ...){
    
    # 
    # Split dataset into train and test based on samples 
    # not on sample/gene pairs to not leak information
    # 
    set.seed(seed)
    sample_chunks <- chunk(unique(dt[,sampleID]), n.chunks=niter, shuffle=TRUE)
    train_set <- dt[sampleID %in% unlist(sample_chunks[-iter]),]
    test_set  <- dt[sampleID %in% unlist(sample_chunks[ iter]),]
    
    # train XGBoost model
    model <- train_xgboost(data=train_set, params=params, ...)
    
    # predict data
    res <- predict_xgboost(model=model, data=test_set)
    res[,iteration:=iter]
    
    # assemble results
    ans <- list(model=model, predictions=res, iteration=iter, niterations=niter,
            train_samples=unlist(sample_chunks[-iter]), 
            test_samples=unlist(sample_chunks[iter]))
    ans
}


train_xgboost <- function(data, params, verbose=1, nrounds=500, ...){
    xgb_data <- xgb.DMatrix(
            data=as.matrix(data[,-c("CAUSAL_GENE", "sampleID", "hgncSymbol")]),
            label=as.integer(data$CAUSAL_GENE))
    xgb_model <- xgb.train(params=params, data=xgb_data, 
            nrounds=nrounds, verbose=verbose, ...)
    xgb_model
}


predict_xgboost <- function(model, data){
    preds <- data.table(model=predict(model, 
            as.matrix(data[, -c("CAUSAL_GENE", "sampleID", "hgncSymbol")])))
    cbind(data, preds)
}
