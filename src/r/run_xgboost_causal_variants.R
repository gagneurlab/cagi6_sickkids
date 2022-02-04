# 
# 
# Runs the XGBoost model on the CAGI dataset
# using manually collected known causal variants 
# from publications of the Prokisch lab
# 
# 
library(data.table)
library(magrittr)
library(xgboost)
library(caTools)
library(caret)
library(plotROC)
source("./src/r/xgsboost_cv.R")


# Process the cagi data
exp_genes_cagi <- fread('Data/project_data/processed_results/aberrant_expression/exp_genes_cagi6.tsv')
browse_cagi <- fread('Data/project_data/processed_results/combined_rna_results_cagi6_with_variants_browseable.tsv.gz')
pt_cagi <- add_features_variants(browse_cagi, exp_genes_cagi)
pt_cagi[, CAUSAL_GENE := F]
pt_cagi[is.na(MAE), MAE := .5] # needed to overcome NAs
pt_cagi <- pt_cagi[!grepl('^hsa', hgncSymbol)] # remove one miRNA annotated as protein coding

# Add Prokisch causal variants
pt_causal_prokisch <- fread('Data/project_data/processed_results/model/variants_prokisch.tsv')
pt_full <- rbind(pt_causal_prokisch, pt_cagi, fill = T)

# Get training dataset for the model
traindata <- getModelData(pt_full)
saveRDS(traindata, 'Data/project_data/processed_results/model/traindata_prokisch_solved_cagi.Rds')


# Define optimization parameters
params_grid <- expand.grid(
  booster="gbtree", # default
  eta=c(.1), # the lower, the more it prevents overfitting
  max_depth=2, 
  subsample=c(.75, 1), 
  colsample_bytree=1,
  gamma=2*2:3,
  objective="binary:logistic",
  eval_metric=c("logloss", "auc"),
  nthread=15) # 0:3
params_grid

# Run crossvalidation
cv_all_res <- lapply(seq_row(params_grid), traindata=traindata, params_grid=params_grid, 
                     FUN=function(i, traindata, params_grid){
                       xgb_params <- lapply(as.list(params_grid[i,]), function(x) ifelse(is.factor(x), as.character(x), I(x)))
                       cv_res <- cv(n=5, dt=traindata, params=xgb_params, 
                                    BPPARAM=SerialParam(progressbar=TRUE))
                       cv_res$all[,param_id:=i]
                       list(all=cv_res$all, i=i, params=xgb_params)
                     })
all_all <- rbindlist(lapply(cv_all_res, "[[", "all"))

# Evaluate the different models
g <- ggplot(all_all, aes(d=as.numeric(CAUSAL_GENE), m=model)) + 
  geom_roc() + facet_wrap(~param_id) + geom_abline() + theme_bw()
g
melt_pred <- melt(all_all[,.(sampleID, hgncSymbol, CAUSAL_GENE, param_id, model)], measure.vars = 'model')

g <- ggplot(melt_pred, aes(d=as.numeric(CAUSAL_GENE), m=value, col = as.factor(param_id))) + 
  geom_roc() + geom_abline() + theme_bw()
g
aucs <- as.data.table(calc_auc(g))
aucs

# Run the model using the best hyperparameters from the parameters search with CV
optim_params <- cv_all_res[[1]]$params # The first combination gave the highest AUC

xgb_model <- xgb.train(params = optim_params, 
                       data = xgb.DMatrix(as.matrix(traindata[, -c("CAUSAL_GENE", "sampleID", "hgncSymbol")]),
                                          label=as.integer(traindata$CAUSAL_GENE)),
                       nrounds = 500, verbose = 1)
cagi_model_data <- getModelData(pt_cagi)

px <- predict_xgboost(model = xgb_model, data = cagi_model_data)
setorder(px, -model)
saveRDS(px, 'Data/project_data/processed_results/model/model_prokisch_causal_vars.Rds')

#####
## Bootstrapping with replacement to compute the standard error
#####

m <- as.matrix(traindata[, -c("CAUSAL_GENE", "sampleID", "hgncSymbol")])
l <- as.integer(traindata$CAUSAL_GENE) # labels
cagi_model_data <- getModelData(pt_cagi)

boot_dt <- lapply(1:10, function(j){
  set.seed(j)
  is <- sample(1:nrow(m), replace = T)
  xgb_model <- xgb.train(params = optim_params, 
                         data = xgb.DMatrix(m[is,], label = l[is]),
                         nrounds = 500, verbose = 1)
  pb <- predict_xgboost(model = xgb_model, data = cagi_model_data)
  pb <- pb[, .(sampleID, hgncSymbol, model)]
  setorder(pb, -model)
  pb[, boot_round := paste0('b', j)]
  return(pb)
}) %>% rbindlist()

# Compute standard error across bootstraps
boot_dt[, se := sd(model)/sqrt(.N), by = .(sampleID, hgncSymbol)]
boot_cast <- dcast.data.table(boot_dt, ... ~ boot_round, value.var = 'model')
saveRDS(boot_cast, 'Data/project_data/processed_results/model/bootstrap_prokisch_causal_vars.Rds')

