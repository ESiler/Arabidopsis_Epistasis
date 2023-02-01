#1. Load processed data, formulas, libraries, and other objects:
source('scripts/load_required_packages.R')
load("rdata/02_workspace.RData")


#1. Get_epi_stats. Gets epistasis values for plant set. 
# Returns a vector of relevant values

get_epi_stats <- function(plantset, formula, df_pred, df){
  dfsubset <- subset(df, Set == plantset)
  model <- lm(formula, data=dfsubset, na.action = na.exclude)
  #Extract stats
  e_est <- as.numeric(coef(model)['DM'])
  lowerCI <- confint(model)[4,1]
  upperCI <- confint(model)[4,2]
  rsquared <- summary(model)$adj.r.squared
  pval_e <- as.numeric(summary(model)$coefficients[,4]['DM'])
  
  pred_fit = predict(model, df_pred, interval="confidence") 
  pred_rel_fit = pred_fit/pred_fit[1,1]
  
  pvals <- summary(model)$coefficients[,4][1:4]
  #Return stats
  result <- c(plantset, e_est, lowerCI, upperCI, rsquared, pval_e, pred_rel_fit, pvals)
  
  return(result)
}

tmp <- get_epi_stats(5, f.tsc.f, df_pred_dummy.f, data) #Need to make dummy frame flexible based on model
tmp
length(tmp)

for (i in sets_with_flats){
  v <- get_epi_stats(i, formula=f.tsc, df=data, df_pred=df_pred_dummy)
  resultlist[[i]] <- as.numeric(v)
}

#2. Gets all epistasis values for trait. 
# Returns a dataframe w epistasis values that you can then plot with plot_epi_forest
get_epistasis_for_formula <- function(plantsets, formula, df_pred, df) {
  
  resultlist <- list()
  
  for (i in plantsets){
    v <- get_epi_stats(i, formula=formula, df=df, df_pred=df_pred)
    resultlist[[i]] <- as.numeric(v)
  }
  these_rownames =  c('Set', 'e_est', 'lowerCI', 'upperCI', 'rsquared', 'pval_e', 
                      'WT_w', 'MA_w', 'MB_w', 'DM_w', 
                      'WT_w_lci', 'MA_w_lci', 'MB_w_lci', 'DM_w_lci', 
                      'WT_w_uci', 'MA_w_uci', 'MB_w_uci', 'DM_w_uci',
                      'pval_int', 'pval_MA', 'pval_MB', 'pval_DM')
  #QC: All sets must actually have enough data to have the full amount of results
  resultlist = resultlist[lengths(resultlist) == length(these_rownames)]
  
  
  df_result <- as.data.frame(t(as.data.frame(resultlist, 
                                             row.names = these_rownames,
                                             byrow=FALSE)))
  #Order by epi value
  df_results <- arrange(df_result, e_est)
  df_results$row <- c(1:1:dim(df_results)[1])
  df_results <- df_results %>% mutate(Epistasis_Direction = case_when(lowerCI < 0 & upperCI < 0  ~ 'Negative',
                                                                      lowerCI > 0 & upperCI > 0 ~ 'Positive',
                                                                      lowerCI <=  0 & upperCI >= 0 ~ 'Not Detected'))
  #Add gene names
  set_key = unique(df[c("Set","mutant_name")])
  df_results = merge(df_results, set_key, sort=FALSE)
  #Return results
  return(df_results)
}

r.df.tsc.f <- get_epistasis_for_formula(sets_with_flats, f.tsc.f, df_pred_dummy.f, data)
r.df.tsc.f