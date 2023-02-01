#1. Load processed data, formulas, libraries, and other objects:
source('scripts/load_required_packages.R')
load("rdata/02_workspace.RData")


## 5. Model Comparison and Selection: CAN SKIP TO PART 7 HERE----

#takes experiment number, dataframe, and two fomulas
#returns model 1 adjusted r-squared,
#model 2 adjusted r-squared, and
#anova p-value
mod_comp_fit <- function(exp, df, f1, f2){
  dfsubset <- subset(df, Set == exp)
  m1 <- lm(f1, dfsubset) 
  m2 <- lm(f2, dfsubset)
  sm1 <- summary(m1)
  sm2 <- summary(m2)
  a <- anova(m1,m2)
  m1ar2 <- sm1$adj.r.squared
  m2ar2 <- sm2$adj.r.squared
  apval <- a$`Pr(>F)`[2]
  results <- data.frame(exp, m1ar2, m2ar2, apval)
  #print(results)
  return(results)
}

#mod_comp_fit(11, data, f.tsc, f.tsc.e) #test data

#Input data, list of experiments, and formulas
#Output a dataframe of model comp stats
mod_comp_fit_loop <- function(explist, df, f1, f2){
  dfr = NULL
  for(i in explist){
    r <- mod_comp_fit(i,df=df,f1=f1,f2=f2)
    dfr=rbind(dfr, r)
  }
  return(dfr)
}

## 6. Generate model comparison heat maps----
make_r2_heatmap <- function(df, main='Title'){
  #number with significant pvals from anova
  edge_exps <- df$exp[(df$apval < .05)]
  no_edge_exps <- df$exp[(df$apval >= .05)]
  nume <- length(edge_exps)
  numnoe <- length(no_edge_exps)
  #number where adj r2 is higher for edge model
  nume2 <- sum(df$m2ar2 > df$m1ar2)
  numnoe2 <- sum(df$m2ar2 <= df$m1ar2)
  #SUMMARY DATA
  r <- paste('ANOVA: m2 wins ', nume, 
             ', m1 wins ', numnoe)
  r2 <- paste('Adj r^2: m2 wins', nume2,
              ', m1 wins ', numnoe2)
  #convert to long for heatmap I guess
  df_m <- df %>% pivot_longer(cols= c('m1ar2', 'm2ar2'),
                              names_to='model',
                              values_to='adj_r_squared')
  #CREATE PLOT
  ggplot(df_m, aes(model, exp, fill=adj_r_squared)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="slateblue4") +
    labs(title = r2, subtitle = r)
}

#### 6a: Edge vs. no edge----
#make dfs for model comps
df_comp_tsc_edge <- mod_comp_fit_loop(setlist, data, f.tsc, f.tsc.e)
df_comp_ln_edge <- mod_comp_fit_loop(setlist, data, f.ln, f.ln.e)
df_comp_dtb_edge <- mod_comp_fit_loop(setlist, data, f.dtb, f.dtb.e)
df_comp_spf_edge <- mod_comp_fit_loop(setlist, data, f.spf, f.spf.e)

#Generates model comp heat maps and results
#for edge effects with different y variables
p9 <- make_r2_heatmap(df_comp_tsc_edge)
p10 <- make_r2_heatmap(df_comp_ln_edge)
p11 <- make_r2_heatmap(df_comp_dtb_edge)
p12 <- make_r2_heatmap(df_comp_spf_edge)

grid.arrange(p9, p10, p11, p12, nrow=1)
#conclusion: inconclusive

#### 6b: Flat vs. No Flat ----
df_comp_tsc_flat <- mod_comp_fit_loop(sets_with_flats, data, f.tsc, f.tsc.f)
df_comp_ln_flat <- mod_comp_fit_loop(sets_with_flats, data, f.ln, f.ln.f)
df_comp_dtb_flat <- mod_comp_fit_loop(sets_with_flats, data, f.dtb, f.dtb.f)
df_comp_spf_flat <- mod_comp_fit_loop(sets_with_flats, data, f.spf, f.spf.f)

p1 <- make_r2_heatmap(df_comp_tsc_flat)
p2 <- make_r2_heatmap(df_comp_ln_flat)
p3 <- make_r2_heatmap(df_comp_dtb_flat)
p4 <- make_r2_heatmap(df_comp_spf_flat)

grid.arrange(p1,p2,p3,p4,nrow=1)
#Conclusion: include flat

#### 6c: Flat with and without edge effects----
df_comp_tsc_flat_e <- mod_comp_fit_loop(sets_with_flats, data, f.tsc.f, f.tsc.ef)
df_comp_ln_flat_e <- mod_comp_fit_loop(sets_with_flats, data, f.ln.f, f.ln.ef)
df_comp_dtb_flat_e <- mod_comp_fit_loop(sets_with_flats, data, f.dtb.f, f.dtb.ef)
df_comp_spf_flat_e <- mod_comp_fit_loop(sets_with_flats, data, f.spf.f, f.spf.ef)

p5 <- make_r2_heatmap(df_comp_tsc_flat)
p6 <- make_r2_heatmap(df_comp_ln_flat)
p7 <- make_r2_heatmap(df_comp_dtb_flat)
p8 <- make_r2_heatmap(df_comp_spf_flat)

grid.arrange(p5,p6,p7,p8,nrow=1)
#Conclusion: Always separate out flat experiments and include flat. 


#Run analysis with and without edge effects. 
#### 6d: EDA Model Comparison ----

#subset set 1
datas1 <- subset(data, data$Set == 1)

#PLot data -- flat doesn't look that important
hist(datas1$logTSC)
plot(logTSC ~ Flat, data = datas1)

#No significant difference caused by flat show in this model
mod1 <- lm(logTSC ~ Flat, data=datas1)
summary(mod1)

#Compare model w and w/out flat
mod2 <- lm(logTSC ~ MA + MB + DM + Flat, data=datas1, na.action = na.exclude)
mod3 <- lmer(logTSC ~ MA + MB + DM + (1|Flat), data=datas1, na.action = na.exclude)

summary(mod2)
summary(mod3)

anova(mod1, mod2)

##Function to plot effects of flat

flatfunclogTSC <- function(x) {
  datan <- subset(data, data$Set == x)
  modn <- lm(logTSC ~ Flat, data=datan)
  modn <- summary(modn)
  r_squared <- modn$adj.r.squared
  P_val <- pf(modn$fstatistic[1],              # Applying pf() function
              modn$fstatistic[2],
              modn$fstatistic[3],
              lower.tail = FALSE)
  stats <- paste("R_squared:", signif(r_squared, digits=2),";", "p-value:", signif(P_val,digits=2))
  plot(logTSC ~ Flat, data = datan, main = paste("Set",x), xlab=stats)
}

#flatfunclogTSC(20)


par(mfrow = c(3, 4))
for (i in sets_with_flats) {
  flatfunclogTSC(i)
}

data_from_flats <- subset(data, data$Set %in% sets_with_flats)

plot(LN ~ Flat, data = data_from_flats)
m_allflats <- lm(LN ~ Flat, data = data_from_flats)
summary(m_allflats)
#add gene names
#Change sig_score to something else

#do this analysis for the mapk data
#(and make cool heat map for same)

plot(TSC ~ Type, data=data)
plot(logTSC ~ Type, data=data)

m_edgemiddle <- lm(TSC ~ Type, data=data)
summ <- summary(m_edgemiddle)
summ$adj.r.squared
p <- pf(summ$fstatistic[1],              # Applying pf() function
        summ$fstatistic[2],
        summ$fstatistic[3],
        lower.tail = FALSE)
