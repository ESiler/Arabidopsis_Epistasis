# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

## 1. Load libraries ----
library(lme4)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(gridExtra)

## 2. Load data ----
data = read.delim("data/Fitness_data_83_sets_24Jan2023.txt", sep = "\t", header = T)


## 3. Format data for analysis ----
numeric_cols <- c('LN', 'WO', "SPF", 'TSC', 'SH')
factor_cols <- c('Set', 'Flat', 'Row', 'Type', 'Genotype', 'MA', 'MB', 'Subline')

data <- data %>% mutate_at(numeric_cols, as.numeric)
data <- data %>% mutate_at(factor_cols, as.factor)

data$MA <- as.integer(data$MA == 'MUT')
data$MB <- as.integer(data$MB == 'MUT')
data$DM <- as.integer(data$Genotype == 'DM')

#I want INSIDE to be the default type (vs BORDER)
levels(data$Type) <- c('INSIDE', 'BORDER')

# create log total seed count 
data$logTSC <- log10(data$TSC)
data <- subset(data, is.finite(data$logTSC)) #remove nans

#In Sets, change 2r to 2222 (it is a repeat of exp 2)
level_2r <- which(levels(data$Set) == '2r')
levels(data$Set)[level_2r] <- 2222

#the setlist is the list of plant sets I'm looking at. 
#Should review this -- supposedly contain duplicate genes. 
setlist <- as.character(sort(as.integer(levels(data$Set))))
sets_with_flats <- as.character(sort(as.integer(as.character(unique(data$Set[data$Flat == '2'])))))
sets_without_flats <- setlist[!(setlist %in% sets_with_flats)]


#Load in manual gene name data.
genenamekey <- read.csv('data/genenamekey_cleaned.csv', header = T)
genenamekey <- select(genenamekey, -1)
tail(genenamekey)

#Convert sets to integer to merge w genenamekey dataset
data$Set_int <- as.integer(as.character(data$Set))

#Merge data with gene name dataset
dataj <- left_join(data, genenamekey, by=c('Set_int' = 'Set'))

# remove extra columns. rename columns. 
data <- select(dataj, -Set_int)
data <- rename(data,locusA = MA.y)
data <- rename(data,locusB = MB.y)
data <- rename(data, MA = MA.x, MB = MB.x)

factor_cols_2 <- c('locusA', 'locusB', 'ma', 'mb', 'ma2', 'mb2', 'mutant_name')
data <- data %>% mutate_at(factor_cols_2, as.factor)

rm(dataj)
tail(data)
str(data)

## 4. Generate Models ----
f.tsc <- formula(logTSC ~ MA + MB + DM )
f.tsc.e <- formula(logTSC ~ MA + MB + DM + Type)
f.tsc.f <- formula(logTSC ~ MA + MB + DM + Flat)
f.tsc.ef <- formula(logTSC ~ MA + MB + DM + Type + Flat)
tsc_formulas <- c(f.tsc, f.tsc.e, f.tsc.f, f.tsc.ef)

f.ln <- formula(log10(LN) ~ MA + MB + DM)
f.ln.e <- formula(log10(LN) ~ MA + MB + DM + Type)
f.ln.f <- formula(log10(LN) ~ MA + MB + DM + Flat)
f.ln.ef <- formula(log10(LN) ~ MA + MB + DM + Type + Flat)

f.dtb <- formula(log10(DTB) ~ MA + MB + DM)
f.dtb.e <- formula(log10(DTB) ~ MA + MB + DM + Type)
f.dtb.f <- formula(log10(DTB) ~ MA + MB + DM + Flat)
f.dtb.ef <- formula(log10(DTB) ~ MA + MB + DM + Type + Flat)

f.spf <- formula(log10(SPF) ~ MA + MB + DM)
f.spf.e <- formula(log10(SPF) ~ MA + MB + DM + Type)
f.spf.f <- formula(log10(SPF) ~ MA + MB + DM + Flat)
f.spf.ef <- formula(log10(SPF) ~ MA + MB + DM + Type + Flat)

#Make a collection of dummy df for model prediction
#Dummy matrix for predictions
df_pred_dummy <- data.frame(MA  = c(0,1,0,1),
                            MB = c(0,0,1,1),
                            DM = c(0,0,0,1)
)
rownames(df_pred_dummy) = c("WT", "MA", "MB", "DM")

df_pred_dummy.e <- data.frame(MA  = c(0,1,0,1),
                            MB = c(0,0,1,1),
                            DM = c(0,0,0,1),
                            Type = c("INSIDE","INSIDE","INSIDE","INSIDE")
)
rownames(df_pred_dummy.e) = c("WT", "MA", "MB", "DM")
df_pred_dummy.e

df_pred_dummy.f <- data.frame(MA  = c(0,1,0,1),
                              MB = c(0,0,1,1),
                              DM = c(0,0,0,1),
                              Flat = c('1','1','1','1')
)

rownames(df_pred_dummy.f) = c("WT", "MA", "MB", "DM")
df_pred_dummy.f

rownames(df_pred_dummy) = c("WT", "MA", "MB", "DM")

df_pred_dummy.ef <- data.frame(MA  = c(0,1,0,1),
                              MB = c(0,0,1,1),
                              DM = c(0,0,0,1),
                              Type = c("INSIDE","INSIDE","INSIDE","INSIDE"),
                              Flat = c('1','1','1','1')
)
rownames(df_pred_dummy.ef) = c("WT", "MA", "MB", "DM")
df_pred_dummy.ef



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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 7. Get epistasis results and crank out graphs ----

# should be a total of 8 figs as currently requested
#tsc, ln, dtb, and spf with and without edge effects. 
#But exps with flat effect need to be modeled separately 

### Get epistasis stats and generate results df ----


#1. Get_epi_stats. Gets epistasis values for plant set. 
# Returns a vector of relevant values

#Dummy matrix for predictions

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
  
  pvals <- summary(model)$coefficients[,4]
  #Return stats
  result <- c(plantset, e_est, lowerCI, upperCI, rsquared, pval_e, pred_rel_fit, pvals)

  return(result)
}

tmp <- get_epi_stats(845, f.tsc, df_pred_dummy, data) #Need to make dummy frame flexible based on model
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
#test
test_results2 <- get_epistasis_for_formula(sets_without_flats, f.dtb, df_pred_dummy, data)
tail(test_results2)

#


# Make forest plots ----

# plot_epi_forest: This function makes a rad epistasis forest plot for all the genes. Woohoo!
# Returns a ggplot plot
plot_epi_forest <- function(epi_data, main="Title") {
  print(str(epi_data))
  plot <- ggplot(epi_data, aes(y = row, x = e_est, color=Epistasis_Direction, ymin=1, ymax=dim(epi_data)[1])) +
    geom_point(shape = 18, size = 3) +  
    geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.5) +
    geom_vline(xintercept = 0, color = "black", cex = .5) +
    scale_y_continuous(name = "Gene Pair", 
                       breaks=1:dim(epi_data)[1], 
                       labels = epi_data$Set, 
                       trans = "reverse", 
                       expand = c(0,0.5)) +
    ggtitle(main) +
    xlab("Epistasis Value [95% CI]") +
    scale_color_manual('Epistasis\nDirection',values = c("#D9027D", 'black', 'blue')) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y.left = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x.bottom = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
  plot
}

plot_epi_forest(test_results2, main="Test/Sample Plot")


#### Figures ----

#Results Plotting Function DELUX EDITION!!!
forest_plot_delux <- function(f1_flats, f2_notflats, main="Title"){
  df_f <- get_epistasis_for_formula(sets_with_flats, f1_flats)
  df_wof <- get_epistasis_for_formula(sets_without_flats, f2_notflats)
  df2 <- rbind(df_f, df_wof)
  df2 <- arrange(df2, e_est)
  df2['row'] <- (1:(dim(df2)[1]))
  fplot <- plot_epi_forest(df2, main = main)           
  return(fplot)
}

forest_plot_delux(f.tsc.ef, f.tsc.e, "Epistasis | Total Seed Count | With Edge Effects")


#TSC results not edgy
efp_tsc <- forest_plot_delux(f.tsc.f, f.tsc, "Epistasis | Total Seed Count | No Edge Effects")

#TSC results edgy
efp_tsc_e <- forest_plot_delux(f.tsc.ef, f.tsc.e, "Epistasis | Total Seed Count | With Edge Effects")

#DTB results not edgy
efp_dtb <- forest_plot_delux(f.dtb.f, f.dtb, "Epistasis | Days To Bolt | No Edge Effects")

#DTB results edgy
efp_dtb_e <- forest_plot_delux(f.dtb.ef, f.dtb.e, "Epistasis | Days To Bolt | With Edge Effects")

#LN results not edgy
efp_ln <- forest_plot_delux(f.ln.f, f.ln, "Epistasis | Leaf Number | No Edge Effects")

#LN results edgy
efp_ln_e <- forest_plot_delux(f.ln.ef, f.ln.e, "Epistasis | Leaf Number | With Edge Effects")

#SPF results not edgy
efp_spf <- forest_plot_delux(f.spf.f, f.spf, "Epistasis | Seeds Per Fruit | No Edge Effects")

#SPF results edgy
efp_spf_e <- forest_plot_delux(f.spf.ef, f.spf.e, "Epistasis | Seeds Per Fruit | With Edge Effects")

## Multiplots THE BIG PIG!~~~~

#ggarrange(efp_tsc, efp_tsc_e, nrow=1, legend='right',common.legend = TRUE)
#ggarrange(efp_dtb, efp_dtb_e, nrow=1, legend='right',common.legend = TRUE)
#ggarrange(efp_ln, efp_ln_e, nrow=1, legend='right',common.legend = TRUE)
#ggarrange(efp_spf, efp_spf_e, nrow=1, legend='right',common.legend = TRUE)
# Use models w no edge effects. Edge effects -> supplement

## Create Forest Plots W Gene Names That Compare Different Fitness Values ----
ggarrange(efp_tsc_e, efp_dtb_e, efp_ln_e, efp_spf_e,nrow=1,legend='right',common.legend=T)

#Step 1: Get Y scale from efp_tcs_e
#hack to extract y axis order
forest_plot_delux_axishack <- function(f1_flats, f2_notflats, main="Title"){
  df_f <- get_epistasis_for_formula(sets_with_flats, f1_flats)
  df_wof <- get_epistasis_for_formula(sets_without_flats, f2_notflats)
  df2 <- rbind(df_f, df_wof)
  df2 <- arrange(df2, e_est)
  df2['row'] <- (1:(dim(df2)[1]))
  return(df2$Set)
}
common_y_axis <- (forest_plot_delux_axishack(f.tsc.ef, f.tsc.e, "Epistasis | Total Seed Count | With Edge Effects"))

common_y_axisf <- as.factor(common_y_axis)
#only changes lables, doesn't change y order
#Does not work


efp_dtb_eo <- efp_dtb_e + scale_y_continuous(name="help!",
                                             breaks=1:length(common_y_axis),
                                             limits = common_y_axisf,
                                             labels = common_y_axis,
                                             )

ggarrange(efp_dtb_e, efp_dtb_eo, efp_tsc_e, nrow=1)

scale_y_continuous(name = "Gene Pair", 
                   breaks=1:dim(epi_data)[1], 
                   labels = epi_data$Set, 
                   trans = "reverse", 
                   expand = c(0,0.5)) 
  
  
  
# Is the coefficient for A and B fitness in the linear model  their relative fitness? ----

# 1. Get mean relative fitness (for set 38 as an example)
dataset38 = data[which(data$Set==38), ]

means_set38 = dataset38 %>%
  group_by(Genotype) %>%
  summarise_at(vars(TSC, logTSC, DTB, SN, SPF), list(mean = mean))

means_set38

rel_fitness_38 = means_set38 %>% mutate(across(where(is.numeric), ~./.[Genotype == "WT"]))

# 2. Get lm coefficients for set 38
model38 <- lm(f.tsc, data=dataset38, na.action = na.exclude)

summary(model38)$coefficients[,4]

df_pred_dummy <- data.frame(MA  = c(0,1,0,1),
                  MB = c(0,0,1,1),
                  DM = c(0,0,0,1)
)
rownames(df_pred_dummy) = c("WT", "MA", "MB", "DM")
df_pred_dummy

pred_38 = predict(model38, df_pred_dummy, interval="confidence") ## Value estimates w confints!!!

wt_38_logtsc = pred_38[1,1]

correct_answer_38 = pred_38/wt_38_logtsc #Normalized! But is it legit to normalized confidence intervals? Seems sketchy.
correct_answer_38 ##GETS CORRECT REL FITNESS...BUT ARE CIs CORRECT? (result: Yes.)
rel_fitness_38

#PRE-NORMALIZE RESULTS TO WT AND COMPARE CIs

#normalize logTSC to WT
dataset38['logTSC.norm'] = dataset38['logTSC']/wt_38_logtsc
str(dataset38)
#Make model
f.tscnorm <- formula(logTSC.norm ~ MA + MB + DM )
model38norm <- lm(f.tscnorm, data=dataset38, na.action = na.exclude)
summary(model38norm)
pred_38_norm = predict(model38norm, df_pred_dummy, interval="confidence")
pred_38_norm
correct_answer_38 #Phew CIs are correct.
all.equal(pred_38_norm, correct_answer_38)

set_38_means_df = as.data.frame.array(correct_answer_38)
set_38_means_df['Genotype'] = rownames(set_38_means_df)

level_order <- c("WT", "MA", "MB", "DM")

#Relative Fitness Bargraph w error bars sample! 
ggplot(set_38_means_df, aes(x = factor(Genotype, level = level_order), y = fit)) +
  geom_col() + 
  geom_errorbar(aes(x=Genotype, ymin=lwr, ymax=upr), width=0.3, color="grey", size=1) + 
  xlab("Genotype") +
  ylab("Relative Fitness") +
  ggtitle("Set 38: Rglg1 and Rglg2") +
  theme_classic()
  
#Works! Hooray!!
#Get variables to make dummy table for predictions
summary(model38norm)

str(coef(model38norm))

length(f.dtb.ef)
f.dtb.ef[[1]]
f.dtb.ef[[2]]
pred_vars <- f.dtb.ef[[3]] #WOOHOOOOOOO!

pred_vars_c = deparse(pred_vars)

pred_vars = strsplit(pred_vars_c, " . ")[[1]]


#Okay now get string in a format that can be used to fill in a table ----


j#Fig 2: Multiplot of Fitness Values by Genotype ----
#1. Fix NAs (2222s)
#2. Normalized by WT fitness
#3. Code by epistasis value
#4. Add CIs/error bars
#5. Add titles
facetplotprelim <- function(Fitness_Value){
  p <- ggplot(data, aes(x=Genotype, y=Fitness_Value)) + 
    geom_bar(position = "dodge",
              stat = "summary",
              fun = "mean") +
    facet_wrap(~mutant_name) +
    scale_x_discrete(limits = c("WT", "MA", "MB", "DM")) +
    theme_bw()
  
  return(p)
  }

facetplotprelim(data$TSC)
facetplotprelim(data$logTSC)
facetplotprelim(data$DTB)
facetplotprelim(data$LN)
facetplotprelim(data$SPF)
  
str(data)
### A. Get results for non-edge models; B. make rel fit plots. ----


