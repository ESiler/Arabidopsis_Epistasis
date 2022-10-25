# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

## Load libraries and data ----
library(lme4)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(gridExtra)

data = read.delim("data/compiled_fitness_data_092922.txt", sep = "\t", header = T)
genenamekey = read.delim("data/genes_in_each_set.txt", sep = "\t", header = T)
## Format data for analysis ----

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
levels(data$Set)[10] <- 2222

#the setlist is the list of plant sets I'm looking at. 
#Should review this -- supposedly contain duplicate genes. 
setlist <- as.character(sort(as.integer(levels(data$Set))))
sets_with_flats <- as.character(sort(as.integer(as.character(unique(data$Set[data$Flat == '2'])))))
sets_without_flats <- setlist[!(setlist %in% sets_with_flats)]

#Add gene pair names to data
genenamekey['gene.pair'] <- paste(genenamekey$MA, '.', genenamekey$MB, sep='')
data['Genes'] <- genenamekey$gene.pair[data$Set]


### Conjunction junction, what's your function?
# Define Functions ----

### Models and model comparison -----
# List of formulas/models
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

## Model Comparison and Selection----

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

#make model comparison heatmaps----
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

####Edge vs. no edge----
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

####Flat vs. No Flat ----
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

####Flat with and without edge effects----
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
## EDA ----

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
setswflats <- c(1, 2, 3, 5, 7, 8, 11, 12, 14, 15, 18, 20)

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
for (i in setswflats) {
  flatfunclogTSC(i)
}

data_from_flats <- subset(data, data$Set %in% setswflats)

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

### Get epistasis results and crank out graphs ----

# should be a total of 8 figs as currently requested
#tsc, ln, dtb, and spf with and without edge effects. 
#But exps with flat effect will need to be modeled separately 

## Functions for this section: ----
#1. Get_epi_stats. Gets epistasis values for plant set. 
# Returns a vector of relevant values
get_epi_stats <- function(plantset, formula, df){
  dfsubset <- subset(df, Set == plantset)
  model <- lm(formula, data=dfsubset, na.action = na.exclude)
  #Extract stats
  e_est <- as.numeric(coef(model)['DM'])
  lowerCI <- confint(model)[4,1]
  upperCI <- confint(model)[4,2]
  rsquared <- summary(model)$adj.r.squared
  pval_e <- as.numeric(summary(model)$coefficients[,4]['DM'])
  #Return stats
  result <- c(plantset, e_est, lowerCI, upperCI, rsquared, pval_e)
  return(result)
}
tmp <- get_epi_stats(2222, f.tsc.e, data) #Works :)

#2. Gets all epistasis values for trait. 
# Returns a dataframe w epistasis values that you can then plot with plot_epi_forest
get_epistasis_for_formula <- function(plantsets, formula, df=data) {
  
  resultlist <- list()
  
  for (i in plantsets){
    v <- get_epi_stats(i, formula=formula, df=df)
    resultlist[[i]] <- as.numeric(v)
  }
  df_result <- as.data.frame(t(as.data.frame(resultlist, 
                                             row.names = c('Set', 'e_est',
                                                           'lowerCI', 'upperCI', 
                                                           'rsquared', 'pval_e'),
                                             byrow=FALSE)))
  df_results <- arrange(df_result, e_est)
  df_results$row <- c(1:1:dim(df_results)[1])
  df_results <- df_results %>% mutate(Epistasis_Direction = case_when(lowerCI < 0 & upperCI < 0  ~ 'Negative',
                                                                      lowerCI > 0 & upperCI > 0 ~ 'Positive',
                                                                      lowerCI <=  0 & upperCI >= 0 ~ 'Not Detected'))
}
#test
test_results <- get_epistasis_for_formula(as.character(sets_with_flats), f.dtb.f)

# plot_epi_forest: This function makes a rad epistasis forest plot for all the genes. Woohoo!
# Returns a ggplot plot
plot_epi_forest <- function(epi_data, main="Title") {
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

plot_epi_forest(test_results, main="whyyyyy")


#### Results! ----

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

ggarrange(efp_tsc, efp_tsc_e, nrow=1, legend='right',common.legend = TRUE)
ggarrange(efp_dtb, efp_dtb_e, nrow=1, legend='right',common.legend = TRUE)
ggarrange(efp_ln, efp_ln_e, nrow=1, legend='right',common.legend = TRUE)
ggarrange(efp_spf, efp_spf_e, nrow=1, legend='right',common.legend = TRUE)

##Comp plot for different Ys ----

ggarrange(efp_tsc_e, efp_dtb_e, efp_ln_e, efp_spf_e,nrow=1,legend='right',common.legend=T)
#Goal: Multiplot w common Y axis showing all 4 epi estimates.

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
  
  
  
  
  
#Fig 2: Multiplot of Fitness Values by Genotype ----

facetplotprelim <- function(Fitness_Value){
  p <- ggplot(data, aes(x=Genotype, y=Fitness_Value)) + 
    geom_bar(position = "dodge",
              stat = "summary",
              fun = "mean") +
    facet_wrap(~Set) +
    scale_x_discrete(limits = c("WT", "MA", "MB", "DM")) +
    theme_bw()
  
  return(p)
  }

facetplotprelim(data$TSC)
facetplotprelim(data$logTSC)  
facetplotprelim(data$DTB)
facetplotprelim(data$LN)
facetplotprelim(data$SPF)
  
