# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

## Load libraries and data ----
library(lme4)
library(ggplot2)
library(dplyr)
library(ggpubr)

data = read.delim("data/compiled_fitness_data_092922.txt", sep = "\t", header = T)

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
data = subset(data, is.finite(data$logTSC)) #remove nans

#the setlist is the list of plant sets I'm looking at. 
#Should review this -- supposedly some are duplicate genes. 
setlist <- as.integer(levels(data$Set))
#2r is a repeat of set 2. Here I change it to 2000. 
setlist[is.na(setlist)] <- '2r'
setlist <- as.character(setlist)

### Conjunction junction, what's your function?
# Define Functions ----

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

#2. Gets all epistasis values for trait. 
# Returns a dataframe w epistasis values that you can then plot with plot_epi_forest

#This is gonna be broken now. 
get_epistasis_fortrait <- function(plantsets, Y, df=data) {
 
  resultlist <- list()
  
  for (i in plantsets){
    v <- get_epi_stats(i, df=df, Y=Y)
    resultlist[[i]] <- as.numeric(v)
  }
  df_result <- as.data.frame(t(as.data.frame(resultlist, 
                                             row.names = c('Set', 'e_est', 'lowerCI', 'upperCI', 'rsquared', 'pval_e'),
                                             byrow=FALSE)))
  
  df_results <- arrange(df_result, e_est)
  df_results$row <- c(1:1:dim(df_results)[1])
  
  df_results <- df_results %>% mutate(Epistasis_Direction = case_when(lowerCI < 0 & upperCI < 0  ~ 'Negative',
                                                                      lowerCI > 0 & upperCI > 0 ~ 'Positive',
                                                                      lowerCI <=  0 & upperCI >= 0 ~ 'Not Detected'))
}

# plot_epi_forest: This function makes a rad epistasis forest plot for all the genes. Woohoo!
  # Returns a ggplot plot
plot_epi_forest <- function(epi_data, main="Title") {
  plot <- ggplot(epi_data, aes(y = row, x = e_est, color=Epistasis_Direction, ymin=1, ymax=dim(epi_data)[1])) +
    geom_point(shape = 18, size = 3) +  
    geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.5) +
    geom_vline(xintercept = 0, color = "black", cex = .5) +
    scale_y_continuous(name = "Gene Pair", breaks=1:dim(epi_data)[1], labels = epi_data$Set, trans = "reverse", expand = c(0,0.5)) +
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

## Step 1: Edge vs. No Edge ----

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
mod_comp_fit(11, data, f.tsc, f.tsc.e)

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

df_comp_tsc_edge <- mod_comp_fit_loop(setlist, data, f.tsc, f.tsc.e)






## I think I'm just going to make an adjusted r2 heatmap for model comparison

#Function to return adjusted r-squared valued from list of formulas ----
mod_comp_ar2 <- function(exp, df, formulalist){
  #Get subset from specific experiment
  dfsubset <- subset(df, Set == exp)
  #initialize results vector (should be 4 adjr2s)
  v <- vector(length = length(formulalist))
  i <- 0
  #Get adjr2 values and add to 
  for(item in formulalist){
    m <- lm(item, dfsubset) 
    i <- i + 1
    summ <- summary(m)
    adjr2 <- summ$adj.r.squared
    v[i] <- adjr2
  }
  return(v)
} 
#mod_comp_ar2(1,data,tsc_formulas) works
m <- matrix(nrow = 4)
for(item in sets_with_flats){
  stuff <- mod_comp_ar2(item, data, tsc_formulas)
  print(stuff)
}

sets_with_flats <- setlist[1:8]
mod_comp_ar2(25, data, tsc_formulas)



# TSC  ----




# LN ----



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




### Getresults and crank out graphs -- should be a total of 16 figs as currently requested. ####

#TSC epistasis results no covariates
get_epi_stats(11, f.tsc.e, data)
## This will be broken til I fix the other function

df_logTSC_results <- get_epistasis_fortrait(setlist, data$logTSC)
plot1 <- plot_epi_forest(df_logTSC_results, "Trait: Total Seed Count")
plot1

#Days to Bolt epistasis results no covariates
df_DTB_results <- get_epistasis_fortrait(setlist, log10(data$DTB))
plot2 <- plot_epi_forest(df_DTB_results, "Trait: Days to Bolt")
plot2

#leaf number epistasis results no covariates
df_LN_results <- get_epistasis_fortrait(setlist, log10(data$LN))
plot3 <- plot_epi_forest(df_LN_results, "Trait: Leaf Number")
plot3

#Seed per Fruit epistasis results no covariates
df_SPF_results <- get_epistasis_fortrait(setlist, log10(data$SPF))
plot4 <- plot_epi_forest(df_SPF_results, "Trait: Seeds per Fruit")
plot4



