# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

## Load libraries ----
library(lme4)
library(ggplot2)
library(dplyr)
library(ggpubr)

## Load data ----

data = read.delim("data/compiled_fitness_data_092922.txt", sep = "\t", header = T)

## Format data for analysis ----

numeric_cols <- c('LN', 'WO', "SPF", 'TSC', 'SH')
factor_cols <- c('Set', 'Flat', 'Row', 'Type', 'Genotype', 'MA', 'MB', 'Subline')

data <- data %>% mutate_at(numeric_cols, as.numeric)
data <- data %>% mutate_at(factor_cols, as.factor)

data$MA <- as.integer(data$MA == 'MUT')
data$MB <- as.integer(data$MB == 'MUT')
data$DM <- as.integer(data$Genotype == 'DM')

# create log total seed count 
data$logTSC <- log10(data$TSC)
data = subset(data, is.finite(data$logTSC)) #remove nans

#the setlist is the list of plant sets I'm looking at. 
#Should review this -- supposedly some are duplicate genes. 
setlist <- as.integer(levels(data$Set))
#2r is a repeat of set 2. Here I change it to 2000. 
setlist[is.na(setlist)] <- 2000
setlist <- as.character(setlist)

### Conjunction junction, what's your function?
# Define Functions ----

#1. Get_epi_stats. Gets epistasis values for plant set. 
    # Returns a vector of relevant values


#2. Gets all epistasis values for trait. 
# Returns a dataframe w epistasis values that you can then plot with plot_epi_forest
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



### Functions being worked on ----
get_epi_stats <- function(plantset, df=df, Y=Y) {
  model = lm(Y ~ MA + MB + DM, data=df, subset = (df$Set == plantset), na.action = na.exclude)
  e_est <- as.numeric(coef(model)['DM'])
  lowerCI <- confint(model)[4]
  upperCI <- confint(model)[8]
  rsquared <- summary(model)$r.squared
  pval_e <- as.numeric(summary(model)$coefficients[,4]['DM'])
  
  result <- c(plantset, e_est, lowerCI, upperCI, rsquared, pval_e)
  return(result)
}


### List of formulas to use ----
#should be 16
f.tsc <- formula(logTSC ~ GeneA + GeneB + DM + Experiment)
f.tsc.e
f.tsc.f <- formula(logTSC ~ GeneA + GeneB + DM + Experiment + Experiment/Flat)
f.tsc.ef

f.ln
f.ln.e
f.ln.f
f.ln.ef

f.dtb
f.dtb.e
f.dtb.f
f.dtb.ef

f.spf
f.spf.e
f.spf.f
f.spf.ef

## Getting results and cranking out graphs -- should be a total of 16 figs as currently requested. ####
  
#TSC epistasis results no covariates
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


################################################################################

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
m_edgemiddle <- lm(TSC ~ Type, data=data)
summ <- summary(m_edgemiddle)
summ$adj.r.squared
p <- pf(summ$fstatistic[1],              # Applying pf() function
   summ$fstatistic[2],
   summ$fstatistic[3],
   lower.tail = FALSE)

