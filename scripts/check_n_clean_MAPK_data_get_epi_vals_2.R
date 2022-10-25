#Load Libraries ----
library(dplyr)
require(lme4)
require(ggplot2)
require(ggpubr)
source("~/Documents/Code/ElliesFavFunctions.R")
library(igraph)
library(ggraph)

### Conjunction Junction, what's your function? ----

#This function turns WT to 0 and MUT to 1. 
func_t1 <- function(x){
  x = ifelse(x=='WT', 0, 1)
  return(x)
}

#This function gets wt/genea/geneb/dm subset of dataset for given gene pair
subset_genepair <- function(p=genpair, df=df){
  genea <- sub('_.*', "", p)
  geneb <- sub('[0-9]*_', "", p)
  df_sub <- subset(df, Genotype == ('Col') | 
                     Genotype == genea | 
                     Genotype == geneb|
                     Genotype == p)
  # Add GeneA, GeneB, and DM columns. All should be 1 or 0. 
  df_sub <- df_sub %>% mutate(DM = ifelse(mutnum==2, 1, 0))
  df_sub <- df_sub %>% mutate(GeneA = ifelse((Genotype==genea|mutnum==2), 1, 0))
  df_sub <- df_sub %>% mutate(GeneB = ifelse((Genotype==geneb|mutnum==2), 1, 0))
  return(df_sub)
}

##This function takes in a list of gene pairs and outputs anova model comparisons:
mod_comp_mapk <- function(genpair, df, f1, f2){
  dfsubset <- subset_genepair(p=genpair, df=df)
  m1 <- lm(f1, dfsubset) #m1
  m2 <- lm(f2, dfsubset)
  aov <- anova(m1, m2)
  pval <- aov[6][2,1]
  if (pval < 0.05) signif= "TRUE" else signif="FALSE"
  print(paste("The p value for gene pair ", genpair,  " is: ", pval, "p<.05 = ", signif))
  return(aov)
}

# A function to collect epistasis estimates.
# Takes gene pair, model formula, and data and return epistasis stats of interest:
get_epi_mapk_stats <- function(genpair, formula, df){
  dfsubset <- subset_genepair(p=genpair, df=df)
  model <- lm(formula, dfsubset) #m1
  e_est <- as.numeric(coef(model)['DM'])
  lowerCI <- confint(model)[4,1]
  upperCI <- confint(model)[4,2]
  rsquared <- summary(model)$adj.r.squared
  pval_e <- as.numeric(summary(model)$coefficients[,4]['DM'])
  
  result <- c(genpair, e_est, lowerCI, upperCI, rsquared, pval_e)
  return(result)
}

#This function makes a df of epistasis stats for all gene pairs:
get_epistasis_all_mapk <- function(plantsets, formula, df=data) {
  #create empty results list
  resultlist <- list()
  
  #Loop through gene pairs adding epistasis stats for all
  for (i in plantsets){
    v <- get_epi_mapk_stats(i, df=df, formula=formula)
    resultlist[[i]] <- v
  }
  
  
  #organize/format data
  df_result <- as.data.frame(t(as.data.frame(resultlist, 
                                             row.names = c('Genepair', 'e_est', 'lowerCI', 'upperCI', 'rsquared', 'pval_e'),
                                             byrow=FALSE)))
  #Make numeric values numeric
  df_result$e_est <- as.numeric(df_result$e_est)
  df_result$lowerCI <- as.numeric(df_result$lowerCI)
  df_result$upperCI <- as.numeric(df_result$upperCI)
  df_result$rsquared <- as.numeric(df_result$rsquared)
  df_result$pval_e <- as.numeric(df_result$pval_e)
  
  #sort by epistasis estimate
  df_results <- arrange(df_result, e_est)
  #Add result column to color code.
  df_results <- df_results %>% mutate(Epistasis_Direction = case_when(lowerCI < 0 & upperCI < 0  ~ 'Negative',
                                                                      lowerCI > 0 & upperCI > 0 ~ 'Positive',
                                                                      lowerCI <=  0 & upperCI >= 0 ~ 'Not Detected'))
  #Add 'row' column for plotting in correct order
  df_results$row <- c(1:1:dim(df_results)[1])
  return(df_results)
}


### Load + Process Data ----
setwd("~/Documents/Projects/Arabidopsis_fitness_Melissa")
df = read.delim("MAPK_DEPI_data_080522.txt", sep = "\t", header = T)
dfOG <- df

#good news someone already got rid of mapk 6_20
factorcol <- c(1, 2, 4, 6:8,14:25)
df[,factorcol] <- lapply(df[,factorcol], factor)

#Turn WT to 0 and MUT to 1
df <- df %>% mutate_at(vars(matches("MPK")), func_t1)

#Add column for mut_num
df <- df %>% mutate(mutnum = MPK1 + MPK3 + MPK5 + MPK6 + MPK8 + MPK9 + MPK13 + 
                MPK14 + MPK16 + MPK17 + MPK18 + MPK20)

#Change "Flat" so that the flat 1's don't match each other. 
df <- df %>% mutate(Flat1 = paste(Experiment, Flat, sep="_"))

#Some TSCs are 0 and cause the model to break
df_tsc <- df %>% mutate(logTSC = log10(TSC))
df_tsc <- subset(df_tsc, is.finite(df_tsc$logTSC)) #Is it bad to throw these out?

df_SN <- df %>% mutate(logSN = log10(SN))
df_SN <- subset(df_SN, is.finite(df_SN$logSN))

df_SPF <- df %>% mutate(logSPF = log10(SPF))
df_SPF <- subset(df_SPF, is.finite(df_SPF$logSPF))


### Model generation and model comparison ----

#List of double mutants (for looping)
double_mutants <- levels(df$Genotype)[grep("_", levels(df$Genotype))]

# Model formulas:
f1 <- formula(logTSC ~ GeneA + GeneB + DM + Experiment)
f2 <- formula(logTSC ~ GeneA + GeneB + DM + Experiment + Experiment/Flat)
f3 <- formula(log(SN) ~ GeneA + GeneB + DM + Experiment)
f4 <- formula(log(SN) ~ GeneA + GeneB + DM + Experiment + Experiment/Flat)
f5 <- formula(log(SPF) ~ GeneA + GeneB + DM + Experiment)
f6 <- formula(log(SPF) ~ GeneA + GeneB + DM + Experiment + Experiment/Flat)


#Model comparison, include flat vs. no flat, total seed count
for (item in double_mutants){
  mod_comp_mapk(item, df_tsc, f1, f2)
} 
#Flat: 18 | No Flat: 5

#Model comparison, include flat vs. no flat, Silique Number
for (item in double_mutants){
  mod_comp_mapk(item, df_SN, f3, f4)
}
# Flat: 20 | No Flat: 3

#Model comparison, include flat vs. no flat, Seeds per Fruit
for (item in double_mutants){
  mod_comp_mapk(item, df_SPF, f5, f6)
}

# Flat: 6 | No Flat: 17


### Get epistasis values:

#get epistasis df for TSC
f1_logTSC_epivals <- get_epistasis_all_mapk(double_mutants, f1, df=df_tsc)
f2_logTSC_epivals <- get_epistasis_all_mapk(double_mutants, f2, df=df_tsc)

#get epistasis df for silique number bc why not
f3_logSN_epivals <- get_epistasis_all_mapk(double_mutants, f3, df=df_SN)
f4_logSN_epivals <- get_epistasis_all_mapk(double_mutants, f4, df=df_SN)

#epistasis for seeds per fruit bc why not
f5_logSPF_epivals <- get_epistasis_all_mapk(double_mutants, f5, df=df_SPF)
f6_logSPF_epivals <- get_epistasis_all_mapk(double_mutants, f6, df=df_SPF)

#odd functions are w/out flat; even functions are with
all_epivals_f <- bind_rows("TCS" = f2_logTSC_epivals, 
                           "SN" = f4_logSN_epivals, 
                           "SPF" = f6_logSPF_epivals,
                           .id="metric")
all_epivals_f$Genepair <- as.factor(all_epivals_f$Genepair)
#Make Figures ----

#Fig 1: Multiplot of Forest Plots ----
#OG forest plot funtion; best for single plots. See plot_epi_forest2
plot_epi_forest <- function(epi_data, main="Title") {
  plot <- ggplot(epi_data, aes(y = row, x = e_est, color=Epistasis_Direction, ymin=1, ymax=dim(epi_data)[1])) +
    geom_point(shape = 18, size = 3) +  
    geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.5) +
    geom_vline(xintercept = 0, color = "black", cex = .5) +
    scale_y_continuous(name = "Gene Pair", breaks=1:dim(epi_data)[1], labels = epi_data$Genepair, trans = "reverse", expand = c(0,0.5)) +
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



plotlogTSC <- plot_epi_forest(f1_logTSC_epivals, main="Epistasis Values\nTotal Seed Count\nw/out flat")
plotlogTSCf <- plot_epi_forest(f2_logTSC_epivals, main="Epistasis Values\nTotal Seed Count\nwith flat")

plotlogSN <- plot_epi_forest(f3_logSN_epivals, main="Epistasis Values\nSilique Number\nw/out flat")
plotlogSNf <- plot_epi_forest(f4_logSN_epivals, main="Epistasis Values\nSilique Number\nwith flat")

plotlogSPF <- plot_epi_forest(f5_logSPF_epivals, main="Epistasis Values\nSeeds Per Fruit\nw/out flat")
plotlogSPFf <- plot_epi_forest(f6_logSPF_epivals, main="Epistasis Values\nSeeds Per Fruit\nwith flat")

#Three plots
threeepiplots <- ggarrange(plotlogTSCf, plotlogSNf, plotlogSPFf, 
          ncol=3,
          common.legend=TRUE,
          legend='right')
threeepiplots

#annotate_figure(threeepiplots, 
#                top= text_grob('Epistasis Values', face = "bold", size = 15))

#Six plots
ggarrange(plotlogTSC, plotlogSN, plotlogSPF, 
          plotlogTSCf, plotlogSNf, plotlogSPFf) #Woot woot!


#Publication quality plot: 

# X+include flat effects for all 
# +Include/Add real gene names (NA for this one??)
# +Use same gene list for Y axis for easy comparison
# +Use same scale for all epi vals on X axis
#Pub quality plot
y_order = factor(rownames(f2_logTSC_epivals), levels = rownames(f2_logTSC_epivals))
y_order

#New and improved plotting function Best for multi-plots
plot_epi_forest2 <- function(epi_data, main="Title") {
  plot <- ggplot(epi_data, aes(y = Genepair, x = e_est, color=Epistasis_Direction, ymin=1, ymax=dim(epi_data)[1])) +
    geom_point(shape = 18, size = 3) +  
    geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.5) +
    geom_vline(xintercept = 0, color = "black", cex = .5) +
    scale_y_discrete(name = "Gene Pair", limits=(rev(y_order))) +
    ggtitle(main) +
    xlab("Epistasis Value [95% CI]") +
    scale_color_manual('Epistasis\nDirection',values = c("#D9027D", 'black', 'blue')) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_rect(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y.left = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x.bottom = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
  plot <- plot + coord_cartesian(xlim = c(-.6, .6))
  return(plot)
}
plotlogTSCf2 <- plot_epi_forest2(f2_logTSC_epivals, main="Total Seed Count")
plotlogSNf2 <- plot_epi_forest2(f4_logSN_epivals, main="Silique Number")
plotlogSPFf2 <- plot_epi_forest2(f6_logSPF_epivals, main="Seeds Per Fruit")


ggthemeb <- theme(axis.text.y = element_blank(),
                  axis.text.y.left = element_blank(),
                  axis.title.y=element_blank(), 
                  axis.ticks.y=element_blank(),
                  axis.line.y=element_blank()
                  )

Figure1 <- ggarrange((plotlogTSCf2 + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0.5), "pt"))), 
                     (plotlogSNf2 + ggthemeb + theme(plot.margin = unit(c(5.5, 0.5, 5.5, 0.5), "pt"))), 
                     (plotlogSPFf2 + ggthemeb + theme(plot.margin = unit(c(5.5, 0.5, 5.5, 0.5), "pt"))), 
                     ncol=3,
                     common.legend=TRUE,
                     legend='right',
                     widths=c(3.9,3,3))

mapkfig1 <- annotate_figure(Figure1, 
                            top = text_grob("Fig 1: Epistasis Values for Map Kinase Double Mutants", size=16))
mapkfig1

#Fig S1: ----
plotlogTSC2 <- plot_epi_forest2(f1_logTSC_epivals, main="Total Seed Count")
plotlogSN2 <- plot_epi_forest2(f3_logSN_epivals, main="Silique Number")
plotlogSPF2 <- plot_epi_forest2(f5_logSPF_epivals, main="Seeds Per Fruit")

FigureS1 <- ggarrange((plotlogTSC2 + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0.5), "pt"))), 
                     (plotlogSN2 + ggthemeb + theme(plot.margin = unit(c(5.5, 0.5, 5.5, 0.5), "pt"))), 
                     (plotlogSPF2 + ggthemeb + theme(plot.margin = unit(c(5.5, 0.5, 5.5, 0.5), "pt"))), 
                     ncol=3,
                     common.legend=TRUE,
                     legend='right',
                     widths=c(3.9,3,3))

mapkfigS1 <- annotate_figure(FigureS1, 
                            top = text_grob("Fig S1: Epistasis Values for Map Kinase Double Mutants (excluding flat effect)", size=16))
mapkfigS1



#Fig 2 and Fig S2: Heatmaps ----

## Function that converts epistasis dataframes into heatmap form:
format_data_for_heatmap <- function(df, genes){
  df['GeneA'] <- gsub("_\\d*", '', df$Genepair)
  df['GeneB'] <- gsub("\\d*_", '', df$Genepair)
  df['GeneA'] <- factor(df$GeneA, levels = genes)
  df['GeneB'] <- factor(df$GeneB, levels = genes)
  dff <- df
  dff['GeneA'] <- df['GeneB']
  dff['GeneB'] <- df['GeneA']
  resultdf <- rbind(df, dff)
  return(resultdf)
}
##Functions to make heat map:
mapepiheat <- function(df, main="Title"){
  ggplot(df, aes(x = GeneA, y = GeneB, fill = Epistasis_Direction)) +
    geom_tile() +
    coord_fixed() + 
    ggtitle(main) +
    scale_fill_manual(values = c("Negative" = "#D9027D", "Not Detected" = "#FFFFCC", "Positive" = "blue")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill='grey'), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
}

# This function does a gradient (which seems sketchy on ns data but idk. Should combine these two somehow)
mapepiheatg <- function(df, main="Title"){
  ggplot(df, aes(x = GeneA, y = GeneB, fill = e_est)) +
    geom_tile() +
    scale_fill_gradient2(low = "#D9027D", mid = "#FFFFCC", high = "blue", na.value='purple') +
    coord_fixed() +
    ggtitle(main) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill='grey'),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
}
heatplotTSCg <- mapepiheatg(f2_logTSC_epivals, main="Epistasis: Total Seed Count")

##Format data:
genelevels <- c("mpk1", "mpk3","mpk5", "mpk6", "mpk8", "mpk9", "mpk13", "mpk14", "mpk16", "mpk17", "mpk18", "mpk20")

df_heatmap_TSC <- format_data_for_heatmap(df=f2_logTSC_epivals, genes=genelevels)
df_heatmap_SN <- format_data_for_heatmap(df=f4_logSN_epivals, genes=genelevels)
df_heatmap_SPF <- format_data_for_heatmap(df=f6_logSPF_epivals, genes=genelevels)

## Generate plots
#heat map positive/negative/nd | TSC
heatplotTSC <- mapepiheat(df_heatmap_TSC, main="Epistasis: Total Seed Count")
heatplotTSCg <- mapepiheatg(df_heatmap_TSC, main="Epistasis: Total Seed Count")
##heat maps positive/negative/nd | SN
heatplotSN <- mapepiheat(df_heatmap_SN, main="Epistasis: Silique Number")
heatplotSNg <- mapepiheatg(df_heatmap_SN, main="Epistasis: Silique Number")
##heat maps positive/negative/nd | SPF
heatplotSPF <- mapepiheat(df_heatmap_SPF, main="Epistasis: Seeds Per Fruit")
heatplotSPFg <- mapepiheatg(df_heatmap_SPF, main="Epistasis: Seeds Per Fruit")

##Display plots :)
ggarrange(heatplotTSCg, heatplotSNg, heatplotSPFg,
          heatplotTSC, heatplotSN, heatplotSPF, 
          ncol=3, nrow=2,
          common.legend=TRUE,
          legend='right')

gradientplots <- ggarrange(heatplotTSCg, heatplotSNg, heatplotSPFg,
                           ncol=3, nrow=1,
                           common.legend=TRUE,
                           legend='none')

sigplots <- ggarrange(heatplotTSC, heatplotSN, heatplotSPF, 
                      ncol=3, nrow=1,
                      common.legend=TRUE,
                      legend='none')

ggarrange(gradientplots, sigplots, ncol=1, nrow=2, common.legend = TRUE)

#Add labels (gradient vs non-gradient) or combine info. 
#Black box for sig? Red/Blue box for sig? Would I need a beige box for ns then??
#Add legend



#Fig 3 and Fig 3S: Network Diagrams----

#each mapk is a node -- check
#gene interactions are edges -- check
#Line width also coded by abs val of epi estimate -- check
#p <05 solid line, p > .05 dotted line -- check
#epi estimate is hue -- magenta-beige-blue (diverging color scheme) -- check

nodes_mapk <- matrix(genelevels)

##Function to make plots :)
plot_gene_network_epi <- function(datadf, main='Graph Title', nodes=nodes_mapk){
  #Create edged dataframe
  edges_mapk_TSC <-datadf[1:((dim(datadf)/2)[1]),
                          c('GeneA', 'GeneB', 'e_est', 'Epistasis_Direction', 'pval_e')]
  #Create graph object
  mapk_graph_data <- graph_from_data_frame(edges_mapk_TSC, directed = FALSE, vertices = nodes)
  #Edge width corresponds to absolute value of epistasis
  edge_widths <- abs(edges_mapk_TSC$e_est)*50 + 1
  #Add line type to df, with solid for statistically significant and dotted for ns
  edges_mapk_TSC <- mutate(edges_mapk_TSC, lty = case_when(
    pval_e < .05 ~ 1,
    pval_e >= .05 ~ 3
  ))
  edge_lty <- edges_mapk_TSC$lty
  #epistasis estimates:
  edge_e_vals <- edges_mapk_TSC$e_est
  #Epistasis direction...p sure this is not used
  edge_result <- factor(edges_mapk_TSC$Epistasis_Direction)
  
  #And now the Big Pig: Color by number
  #Maybe I should move this
  epi_colors <- colorRampPalette(c('#D9027D', 'grey', 'blue'))
  my_colors <- epi_colors(100)
  values <- seq(-.3, .3, length.out=100)
  color_df <- data.frame(values, my_colors)
  colors_eTSC<-my_colors[findInterval(edge_e_vals, values)]
  
  #Generate plot, yay!
  plot(mapk_graph_data, vertex.color='lightgrey', vertex.frame.color='lightgrey', vertex.size=25, 
       edge.width=edge_widths, edge.lty=edge_lty, edge.color = colors_eTSC,
       layout=layout_in_circle, 
       main=main)
}

#Generate Plots Network :) ---=-
plot_gene_network_epi(df_heatmap_TSC, main = "Epistasis in MapK Genes | Total Seed Count")
plot_gene_network_epi(df_heatmap_SN, main = "Epistasis in MapK Genes | Silique Number")
plot_gene_network_epi(df_heatmap_SPF, main = "Epistasis in MapK Genes | Seeds Per Fruit")


#Fig 4 and Fig 4S: Barplot multiplot
