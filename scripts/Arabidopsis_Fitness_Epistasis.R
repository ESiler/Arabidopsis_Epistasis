# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

#1. Load processed data, formulas, libraries, and other objects:
source('scripts/load_required_packages.R')
load("rdata/02_workspace.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 7. Get epistasis results and crank out graphs ----

# should be a total of 8 figs as currently requested
#tsc, ln, dtb, and spf with and without edge effects. 
#But exps with flat effect need to be modeled separately 


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
plot(model38norm)
str(coef(model38norm))

length(f.dtb.ef)
f.dtb.ef[[1]]
f.dtb.ef[[2]]
pred_vars <- f.dtb.ef[[3]] #WOOHOOOOOOO!

pred_vars_c = deparse(pred_vars)

pred_vars = strsplit(pred_vars_c, " . ")[[1]]


#Okay now get string in a format that can be used to fill in a table ----


#Fig 2: Multiplot of Fitness Values by Genotype ----
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

#A. Create results DFs for non-edge models

#Make epistasis and rel. fitness results dataframe combining sets with and without flats
get_rfit_results_df <- function(woflatfunc, wflatfunc, data, d1=df_pred_dummy, d2=df_pred_dummy.f){
  r.df <- get_epistasis_for_formula(sets_without_flats, woflatfunc, d1, data)
  r.df.f <- get_epistasis_for_formula(sets_with_flats, wflatfunc, d2, data)
  r.df.r <- rbind(r.df, r.df.f)
  r.df.r <- arrange(r.df.r, e_est)
  r.df.r$row <- c(1:1:dim(r.df.r)[1])
  return(r.df.r)
}

results.df.tsc <- get_rfit_results_df(woflatfunc = f.tsc, wflatfunc = f.tsc.f, data = data)
results.df.dtb <- get_rfit_results_df(woflatfunc = f.dtb, wflatfunc = f.dtb.f, data = data)
results.df.ln <- get_rfit_results_df(woflatfunc = f.ln, wflatfunc = f.ln.f, data = data)
results.df.spf <- get_rfit_results_df(woflatfunc = f.dtb, wflatfunc = f.ln.f, data = data)

#GET PLANT SET ORDER -- IMPORTANT
set.order <- r.df.tsc$Set
names.order <- r.df.tsc$mutant_name


