# Script to get epistasis results from Arabidopsis fitness data
# on duplicated homologous genes. 
# By E Siler

#1. Load processed data, formulas, libraries, and other objects:
source('scripts/load_required_packages.R')
load("rdata/02_workspace.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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



#GET PLANT SET ORDER -- IMPORTANT
set.order <- r.df.tsc$Set
names.order <- r.df.tsc$mutant_name


