#Make relative fitness heatmaps

#get packages and data
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")


#1: Relative Fitness Heatmap
head(r.df.tsc)
#Get results in long format
data_hm_tmp <- cbind(r.df.tsc[c(1,2,24,25)], stack(r.df.tsc[7:10]))
data_hm_tmp$Set <- as.factor(data_hm_tmp$Set)
data_hm_tmp <- data_hm_tmp %>% filter(ind != 'WT_w')
head(data_hm_tmp)

ggplot(data_hm_tmp, aes(ind, mutant_name, fill=values)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "#D9027D", midpoint = 1) +
  geom_tile() +
  scale_x_discrete(labels=c('geneA', 'geneB', 'double.mutant')) +
  scale_y_discrete(limits = rev(names.order)) +
  xlab('Mutation') +
  ylab('Gene Pair') +
  labs(fill = "Relative Fitness") +
  ggtitle('Relative Fitness', subtitle='Total Seed Count')

#add significance (boxes?)
#set version and gene name version
#add fruit 
#values = c("#D9027D", 'black', 'blue'))