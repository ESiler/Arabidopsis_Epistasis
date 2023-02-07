#Make relative fitness heatmaps
#T suggests using complex heatmaps package
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

#function that converts desired data to long format:
#YOU ARE HERE ----
get_heat_data <- function(rdf){
  data_hm <- cbind(rdf[c(1,2,24,25)], stack(rdf[7:10]))
  data_hm$Set <- as.factor(data_hm$Set)
  data_hm <- data_hm %>% filter(ind != 'WT_w') 
  return(data_hm)
}



#Function that makes a ggplot of given data:
plot_heatmap <- function(data, trait="trait"){
  #ggplot blah blah blah
  ggplot(data, aes(ind, Set, fill=values)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "#D9027D", midpoint = 1) +
    geom_tile() +
    scale_x_discrete(labels=c('geneA', 'geneB', 'double.mutant')) +
    scale_y_discrete(limits = rev(factor(set.order))) +
    xlab('Mutation') +
    ylab('Set') +
    labs(fill = "Relative Fitness") +
    ggtitle('Relative Fitness', subtitle=trait)
}


p.heat.tsc <- r.df.tsc %>% get_heat_data() %>% plot_heatmap(trait="Total Seed Number")
p.heat.sn <- r.df.sn %>% get_heat_data() %>% plot_heatmap(trait="Silique Number")
p.heat.spf <- r.df.spf %>% get_heat_data() %>% plot_heatmap(trait="Seeds per Fruit")
p.heat.dtb <- r.df.dtb %>% get_heat_data() %>% plot_heatmap(trait="Days to Bolt")
p.heat.ln <- r.df.ln %>% get_heat_data() %>% plot_heatmap(trait="Leaf Number")


p.heat.tsc
p.heat.sn
p.heat.spf
p.heat.dtb
p.heat.ln

comheatlegend = get_legend(p.heat.tsc, position='right')

p.heat.5traits <- ggarrange(p.heat.tsc, p.heat.sn, p.heat.spf, p.heat.ln, p.heat.dtb,
                         nrow=1, 
                         legend='right',
                         legend.grob=comheatlegend)


p.heat.5traits
#TODO
#add significance (boxes?)
#gene name version

#NOTES
#values = c("#D9027D", 'black', 'blue'))