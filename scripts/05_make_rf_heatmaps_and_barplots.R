#Make relative fitness heatmaps

#NOTES
#values = c("#D9027D", 'black', 'blue'))
#T suggests using complex heatmaps package
#get packages and data
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")


#1: Relative Fitness Heatmap
head(r.df.tsc)
#Get results in long format

data_hm_tmp <- cbind(r.df.tsc[c(1,2,23,24,25)], stack(r.df.tsc[c(8:10)]))
data_hm_tmp <- cbind(data_hm_tmp, stack(r.df.tsc[20:22]))
data_hm_tmp$Set <- as.factor(data_hm_tmp$Set)
colnames(data_hm_tmp)[8] <- 'p-val'
colnames(data_hm_tmp)[6] <- 'rel_fitness'

head(data_hm_tmp)

#function that converts desired data to long format:
#YOU ARE HERE ----
get_heat_data <- function(rdf){
  data_hm <- cbind(rdf[c(1,2,23,24,25)], stack(rdf[8:10]))
  data_hm <- cbind(data_hm, stack(rdf[20:22]))
  data_hm$Set <- as.factor(data_hm$Set)
  colnames(data_hm)[8] <- 'pval'
  colnames(data_hm)[9] <- 'sig'
  colnames(data_hm)[6] <- 'rel_fitness'
  data_hm$sig <- (data_hm$pval < .05)
  return(data_hm)
}
heat.data.tsc <- get_heat_data(r.df.tsc)
tail(heat.data.tsc)

#Function that makes a ggplot of given data:
plot_heatmap <- function(data, trait="trait"){
  #ggplot blah blah blah
  ggplot(data, aes(ind, Set, fill=rel_fitness)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "#D9027D", midpoint = 1) +
    geom_tile() +
    # Add outline via a geom_rect
    #geom_tile(
    #  data = subset(data, sig == TRUE),
    #  aes(), fill = NA, color = "#707070", size = .5
    #) +
    scale_x_discrete(labels=c('geneA', 'geneB', 'double.mutant')) +
    scale_y_discrete(limits = rev(factor(set.order))) +
    xlab('Mutation') +
    ylab('Set') +
    labs(fill = "Relative Fitness") +
    ggtitle('Relative Fitness', subtitle=trait)
}
p.heat.tsc <- r.df.tsc %>% get_heat_data() %>% plot_heatmap(trait="Total Seed Number")
p.heat.tsc


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
#make boxes prettier
#gene name version

