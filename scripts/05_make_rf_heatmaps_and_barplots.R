#Make relative fitness heatmaps

#get packages and data
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")


#1: Relative Fitness Heatmap
make_heatmap <- function(r.df, trait = 'title') {
  df <- pivot_longer(r.df, cols = c(MA_w, MB_w, DM_w), names_to = "mutant_type", values_to = 'rel_fitness')
  df$mutant_rank <- match(df$mutant_name, mutant_order)
  
  p <- ggplot(df, aes(mutant_type, mutant_rank, fill = rel_fitness)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "#D9027D", midpoint = 1) +
    geom_tile() +
    scale_x_discrete(labels = c('geneA', 'geneB', 'double.mutant')) +
    xlab('Mutant Type') +
    labs(fill = "Relative Fitness") +
    scale_y_continuous(name = "Gene Pair",
                       breaks = 1:length(mutant_order), 
                       labels = rev(mutant_order),
                       expand = c(.01, .01)) +
    ggtitle(trait) +
    theme_bw() +
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(colour = "black"),
            axis.text.y.left = element_text(size = 6, colour = "black"),
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.text.x.bottom = element_text(size = 10, colour = "black"),
            axis.title.x = element_text(size = 12, colour = "black"))
  
  return(p)
}


p.heat.tsc <- make_heatmap(r.df.tsc, trait="Total Seed Number")
p.heat.sn <- make_heatmap(r.df.sn, trait="Silique Number")
p.heat.spf <- make_heatmap(r.df.spf, trait="Seeds per Fruit")
p.heat.dtb <- make_heatmap(r.df.dtb, trait="Days to Bolt")
p.heat.ln <- make_heatmap(r.df.ln, trait="Leaf Number")


p.heat.tsc
p.heat.sn
p.heat.spf
p.heat.dtb
p.heat.ln
