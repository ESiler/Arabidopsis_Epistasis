#load libraries and data
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")

# should be a total of 8 figs as currently requested
#tsc, ln, dtb, and spf with and without edge effects. 
#But exps with flat effect need to be modeled separately 
max <- max(r.df.all_results$upperCI)
min <- min(r.df.all_results$lowerCI)
nsets <- nrow(r.df.tsc)
forbreaks <- 1:nsets

no_y_axis_theme <- function() {
  theme(axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank() 
    
  )
}
# Make forest plots ----

# plot_epi_forest: This function makes a epistasis forest plots for all the genes. 
# Returns a ggplot plot
plot_epi_forest <- function(epi_data, main="Title") {
  
  epi_data$mutant_rank = match(epi_data$mutant_name, mutant_order)
  
  plot <- ggplot(epi_data, aes(y = mutant_rank, x = e_est, color=Epistasis_Direction, ymin=1, ymax=dim(epi_data)[1])) +
    geom_point(shape = 18, size = 3) +  
    geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.5) +
    geom_vline(xintercept = 0, color = "black", cex = .5) +
    ggtitle(main) +
    xlab("Epistasis Value [95% CI]") +
    xlim(-1, 0.6) +
    scale_y_continuous(name="Gene Pair",
                       breaks = 1:length(mutant_order), 
                       labels = rev(mutant_order),
                       expand=c(.01,.01))+
    scale_color_manual('Epistasis\nDirection',values = c('blue', 'black', "#D9027D")) +
    
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),
          axis.text.y.left = element_text(size = 5, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x.bottom = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
  plot
}

fp.tsc <- plot_epi_forest(r.df.tsc, main="Total Seed Count")
fp.dtb <- plot_epi_forest(r.df.dtb, main="Days to Bolt")
fp.ln <- plot_epi_forest(r.df.ln, main="Leaf Number")
fp.sn <- plot_epi_forest(r.df.sn, main="Silique Number")
fp.spf <- plot_epi_forest(r.df.spf, main="Seed per Fruit")

big_forest_plot <- ggarrange(fp.tsc, 
          fp.sn + no_y_axis_theme(),
          fp.spf + no_y_axis_theme(), 
          fp.dtb + no_y_axis_theme(), 
          fp.ln + no_y_axis_theme(), 
          widths = c(1.7,1,1,1,1),
          nrow = 1, 
          legend='right',
          common.legend=T)

annotate_figure(big_forest_plot, 
                top = text_grob("Epistasis", size=16))
fp.tsc
fp.dtb
fp.ln
fp.sn
fp.spf
big_forest_plot
