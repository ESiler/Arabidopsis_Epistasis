#load libraries and data
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")

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



#plot_epi_forest(test_results2, main="Test/Sample Plot")

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



   