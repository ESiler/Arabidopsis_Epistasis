
#Make a scatterplot comparing the what the additive fitness values would be to the actual values w epistasis

source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")

r.df.3seedtraits <- filter(r.df.all_results, trait == 'total_seed_count' | trait == 'seeds_per_fruit' | trait == 'silique_number')
r.df.2vegtraits <- filter(r.df.all_results, trait == 'days_to_bolt' | trait == 'leaf_number')


ggplot(data=r.df.all_results, 
       aes(x=AF_w, y=DM_w, color=trait, shape = Epistasis_Direction, alpha = Epistasis_Direction)) + 
  geom_point(alpha=.33, size=7) + 
  xlim(.5, 1.7) +
  ylim(.5, 1.7) +
  scale_shape_manual(values=c(17,1,15)) +
  geom_abline(linetype=2, alpha=.33) +
  geom_vline(xintercept = 1, linetype=3, alpha = .33) +
  geom_hline(yintercept = 1, linetype=3, alpha = .33) +
  xlab("Predicted Double Mutant Fitness Relative to Wildtype\nIf there was No Epistasis") +
  ylab("Actual Double Mutant Fitness Relative to Wildtype") +
  theme_classic()

ggplot(data=r.df.3seedtraits, 
       aes(x=AF_w, y=DM_w, color=trait, shape = Epistasis_Direction, alpha = Epistasis_Direction)) + 
  geom_point(alpha=.33, size=7) + 
  xlim(.5, 1.7) +
  ylim(.5, 1.7) +
  scale_color_manual(values=c("#00BF7D", "#00B0F6", "#E76BF3")) +
  scale_shape_manual(values=c(17,1,15)) +
  geom_abline(linetype=2, alpha=.33) +
  geom_vline(xintercept = 1, linetype=3, alpha = .33) +
  geom_hline(yintercept = 1, linetype=3, alpha = .33) +
  xlab("Predicted Double Mutant Fitness Relative to Wildtype\nIf there was No Epistasis") +
  ylab("Actual Double Mutant Fitness Relative to Wildtype") +
  theme_classic()

ggplot(data=r.df.2vegtraits, 
       aes(x=AF_w, y=DM_w, color=trait, shape = Epistasis_Direction, alpha = Epistasis_Direction)) + 
  geom_point(alpha=.5, size=7) + 
  xlim(.9, 1.25) +
  ylim(.9, 1.25) +
  scale_color_manual(values=c("#F8766D", "#A3A500")) +
  scale_shape_manual(values=c(17,1,15)) +
  geom_abline(linetype=2, alpha=.33) +
  geom_vline(xintercept = 1, linetype=3, alpha = .33) +
  geom_hline(yintercept = 1, linetype=3, alpha = .33) +
  xlab("Predicted Double Mutant Fitness Relative to Wildtype\nIf there was No Epistasis") +
  ylab("Actual Double Mutant Fitness Relative to Wildtype") +
  theme_classic()




