#Epistasis_VS_gene_traits
source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")

gene_trait_df = read.delim("data/set_gene_info_df.txt", sep = ",", header = T)

str(gene_trait_df)
str(r.df.all_results)

df = merge(r.df.all_results, gene_trait_df, by.x = 'Set', by.y = 'Set_2_A', all.x = TRUE)

#str(df)

epi_vs_famsize_plt = ggplot(data=df, aes(x=log10(gene_fam_size_C), y=e_est, color=trait)) + 
  geom_point(alpha=.33, size=7) + 
  xlab('log Gene family size') +
  theme_classic() +
  ylab("Epistasis Value (e)")


epi_vs_medexp_plt = ggplot(data=df, aes(x=log10(med_exp_ratio_A.B), y=e_est, color=trait)) + 
  geom_point(alpha=.33, size=7) + 
  xlab('Gene expression ratio (log of median expression level)') +
  theme_classic() +
  ylab("Epistasis Value (e)") 


epi_vs_maxexp_plt = ggplot(data=df, aes(x=log10(max_exp_ratio_A.B), y=e_est, color=trait)) + 
  geom_point(alpha=.33, size=7) + 
  xlab('Gene expression ratio (log of maximum expression level)') +
  ylab("Epistasis Value (e)") +
  theme_classic()

epi_vs_kssim_plt = ggplot(data=df, aes(x=sim1_paralog_ks_C, y=e_est, color=trait)) + 
  geom_point(alpha=.33, size=7) + 
  xlab('Paralog similarity (ks)') +
  ylab("Epistasis Value (e)") +
  theme_classic()

epi_vs_kaks_plt = ggplot(data=df, aes(x=ka_ks_evolve_rate_B, y=e_est, color=trait)) + 
  geom_point(alpha=.33, size=7) + 
  xlab('Selection pressure (ka/ks)') +
  ylab("Epistasis Value (e)") +
  theme_classic()


epi_vs_famsize_plt
epi_vs_medexp_plt
epi_vs_maxexp_plt
epi_vs_kssim_plt
epi_vs_kaks_plt
