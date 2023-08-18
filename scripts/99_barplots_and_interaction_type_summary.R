source('scripts/load_required_packages.R')
load("rdata/03_workspace.RData")

df = r.df.all_results
#Add six columns: AvsWT, BvsWT, ABvsWT, ABvsA, ABvsB, BvsA 
#For each column and each row it should tell you if the fitness is higher (1), lower (-1), or nonsignificant(0)
#Define function for getting fitness comparison value:
get_effect = function(X_lci, X_uci, Y_lci, Y_uci){
  if (X_uci < Y_lci) {
    dif=-1}
  else if (X_lci > Y_uci) {
    dif=1}
  else if ((X_uci >= Y_lci) & (X_lci <= Y_uci)){
    dif = 0}
  else{dif=NA
  print("invalid input?")}
  return(dif)
}
##It turns out you need to "vectorize" udfs for them to work with mutate. Who knew?
get_effect <- Vectorize(get_effect)

#Add comparison columns
df <- df %>% mutate(A_vs_WT = get_effect(MA_w_lci, MA_w_uci, WT_w_lci, WT_w_uci))
df <- df %>% mutate(B_vs_WT = get_effect(MB_w_lci, MB_w_uci, WT_w_lci, WT_w_uci))
df <- df %>% mutate(AB_vs_WT = get_effect(DM_w_lci, DM_w_uci, WT_w_lci, WT_w_uci))
df <- df %>% mutate(AB_vs_A = get_effect(DM_w_lci, DM_w_uci, MA_w_lci, MA_w_uci))
df <- df %>% mutate(AB_vs_B = get_effect(DM_w_lci, DM_w_uci, MB_w_lci, MB_w_uci))
df <- df %>% mutate(A_vs_B = get_effect(MA_w_lci, MA_w_uci, MB_w_lci, MB_w_uci))

#Make effect comparison dataframe
effect_comp_df <- df[, c(1, (ncol(df) - 9):ncol(df))]

effect_comp_df <- effect_comp_df %>% arrange(trait, A_vs_WT, B_vs_WT, AB_vs_WT, AB_vs_A, AB_vs_B, A_vs_B)
head(effect_comp_df)

#Shows number of effect combinations
effect_comp_df %>%
  select(3:8) %>%
  distinct() %>%
  nrow()

##Top groups by number rows
top_combinations <- effect_comp_df %>%
  group_by_at(vars(3:8)) %>%
  summarise(num_rows = n()) %>%
  arrange(desc(num_rows)) %>%
  head(20)


print(top_combinations)

##Top groups by proportion
total_rows <- nrow(effect_comp_df)

top_combinations <- effect_comp_df %>%
  group_by_at(vars(3:8)) %>%
  summarise(num_rows = n()) %>%
  arrange(desc(num_rows)) %>%
  head(20) %>%
  mutate(proportion = num_rows / total_rows)

print(top_combinations)
#!!! More than 2/3 showed no effects!

traitlist = c('total_seed_count', 'silique_number', 'seeds_per_fruit', 'leaf_number', 'days_to_bolt')
comp_order = c('AB_vs_WT', 'A_vs_WT', 'B_vs_WT', 'AB_vs_A', 'AB_vs_B', 'A_vs_B')
#Visualize data so we can see what the patterns are

make_heatmap2 <- function(df){
  hm_df = df %>% pivot_longer(
    cols = 6:11,
    names_to = "comparison",
    values_to = "effect",
  )

  plot = ggplot(hm_df, aes(x=factor(comparison, levels = comp_order), 
                           y=factor(Set, levels = rev(set.order)), 
                           fill=factor(effect),
                           alpha = Epistasis_Direction,
                           color = )) +
    
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust=.6)) + 
    scale_fill_discrete(type = c('magenta', 'grey', 'blue'))
  plot = plot + facet_wrap(. ~ factor(trait, levels=traitlist), nrow=1)
  return(plot)
}

make_heatmap2(effect_comp_df)

