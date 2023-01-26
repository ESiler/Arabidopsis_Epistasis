library(dplyr, tidyr)

#Get and format common names for gene loci in Huan's list
gene_names_from_huan = read.delim("data/Arabidopsis_locus-namelist.txt", sep = "|", header=T)

gene_names_from_huan$gene1 <- gsub(",.*$", "", gene_names_from_huan$CommonName)
gene_names_from_huan$genes <- tolower(gene_names_from_huan$CommonName)

for (i in 1:(nrow(gene_names_from_huan))){
  if (gene_names_from_huan$locus[i] != gene_names_from_huan$CommonName[i]) {
    gene_names_from_huan$gene1[i] <- tolower(gene_names_from_huan$gene1[i])
  }
}

#genenamekey: load data and change loci to uppercase
genenamekey = read.delim("data/genes_in_each_set.txt", sep = "\t", header = T)

genenamekey$MA <- toupper(genenamekey$MA)
genenamekey$MB <- toupper(genenamekey$MB)

#merge common names into genenamekey
genenamekey['ma'] <- left_join(genenamekey, gene_names_from_huan, by=c('MA' = 'locus'))$gene1
genenamekey['mb'] <- left_join(genenamekey, gene_names_from_huan, by=c('MB' = 'locus'))$gene1

genenamekey['ma2'] <- left_join(genenamekey, gene_names_from_huan, by=c('MA' = 'locus'))$genes
genenamekey['mb2'] <- left_join(genenamekey, gene_names_from_huan, by=c('MB' = 'locus'))$genes

genenamekey$mutant_name <- paste(genenamekey$ma, genenamekey$mb, sep='/')

head(genenamekey)
#write.csv(genenamekey, 'data/genenamekey.csv') don't accidentally overwrite