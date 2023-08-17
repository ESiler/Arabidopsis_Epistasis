## 1. Load libraries ----
source("scripts/load_required_packages.R")

## 2. Load data ----
data = read.delim("data/FITNESS_DATA_06122023_100sets.csv", sep = ",", header = T)

## 3. Format data for analysis ----
numeric_cols <- c('LN', 'WO', "SPF", 'TSC', 'SH')
factor_cols <- c('Set', 'Flat', 'Row', 'Type', 'Genotype', 'MA', 'MB', 'Subline')

data <- data %>% mutate_at(numeric_cols, as.numeric)
data <- data %>% mutate_at(factor_cols, as.factor)

data$MA <- as.integer(data$MA == 'MUT')
data$MB <- as.integer(data$MB == 'MUT')
data$DM <- as.integer(data$Genotype == 'DM')

#I want INSIDE to be the default type (vs BORDER)
levels(data$Type) <- c('INSIDE', 'BORDER')

# create log total seed count 
data$logTSC <- log10(data$TSC)
data <- subset(data, is.finite(data$logTSC)) #remove nans

#Load in manual gene name data
##Label Duplicates
genenamekey <- read.csv('data/gene_name_key.csv', header = T)
genenamekey <- select(genenamekey, -1)
genenamekey$duplicate <- duplicated(genenamekey[, c('MA', 'MB')])
# Duplicate sets: 704  705  827 2222

#In Sets, change 2r to 2222 (it is a repeat of exp 2)
level_2r <- which(levels(data$Set) == '2r')
levels(data$Set)[level_2r] <- 2222

#

#Convert sets to integer to merge w genenamekey dataset
data$Set_int <- as.integer(as.character(data$Set))

#Merge data with gene name dataset
dataj <- left_join(data, genenamekey, by=c('Set_int' = 'Set'))

# remove extra columns. rename columns. 
data <- select(dataj, -Set_int)
data <- rename(data,locusA = MA.y)
data <- rename(data,locusB = MB.y)
data <- rename(data, MA = MA.x, MB = MB.x)

#more factors
factor_cols_2 <- c('locusA', 'locusB', 'ma', 'mb', 'ma2', 'mb2', 'mutant_name')
data <- data %>% mutate_at(factor_cols_2, as.factor)

#separate duplicate gene pairs
data = subset(data, data$duplicate==FALSE)
duplicates_data = subset(data, data$duplicate==TRUE)

data <- data %>% select(-duplicate)
data$Set <- droplevels(data$Set)
duplicates_data <- duplicates_data %>% select(-duplicate)
duplicates_data$Set <- droplevels(duplicates_data$Set)
#str(data)
rm(dataj)
#Save processed data :-) 

#uncomment to overwrite
#saveRDS(data, file = "rdata/01_data.rds")

##dev area

