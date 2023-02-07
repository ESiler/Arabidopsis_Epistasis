## 1. Load libraries ----
source("scripts/load_required_packages.R")

## 2. Load data ----
data = read.delim("data/Fitness_data_83_sets_24Jan2023.txt", sep = "\t", header = T)


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

#In Sets, change 2r to 2222 (it is a repeat of exp 2)
level_2r <- which(levels(data$Set) == '2r')
levels(data$Set)[level_2r] <- 2222

#Load in manual gene name data.
genenamekey <- read.csv('data/genenamekey_cleaned.csv', header = T)
genenamekey <- select(genenamekey, -1)
tail(genenamekey)

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

#Save processed data :-)
rm(dataj)
saveRDS(data, file = "rdata/01_data.rds")