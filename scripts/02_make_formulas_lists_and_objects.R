source("scripts/load_required_packages.R")

#Load data from 
data <- readRDS("rdata/01_data.rds")

#the setlist is the list of plant sets I'm looking at. 
#Should review this -- supposedly contain duplicate genes. 
setlist <- as.character(sort(as.integer(levels(data$Set))))
sets_with_flats <- as.character(sort(as.integer(as.character(unique(data$Set[data$Flat == '2'])))))
sets_without_flats <- setlist[!(setlist %in% sets_with_flats)]


## 4. Generate Model Formulas ----

#Total seed count
f.tsc <- formula(log10(TSC+1) ~ MA + MB + DM )
f.tsc.e <- formula(log10(TSC+1) ~ MA + MB + DM + Type)
f.tsc.f <- formula(log10(TSC+1) ~ MA + MB + DM + Flat)
f.tsc.ef <- formula(log10(TSC+1) ~ MA + MB + DM + Type + Flat)
tsc_formulas <- c(f.tsc, f.tsc.e, f.tsc.f, f.tsc.ef)

#Leaf number
f.ln <- formula(log10(LN+1) ~ MA + MB + DM)
f.ln.e <- formula(log10(LN+1) ~ MA + MB + DM + Type)
f.ln.f <- formula(log10(LN+1) ~ MA + MB + DM + Flat)
f.ln.ef <- formula(log10(LN+1) ~ MA + MB + DM + Type + Flat)

#Days to bolt
f.dtb <- formula(log10(DTB+1) ~ MA + MB + DM)
f.dtb.e <- formula(log10(DTB+1) ~ MA + MB + DM + Type)
f.dtb.f <- formula(log10(DTB+1) ~ MA + MB + DM + Flat)
f.dtb.ef <- formula(log10(DTB+1) ~ MA + MB + DM + Type + Flat)

#Seeds per fruit
f.spf <- formula(log10(SPF+1) ~ MA + MB + DM)
f.spf.e <- formula(log10(SPF+1) ~ MA + MB + DM + Type)
f.spf.f <- formula(log10(SPF+1) ~ MA + MB + DM + Flat)
f.spf.ef <- formula(log10(SPF+1) ~ MA + MB + DM + Type + Flat)

#Silique number
f.sn <- formula(log10(SN+1) ~ MA + MB + DM)
f.sn.e <- formula(log10(SN+1) ~ MA + MB + DM + Type)
f.sn.f <- formula(log10(SN+1) ~ MA + MB + DM + Flat)
f.sn.ef <- formula(log10(SN+1) ~ MA + MB + DM + Type + Flat)

#Make a collection of dummy dataframes to make prediction w models
## Adding "additive fitness" with effects of A and B but not epistasis
df_pred_dummy <- data.frame(MA  = c(0,1,0,1,1),
                            MB = c(0,0,1,1,1),
                            DM = c(0,0,0,1,0)
)
rownames(df_pred_dummy) = c("WT", "MA", "MB", "DM","AF")

df_pred_dummy.e <- data.frame(MA  = c(0,1,0,1,1),
                              MB = c(0,0,1,1,1),
                              DM = c(0,0,0,1,0),
                              Type = c("INSIDE","INSIDE","INSIDE","INSIDE","INSIDE")
)
rownames(df_pred_dummy.e) = c("WT", "MA", "MB", "DM", "AF")

df_pred_dummy.f <- data.frame(MA  = c(0,1,0,1,1),
                              MB = c(0,0,1,1,1),
                              DM = c(0,0,0,1,0),
                              Flat = c('1','1','1','1','1')
)

rownames(df_pred_dummy.f) = c("WT", "MA", "MB", "DM","AF")


df_pred_dummy.ef <- data.frame(MA  = c(0,1,0,1,1),
                               MB = c(0,0,1,1,1),
                               DM = c(0,0,0,1,0),
                               Type = c("INSIDE","INSIDE","INSIDE","INSIDE","INSIDE"),
                               Flat = c('1','1','1','1','1')
)
rownames(df_pred_dummy.ef) = c("WT", "MA", "MB", "DM", "AF")

save.image(file = "rdata/02_workspace.RData")

