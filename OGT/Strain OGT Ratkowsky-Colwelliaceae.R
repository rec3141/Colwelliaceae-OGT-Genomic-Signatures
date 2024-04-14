library(readxl) # read XLS files
library(dplyr)

# ----------------------------------------------------load files -----------------------------------------------------------------------
lit_col = readxl::read_excel("./Refs& Raw Input/Literature OGT.xlsx", sheet = "Literature OGT") # other OGT data from the literature
init.save.topt = read.table("./OGT/Strain Ratk OGT Outputs Raw.tsv",fill=T,head=T,sep="\t")
biospec = read.csv("./Refs& Raw Input/Literature Colwellia Growthrates.csv")

genus = "Colwellia"

#----------------------------------------------- format final OGT table --------------------------------------------------------------
genome.list = as.character(unique(biospec$binomial.name[grepl(genus,biospec$binomial.name)]))

#format data table
OGT = init.save.topt
OGT$j = NULL
cols.num_OGT = colnames(OGT)[2:8]
OGT[cols.num_OGT] = sapply(OGT[cols.num_OGT],as.character)
OGT[cols.num_OGT] = sapply(OGT[cols.num_OGT],as.numeric)

#subset out literature strains and take the weighted OGT of strain as their final OGT
genome.list = paste(genome.list,collapse = "|")
ex_topt = subset.data.frame(OGT,grepl(genome.list,strn))
ex_topt = ex_topt %>% group_by(strn) %>% dplyr::summarize(OGT = median(Weight_OGT))
ex_topt = as.data.frame(ex_topt)
ex_topt$Type = "lit_ratk"

OGT = OGT %>% filter(!grepl(genome.list,strn))
OGT = as.data.frame(OGT)

# filter monte carlo outputs by RSME and take the median of the OGT by strain
OGT$Weight_OGT = NULL
OGT_params = subset(OGT,OGT$`MC.RMSE` < exp(-30))
Final_OGT = OGT_params %>% group_by(strn) %>% dplyr::summarize(OGT = median(MC_OGT))
MAD_OGT = OGT_params %>% group_by(strn) %>% dplyr::summarize(MAD = mad(MC_OGT))

Final_OGT$Type = "ratk"
Final_OGT = rbind(Final_OGT,ex_topt)

#------------------------------------------------- Format strain names-------------------------------------------------------------------
# format dataframe
Final_OGT = as.data.frame(Final_OGT)
Final_OGT2 = Final_OGT
Final_OGT2$strn = as.character(Final_OGT2$strn)

# set strain names of literature based growth rate Colwellia to actually strain names (PATRIC format)
Final_OGT2$strn = gsub("Colwellia psychrerythraea","Colwellia psychrerythraea ACAM 605",Final_OGT2$strn) 
Final_OGT2$strn = gsub("Colwellia hornerae","Colwellia hornerae strain ACAM 607",Final_OGT2$strn)
Final_OGT2$strn = gsub("Colwellia piezophila","Colwellia piezophila ATCC BAA-637",Final_OGT2$strn) 
Final_OGT2$strn = gsub("Colwellia demingiae","Colwellia demingiae strain ACAM 459",Final_OGT2$strn) 

#-------------------------- Merge dataframe with OGTs of other Colwelliaceae found in the literature ------------------------------------
lit_col = lit_col[,c(1:2)]
colnames(lit_col) = c("strn","OGT")
lit_col$Type = "lit"

Final_OGT2 = rbind(Final_OGT2,lit_col)

#----------------------------------------write table of final Colwelliaceae OGT----------------------------------------------------------
write.table(Final_OGT2,"./OGT Figure/Strain Ratk OGT Final.tsv",sep = "\t",quote = F,row.names = F)
write.table(MAD_OGT,"./Supplementary Files&Figures/OGT MAD.tsv",sep = ".\t",quote = F,row.names = T)
