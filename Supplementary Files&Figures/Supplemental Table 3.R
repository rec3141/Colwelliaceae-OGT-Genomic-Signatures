# Combining to Make Supplementary Table 1

OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table
MAD = read.table("./Supplementary Files&Figures/OGT MAD.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table
growth = read.table("./Supplementary Files&Figures/TempGrowthrates.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table

# combine all three files 
MAD$strn. = gsub('.{1}$','',MAD$strn.)
colnames(MAD)[1] = "strn"
colnames(growth)[1] = "strn"
OGT$Type = NULL

df_list = list(growth,OGT,MAD)
total = df_list %>% reduce(full_join, by='strn')
