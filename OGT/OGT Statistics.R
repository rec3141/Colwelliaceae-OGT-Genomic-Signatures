#OGT Statistics

# load growth rates, OGT, and OGT parameters of Colwelliaceae strains (my strains)
GR = read.table("./Refs& Raw Input/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character") # my growthrates
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character") # OGT of my strains and those that I literature have growthrates

#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae"),clustCol = c("seagreen3","olivedrab4","darkorange2","purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)

#set up cluster name & IDS
tree.DF.main = read.table("./Phylogenetics/Genera Clade ID.tsv", sep = "\t", header = T)
tree.DF.temp = read.table("./Phylogenetics/Sub-Clade ID.tsv",sep = "\t", header = T)

#set up column names and format to numeric 
colnames(OGT)[1] = "strain"
OGT$OGT = as.numeric(OGT$OGT)


# merge OGT data, and cluster IDs
OGT_Clade = OGT
OGT_Clade$`Temp ClusterName` = tree.DF.temp$ClusterName[match(OGT_Clade$strain,tree.DF.temp$label)]
OGT_Clade$`Main ClusterName` = tree.DF.main$ClusterName[match(OGT_Clade$strain,tree.DF.main$label)]
OGT_Clade$Type = NULL
OGT_Clade = OGT_Clade[!is.na(OGT_Clade$`Temp ClusterName`),]

# add Clade marking in GR
GR$`Temp ClusterName` = tree.DF.temp$ClusterName[match(GR$strain,tree.DF.temp$label)]

# calculate by Temperature Clade the OGT
OGT_TC = OGT_Clade %>% dplyr::group_by(`Temp ClusterName`) %>% dplyr::summarize(`Average OGT` = paste0(signif(mean(OGT),2)," ? ",signif(sd(OGT),2)),
                                                                                `Min OGT` = signif(min(OGT)),`Max OGT` = signif(max(OGT)), N = length(unique(strain)))
# calculate by Main Clade the OGT
OGT_MC = OGT_Clade %>% dplyr::group_by(`Main ClusterName`) %>% dplyr::summarize(`Average OGT` = paste0(signif(mean(OGT),2)," ? ",signif(sd(OGT),2)),
                                                                                `Min OGT` = signif(min(OGT)),`Max OGT` = signif(max(OGT)), N = length(unique(strain)))

# calculate by Temperature Clade the GR
GR$mumax = as.numeric(GR$mumax)
GR_TC = GR %>% dplyr::group_by(`Temp ClusterName`) %>% dplyr::summarize(`Average GR` = paste0(signif(mean(mumax),2)," ? ",signif(sd(mumax),2)),
                                                                                `Min GR` = signif(min(mumax)),`Max GR` = signif(max(mumax)), N = length(unique(strain)))


# ANOVA of OGT between Main Clades
OGT_MC_aov = aov(OGT~`Main ClusterName`,data = OGT_Clade)
summary(OGT_MC_aov)

# ANOVA of OGT between Sub Clades
OGT_TC_aov = aov(OGT~`Temp ClusterName`,data = OGT_Clade)
summary(OGT_TC_aov)

model = lm(OGT_Clade$OGT~OGT_Clade$`Temp ClusterName`)
ANOVA =aov(model) # anova with all groups
summary(ANOVA)
TUKEY = TukeyHSD(x =ANOVA,"OGT_Clade$`Temp ClusterName`", conf.level = 0.95)
plot(TUKEY,las =1,col="brown")









