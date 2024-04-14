# Plot Figure 3

# load growth rates, OGT, and OGT parameters of Colwelliaceae strains (my strains)
GR = read.table("./Refs& Raw Input/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character") # my growthrates
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character") # OGT of my strains and those that I literature have growthrates
colnames(OGT)[1] = "strain"

#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae"),clustCol = c("seagreen3","olivedrab4","darkorange2","purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)

#set up cluster name & IDS
tree.DF.main = read.table("./Phylogenetics/Genera Clade ID.tsv", sep = "\t", header = T)
tree.DF.temp = read.table("./Phylogenetics/Sub-Clade ID.tsv",sep = "\t", header = T)

# merge temp ID with OGT and clade colors and make the strains rownames
TC_df = merge(OGT, tree.DF.temp, by.x = "strain", by.y = "label", all.y= TRUE)
TC_df$`SubColor` = Temp.Clust.colors$clustCol[match(TC_df$ClusterName,Temp.Clust.colors$Clade)]
rownames(TC_df) = TC_df$strain
TC_df = as.matrix(TC_df)

TC_df = as.data.frame(TC_df)
TC_df$SubColor = as.character(TC_df$SubColor)

# edit the "Type" in the dataframe so that its prettier as a legend
TC_df$Type = gsub("lit_ratk","OGT derived from Ratkowsky fitting\n of reported temperature growth data",TC_df$Type)
TC_df$Type = gsub("lit","OGT in Reported Literature",TC_df$Type)
TC_df$Type = gsub("ratk","OGT from this Study",TC_df$Type)
colnames(TC_df)[3] = "OGT Source"


#-------------------------plot and save the figure - my strains and literature-----------------------------------------------
tiff("./Figures/Figure 3.tiff",width = 8,height = 11, units = "in", res = 1200)
ggplot(TC_df,aes(x = `ClusterName`,y = as.numeric(OGT))) + 
  geom_boxplot(lwd = 1, color = "grey") +
  geom_jitter(color = TC_df$SubColor,alpha = 0.9, size = 2.5, aes(shape = `OGT Source`))+
  labs(x = "Taxonomic Clades", y = "Optimal Growth Temperature (°C)") +
  theme_classic(base_line_size = 1, base_size = 16) + 
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1),axis.text = element_text(size = 10),
        axis.title = element_text(size =  14),legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 13))
dev.off()
