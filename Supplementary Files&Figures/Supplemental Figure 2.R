#Supplementary Figure 2
# Seperate family tree into main genera and into sub-genera clusters based on OGT calculations
library(ape) # read.tree()
library(phytools) # midpoint.root()
library(BiocManager) # to install the package "treeio", BiocManager::install("treeio)
library(treeio) # rename_taxa()
library(DECIPHER) # ReadDendrogram() , BiocManager::install("DECIPHER)
library(ggtree) # BiocManager::install("ggtree", force = T)
library(ggpubr) # stat_cor

#-----------------------------------------------------------------------load files------------------------------------------------------------------------

# load growth rates, OGT, and OGT parameters of Colwelliaceae strains (my strains)
GR = read.table("./Refs& Raw Input/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character") # my growthrates
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character") # OGT of my strains and those that I literature have growthrates

tree.DF.temp = read.table("./Phylogenetics/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")


#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae"),clustCol = c("seagreen3","olivedrab4","darkorange2","purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)


#-----------------------------------------------------------------------set up files----------------------------------------------------------------------- 
#set up column names and format to numeric 
colnames(OGT)[1] = "strain"
OGT$OGT = as.numeric(OGT$OGT)
GR = GR[,c(1,3,4)]
GR$mumax = as.numeric(GR$mumax)
GR$strain<-sub("^","Colwellia ",GR$strain)

#--------------------------------Calculate average and standard devation of growth rate and OGT by Main Clade and Sub-Clade---------------------------
# merge growth rate data, OGT data, temperature, and main cluster ID
GR_OGT_Clade = merge.data.frame(GR,OGT,by = "strain", all = F)
GR_OGT_Clade$`Temp ClusterName` = tree.DF.temp$ClusterName[match(GR_OGT_Clade$strain,tree.DF.temp$label)]
#GR_OGT_Clade$`Main ClusterName` = tree.DF.main$ClusterName[match(GR_OGT_Clade$strain,tree.DF.main$label)]
GR_OGT_Clade$Type = NULL
GR_OGT_Clade = GR_OGT_Clade[!is.na(GR_OGT_Clade$`Temp ClusterName`),]

# add the colors to the dataframe and filter it by OGTs under 25 degrees C and set up facet names
GR_OGT_Clade$clustCol = Temp.Clust.colors$clustCol[match(GR_OGT_Clade$`Temp ClusterName`,Temp.Clust.colors$Clade)]
#GR_OGT_Clade$temp = as.numeric(GR_OGT_Clade$temp)
Temp.labs = c("(A) Growth at -1?C","(B) Growth at 4?C","(C) Growth at 11?C","(D) Growth at 17?C")
names(Temp.labs) = c("-1","4","11","17")


ggplot(GR_OGT_Clade,aes(x = OGT,y = mumax)) +
  geom_point(aes(color = `Temp ClusterName`, shape = `Temp ClusterName`)) +
  stat_smooth(method = "lm", se = F,aes(group = 1), colour = "black")+
  #stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")),label.x = 6,size = 4) +
  stat_cor(method = "pearson",aes(label = paste(..rr.label..)),label.x = 6,size = 4, label.y = 3.49) +
  stat_cor(method = "pearson",aes(label = paste(..r.label..)),label.x = 6,size = 4, label.y = 3.34) +
  stat_cor(method = "pearson",aes(label = paste(..p.label..)),label.x = 6,size = 4, label.y = 3.19) +
  labs(x = "Optimal Growth Temperature (?C)", y = expression("Growth Rate" ~ (d^{-1})))+
  scale_color_manual("Taxonomic Clades",values = c("seagreen3","olivedrab4","darkorange2","purple","black"), labels = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae")) +
  scale_shape_manual("Taxonomic Clades",values = c(15,16,17,18,8,9),labels = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae")) +
  facet_grid(.~temp,scales = "free",labeller = labeller(temp = Temp.labs)) +
  theme_classic(base_line_size = 1, base_size = 16, base_rect_size = 1) +
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1),
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 13),
        axis.text = element_text(size = 10), axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10)) +
  ylim(-0.5,3.5)

