library(ggpubr) # stat_cor() & ggplot()
library(reshape2) # melt()
library(RColorBrewer) # brewer.pal()
library(readr)
library(dplyr)

# ------------------------------load files and set them up -------------------------------------------------------------
# read amino acid data and phylogenetic data
genomes = read.table("./Patric Input/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
WGAA = read.table("./Amino Acid Indices/WG Indice Data.tsv",fill = T,head = T, sep = "\t",colClasses = 'character', row.names = 1)
Main.Clust.ID = read.table("./Phylogenetics/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

#read optimal growth temperature data
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table

# set up matching column names
colnames(Main.Clust.ID)[1] = "strain"
colnames(OGT)[1] = "strain"

#-------------------------------------set up color schemes ------------------------------------------------------------
# set colors for genera clusters
#Main.Colors = c("Clade A" = "darkgreen","Clade B" = "orange","Clade C" = "purple","Other Colwelliaceae" = "black")
Main.Colors = c("Clade A1" = "seagreen3","Clade A2" = "olivedrab4","Clade B" = "darkorange2","Clade C" = "purple","Other Colwelliaceae" = "black")

#---------------------------------------set up main dataframe --------------------------------------------------------
#merge to have the OGT with Genera CLuster ID and add PATRIC ID
OGT.Clust = merge.data.frame(Main.Clust.ID,OGT, by = "strain", all = F)
OGT.Clust$ID = genomes$genome.genome_id[match(OGT.Clust$strain,genomes$genome.genome_name)]
rownames(OGT.Clust) = OGT.Clust$ID
OGT.Clust$ID = NULL

# merge dataframe with genome amino acid indices and make column numeric
WG_Data = merge.data.frame(OGT.Clust,WGAA,by = "row.names",all = F)
WG_Data$Row.names = NULL
cols.num = colnames(WG_Data)[c(3,5,6,7,8,9,10)]
WG_Data[cols.num] = sapply(WG_Data[cols.num],as.numeric)


#-----------------------------------------Plot OGT vs. Genome Amino Acid indices------------------------------------------------------------------
#melt dataframes
WG_Data_melt = melt(WG_Data,id.vars = c("strain","Type","ClusterName","OGT"))
WG_Data_melt$variable = as.character(WG_Data_melt$variable)

#convert the variable names to proper format
WG_Data_melt$variable[WG_Data_melt$variable == "Arg.Lys.Ratio"] = "Arginine-Lysine Ratio"
WG_Data_melt$variable[WG_Data_melt$variable == "Aromaticity.Index"] = "Aromaticity Index"
WG_Data_melt$variable[WG_Data_melt$variable == "Aliphatic.Index"] = "Aliphatic Index"
WG_Data_melt$variable[WG_Data_melt$variable == "Proline.Residue.Proportion"] = "Proline Residue Proportion"
WG_Data_melt$variable[WG_Data_melt$variable == "Acidic.Residue.Proportion"] = "Acidic Residue Proportion"

WG_Data_melt$variable = factor(WG_Data_melt$variable, levels = c("Acidic Residue Proportion", "Aliphatic Index", "Arginine-Lysine Ratio", "Aromaticity Index", "GRAVY", "Proline Residue Proportion"), 
                               labels = c("(A)~~~Acidic~Residue~Proportion", "(B)~~~Aliphatic~Index", "(C)~~~Arginine-Lysine~Ratio", "(D)~~~Aromaticity~Index", "(E)~~~GRAVY", "(F)~~~Proline~Residue~Proportion"))
WG_Data_melt$value = as.numeric(WG_Data_melt$value)
WG_Data_melt$ClusterName[WG_Data_melt$ClusterName == "Clade A1"] = "Clade A"
WG_Data_melt$ClusterName[WG_Data_melt$ClusterName == "Clade A2"] = "Clade A"

# edit the "Type" in the dataframe so that its prettier as a legend
WG_Data_melt$Type = gsub("lit_ratk","OGT derived from Ratkowsky fitting of reported temperature growth data",WG_Data_melt$Type)
WG_Data_melt$Type = gsub("lit","OGT in Reported Literature",WG_Data_melt$Type)
WG_Data_melt$Type = gsub("ratk","OGT from this Study",WG_Data_melt$Type)
colnames(WG_Data_melt)[2] = "OGT Source"



tiff("./Figures/Figure 4.tiff",width = 22,height = 10, units = "in", res = 500)
ggplot(WG_Data_melt,aes(x = OGT,y = value)) + 
  geom_point(aes(color = ClusterName,shape = `OGT Source`),size = 2)+ 
  geom_smooth(method = "lm", se = F,aes(group = 1)) + 
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,if_else(readr::parse_number(..p.label..) <0.01,"p<0.01","p<0.05"), sep = "~`,`~")),label.x = 15)+
  labs(x = "Optimal Growth Temperature (°C)", y = "Amino Acid Index")+
  scale_color_manual("Taxanomic Clades",values = c("darkgreen","orange","purple","black"),labels = c("Clade A","Clade B","Clade C","Other Colwelliaceae")) + 
  scale_shape_manual("OGT Sources",values = c(16,17,15),labels = c("OGT derived from Ratkowsky fitting of reported temperature growth data","OGT in Reported Literature","OGT from this Study")) +
  facet_wrap(~variable,scales = "free" ,labeller = label_parsed) +
  scale_y_continuous() + 
  theme_minimal(base_line_size = 1, base_size = 16, base_rect_size = 1) + 
  theme(panel.border = element_rect(colour = "black", size = 1, fill = NA),
        legend.position = "bottom",
        legend.text = element_text(size = 13),legend.title = element_text(size = 15),
        axis.text = element_text(size = 15), axis.title = element_text(size = 17),
        strip.text.x = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )
dev.off()

