# Seperate family tree into main genera and into sub-genera clusters based on OGT calculations & make phylogenetic Tree( Figure 1)
library(ape) # read.tree()
library(phytools) # midpoint.root()
library(BiocManager) # to install the package "treeio", BiocManager::install("treeio)
library(treeio) # rename_taxa()
library(DECIPHER) # ReadDendrogram() , BiocManager::install("DECIPHER)
library(ggtree) # BiocManager::install("ggtree)
library(ggpubr) # stat_cor

#-----------------------------------------------------------------------load files------------------------------------------------------------------------
# load PATRIC phylogenetic files
treefile="./Patric Input/Colwelliaceae Rooted Tree.nwk"

# load growth rates, OGT, and OGT parameters of Colwelliaceae strains (my strains)
GR = read.table("./Refs& Raw Input/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character") # my growthrates
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character") # OGT of my strains and those that I literature have growthrates

#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae"),clustCol = c("seagreen3","olivedrab4","darkorange2","purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)

# set up genome Patric Data and edit names
genomes = read.table("./Patric Input/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")

main_nodes = c(117,145,142,169) # order of shallow to deep branched cladee
temp_nodes = c(117,145,172,179)


#-----------------------------------------------------------------------set up files----------------------------------------------------------------------- 
#read tree
tre.phylo = read.tree(treefile)
tre.phylo = ape::root(tre.phylo,outgroup = "1938339.4",resolve.root = T)

# remove any genomes that had poor assembly
tre.phylo = drop.tip(tre.phylo,"1380381.3")

# subset genomes that are only used in the phylogenetic tree 
common_ids = intersect(tre.phylo$tip.label,genomes$genome.genome_id)
tre.phylo = keep.tip(tre.phylo, common_ids) # keeping the tips that are in the dataframe 

#set up tree and rename tips with strain names 
tre.phylo = midpoint.root(tre.phylo)
tre.phylo = ladderize(tre.phylo)
tre.phylo = rename_taxa(tre.phylo,genomes,genome.genome_id,genome.genome_name) # rename the ends. 

# convert to dendrogram
dend = tre.phylo
dend$node.label = NULL
dend <- ReadDendrogram(textConnection(write.tree(dend)))

#set up column names and format to numeric 
colnames(OGT)[1] = "strain"
OGT$OGT = as.numeric(OGT$OGT)
GR = GR[,c(1,3,4)]
GR$mumax = as.numeric(GR$mumax)

#----------------------------------------------------------------------set up tree splits and plot--------------------------------------------------------
#set up clusters for the main groups
tre.phylo.main = groupClade(tre.phylo,main_nodes)
tree.DF.main = as.data.frame(as_tibble(tre.phylo.main))
tree.DF.main = tree.DF.main[1:length(tre.phylo.main$tip.label),] # get just the groups associated with the tips

#set up clusters for the temperature phenotype groups
tre.phylo.temp = groupClade(tre.phylo,temp_nodes)
tree.DF.temp = as.data.frame(as_tibble(tre.phylo.temp))
tree.DF.temp = tree.DF.temp[1:length(tre.phylo.temp$tip.label),]

#--------------------------------------------------------------------set up ID files for each cluster type ------------------------------------------------
# set up the names of main clusters based on the group number
Main_ClusterID = data.frame(group = 0:4,ClusterName = c("Other Colwelliaceae","Clade A","Clade B","Other Colwelliaceae","Clade C"))
tree.DF.main = merge.data.frame(tree.DF.main,Main_ClusterID, by = "group", all = T)
tree.DF.main = tree.DF.main[,5:6]

# set up the names of temperature clusters based on the group number
Temp_ClusterID = data.frame(group = 0:4,ClusterName = c("Other Colwelliaceae","Clade C","Clade B","Clade A2","Clade A1"))
tree.DF.temp = merge.data.frame(tree.DF.temp,Temp_ClusterID, by = "group", all = T)
tree.DF.temp = tree.DF.temp[,5:6]

#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade B","Clade C","Other Colwelliaceae"),
                               clustCol = c("seagreen3","olivedrab4","darkorange2","purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)

# --------------------------------------------- write files -----------------------------------------------------------------------------
write.table(tree.DF.main,"./Phylogenetics/Genera Clade ID.tsv",sep = "\t",quote = F,row.names = F)
write.table(tree.DF.temp,"./Phylogenetics/Sub-Clade ID.tsv",sep = "\t",quote = F,row.names = F)

#---------------------------------------------------------------- Plot Phylogenetic Tree -------------------------------------------------------------
# merge temp ID with OGT and clade colors and make the strains rownames
TC_df = merge(OGT, tree.DF.temp, by.x = "strain", by.y = "label", all.y= TRUE)
TC_df$`SubColor` = Temp.Clust.colors$clustCol[match(TC_df$ClusterName,Temp.Clust.colors$Clade)]
rownames(TC_df) = TC_df$strain
TC_df = as.matrix(TC_df)

#set the strain names in the order of the phylogenetic tree
TC_df = TC_df[tre.phylo$tip.label,]
TC_df = as.data.frame(TC_df)
TC_df$SubColor = as.character(TC_df$SubColor)

#set Sub-Clade colors as a factor list
SubClade_color = setNames(TC_df$SubColor,TC_df$ClusterName)

# create a dataframe for colors based on bootstrap values 
nodeColor = data.frame(node = tre.phylo$node.label,color = NA)
nodeColor$node = as.character(nodeColor$node)
nodeColor$nodenum = as.numeric(nodeColor$node)
for (i in 1:nrow(nodeColor)){
  if((is.na(nodeColor$nodenum[i]) & nodeColor$node[i] == "Root") == T){
    nodeColor$color[i] = "white"
  }else if(nodeColor$node[i] == "100"){
    nodeColor$color[i] = "dark grey"
  } else if(is.na(nodeColor$nodenum[i]) & nodeColor$node[i] == ""){
    nodeColor$color[i] = "brown"
  }else if(nodeColor$nodenum[i] < 100 && nodeColor$nodenum[i] > 90){
    nodeColor$color[i] = "dark grey"
  }else {
    nodeColor$color[i] = "brown"
  }
}

# edit the "Type" in the dataframe so that its prettier as a legend
TC_df$Type = gsub("lit_ratk","OGT derived from Ratkowsky fitting\n of reported temperature growth data",TC_df$Type)
TC_df$Type = gsub("lit","OGT in Reported Literature",TC_df$Type)
TC_df$Type = gsub("ratk","OGT from this Study",TC_df$Type)
colnames(TC_df)[3] = "OGT Source"

#plain tree with strains as tiplabels
p = ggtree(tre.phylo, size = 1) + # size used to be 2
  theme_tree2() +
  #geom_text(aes(label=node), hjust=-.3) +
  geom_nodepoint(color = nodeColor$color,size = 2) # size used to be 4

# add an X where the main clades are, and their main cluster names
p1 = p + geom_nodepoint(aes(subset = node == 117, x = x - branch.length * 0.5),shape = 4,size = 3) + # size used to be 7
  geom_text2(aes(subset = node == 117,x = x - branch.length * 0.5, label = "Clade C \n(Colwellia)", fontface = "italic"),vjust = 0.4,size = 4,hjust = -0.05 ,color = "purple",angle = 90) +
  geom_nodepoint(aes(subset = node == 145, x = x - branch.length * 0.5),shape = 4,size = 3) +
  geom_text2(aes(subset = node == 145,x = x - branch.length * 0.5, label = "Clade B \n(Cognaticolwellia)",fontface = "italic"),vjust = 0.4,size = 4,hjust = -0.05,color = "orange", angle = 90) +
  geom_nodepoint(aes(subset = node == 169, x = x - branch.length * 0.5),shape = 4,size = 3) +
  geom_text2(aes(subset = node == 169,x = x - branch.length * 0.5, label = "Clade A \n(Colwellia_A)",fontface = "italic"),vjust = 0.4,size = 4,hjust = -0.05,color = "darkgreen", angle = 90) # all sizes used to be 8


#have colored strain names according to sub-clade and all types of OGT source 
p2 = p1 %<+% TC_df + 
  geom_tiplab(aes(color = ClusterName),cex = 3.5,offset = .05, size = 2.5,show.legend = F) + # cex used to be 3.5
  geom_tippoint(aes(subset = !is.na(`OGT Source`),shape = `OGT Source`,x = x + 0.03),size = 2) + # size use to be 2.5 
  scale_color_manual(values = SubClade_color) + 
  theme(legend.justification = "top") + 
  #guides(fill = guide_legend("OGT Source"))+
  xlim(0,3.2) 

# add bar to designate the subclades & legend from source of OGTs
p3 = p2 +
  geom_strip("Colwellia BRX8.6","Colwellia hornerae strain IC036",barsize = 2,label = "Clade A1",offset = 1.1,fontsize = 3,hjust = -0.1) + 
  geom_strip("Colwellia MB3u-64","Colwellia Bg11-12",barsize = 2,label = "Clade A2",offset = 1.1,fontsize = 3,hjust = -0.1)

p3

#save the figure 
tiff("./Figures/Figure 1.tiff",width = 8,height = 11, units = "in", res = 1200) # 15 by 17 before
p3
dev.off()
