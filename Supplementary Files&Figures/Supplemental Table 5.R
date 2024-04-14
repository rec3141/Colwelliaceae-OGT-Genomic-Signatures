library(openxlsx)
library(ape)
library(phytools)
library(phylogram)
library(dplyr)

#---------------------------------------------------------load files---------------------------------------------------------------------
# load PATRIC phylogenetic and protein family files
treefile="./Patric Input/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Patric Input/Colwelliaceae UnRooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
feat = read.table("./Patric Input/Colwelliaceae PF Features.tsv",fill = T,head = T,sep = "\t",colClasses = "character")
cfg = read.table("./Patric Input/Colwelliaceae_families_global.tsv",fill=T,head=F,sep="\t",colClasses="character",quote="")
CladeID = read.table("./Phylogenetics/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

#load presence/absence dataframe calculated from Protein Family Indices.R
Gene_Prescence = read.csv("./Amino Acid Indices/Gene Prescence.csv",fill = T,header = T, check.names = F)

#load OGT dataframe calculated from Strain OGT.R
ratk = read.table("./OGT/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character")

#blast output of identified cold adaptation genes against strains of interest ( all strains in my phylogenetic tree)
blast_output = read.table("./Supplementary Files&Figures/blast_CA_Colwelliaceae.txt",fill = T, header = T, sep = '\t',quote="")

#---------------------------------------------------------format files ------------------------------------------------------------------
colnames(CladeID)[1] = "strn"

rownames(genomes) = genomes$genome.genome_id
colnames(feat)[1] = "genome.genome_id"
rownames(Gene_Prescence) = Gene_Prescence[,1]; Gene_Prescence[,1] = NULL
colnames(ratk) = c("strn","topt","Type")

rownames(cfg)=cfg$V1
cfg$V1=NULL
colnames(cfg) = c("feature.pgfam_id","type","description")

#set up tree
tre.phylo = read.tree(treefile)
tre.phylo = drop.tip(tre.phylo,"1380381.3") #drop long tip (poor genome assembly)
common_ids = intersect(tre.phylo$tip.label,rownames(genomes))
tre.phylo = keep.tip(tre.phylo, common_ids) # keeping the tips that are in the dataframe 
tre.phylo$node.label = NULL
tre.phylo = midpoint.root(tre.phylo)
tre.phylo = ladderize(tre.phylo)

dend <- read.dendrogram(textConnection(write.tree(tre.phylo))) # convert phylo to dend
dend = as.dendrogram(dend)

#set up whole dataframe to analyze by cluster
#CladeID$patric.id = genome_ID$genome.genome_id[match(CladeID$strn,genome_ID$genome.genome_name)]

#----------------------------------------create dataframe of cold adaptation proteins take from the family of interest-------------------------------------------------------
# set up blast information dataframe 
blast_output$sid = gsub(":","|",blast_output$sid)
colnames(blast_output)[4] = "feature.patric_id"
blast_output_feat = merge.data.frame(blast_output,feat,by = "feature.patric_id", all = F)
blast_output_feat = blast_output_feat[,-c(2,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21)]

# add the patric ids to the OGT strains
genome_ID = genomes[,c(2,3)]
genome_ID = subset(genome_ID,genome.genome_id %in% tre.phylo$tip.label)
ratk$id.patric = genome_ID$genome.genome_id[match(ratk$strn,genome_ID$genome.genome_name)]
ratk$topt = as.numeric(ratk$topt)
ratk = na.omit(ratk)
rownames(ratk) = ratk$id.patric

#a list of the all of the protein families that the result of the U.S. blast 
# against all Colwelliaceae in my tree belong to
blast_match = as.character(unique(blast_output_feat$feature.pgfam_id))

#set up the Gene Presence dataframe
Gene_Prescence = Gene_Prescence[,labels(dend)] # reorders the matrix according to the dendrogram tip labels
GeneP1 = as.matrix(Gene_Prescence)
GeneP1[is.infinite(GeneP1)] = NA
GeneP1[is.nan(GeneP1)] = NA
GeneP1 = na.omit(GeneP1)

# filter out dataframe of presence/abscence of all protein families in all strains by those
# that were a match to those where the cold adaptation proteins belong to
GeneP0 = GeneP1
GeneP1 = as.matrix(GeneP1[rownames(GeneP1) %in% blast_match,])

#set up the strain names instead of the Patric ID numbers 
colnames(GeneP1) = genomes[colnames(GeneP1),"genome.genome_name"]


#--------------------------------------------------Analyze cold and warm clusters of hornerea ----------------------------------------------
# isolate the cold clade : Clade A2, without Bg11-12 - creates a dataframe of the protein families found and absent in each strain of this clade
CH_cold_id = subset.data.frame(CladeID, ClusterName == "Clade A2");rownames(CH_cold_id)=CH_cold_id$patric.id
CH_cold_id = filter(CH_cold_id,CH_cold_id$strn != "Colwellia Bg11-12")
CH_cold = GeneP1[,colnames(GeneP1) %in% CH_cold_id$strn]


CH_cold = t(CH_cold>0)+0

CH_cold_PFG = CH_cold[,colSums(CH_cold)>0] # remove protein families that are not present in any strain in the cluster
CH_cold_commonPFG = CH_cold_PFG[,colSums(CH_cold_PFG) == nrow(CH_cold_PFG)] # keep protein families that are only present in all strains of that cluster
CH_cold_commonPFG = colnames(CH_cold_commonPFG) # all cold adaptation associated protein families that are present in at least one strain of the cluster

#isolate Colwellia Bg11-12 - creates a dataframe of the protein families found and absent in each strain of this clade
CH_warm_1_id = filter(CladeID,CladeID$strn == "Colwellia Bg11-12");rownames(CH_warm_1_id)=CH_warm_1_id$patric.id
CH_warm_1 = GeneP1[,colnames(GeneP1) %in% CH_warm_1_id$strn,drop = F]

CH_warm_1 = t(CH_warm_1>0)+0

CH_warm_1_PFG = CH_warm_1[,colSums(CH_warm_1)>0,drop = F] # remove protein families that are not present in any strain in the cluster
CH_warm_1_commonPFG = colnames(CH_warm_1_PFG) # all cold adaptation associated protein families that are present in this strain

#isolate the warm clade: Clade A1 - creates a dataframe of the protein families found and absent in each strain of this clade
CH_warm_id = filter(CladeID,ClusterName == "Clade A1");rownames(CH_warm_id)=CH_warm_id$patric.id
CH_warm = GeneP1[,colnames(GeneP1) %in% CH_warm_id$strn,drop = F]

CH_warm = t(CH_warm>0)+0

CH_warm_PFG = CH_warm[,colSums(CH_warm)>0, drop = F] # remove protein families that are not present in any strain in the cluster
CH_warm_commonPFG = CH_warm_PFG[,colSums(CH_warm_PFG) == nrow(CH_warm_PFG)] # keep protein families that are only present in all strains of that cluster
CH_warm_commonPFG = colnames(CH_warm_PFG) # all cold adaptation associated protein families that are present in this strain


# set diff of the the different unique global protein families 
diffcw1 = setdiff(CH_cold_commonPFG,CH_warm_1_commonPFG) # present in Clade A2 but not in Bg11-12 *
diffw1c = setdiff(CH_warm_1_commonPFG,CH_cold_commonPFG) # present in Bg11-12 but not in Clade A2
diffwC = setdiff(CH_warm_commonPFG,CH_cold_commonPFG) # present in Clade A1 but not Clade A2
diffcw = setdiff(CH_cold_commonPFG,CH_warm_commonPFG) # present in Clade A2 but not in Clade A1*

diffall_warm = intersect(diffcw,diffcw1) # find the protein families that are present in A2 but absent in both Bg11-12 & Clade A1

# this removes the absent protein family in A1 and Bg11-12 but present in A2 from the set of proteins previously
unq_diffcw1 = diffcw1[! diffcw1 %in% diffall_warm]
unq_cw = diffcw[! diffcw %in% diffall_warm]

#setup dataframe of proteins present
unq_diffcw1 = data.frame(`Global Protein Family` = diffcw1,"Clade Absence" = "Colwellia Bg11-12")
unq_cw = data.frame(`Global Protein Family` = diffcw,"Clade Absence" = "Clade A1")
diffall_warm = data.frame(`Global Protein Family` = diffall_warm,"Clade Absence" = "Both")

# add rows together
P_A_df = rbind(unq_diffcw1,unq_cw,diffall_warm)
P_A_df$`GPF Description` = cfg$description[match(P_A_df$Global.Protein.Family,cfg$feature.pgfam_id)]

#write.csv(P_A_df,"./Figures/Supplemental Table 5.csv")

