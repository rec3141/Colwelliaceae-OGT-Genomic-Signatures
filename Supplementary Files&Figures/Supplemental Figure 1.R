# create pangenome vs core genome plot, volcano plot, and pairwise differences of protein family and OGT
source("./Supplementary Files&Figures/f-pangenome.R")
library(DECIPHER) # ReadDendrogram()
library(ape) # read.tree()
library(reshape2) # melt()
library(RColorBrewer) # brewer.pal()
library(phytools) #midpoint.root()
library(dendextend)  # set()
library(colorspace) # rainbow_hcl()
library(phylogram) # read.dendrogram()
library(rlang) # is_empty()
library(dplyr) # filter()
library(ggplot2) # ggplot()
library(ggridges) # geom_density_ridges()
library(gridExtra) # grid.arrange
theme_set(theme_ridges())

# ------------------------------------------------------------load files -------------------------------------------------------------
# load PATRIC phylogenetic and protein family files
treefile="./Patric Input/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Patric Input/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
colfeat = read.table("./Patric Input/Colwelliaceae PF Features.tsv",fill = T,head = T,sep = "\t",colClasses = "character") # patric protein features
plfdat = read.table("./Patric Input/Colwelliaceae_families_local.tsv",sep="\t",row.names = 1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]
pgfdat = read.table("./Patric Input/Colwelliaceae_families_global.tsv",sep="\t",row.names=1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]

# load Protein family amino acid indice calculations
colind = read.table("./Supplementary Files&Figures/Protein Indices Data.tsv",row.names = 1, sep="\t",stringsAsFactors = F) # Protein Amino Acid Index
colsave = readRDS("./Supplementary Files&Figures/Colwelliaceae PF AA Calc.RDS")

# load OGT by strain and  Sub-Clade IDs
ratk = read.table("./OGT/Strain Ratk OGT Final.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table
Temp.Clust.ID = read.table("./Phylogenetics/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

# set up list of literature based Colwellia growth rates 
ex_col = paste(c("Colwellia psychrerythraea ACAM 605","Colwellia hornerae strain ACAM 607","Colwellia piezophila ATCC BAA-637","Colwellia demingiae strain ACAM 459",
                 "Colwellia psychotropica"),collapse = "|")
taxaname = "Colwellia (PGF)"
#---------------------------------------------------------- Format Files -----------------------------------------------------------------
#read whole-genome tree from PATRIC
tre.dend = ReadDendrogram(treefile)
tre.phylo = read.tree(treefile)

#set up patric metatable
rownames(genomes) = genomes$genome.genome_id

#set up protein features
rownames(colfeat) = colfeat$feature.patric_id
colfeat = colfeat[,-2]
colnames(colfeat)[1] = "genome.genome_id"

# limit datasets to common genomes and merge protein features and amino acid indices into dataframe
common_ids = intersect(rownames(colfeat), rownames(colind))
colfeat = colfeat[common_ids, ]
colind = colind[common_ids, ]
colout = data.frame(colfeat, colind)

# add the patric ids to the OGT strains
rownames(ratk) = ratk$strn
genome_ID = genomes[,c(2,3)]
genome_ID = subset(genome_ID,genome.genome_id %in% tre.phylo$tip.label)
ratk$id.patric = genome_ID$genome.genome_id[match(ratk$strn,genome_ID$genome.genome_name)]

#-------------------------------------------------------set up phylogeny -----------------------------------------------------------------
#use the global protein family
mat = as.matrix(colsave[[1]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

# keep only genomes that are in both the tree and the amino acid indices calculations
mattmp = mat[,intersect(tre.phylo$tip.label, colnames(mat))]
dend = keep.tip(tre.phylo, intersect(tre.phylo$tip.label,colnames(mattmp)))

# clean up and reorder tree
dend$node.label = NULL # Need to remove node labels
dend = midpoint.root(dend)
dend = ladderize(dend)
dend = drop.tip(dend,"1380381.3") #drop long tip (poor genome assembly)

# convert from 'phylo' to 'dendrogram'
dend1 = ReadDendrogram(textConnection(write.tree(dend)))
#dend1 = set(dend1, "branches_lwd", 5)

#--------------------------------------------------------------set up color scheme----------------------------------------------------------
#set up OGT colors
cols = rev(brewer.pal(11, "RdBu"))
pal = colorRampPalette(c("blue", "red")) # Define colour pallete
pal = colorRampPalette(cols) # Use the following line with RColorBrewer

#Rank variable for colour assignment
rowcols = ratk$OGT[match(labels(dend1),ratk$id.patric)]
rc.order = findInterval(rowcols, sort(rowcols))
rowcols = pal(length(sort(rowcols)))[rc.order]

#set up Sub-Clade Colors
Temp.Clust.colors = c("Clade A1" = "seagreen3","Clade A2" = "olivedrab4","Clade A3" = "chartreuse3", 
                      "Clade B1" = "darkorange2", "Clade B2" = "darkgoldenrod3", 
                      "Clade C" = "purple")

#------------------------------------------------------------------ pangenome vs. core genome plot -----------------------------------------------------
# use global protein families dataframe
mat = as.matrix(colsave[[1]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

#subset data
matsmall = mat[,colnames(mattmp)]
matcols = colnames(matsmall)

#calculate gene family frequency spectrum
mat = !is.na(as.matrix(colsave[[1]])) # use global families instead of local families
mat = mat + 0
mat = mat[,matcols]

Gk <- f.getspectrum(mat)

genomesize <- median(colSums(mat>0)) # median genome size measured in gene families
ng <- dim(mat)[2] #number of genomes

# Calculate 100 permutations each of the pangenome and core genome
perm.pangenome <- f.pangenome(mat,100)
perm.core <- f.core(mat,100)

# Calculate the exact mean pan and core genome curves
# from the gene frequency spectrum G(k)
mean.pangenome <- f.meanpancore(Gk)$pan
mean.core <- f.meanpancore(Gk)$core
pancore <- c(mean.pangenome,mean.core)

# Calculate the RMS value for the permutations
rms.perm <- mean(f.rms(c(mean.pangenome,mean.core),rbind(perm.pangenome,perm.core)))

tiff(file="./Figures/Test figures/pancore.tiff",width = 11,height = 8, units = "in", res = 1200)
par(mar =c(5,5,1,1))
# Prepare a new plot window
plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log="",
     xlab="Number of Genomes Added", ylab="PATRIC Global Protein Family",pch='',cex.lab = 1.2,cex.axis =0.83)

# Plot polygons outlining permutations
polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

# Add the mean pan and core genome curves to the plot
points(1:ng,mean.pangenome,type='l',col = "red")
points(1:ng,mean.core,type='l')

dev.off()

