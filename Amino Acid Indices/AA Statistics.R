# Amino Acid Statistics
library(ggpubr) # stat_cor() & ggplot()
library(reshape2) # melt()
library(RColorBrewer) # brewer.pal()

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

#----------calculate  the linear correlation between AA and Clade A---------------------------------------------------
CladeA = filter(WG_Data,ClusterName == "Clade A1" | ClusterName == 'Clade A2')
  summary(lm(OGT ~ Arg.Lys.Ratio, data = CladeA )) # CladeA_ArgLys
  ggplot(CladeA,aes(OGT,Arg.Lys.Ratio)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ GRAVY, data = CladeA )) # CladeA_GRAVY
  ggplot(CladeA,aes(OGT,GRAVY)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aromaticity.Index, data = CladeA )) # CladeA_Aromaticity
  ggplot(CladeA,aes(OGT,Aromaticity.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Proline.Residue.Proportion, data = CladeA )) # CladeA_Proline
  ggplot(CladeA,aes(OGT,Proline.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aliphatic.Index, data = CladeA )) # CladeA_aliphatic
  ggplot(CladeA,aes(OGT,Aliphatic.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Acidic.Residue.Proportion, data = CladeA )) # CladeA_acidic
  ggplot(CladeA,aes(OGT,Acidic.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))

#----------calculate  the linear correlation between AA and Clade B---------------------------------------------------
  Clade = filter(WG_Data,ClusterName == "Clade B")
  summary(lm(OGT ~ Arg.Lys.Ratio, data = Clade )) # CladeA_ArgLys
  ggplot(Clade,aes(OGT,Arg.Lys.Ratio)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ GRAVY, data = Clade )) # CladeA_GRAVY
  ggplot(Clade,aes(OGT,GRAVY)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aromaticity.Index, data = Clade )) # CladeA_Aromaticity
  ggplot(Clade,aes(OGT,Aromaticity.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Proline.Residue.Proportion, data = Clade )) # CladeA_Aromaticity
  ggplot(Clade,aes(OGT,Proline.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aliphatic.Index, data = Clade )) # CladeA_aliphatic
  ggplot(Clade,aes(OGT,Aliphatic.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Acidic.Residue.Proportion, data = Clade )) # CladeA_aliphatic
  ggplot(Clade,aes(OGT,Acidic.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))

#----------calculate  the linear correlation between AA and Clade C---------------------------------------------------
Clade = filter(WG_Data,ClusterName == "Clade C")
  summary(lm(OGT ~ Arg.Lys.Ratio, data = Clade )) # CladeA_ArgLys
  ggplot(Clade,aes(OGT,Arg.Lys.Ratio)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ GRAVY, data = Clade )) # CladeA_GRAVY
  ggplot(Clade,aes(OGT,GRAVY)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aromaticity.Index, data = Clade )) # CladeA_Aromaticity
  ggplot(Clade,aes(OGT,Aromaticity.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Proline.Residue.Proportion, data = Clade )) # CladeA_Aromaticity
  ggplot(Clade,aes(OGT,Proline.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Aliphatic.Index, data = Clade )) # CladeA_aliphatic
  ggplot(Clade,aes(OGT,Aliphatic.Index)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
  summary(lm(OGT ~ Acidic.Residue.Proportion, data = Clade )) # CladeA_aliphatic
  ggplot(Clade,aes(OGT,Acidic.Residue.Proportion)) + 
    geom_point()+ 
    geom_smooth(method = "lm", color = "black", size = 0.3) + 
    stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")))
  
# ----------------- Finding the best combination of amino acid indices for a predictive model--------------------------
  # using a step function to find the combinaiton with the lowest AIC value
  step(lm(OGT~ Arg.Lys.Ratio + Aliphatic.Index + Aromaticity.Index + Acidic.Residue.Proportion + GRAVY + Proline.Residue.Proportion,data = WG_Data),direction = "both")
  
  # statistics of linear model with all amino acid indices
  mlm_all = with(WG_Data,lm(OGT ~ Proline.Residue.Proportion + Arg.Lys.Ratio + Aliphatic.Index + Aromaticity.Index + Acidic.Residue.Proportion + GRAVY))
  
  # statistics of linear model with individual amino acid indices
  mlm_Pro = with(WG_Data,lm(OGT ~ Proline.Residue.Proportion))
  mlm_Arg = with(WG_Data,lm(OGT ~ Arg.Lys.Ratio))
  mlm_Ali = with(WG_Data,lm(OGT ~ Aliphatic.Index))
  mlm_Aro = with(WG_Data,lm(OGT ~ Aromaticity.Index))
  mlm_Aci = with(WG_Data,lm(OGT ~ Acidic.Residue.Proportion))
  mlm_Gra = with(WG_Data,lm(OGT ~ GRAVY))
  
  #statistics of linear model of the best combination from original step function
  mlm_best = lm(OGT ~ Aliphatic.Index + GRAVY, data = WG_Data )
  
  # compares AIC between all linear models
  AIC(mlm_all,mlm_Aci,mlm_Aro,mlm_Ali,mlm_Arg,mlm_Pro,mlm_Gra,mlm_best)
  
  # calculate the predicted OGT for each genome based on every possible linear model (all, individual, and best)
  # and find the OGT difference with the observed OGT
  Calc_OGT = matrix(ncol = 11, nrow = nrow(WG_Data))
  OGT_diff = matrix(ncol = 10, nrow = nrow(WG_Data))
  for(i in 1:nrow(WG_Data)){
    
    #save strain name
    strn = paste(WG_Data$strain[i])
    Calc_OGT[i,1] = strn
    OGT_diff[i,1] = strn
    
    #save clustername
    clade = paste(WG_Data$ClusterName[i])
    Calc_OGT[i,2] = clade
    OGT_diff[i,2] = clade
    
    #save observed OGT
    Meas_OGT = WG_Data$OGT[i]
    Calc_OGT[i,3] = Meas_OGT
    
    #save each amino acid index
    a = WG_Data$Proline.Residue.Proportion[i]
    b = WG_Data$Arg.Lys.Ratio[i]
    c = WG_Data$Aliphatic.Index[i]
    d = WG_Data$Aromaticity.Index[i]
    e = WG_Data$Acidic.Residue.Proportion[i]
    f = WG_Data$GRAVY[i]
    
    # calculate predicted OGT ( a = proline residue, b = arginine lysine, c = aliphatic index, d = aromaticity, e = acidic residue, f = GRAVY)
    Calc_OGT[i,4] = coef(mlm_all)[[1]] + coef(mlm_all)[[2]]*a + coef(mlm_all)[[3]]*b + coef(mlm_all)[[4]]*c + coef(mlm_all)[[5]]*d + coef(mlm_all)[[6]]*e + coef(mlm_all)[[7]]*f
    Calc_OGT[i,5] = coef(mlm_Pro)[[1]] + coef(mlm_Pro)[[2]]*a
    Calc_OGT[i,6] = coef(mlm_Arg)[[1]] + coef(mlm_Arg)[[2]]*b
    Calc_OGT[i,7] = coef(mlm_Ali)[[1]] + coef(mlm_Ali)[[2]]*c
    Calc_OGT[i,8] = coef(mlm_Aro)[[1]] + coef(mlm_Aro)[[2]]*d
    Calc_OGT[i,9] = coef(mlm_Aci)[[1]] + coef(mlm_Aci)[[2]]*e
    Calc_OGT[i,10] = coef(mlm_Gra)[[1]] + coef(mlm_Gra)[[2]]*f
    Calc_OGT[i,11] = coef(mlm_best)[[1]] + coef(mlm_best)[[2]]*c + coef(mlm_best)[[3]]*f
    
    # calculate absolute OGT difference between observed and predicted
    OGT_diff[i,3] = abs(Meas_OGT - (coef(mlm_all)[[1]] + coef(mlm_all)[[2]]*a + coef(mlm_all)[[3]]*b + coef(mlm_all)[[4]]*c + coef(mlm_all)[[5]]*d + coef(mlm_all)[[6]]*e + coef(mlm_all)[[7]]*f))
    OGT_diff[i,4] = abs(Meas_OGT - (coef(mlm_Pro)[[1]] + coef(mlm_Pro)[[2]]*a))
    OGT_diff[i,5] = abs(Meas_OGT - (coef(mlm_Arg)[[1]] + coef(mlm_Arg)[[2]]*b))
    OGT_diff[i,6] = abs(Meas_OGT - (coef(mlm_Ali)[[1]] + coef(mlm_Ali)[[2]]*c))
    OGT_diff[i,7] = abs(Meas_OGT - (coef(mlm_Aro)[[1]] + coef(mlm_Aro)[[2]]*d))
    OGT_diff[i,8] = abs(Meas_OGT - (coef(mlm_Aci)[[1]] + coef(mlm_Aci)[[2]]*e))
    OGT_diff[i,9] = abs(Meas_OGT - (coef(mlm_Gra)[[1]] + coef(mlm_Gra)[[2]]*f))
    OGT_diff[i,10] = abs(Meas_OGT - (coef(mlm_best)[[1]] + coef(mlm_best)[[2]]*c + coef(mlm_best)[[3]]*f))
  }
  
  colnames(Calc_OGT) = c("strain","clade","Meas_OGT", "mlm_all","mlm_Pro","mlm_Arg","mlm_Ali","mlm_Aro","mlm_Aci","mlm_Gra","mlm_best")
  colnames(OGT_diff) = c("strain","clade","mlm_all","mlm_Pro","mlm_Arg","mlm_Ali","mlm_Aro","mlm_Aci","mlm_Gra","mlm_best")
  Calc_OGT = as.data.frame(Calc_OGT)
  OGT_diff = as.data.frame(OGT_diff)
  
  write.table(x= Calc_OGT,file="./Amino Acid Indices/Genome AA model.tsv",sep="\t",quote=F)
  write.table(x= OGT_diff,file="./Amino Acid Indices/Genome AA model diff.tsv",sep="\t",quote=F)
