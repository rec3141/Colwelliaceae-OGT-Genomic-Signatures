library(scales)
library(dplyr)
library(reshape2)
library(DescTools) #medianCI
# ----------------------------------------------------read in files---------------------------------------------------------
ratk_TC = read.table("./Phylogenetics/Clade Ratk Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk_TC_avg = read.table("./Phylogenetics/Clade Ratk Avg Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk_TC_med = read.table("./Phylogenetics/Clade Ratk Med Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk_TC_sd = read.table("./Phylogenetics/Clade Ratk SD Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
Temp.Clust.ID = read.table("./Phylogenetics/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

# ratk = read.table("./Data/Strain Ratk OGT Outputs Clean.tsv",sep = "\t",header = T,colClasses = "character")
GR = read.table("./Refs& Raw Input/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character")
#-------------------------------------set up color schemes ------------------------------------------------------------
Temp.Clust.colors = c("Clade A1" = "seagreen3","Clade A2" = "olivedrab4", 
                      "Clade B" = "darkorange2","Clade C" = "purple")
#----------------------------------------------------set up functions ------------------------------------------------------------------------
# load functions
rat83 = function(x, b, T1, cc, T2) {
  rT = b * (x - T1) * (1 - exp(cc * (x - T2)))
  return(rT)
}

rat82a = function(x, b, T0) {
  r = b*(x - T0)
}

#------------------------------------------------------set up  files -------------------------------------------------------------------------
#add "Colwellia" in front all strains in growthrates file
GR$strain = paste("Colwellia",GR$strain)

colnames(Temp.Clust.ID)[1] = "strain"

# set up growth rate file
cols.num_GR = colnames(GR)[2:5]
GR[cols.num_GR] = sapply(GR[cols.num_GR],as.numeric)

# set up the ratkowsky fit by temperature cluster 
init.plot = as.data.frame(ratk_TC_med)
rownames(init.plot) = init.plot[,1]
cols.num = colnames(init.plot)[2:6]
init.plot[cols.num] = sapply(init.plot[cols.num],as.numeric)
ratk_TC_sd[cols.num] = sapply(ratk_TC_sd[cols.num],as.numeric)

# add color column to dataframe
init.plot$col = Temp.Clust.colors[match(rownames(init.plot),names(Temp.Clust.colors))]
ratk_TC$col = Temp.Clust.colors[match(ratk_TC$Clade,names(Temp.Clust.colors))]
grd = merge.data.frame(Temp.Clust.ID,GR,by = "strain", all = F)

#set up growth rates data to plot
colnames(grd) = c("strain_original","strain","rep","temp","rate","r2") #'strain' is actually 'clade'
grd$strain_original = as.factor(grd$strain_original)
grd$temp = as.numeric(grd$temp) + 273
grd$rep = as.factor(grd$rep)
grd$rate = as.numeric(grd$rate)
grd = na.omit(grd)
grd$sqperday = sqrt(grd$rate) 

#------------------------------------Test figures with unique clades with overlap plot next to them ----------------------------------
tiff("./Figures/Figure 2.tiff",width = 11,height = 6, units = "in",res = 1200) # original width of 6.8

layout.matrix = matrix(c(1,2,3,4,5,5),nrow = 2,ncol = 3, byrow = F)
graphics::layout(mat = layout.matrix)

labels = c("A", "B", "C", "D")
op = par(mar = c(3,3,2,1),oma = c(2,3,0,0))

for(i in 1:length(unique(ratk_TC$Clade))) {
  #i = 1
  #set up clade name and temperature range
  clade = sort(unique(ratk_TC$Clade))[i]
  str_nm = length(unique(grd[grd$strain == clade,c("strain_original")]))
  print(as.character(clade))
  xs = seq(268,298,0.1)
  
  # subset the dataframes of growth rates and ratkowksy parameters by Clade
  sgrd = grd[grd$strain == clade,c("temp","sqperday")]
  Clade_df = ratk_TC[ratk_TC$Clade == clade,c("b","Tmin","c","Tmax","topt","col")]
  cols.num = colnames(Clade_df)[1:5]
  Clade_df[cols.num] = sapply(Clade_df[cols.num],as.numeric)
  
                                                         
  #isolate the parameters for the median line
  Clade_param_quant = ratk_TC_med[ratk_TC_med$Clade == clade,]
  
  # calcualte the fit of the clade with avg but listed as median for convenience
  rt.fit_med = rat83(xs, as.numeric(Clade_param_quant$b), as.numeric(Clade_param_quant$Tmin), as.numeric(Clade_param_quant$c), as.numeric(Clade_param_quant$Tmax))
  
  # plot original growth rates, individual ratkowsky fits, and median ratkowski fits
  plot(sgrd$temp - 273, sgrd$sqperday^2, pch=19, col=Clade_df$col,cex =0.5, # remove the ^2 to look at the shape, not for the paper
       xlim=c(-1,17), 
       ylim=c(0,3),
       cex.axis= 0.83)
  for ( j in 1:nrow(Clade_df)){
    rt.fit = rat83(xs, as.numeric(Clade_df$b[j]), as.numeric(Clade_df$Tmin[j]), as.numeric(Clade_df$c[j]), as.numeric(Clade_df$Tmax[j]))
    rt.fit[rt.fit < 0] = NA
    lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df$col, 0.03))
  }
  #lines(xs-273, rt.fit_med^2, type="l",lwd=3, col = Clade_df_med$col)
  lines(xs-273, rt.fit_med^2, type="l",lwd=3, col = Clade_df$col)
  
  # Find the OGT from the maximum growth rate of the median line
  topt_med = round(xs[which.max(rt.fit_med)],2)
  GR_med = max(rt.fit_med,na.rm = T)^2
  mtext(paste0("(", labels[i], ")"), side = 3, adj = 0.05, line = -1.3)
  title(paste0(clade," (N = ",str_nm,")"),cex.main=1.3)
  
  b_text = paste0("b =  ",round(as.numeric(Clade_param_quant$b),2))
  c_text = paste0("c =  ",round(as.numeric(Clade_param_quant$c),2))
  Tmin_text = paste0("Tmin =  ",round(as.numeric(Clade_param_quant$Tmin),2)-273,"(°C)")
  Tmax_text = paste0("Tmax =  ",round(as.numeric(Clade_param_quant$Tmax),2)-273,"(°C)")
  text(x = -0.8, y = 0,paste(b_text,c_text,Tmin_text,Tmax_text,sep = "\n") ,cex = 1.1,adj = c(0,0))
  
  }

title(xlab = "Temperature (°C)",ylab = expression("Growth Rate" ~ (d^{-1})),outer = T, line = 1, cex.lab = 2) # line = 1, cex.lab = 1.8

plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (°C)", ylab = "",
     cex.axis= 0.83, cex.lab= 1.2) 
mtext(paste0("(E)"), side = 3, adj = 0.05, line = -1.3)

for(i in 1:nrow(init.plot)) {

  strn = rownames(init.plot)[i]
  print(as.character(strn))
  
  xs = seq(0,1000,0.1)
  rt.fit = rat83(xs, init.plot$b[i], init.plot$Tmin[i], init.plot$c[i], init.plot$Tmax[i])
  rt.fit[rt.fit < 0] = NA
  
  xs = xs - 273
  lines(xs, rt.fit^2, type="l", col=init.plot$col[i],lwd=6, lty = i)
}
legend("topright",legend = rownames(init.plot), bty = "n",
       lwd = 4, col = init.plot$col,lty = 1:6 ,cex = 1.3)


par(op)
dev.off()

