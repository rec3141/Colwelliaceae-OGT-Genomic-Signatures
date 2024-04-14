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
colnames(Temp.Clust.ID)[1] = "strain"

# set up growth rate file
cols.num_GR = colnames(GR)[2:5]
GR[cols.num_GR] = sapply(GR[cols.num_GR],as.numeric)
GR$strain<-sub("^","Colwellia ",GR$strain)

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

# Filter by Clade A2, replace Bg11-12 Clade with Bg11-12
grd = filter(grd,strain == "Clade A2")
grd$strain = ifelse(grd$strain_original == "Colwellia Bg11-12","Colwellia Bg11-12","Clade A2")

#filter Ratk TC by Clade A2, replace Bg11-12 Clade with Bg11-12
ratk_TC = filter(ratk_TC,Clade == "Clade A2")
ratk_TC$Clade = ifelse(ratk_TC$Clade == "Colwellia Bg11-12","Colwellia Bg11-12","Clade A2")



#----------------------------------------Plot Clade ratkowsky fits for Clades of Interest (Clade A2 & Clade A3)-------------------------------------------------------
tiff("./Figures/Clade Ratk A2 & A3.tiff",width = 6.8,height = 9.2, units = "in", res = 1200)
par(mfrow=c(1,1),mar = c(4,5,2,1))
#set up empty plot
plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (?C)", ylab = "",
     cex.axis=.83, cex.lab=1.2)
title(ylab = expression("Average Growth Rate" ~ (d^{-1})), line = 2, cex.lab = 1.2)
for(i in 1:2) {
  i = 1
  #set up clade name and temperature range
  clade = sort(unique(ratk_TC$Clade))[i] #i
  print(as.character(clade))
  xs = seq(268,298,0.1)
  
  # subset the dataframes of growth rates and ratkowksy parameters by Clade
  sgrd = grd[grd$strain == clade,c("temp","sqperday")]
  sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday")])
  Clade_df = ratk_TC[ratk_TC$Clade == clade,c("b","Tmin","c","Tmax","topt","col")]
  
  # set up clade parameters to numeric and melt the dataframe
  cols.num = colnames(Clade_df)[1:5]
  Clade_df[cols.num] = sapply(Clade_df[cols.num],as.numeric)
  Clade_df_melt = melt(Clade_df,id.vars = "col")
  
  # identify the median,each clade and make a dataframe of each
  Clade_df_med = Clade_df_melt %>% group_by(col,variable) %>% dplyr::summarize(med = median(value))
  Clade_df_med = recast(Clade_df_med,col~variable)
  
  #Calculate the median fit of the clade 
  rt.fit_med = rat83(xs, Clade_df_med$b, Clade_df_med$Tmin, Clade_df_med$c, Clade_df_med$Tmax)
  rt.fit_med[rt.fit_med < 0] = NA
  
  if ( i == 2){
    for ( j in 1:nrow(Clade_df)){
      rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
      rt.fit[rt.fit < 0] = NA
      lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df_med$col, 0.3))
    }
  } else {
    for ( j in 1:nrow(Clade_df)){
      rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
      rt.fit[rt.fit < 0] = NA
      lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df_med$col, 0.03))
    }
    
  }
  lines(xs-273, rt.fit_med^2, type="l",lwd=5, lty = i) #,col = Clade_df_med$col)
  
}
legend("topleft",legend = c("Clade A2", "Clade A3"), bty = "n",
       lwd = 4, col =Temp.Clust.colors[2:3],lty = 2:3,cex = 1)
dev.off()

#----------------------------------------Plot Overlapped Clade ratkowsky fits -------------------------------------------------------
tiff("./Figures/Clade Ratk All.tiff",width = 6.8,height = 9.2, units = "in", res = 1200)
plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (?C)", ylab = "",
     cex.axis= 0.83, cex.lab= 1.2) 
title(ylab = expression("Average Growth Rate" ~ (d^{-1})), line = 2, cex.lab = 1.2)

for(i in 1:nrow(init.plot)) {
  #i = 1
  strn = rownames(init.plot)[i]
  print(as.character(strn))
  
  xs = seq(0,1000,0.1)
  rt.fit = rat83(xs, init.plot$b[i], init.plot$Tmin[i], init.plot$c[i], init.plot$Tmax[i])
  rt.fit[rt.fit < 0] = NA
  
  xs = xs - 273
  lines(xs, rt.fit^2, type="l", col=init.plot$col[i],lwd=6, lty = i)
}
legend("topleft",legend = rownames(init.plot), bty = "n",
       lwd = 4, col = init.plot$col,lty = 1:6 ,cex = 1)

dev.off()

#------------------------------------Test figures with unique clades with overlap plot next to them ----------------------------------
tiff("./Figures/Test figures/Clade Ratk Individual and overlap.tiff",width = 11,height = 6, units = "in",res = 1200) # original width of 6.8

layout.matrix = matrix(c(1,2,3,4,5,5),nrow = 2,ncol = 3, byrow = F)
graphics::layout(mat = layout.matrix)

labels = c("A", "B", "C", "D")
op = par(mar = c(3,3,2,1),oma = c(2,3,0,0))

for(i in 1:length(unique(ratk_TC$Clade))) {
  #i = 1
  #set up clade name and temperature range
  clade = sort(unique(ratk_TC$Clade))[i] 
  print(as.character(clade))
  xs = seq(268,298,0.1)
  
  # subset the dataframes of growth rates and ratkowksy parameters by Clade
  sgrd = grd[grd$strain == clade,c("temp","sqperday")]
  #sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday")])
  Clade_df = ratk_TC[ratk_TC$Clade == clade,c("b","Tmin","c","Tmax","topt","col")]
  cols.num = colnames(Clade_df)[1:5]
  Clade_df[cols.num] = sapply(Clade_df[cols.num],as.numeric)
  
  # calculate the clade CI to the graph
  # Clade_param_CI = ratk_TC %>% dplyr::group_by(Clade) %>% dplyr::summarize(med_b = median(as.numeric(b)),
  #                             med_c = median(as.numeric(c)),
  #                             med_Tmin = median(as.numeric(Tmin)),
  #                             med_Tmax = median(as.numeric(Tmax)),
  #                             CI_b = paste0("[",round(MedianCI(as.numeric(b))[2],3),",",round(MedianCI(as.numeric(b))[3],3),"]"),
  #                             CI_c = paste0("[",round(MedianCI(as.numeric(c))[2],3),",",round(MedianCI(as.numeric(c))[3],3),"]"),
  #                             CI_Tmin = paste0("[",round(MedianCI(as.numeric(Tmin))[2],3)-273,",",round(MedianCI(as.numeric(Tmin))[3],3)-273,"]"),
  #                             CI_Tmax = paste0("[",round(MedianCI(as.numeric(Tmax))[2],3)-273,",",round(MedianCI(as.numeric(Tmax))[3],3)-273,"]"), 
  #                             col = unique(col))
  # Clade_param_CI = Clade_param_CI[Clade_param_CI$Clade == clade,]
  
  
  
  #calculate the Clade median absolute deviation
  # Clade_param_mad = ratk_TC %>% dplyr::group_by(Clade) %>% dplyr::summarize(med_b = median(as.numeric(b)),
  #                                                                           med_c = median(as.numeric(c)),
  #                                                                           med_Tmin = median(as.numeric(Tmin)),
  #                                                                           med_Tmax = median(as.numeric(Tmax)),
  #                                                                           mad_b = paste(round(mad(as.numeric(b)),3)),
  #                                                                           mad_c = paste(round(mad(as.numeric(c)),3)),
  #                                                                           mad_Tmin = paste(round(mad(as.numeric(Tmin)-273),3)),
  #                                                                           mad_Tmax = paste(round(mad(as.numeric(Tmax)-273),3)))
  # Clade_param_mad = Clade_param_mad[Clade_param_mad$Clade == clade,]
  
  #calculate the clade inter-quantile range 
  Clade_param_quant = ratk_TC %>% dplyr::group_by(Clade) %>% dplyr::summarize(med_b = round(median(as.numeric(b)),2),
                                                                              med_c = round(median(as.numeric(c)),2),
                                                                              med_Tmin = round(median(as.numeric(Tmin)),1),
                                                                              med_Tmax = round(median(as.numeric(Tmax)),1)#,
                                                                              # quant_b = round(IQR(as.numeric(b)),3),
                                                                              # quant_c = round(IQR(as.numeric(c)),3),
                                                                              # quant_Tmin = round(IQR(as.numeric(Tmin)),3),
                                                                              # quant_Tmax = round(IQR(as.numeric(Tmax)),3)
  )
  Clade_param_quant = Clade_param_quant[Clade_param_quant$Clade == clade,]
  
  #Calculate the median fit of the clade
  rt.fit_med = rat83(xs, Clade_param_quant$med_b, Clade_param_quant$med_Tmin, Clade_param_quant$med_c, Clade_param_quant$med_Tmax)
  
  # plot original growth rates, individual ratkowsky fits, and median ratkowski fits
  plot(sgrd$temp - 273, sgrd$sqperday^2, pch=19, col=Clade_df$col,cex =0.5, # remove the ^2 to look at the shape, not for the paper
       xlim=c(-1,17), 
       ylim=c(0,3),
       cex.axis= 0.83)
  for ( j in 1:nrow(Clade_df)){
    rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
    rt.fit[rt.fit < 0] = NA
    lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df$col, 0.03))
  }
  #lines(xs-273, rt.fit_med^2, type="l",lwd=3, col = Clade_df_med$col)
  lines(xs-273, rt.fit_med^2, type="l",lwd=3, col = Clade_df$col)
  
  # Find the OGT from the maximum growth rate of the median line
  topt_med = round(xs[which.max(rt.fit_med)],2)
  GR_med = max(rt.fit_med,na.rm = T)^2
  mtext(paste0("(", labels[i], ")"), side = 3, adj = 0.05, line = -1.3)
  title(paste0(clade),cex.main=1.3)
  
  b_text = paste0("b =  ",Clade_param_quant$med_b)
  c_text = paste0("c =  ",Clade_param_quant$med_c)
  Tmin_text = paste0("Tmin =  ",Clade_param_quant$med_Tmin-273,"(°C)")
  Tmax_text = paste0("Tmax =  ",Clade_param_quant$med_Tmax-273,"(°C)")
  
  
  # b_text = paste0("b =  ",round(Clade_param_CI$med_b,2),", 95% CI ", Clade_param_CI$CI_b )
  # c_text = paste0("c =  ",round(Clade_param_CI$med_c,2),", 95% CI ", Clade_param_CI$CI_c)
  # Tmin_text = paste0("Tmin =  ",round(Clade_param_CI$med_Tmin,2)-273,"?C, 95% CI ", Clade_param_CI$CI_Tmin)
  # Tmax_text = paste0("Tmax =  ",round(Clade_param_CI$med_Tmax,2)-273,"?C, 95% CI ", Clade_param_CI$CI_Tmax)
  
  text(x = -0.8, y = 0,paste(b_text,c_text,Tmin_text,Tmax_text,sep = "\n") ,cex = 1.1,adj = c(0,0))
  
}

title(xlab = "Temperature (°C)",ylab = expression("Growth Rate" ~ (d^{-1})),outer = T, line = 1, cex.lab = 2) # line = 1, cex.lab = 1.8

plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (°C)", ylab = "",
     cex.axis= 0.83, cex.lab= 1.2) 
#title(ylab = expression("Average Growth Rate" ~ (d^{-1})), line = 2, cex.lab = 1.2)
mtext(paste0("(E)"), side = 3, adj = 0.05, line = -1.3)

for(i in 1:nrow(init.plot)) {
  #i = 1
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

