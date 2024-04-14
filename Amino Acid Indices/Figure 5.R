library(ggpubr) # stat_cor() & ggplot()
library(reshape2) # melt()
library(RColorBrewer) # brewer.pal()

# ------------------------------load files and set them up -------------------------------------------------------------
Calc_OGT = read.table("./Amino Acid Indices/Genome AA model.tsv", head = T, sep = "\t", fill = T)
OGT_diff = read.table("./Amino Acid Indices/Genome AA model diff.tsv", head = T, sep = "\t", fill = T)

#read optimal growth temperature data - just to get the source of each OGT
OGT = read.table("./OGT/Strain Ratk OGT Final.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table

Calc_OGT$Type = OGT$Type[match(Calc_OGT$strain,OGT$strn)]
Calc_OGT$Meas_OGT = as.numeric(Calc_OGT$Meas_OGT)
Calc_OGT$mlm_best = as.numeric(Calc_OGT$mlm_best)
RMSE = signif(sqrt(mean((OGT_diff$mlm_best)^2)),2)

# edit the "Type" in the dataframe so that its prettier as a legend
Calc_OGT$Type = gsub("lit_ratk","OGT derived from Ratkowsky fitting of reported temperature growth data",Calc_OGT$Type)
Calc_OGT$Type = gsub("lit","OGT in Reported Literature",Calc_OGT$Type)
Calc_OGT$Type = gsub("ratk","OGT from this Study",Calc_OGT$Type)
colnames(Calc_OGT)[12] = "OGT Source"


#------------------------------------------ Plot predicted OGT (best combination) vs. measured OGT  All------------------------------------

tiff("./Figures/Figure 5.tiff",width = 20,height = 11, units = "in", res = 500)
ggplot(Calc_OGT,aes(y = as.numeric(mlm_best),x = as.numeric(Meas_OGT))) + 
  geom_point(aes(color = clade, shape = `OGT Source`), size = 3.5) + 
  geom_smooth(method = "lm", color = "black", size = 0.3) + 
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")), label.x = 5, label.y = 37, size = 7) +
  scale_color_manual("Taxonomic Clades",values = c("seagreen3","olivedrab4","darkorange2","purple", "black"), labels = c("Clade A1", "Clade A2", "Clade B","Clade C","Other Colwelliaceae")) +
  scale_shape_manual("OGT Sources",values = c(16,17,15),labels = c("OGT derived from Ratkowsky fitting of reported temperature growth data","OGT in Reported Literature","OGT from this Study")) +
  labs(x = "Fitted OGT (°C)", y = "Predicted OGT (°C)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 13),
         axis.text = element_text(size = 17), axis.title = element_text(size = 19),
         strip.text.x = element_text(size = 17)) +
    annotate("text",x = 10.3, y = 35.5,label = paste0("RMSE = ",RMSE,"(°C)"), size = 7)
dev.off()

#------------------------------------------ Plot predicted OGT (best combination) vs. measured OGT  Clade A------------------------------------
Calc_OGTA = filter(Calc_OGT, clade == "Clade A1" | clade == "Clade A2")
OGT_diffA = filter(Calc_OGT, clade == "Clade A1" | clade == "Clade A2")

RMSE = signif(sqrt(mean(OGT_diffA$mlm_best)^2),2)

ggplot(Calc_OGTA,aes(y = as.numeric(mlm_best),x = as.numeric(Meas_OGT))) + 
  geom_point(aes(color = clade, shape = clade), size = 3.5) + 
  geom_smooth(method = "lm", color = "black", size = 0.3) + 
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~"))) +
  scale_color_manual("Taxonomic Clades",values = c("seagreen3","olivedrab4","darkorange2","purple", "black"), labels = c("Clade A1", "Clade A2", "Clade B","Clade C","Other Colwelliaceae")) +
  scale_shape_manual("Taxonomic Clades",values = c(15,16,17,18,8,9),labels = c("Clade A1", "Clade A2","Clade B", "Clade C", "Other Colwelliaceae")) +
  labs(x = "Fitted OGT (°C)", y = "Predicted OGT (°C)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        axis.text = element_text(size = 17), axis.title = element_text(size = 19),
        strip.text.x = element_text(size = 17)) +
  annotate("text",x = 10.3, y = 35.5,label = paste0("RMSE = ",RMSE,"(°C)"), size = 7)

#------------------------------------------ Plot predicted OGT (best combination) vs. measured OGT  Clade B------------------------------------
Calc_OGTB = filter(Calc_OGT, clade == "Clade B")
OGT_diffB = filter(Calc_OGT, clade == "Clade B")

RMSE = signif(sqrt(mean(OGT_diffB$mlm_best)^2),2)

ggplot(Calc_OGTB,aes(y = as.numeric(mlm_best),x = as.numeric(Meas_OGT))) + 
  geom_point(aes(color = clade, shape = clade), size = 3.5) + 
  geom_smooth(method = "lm", color = "black", size = 0.3) + 
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~"))) +
  scale_color_manual("Taxonomic Clades",values = c("seagreen3","olivedrab4","darkorange2","purple", "black"), labels = c("Clade A1", "Clade A2", "Clade B","Clade C","Other Colwelliaceae")) +
  scale_shape_manual("Taxonomic Clades",values = c(15,16,17,18,8,9),labels = c("Clade A1", "Clade A2","Clade B", "Clade C", "Other Colwelliaceae")) +
  labs(x = "Fitted OGT (°C)", y = "Predicted OGT (°C)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        axis.text = element_text(size = 17), axis.title = element_text(size = 19),
        strip.text.x = element_text(size = 17)) +
  annotate("text",label = paste0("RMSE = ",RMSE,"(°C)"), size = 7)

#------------------------------------------ Plot predicted OGT (best combination) vs. measured OGT  Clade B------------------------------------
Calc_OGTC = filter(Calc_OGT, clade == "Clade C")
OGT_diffC = filter(Calc_OGT, clade == "Clade C")

RMSE = signif(sqrt(mean(OGT_diffC$mlm_best)^2),2)

ggplot(Calc_OGTC,aes(y = as.numeric(mlm_best),x = as.numeric(Meas_OGT))) + 
  geom_point(aes(color = clade, shape = clade), size = 3.5) + 
  geom_smooth(method = "lm", color = "black", size = 0.3) + 
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~"))) +
  scale_color_manual("Taxonomic Clades",values = c("seagreen3","olivedrab4","darkorange2","purple", "black"), labels = c("Clade A1", "Clade A2", "Clade B","Clade C","Other Colwelliaceae")) +
  scale_shape_manual("Taxonomic Clades",values = c(15,16,17,18,8,9),labels = c("Clade A1", "Clade A2","Clade B", "Clade C", "Other Colwelliaceae")) +
  labs(x = "Fitted OGT (°C)", y = "Predicted OGT (°C)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        axis.text = element_text(size = 17), axis.title = element_text(size = 19),
        strip.text.x = element_text(size = 17)) +
  annotate("text",label = paste0("RMSE = ",RMSE,"(°C)"), size = 7)

