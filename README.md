# Colwelliaceae-OGT-Genomic-Signatures
Supporting documentation for publication. The documents included in this repository will allow you to recreate the analysis and figures in the paper that were made in R. 

Figure 1. Phylogenetic Tree.
Folder: Phylogenetics
R Script: Figure 1.R
Input Files: 
-	./Patric Input/Colwelliaceae Rooted Tree.nwk
-	./Refs& Raw Input/growthrates.tsv
-	./OGT/Strain Ratk OGT Final.tsv
-	./Patric Input/Colwelliaceae Rooted PATRIC Data.tsv

Figure 2. Growth Curve by Clade
Folder: Phylogenetics
R Script: Figure 2.R
Input Files: 
-	./Phylogenetics/Clade Ratk Params.tsv
-	./Phylogenetics/Clade Ratk Avg Params.tsv
-	./Phylogenetics/Clade Ratk Med Params.tsv
-	./Phylogenetics/Clade Ratk SD Params.tsv
-	./Phylogenetics/Sub-Clade ID.tsv
-	./Refs& Raw Input/growthrates.tsv

Figure 3. Boxplot of OGT by Clade
Folder: Phylogenetics
R Script: Figure 3.R
Input Files: 
-	./Refs& Raw Input/growthrates.tsv
-	./OGT/Strain Ratk OGT Final.tsv
-	./Phylogenetics/Genera Clade ID.tsv
-	./Phylogenetics/Sub-Clade ID.tsv

Figure 4. Amino Acid Indices Correlation to OGT
Folder: Amino Acid Indices
R Script: Figure 4.R
Input Files: 
-	./Patric Input/Colwelliaceae Rooted Patric Data.tsv
-	./Amino Acid Indices/WG Indice Data.tsv
-	./Phylogenetics/Sub-Clade ID.tsv
-	./OGT/Strain Ratk OGT Final.tsv

Figure 5. Fitted OGT vs. Predicted OGT
Folder: Amino Acid Indices
R Script: Figure 4.R
Input Files: 
-	./Amino Acid Indices/Genome AA model.tsv
-	./Amino Acid Indices/Genome AA model diff.tsv
-	./OGT/Strain Ratk OGT Final.tsv

Supplemental Table 1 - Genome sequence summary for Colwelliaceae strains.
-	No Code
Supplemental Table 2 - Phylogenetic cluster summary for novel Colwelliaceae strains
-	No Code
Supplemental Table 3 - Temperature-Dependent Growth of Colwelliaceae strains isolated in this
Folder: Supplementary Files&Figures
R Script: Supplemental Table 3.R
Input: 
-	./OGT/Strain Ratk OGT Final.tsv
-	./Supplementary Files&Figures/OGT MAD.tsv
-	./Supplementary Files&Figures/TempGrowthrates.tsv
Supplemental Table 4 - Summary of genes associated with cold adaptation found in Colwelliaceae through a literature compilation.
-	No Code
Supplemental Table 5 - Differential cold-adaptive protein families present in all Clade A2 strains but absent in among proteins present in Colwellia Bg11-12 or all strains of Clade A1
Folder: Supplementary Files&Figures
R Script: Supplemental Table 5.R
Input: 
-	./Patric Input/Colwelliaceae Rooted Tree.nwk
-	./Patric Input/Colwelliaceae UnRooted PATRIC Data.tsv
-	./Patric Input/Colwelliaceae PF Features.tsv
-	./Patric Input/Colwelliaceae_families_global.tsv
-	./Phylogenetics/Sub-Clade ID.tsv
-	./Amino Acid Indices/Gene Prescence.csv
-	./OGT/Strain Ratk OGT Final.tsv
-	./Supplementary Files&Figures/blast_CA_Colwelliaceae.txt

Supplementary Figure 1 - Number of gene families present in the Colwelliaceae core genome (black line) vs. pangenome (red line)
Folder: Supplementary Files&Figures
R Script: Supplemental Figure 1.R
Input : 
-	./Supplementary Files&Figures/f-pangenome.R
-	./Patric Input/Colwelliaceae Rooted Tree.nwk
-	./Patric Input/Colwelliaceae Rooted PATRIC Data.tsv
-	./Patric Input/Colwelliaceae PF Features.tsv
-	./Patric Input/Colwelliaceae_families_local.tsv
-	./Patric Input/Colwelliaceae_families_global.tsv
-	./Supplementary Files&Figures/Protein Indices Data.tsv
-	./Supplementary Files&Figures/Colwelliaceae PF AA Calc.RDS
-	./OGT/Strain Ratk OGT Final.tsv
-	./Phylogenetics/Sub-Clade ID.tsv
Supplementary Figure 2 - Growth rate versus optimal growth temperature of Colwelliaceae
Folder: Supplementary Files&Figures
R Script: Supplemental Figure 1.R
Input : 
-	./Refs& Raw Input/growthrates.tsv
-	./OGT/Strain Ratk OGT Final.tsv
-	./Phylogenetics/Sub-Clade ID.tsv
