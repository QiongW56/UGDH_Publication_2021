# UGDH_Publication_2021
A description of the figures in this "UGDH,2021" paper.

## Figure 1
Figures were generated in Excel and Powerpoints.

## Figure 2
Data used:
### Fig.2.A&B
“C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults”
"1 NTNU_proteomics_six cell lines_no gene names_28.12.2018.xlsx"

### Fig.2.C&D (CDH1&2)
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH”
“RNA-Expression-CDH1-CDH2-UGDH.xlsx”

The location of the original data:
"C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\4 Experiments documents_QIONG\20 UGDH-Paper\1-WideType-D492-HMLE-PMC42"

### Fig.2.E-G
"C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults"
"2 NTNU_Proteomics_D492 vs. D492M_28.12.2018.xlsx"
"2 NTNU_Proteomics_HMLE vs. HMLE_M_28.12.2018.xlsx"
"2 NTNU_Proteomics_PMC42_LA vs. PMC42_ET_28.12.2018.xlsx"

### Fig.2.H (PKP3)
C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH
“NTNU-LFQ-Expression-CTGF-FDFT1-PKP3.xlsx” 

The location of the original data:
"C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Perseus”
“PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx"

Why did PKP3 and FDFT1 (did not use CTGF) be plotted separately?
Because when using “BarPlotPlotting_NTNU_ThreeEMT_2021.01.09.R” to plot, it uses these data:
"C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults"
"2 NTNU_Proteomics_D492 vs. D492M_28.12.2018.xlsx"
"2 NTNU_Proteomics_HMLE vs. HMLE_M_28.12.2018.xlsx"
"2 NTNU_Proteomics_PMC42_LA vs. PMC42_ET_28.12.2018.xlsx"

There are missing values in these data, so I have to use the Perseus output data.

Figure generated:
### Fig.2.A&B 
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”
PCA:
“PCAplot_NTNU_ThreeEMT_Normalized_2021.01.09.R”
Dendrogram:
“Dendrogram_NTNU_ThreeEMT_Imputation_2021.01.09.R”

### Fig.2.C&D (CDH1&2)
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”
“BarPlot_RNA-Expression-CDH1-CDH2-UGDH-2020.01.09.R”

### Fig.2.E-G
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”
“BarPlotPlotting_NTNU_ThreeEMT_2021.01.09.R”

### Fig.2.H (PKP3)
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”
“BarPlot_NTNU-LFQ-Expression-CTGF-FDFT1-PKP3-2020.01.13.R”

All figures:
PCAplot: C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Figures-PCA-Plot
Dendrogram - “C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH”
Barplots - “C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Figures-EMTmarkers\Fig.1”

Note:
1.	Used “EMT-Markers-NTNU-ThreeEMT.R” to find the EMT targets. (first analysis)
2.	Realized the data in “EMT-Markers-NTNU-ThreeEMT.R” is not good.
3.	Used “EMT-Markers-NTNU-ThreeEMT_V2.R” to find again the EMT targets
4.	Did not added the missing EMT markers from the first analysis because they did not pass the FDR threshold in Perseus. (For the consistent markers)
5.	Added the missing EMT markers for the inconsistent markers since they have passed the FDR cutoff in Perseus. - EGFR, S100A2, NDRG1

Even though all bar plots for proteomic data had p values, the “*” marks were calculated separately.
