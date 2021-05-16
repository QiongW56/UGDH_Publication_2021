# UGDH_Publication_2021
A description of the figures in this "UGDH,2021" paper. <br/>
For details: <br/>
Proteomics_UGDH_2021_Manuscript_FigureLegends_V3.docx

## Figure 1
Study workflow and a summary of the three EMT breast cell models. <br/>
Generated in Excel and powerpoint.

D492: <br/>
Isolation, immortalization, and characterization of a human breast epithelial cell line with stem cell properties.<br/>
D492M: <br/>
Endothelial Induced EMT in Breast Epithelial Cells with Stem Cell Properties<br/>
ECM1 secreted by HER2-overexpressing breast cancer cells promotes formation of a vascular niche accelerating cancer cell migration and invasion. (non-tumorigenic)<br/>

HMLE:<br/>
ISOLATION AND GROWTH OF HUMAN MAMMARY EPITHELIAL CELLS (Stampfer 1985)<br/>
Human breast cancer cells generated by oncogenic transformation of primary mammary epithelial cells.<br/>
HMLEM:<br/>
Protein Kinase C α Is a Central Signaling Node and Therapeutic Target for Breast Cancer Stem Cells.<br/>

PMC42ET:<br/>
A New Human Breast Carcinoma Cell Line (PMC42) With Stem Cell Characteristics. I. Morphologic Characterization.<br/>
PMC42LA:<br/>
PMC42, A Novel Model for the Differentiated Human Breast<br/>

## Figure 2A
PCA plot of the proteomic dataset

## Figure 2B
Dendrogram plot of the proteomic dataset

## Figure 2C-D
RNA expression of CDH1 and CDH2

## Figure 2E-G
Proteomic expression of the consistent changes in EMT<br/>
http://dbemt.bioinfo-minzhao.org/

## Figure 2H
Proteomic expression of the consistent changes in EMT- PKP3

## Figure 3A-B
Analysis of the proteomic datasets. <br/>
Generated in Excel and powerpoint.

## Figure 3C-E
Volcano plots of the proteomic datasets <br/>

## Figure 3F
Heatmap of the consistent proteomic changes.

## Figure 4A-C
GO annotation - BP <br/>
DAVID Bioinformatics Resources 6.8<br/>
Date: 2021.02.09<br/>

Procedures:<br/>
1. Go to Perseus Output data: 
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Perseus\”
“PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx”
2. Find genes with significance "+" for each EMT model.
3. Copy these genes into DAVID for "Functional Annotation Chart" analysis, for GO-BP, GO-CC and GO-MF separately
4. Export the output and copy the output .txt file into an Excel file.
5. Use R file for plotting: <br/>
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”<br/>
“BarPlot-DoxPlot-GO-Annotation-2021.02.08.R”

## Figure 4D-E
Reactome pathway analysis<br/>
Reactome Version:<br/>
Reactome Database Release: 75<br/>
Pathway Browser Version: 3.7<br/>
Date: 2021.02.09<br/>

Procedures:<br/>
1. Go to Perseus Output data: 
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Perseus\”
“PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx”
2. Find genes with significance "+" for each EMT model.
3. Copy these genes into Reactome for pathway enrichment.
4. Export the output .csv file
5. Use R file for plotting: <br/>
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles\”<br/>
“TreePlot_Reactome_Pathway_Enrichment_2021.02.09.R”<br/>

## Figure 5A-D
Consistent metabolic EMT markers - proteomics

## Figure 5E
RNA expression of UGDH

## Figure 5F
UGDH proteomic expression confirmed in another study

## Figure 6A
PCA plot of the metabolomic dataset

## Figure 6B
Dendrogram plot of the metabolomic dataset

## Figure 6C-E
Metabolite levels <br/>
UDP-Glc: BASIC mode - Phenylalanine IS; Lysine IS; Adenine IS <br/>
UDP-GlcA: BASIC mode - Cysteine IS; Glutamine IS; Glutamic acid IS <br/>
Glycerophosphocholine: Acidic negative mode - Lysine IS.

## Figure 6F
Glycerophosphhocholine level in D492, D492M and D492HER2 <br/>
Acidic negative mode - Glutamic acid IS. <br/>
Generated in Excel and powerpoint.

## Figure 6G
IPA analysis of the phosphoproteoimc datasets

## Figure 7A
Patient survival analysis of UGDH <br/>
Generated in https://kmplot.com/analysis/index.php?p=service&cancer=breast.

Go to website: https://kmplot.com/analysis/index.php?p=service <br/>
input gene UGDH<br/>
Auto select best cutoff<br/>
intrinsic subtype: basal<br/>


## Figure 7B-C
Proliferation - siUGDH (only one siRNA)

## Figure 7D-E
Invasion - siUGDH (first siRNA) <br/>
Generated in Excel and powerpoint.

## Figure 7F-G
SNAI1 - siUGDH (first siRNA)

## Figure 8A-B
Generated in Excel and powerpoint.

## Figure 8C
Motif enrichment analysis from Perseus version 1.6.14.0

## Figure 8D-F
siPDGFRB in D492M (only one siRNA)

## Figure 8G-H
siRELA in D492M with the first siRNA

## Supplementary figure 1
Cell images

## Supplementary figure 2
Scatter plot of the datasets consistency <br/>
Dundee vs. NTNU (D492 vs. D492M)

## Supplementary figure 3
Inconsistent EMT markers; http://dbemt.bioinfo-minzhao.org/

## Supplementary figure 4
GO annotation - CC and MF<br/>

DAVID Bioinformatics Resources 6.8<br/>
Date: 2021.02.09<br/>

Procedures: <br/>
1. Go to Perseus Output data: <br/>
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Perseus\”
“PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx”
2. Find genes with significance "+" for each EMT model.
3. Copy these genes into DAVID for "Functional Annotation Chart" analysis, for GO-BP, GO-CC and GO-MF separately
4. Export the output and copy the output .txt file into an Excel file.
5. Use R file for plotting: <br/>
“C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Rfiles”<br/>
“BarPlot-DoxPlot-GO-Annotation-2021.02.08.R”<br/>

## Supplementary figure 5A
KD efficiency of UGDH in the metabolomic study

## Supplementary figure 5B-E
KD efficiency of UGDH in the other experiments

## Supplementary figure 6A-C
Patient survival analysis <br/>
Generated in https://kmplot.com/analysis/index.php?p=service&cancer=breast.

## Supplementary figure 6D
GPC level changes with KD in D492HER2 <br/>
Generated in Excel and powerpoint.

## Supplementary figure 6E-F
Invasion - siUGDH (second siRNA) <br/>
Generated in Excel and powerpoint.

## Supplementary figure 6G-H
SNAI1- siUGDH (second siRNA)

## Supplementary figure 7A
Generated in Excel and powerpoint.

## Supplementary figure 7B-D
siPDGFRB in D492HER2 (only one siRNA)

## Supplementary figure 7E-F
siRELA in D492M with the second siRNA

## Supplementary figure 7G-J
siRELA in D492HER2 (with two siRNAs)






