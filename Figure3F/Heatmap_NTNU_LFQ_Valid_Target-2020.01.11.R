# upload data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Perseus")
datPerZ <- read.delim("PerseusOutputData-NTNU-Z-Score-1.6.14.0-2021.01.11_Grouping.txt")
datPerZ <- data.frame(datPerZ)

# Delete the first two rows: type, group
datPerZ <- datPerZ[-1:-2, ]

datPerZ <- datPerZ[, c(48, 1:18)]

for (i in c(2:19)){
  datPerZ[, i] <- as.numeric(datPerZ[, i])
}

# delete the duplicated genes
datPerZ <- datPerZ[-which(duplicated(datPerZ[, 1])), ]
rownames(datPerZ) <- datPerZ[, 1]

# upload the valid genes for plotting
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
library(readxl)
datVal <- read_excel("NTNU_ThreeEMT_Targets_2021.01.10.xlsx", sheet = 2)
datVal <- data.frame(datVal)

# Select only the valid genes from Perseus output Z-score dataset
GNval <- datVal$Gene.Name

# Added two enzymes which are significant in R-pvalues: FDFT1 and GANAB
GNval <- c(GNval, "FDFT1") # did not use "GANAB"

datZval <- datPerZ[which(datPerZ$Gene.Name %in% GNval), ]

# Change the column names
colnames(datZval) <- gsub("\\.", "", colnames(datZval))
colnames(datZval) <- gsub("42_", "42", colnames(datZval))

# split the data based on EMT models
datD <- datZval[, 1:7]
datH <- datZval[, c(1, 8:13)]
datP <- datZval[, c(1, 14:19)]

datD_matrix <- as.matrix(datD[, -1])
datH_matrix <- as.matrix(datH[, -1])
datP_matrix <- as.matrix(datP[, -1])

# Heatmap
library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "NTNU_ThreeEMT_HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 12, height = 10, res = 300)

col_fun_1 = circlize::colorRamp2(c(-1.88, 0, 1.77), c("steelblue", "white", "red"))
ht_1 <- Heatmap(datD_matrix, 
                col = col_fun_1,
                heatmap_legend_param = list(title = "",
                                            title_gp = gpar(fontsize = 30),
                                            labels_gp = gpar(fontsize = 30),
                                            at = c(-1.88, 0, 1.77),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(5, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 2,
                # combined_name_fun = NULL,
                row_title = NULL,
                row_title_gp = gpar(fontsize = 30),
                column_title = "D492",
                column_title_gp = gpar(fontsize = 30))

col_fun_2 = circlize::colorRamp2(c(-1.15, 0, 1.16), c("steelblue", "white", "red"))
ht_2 <- Heatmap(datH_matrix, 
                col = col_fun_2,
                heatmap_legend_param = list(title = "",
                                            title_gp = gpar(fontsize = 30),
                                            labels_gp = gpar(fontsize = 30),
                                            at = c(-1.15, 0, 1.16),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 2,
                # combined_name_fun = NULL,
                row_title_gp = gpar(fontsize = 30),
                column_title = "HMLE",
                column_title_gp = gpar(fontsize = 30),
                show_heatmap_legend = F,
                show_row_names = F)

col_fun_3 = circlize::colorRamp2(c(-1.15, 0, 1.16), c("steelblue", "white", "red"))
ht_3 <- Heatmap(datP_matrix, 
                col = col_fun_3,
                heatmap_legend_param = list(title = "",
                                            title_gp = gpar(fontsize = 30),
                                            labels_gp = gpar(fontsize = 30),
                                            at = c(-1.15, 0, 1.16),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 2,
                # combined_name_fun = NULL,
                row_title_gp = gpar(fontsize = 30),
                column_title = "PMC42",
                column_title_gp = gpar(fontsize = 30),
                show_heatmap_legend = F,
                show_row_names = F)

ht_list <- ht_1 + ht_2 + ht_3
draw(ht_list, ht_gap = unit(5, "mm"))

dev.off() # export the heatmap plot into .tiff


