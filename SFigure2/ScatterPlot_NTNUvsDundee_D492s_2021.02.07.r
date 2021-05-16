# This file is to compare Perseus with R output in terms of significant targets
# this file is also to compare D492 VS. D492M for NTNU and Dundee

# -------------------------------------------------------------------- # 

# Perseus outputs VS. R outputs

# -------------------------------------------------------------------- #

# load the Perseus output data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Perseus")
library(readxl)
datPerseus <- read_excel("PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx", sheet = 1)
datPerseus <- data.frame(datPerseus) # 873

# 188
datPerD <- datPerseus[which(datPerseus$Student.s.T.test.Significant.D492_D492M == "+"), ]

# load the R output data for D492 vs. D492M
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults")
datRD <- read_excel("2-1 2018-12-28 17-51 Two.groups.statistics.LFQ.xlsx", sheet = "valid")
datRD <- data.frame(datRD) # 866

datRDval <- datRD[which(datRD$q.value < 0.05), ] # 191
  
length(which(datPerD$Gene.Name %in% datRDval$Gene.Name)) # 185

`%!in%` <- Negate(`%in%`)
datPerD$Gene.Name[which(datPerD$Gene.Name %!in% datRDval$Gene.Name)] # 3
# [1] "CRIP1" "CDH6"  "HMGA1"

datRDval$Gene.Name[which(datRDval$Gene.Name %!in% datPerD$Gene.Name)] # 6
# [1] "TPM1"    "PPP2R1A" "LAMC2"   "ASS1"    "RAB2A"   "MYADM" 

# Note: 
# the difference between Perseus and R is due to the imputation methods
# or, it is because, in R output files, the p value is not less than 0.05 even though FDR < 0.05.

# -------------------------------------------------------------------- # 

# NTNU vs. Dundee (D492 vs. D492M), We will use the Perseus output data for comparsion

# -------------------------------------------------------------------- #
# load the Perseus output data (NTNU)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Perseus")
library(readxl)
datPerseus <- read_excel("PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx", sheet = 1)
datPerseus <- data.frame(datPerseus) # 873

# 188 of genes are sigficantly different in NTNU-Perseus
datPerD <- datPerseus[which(datPerseus$Student.s.T.test.Significant.D492_D492M == "+"), ]

# load the Dundee data (LFQ)
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
datDun <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
datDun <- data.frame(datDun)

# load the Dundee data (SILAC)
datSILAC <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx")
datSILAC <- data.frame(datSILAC)

datDunVal <- datDun[which(datDun$EM_significant == "+"), ] # 1237
datSILval <- datSILAC[which(datSILAC$p.value < 0.05), ] # 3177

# number of genes in Perseus (significant) also detected in Dundee
length(which(datPerD$Gene.Name %in% datDun$Gene.Name..GN.)) # 175/2734

# number of genes in Perseus (significant) also in Dundee (significant)
length(which(datPerD$Gene.Name %in% datDunVal$Gene.Name..GN.)) # 158/188

`%!in%` <- Negate(`%in%`)
gnPer <- datPerD$Gene.Name[which(datPerD$Gene.Name %!in% datDunVal$Gene.Name..GN.)]

table(gnPer %in% datSILval$Gene.Name) # 18/30 not in SILAC (Dundee) either

gnPer1 <- gnPer[which(gnPer %in% datSILval$Gene.Name)]
# [1] "CAVIN1" "RBBP7"  "GARS"   "HDAC1"  "PPP3CB" "DLST"   "RPL34"  "DRAP1" 
# [9] "DHRS4"  "ASF1B"  "ARPC5"  "RAB6A" 

# (158+12)/188 = 90.4%

# ------------------------------------------------------# Dundee-LFQ
# Check the direction of the data (same up or downs)

datPerD1 <- datPerD[which(datPerD$Gene.Name %in% datDunVal$Gene.Name..GN.), c(3, 40)] # 158

datDunVal <- datDunVal[which(datDunVal$Gene.Name..GN. %in% datPerD1$Gene.Name), c(3, 18)] # 161
datDunVal <- datDunVal[-which(duplicated(datDunVal$Gene.Name..GN.)), ] # 158

# change to the same colname for merging
colnames(datDunVal)[1] <- colnames(datPerD1)[1]

# merge two datasets (NTNU and Dundee-LFQ)
datMer <- merge(datPerD1, datDunVal) # 158

for (i in 1:nrow(datMer)){
  datMer$Compare[i] <- ifelse(datMer[i, 2] > 0 & datMer[i, 3] > 0, 1, 
                              ifelse(datMer[i, 2] < 0 & datMer[i, 3] < 0, -1, 0))
}

table(datMer$Compare == 0) # 155/158 are consistent
datMer[which(datMer$Compare == 0), ] # which genes are not consistent

# ------------------------------------------------------# Dundee-SILAC
# Check the direction of the data (same up or downs)

datPerD2 <- datPerD[which(datPerD$Gene.Name %in% datSILval$Gene.Name), c(3, 40)] # 114

datSILval <- datSILval[which(datSILval$Gene.Name %in% datPerD2$Gene.Name), c(3, 7)] # 122
datSILval <- datSILval[-which(duplicated(datSILval$Gene.Name)), ] # 114

# change to the same colname for merging
colnames(datSILval)[1] <- colnames(datPerD2)[1]

# merge two datasets (NTNU and Dundee-LFQ)
datMer1 <- merge(datPerD2, datSILval) # 114

for (i in 1:nrow(datMer1)){
  datMer1$Compare[i] <- ifelse(datMer1[i, 2] > 0 & datMer1[i, 3] > 0, 1, 
                              ifelse(datMer1[i, 2] < 0 & datMer1[i, 3] < 0, -1, 0))
}

table(datMer1$Compare == 0) # 108/114
gnSIL <- datMer1$Gene.Name[which(datMer1$Compare == 0)] # 6 of genes are inconsistent
# [1] "DLST"   "HDAC1"  "MACF1"  "MAP2K1" "PPP3CB" "RPL34" 
table(gnSIL %in% gnPer1) # 4 out of the 6 genes are in the 12 genes before
# so, 12-4=8

# mordified after consistency checking, (155+8)/188 = 86.7%

#-------------------------------------------------------------------#

# Scatter Plotting

#-------------------------------------------------------------------#
# Data for plotting
dat <- datMer[, -ncol(datMer)]
colnames(dat)[2:3] <- c("x", "y")

# ploting with pearson cor
cor <- cor.test(dat$x, dat$y, 
                method = "pearson", conf.level = 0.95)
cor1 <- as.numeric(format(cor[4], digits = 3))
cor2 <- paste("Pearson Correlation =", cor1, sep = " ")

xmin <- min(dat[, 2])
xmax <- max(dat[, 2])
xtk <- round((xmax - xmin)/5)

ymin <- min(dat[, 3])
ymax <- max(dat[, 3])
ytk <- round((ymax - ymin)/5)

xtitle <- "(Current Dataset)"
ytitle <- "(Previous Dataset)"

library(ggplot2)
p <- ggplot(dat, aes(x = x, y = y)) +
  geom_point(color = "steelblue", size = 6) +
  theme_classic() +
  scale_x_continuous(breaks = seq(xmin, xmax, xtk), 
                     labels = scales::number_format(accuracy = 0.01, 
                                                    decimal.mark = ".")) +
  scale_y_continuous(breaks = seq(ymin, ymax, ytk), 
                     labels = scales::number_format(accuracy = 0.01, 
                                                    decimal.mark = ".")) +
  labs(title = "",
       x = bquote("Log"[2]*.(xtitle)),
       y = bquote("Log"[2]*.(ytitle))) +
  geom_smooth(method = lm, color = "black", size = 3) +
  theme(legend.position = "none",
        axis.line = element_line(colour = 'black', size = 2),
        axis.title.y.left = element_text(size = 32, face = "bold",
                                         margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size = 32, face = "bold", color = "black"),
        axis.title.x.bottom = element_text(size = 32, face = "bold",
                                           margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.ticks.y.left = element_line(colour = "black", size = 2)) +
  geom_text(x = 0, y = 8, label = cor2, color = "black", fontface = "bold", size = 8)
# geom_text_repel(aes(x = x, y = y, 
# label = ifelse(Grouping == "YES", GeneName, "")), size = 3)
p

# save the plot
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-ScatterPlot")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "_", "DataValidation", "_","Scatter.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 8, height = 6)
