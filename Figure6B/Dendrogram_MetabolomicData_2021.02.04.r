# load the AcidicNeg, AcidicPos and BasicNeg data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/PCAplotting_Data")
library(readxl)
dat <- read_excel("MetabolomicData_For_PCA.xlsx", sheet = 1)
dat <- data.frame(dat)

# factor two columns "Treatment" and "CellType"
dat$Treatment <- factor(dat$Treatment, levels = c("Scramble", "siUGDH_1", "siUGDH_2"))
dat$CellType <- factor(dat$CellType, levels = c("D492M", "HMLEM", "PMC42ET"))

rowScrD <- which(dat$CellType == "D492M" & dat$Treatment == "Scramble")
rowU1D <- which(dat$CellType == "D492M" & dat$Treatment == "siUGDH_1")
rowU2D <- which(dat$CellType == "D492M" & dat$Treatment == "siUGDH_2")

rowScrH <- which(dat$CellType == "HMLEM" & dat$Treatment == "Scramble")
rowU1H <- which(dat$CellType == "HMLEM" & dat$Treatment == "siUGDH_1")
rowU2H <- which(dat$CellType == "HMLEM" & dat$Treatment == "siUGDH_2")

rowScrP <- which(dat$CellType == "PMC42ET" & dat$Treatment == "Scramble")
rowU1P <- which(dat$CellType == "PMC42ET" & dat$Treatment == "siUGDH_1")
rowU2P <- which(dat$CellType == "PMC42ET" & dat$Treatment == "siUGDH_2")

ScrD <- mean(dat$CrystalViolet[rowScrD])
U1D <- mean(dat$CrystalViolet[rowU1D])
U2D <- mean(dat$CrystalViolet[rowU2D])

ScrH <- mean(dat$CrystalViolet[rowScrH])
U1H <- mean(dat$CrystalViolet[rowU1H])
U2H <- mean(dat$CrystalViolet[rowU2H])

ScrP <- mean(dat$CrystalViolet[rowScrP])
U1P <- mean(dat$CrystalViolet[rowU1P])
U2P <- mean(dat$CrystalViolet[rowU2P])

dat[rowScrD, -1:-4] <- dat[rowScrD, -1:-4]/ScrD
dat[rowU1D, -1:-4] <- dat[rowU1D, -1:-4]/U1D
dat[rowU2D, -1:-4] <- dat[rowU2D, -1:-4]/U2D

dat[rowScrH, -1:-4] <- dat[rowScrH, -1:-4]/ScrH
dat[rowU1H, -1:-4] <- dat[rowU1H, -1:-4]/U1H
dat[rowU2H, -1:-4] <- dat[rowU2H, -1:-4]/U2H

dat[rowScrP, -1:-4] <- dat[rowScrP, -1:-4]/ScrP
dat[rowU1P, -1:-4] <- dat[rowU1P, -1:-4]/U1P
dat[rowU2P, -1:-4] <- dat[rowU2P, -1:-4]/U2P

# delete two columns "Treatment" and "CellType"
dat <- dat[, -2:-4]

dat <- data.frame(t(dat))
colnames(dat) <- dat[1, ]
dat <- dat[-1, ]

#colnames(dat) <- c(paste0("D492M", "_", 1:18),
#                   paste0("HMLEM", "_", 1:18),
#                   paste0("PMC42ET", "_", 1:17))
#

colnames(dat) <- c(paste0("D492M", "_", "Scr", 1:6),
                   paste0("D492M", "_", "siUGDH_1_", 1:6),
                   paste0("D492M", "_", "siUGDH_2_", 1:6),
                   paste0("HMLEM", "_", "Scr", 1:6),
                   paste0("HMLEM", "_", "siUGDH_1_", 1:6),
                   paste0("HMLEM", "_", "siUGDH_2_", 1:6),
                   paste0("PMC42ET", "_", "Scr", 1:6),
                   paste0("PMC42ET", "_", "siUGDH_1_", 1:5),
                   paste0("PMC42ET", "_", "siUGDH_2_", 1:6))

dat$Metabolites <- rownames(dat)
dat <- dat[, c(ncol(dat), 1:(ncol(dat)-1))]

for (i in 2:ncol(dat)){
  dat[, i] <- as.numeric(dat[, i])
}

# z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat[, -1] <- t(apply(dat[, -1], 1, cal_z_score))


# function for dendrogram plotting
funDen <- function(file){
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  
  file.cor <- cor(file[, -1], method = "pearson") # a dataset with samples on both x-axis and y-axis
  d.cor <- as.dist(1 - file.cor) # a list of distance between samples in a triangular shape
  file.hclust <- hclust(d.cor, method = "ward.D2") # a list
  
  dend <- as.dendrogram(file.hclust) # not a plot, but a list, already have the hierarchical structure
  dend_data <- dendro_data(dend, type = "rectangle") # a list
  file.cut <- cutree(file.hclust, k = 3) # to color the labels into 4 groups
  file.cut.df <- data.frame(label = names(file.cut), 
                            cluster = factor(file.cut))
  dend_data[["labels"]] <- merge(dend_data[["labels"]], file.cut.df, by="label")
  
  # add extra rectangles
  # dat <- data.frame(x1 = c(7, 18), x2 = c(8, 19), y1 = c(1, 1), y2 = c(2, 2))
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/DendrogramPlotting")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "NTNU_ThreeEMT_dendrogram.plot.tiff", sep = " ")
  tiff(filename = file, units = "in", res = 300, width = 10, height = 5)
  
  col1 <- c("steelblue", "red", "orange1")
  dendrogram <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, 
                     y = y, 
                     xend = xend, 
                     yend = yend), 
                 size = 1.5) + # the boldness of the dendrogram
    geom_text(data = dend_data$labels, 
              aes(x, y, label = label, hjust = 0, color = cluster),
              angle = 270, 
              size = 4,
              fontface = "bold") +
    scale_colour_manual(values = col1) +
    ylim(-4, 6.5) +
    xlim(0, 60) +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5)) +
    # geom_segment(aes(x = 10.5, y = -0.55, xend = 13.2, yend = -0.55), size = 1) +
    # geom_segment(aes(x = 33.5, y = -0.45, xend = 36.2, yend = -0.45), size = 1) +
    ggtitle("Metabolome of The Mesenchymal Cells") +
    theme(plot.title = element_text(size = 20, face = "bold"))
  
  print(dendrogram)
  
  dev.off()
  
  print(dendrogram)
  
}

funDen(dat)
