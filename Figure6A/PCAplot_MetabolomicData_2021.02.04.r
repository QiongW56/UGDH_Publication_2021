# load the AcidicNeg, AcidicPos and BasicNeg data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/PCAplotting_Data")
library(readxl)
dat <- read_excel("MetabolomicData_For_PCA.xlsx", sheet = 1) # sheet 1 is mine and sheet 2 is from Freyr
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

datD <- dat[which(dat$CellType == "D492M"), ]
datH <- dat[which(dat$CellType == "HMLEM"), ]
datP <- dat[which(dat$CellType == "PMC42ET"), ]

# Copy FJfunctions folder to R working directory.
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019")
source("./FJfunctions/metPCA.r")

# For PCA:
p <- metPCA(dat,
            "CellType",
            minvalue = T, 
            log = T, 
            colors = c("steelblue", "red", "orange1"),
            pointsize = 5,
            CI = TRUE)

library(ggplot2)
library(ggrepel)
p1 <- p +
  # geom_text_repel(aes(label = tdata[, 3]), fontface = "bold") +
  labs(title = "") +
  theme_classic() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "right",
        title = element_text(""),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 30, face = "bold", color = "black"),
        axis.text.x.bottom = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y.left = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.line = element_line(size = 2))
p1

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/PCAplotting_Data")
file <- paste0(format(Sys.time(), "%F %H-%M-%S"), " ", "PCAplot.tiff")
ggsave(filename = file, units = "in", dpi = 300, width = 9, height = 6)
