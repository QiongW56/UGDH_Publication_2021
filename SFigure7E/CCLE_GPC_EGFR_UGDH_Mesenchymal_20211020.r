# upload the epithelial and mesenchymal cell lines
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/UGDH-GPC-NAA_PublicData")
library(readxl)
datBC <- data.frame(read_excel("mmc1.xlsx", sheet = 2))
datBC$Cell.line <- gsub("-", "", datBC$Cell.line)

# upload the GPC data
datGPC <- data.frame(read_excel("41591_2019_404_MOESM2_ESM.xlsx", sheet = "1-clean data"))

grep("glycerophosphocholine", names(datGPC))
datGPC <- datGPC[, c(1, 101)]
datGPC <- cbind(datGPC, datGPC[, 1])

colnames(datGPC)[1] <- "Cell.line"
colnames(datGPC)[3] <- "Cell.line_cell.type"

# z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

datGPC$GPClevel <- 10^(datGPC$alpha.glycerophosphocholine)/100000
# datGPC$GPClevel <- cal_z_score(datGPC$GPClevel)

# delete the cell types
datGPC$Cell.line <- sub("_[^_]+$", "", datGPC$Cell.line) # gsub('(.*)_\\w+', '\\1', datGPC$Cell.line)
datGPC$Cell.line <- sub("_[^_]+$", "", datGPC$Cell.line)
datGPC$Cell.line <- sub("_[^_]+$", "", datGPC$Cell.line)
datGPC$Cell.line <- sub("_[^_]+$", "", datGPC$Cell.line)

dat <- merge(datBC, datGPC)

# Only epithelial and mesenchymal groups
emt <- c("Epithelial", "Mesenchymal")
dat <- dat[which(dat$Group %in% emt), ]

dat$Group <- factor(dat$Group)

aggregate(dat[, 7], list(dat$Group), mean)
aggregate(dat[, 7], list(dat$Group), sd)

table(dat$Group) 
# only 89 mesenchymal, 110 were found by manually matching the cell line names
# in "1_BC_18_GPC_NAA_EGFR_UGDH_Oct2021.xlsx"

#--------------------------------------------------------------------------#

# correlation to EGFR or UGDH

#--------------------------------------------------------------------------#
# or can change to "UGDH Expression 21Q3 Public.csv"
EGFR <- read.csv("EGFR Expression 21Q3 Public.csv") 

EGFR$Expression.21Q3.Public <- cal_z_score(EGFR$Expression.21Q3.Public)

names(EGFR)[3] <- "Cell.line"

dat <- merge(dat, EGFR)

cor_temp_p <- cor.test(as.numeric(dat[, 7]), 
                       as.numeric(dat[, 9]),
                       method = 'spearman')[3]
cor_temp_p

cor_temp_val <- cor.test(as.numeric(dat[, 7]),
                         as.numeric(dat[, 9]),
                         method = 'spearman')$estimate
cor_temp_val

#----------------------------------------------------------------------#

# Plot two images (using ggplot) and combine into a panel:

#----------------------------------------------------------------------#
# if plotting EGFR/UGDH - GPC
dat$Expression.21Q3.Public <- cal_z_score(dat$Expression.21Q3.Public)

dat$EGFR_grouped <- ifelse(dat$Expression.21Q3.Public < 0,'EGFR_Low','EGFR_High')
dat$EGFR_grouped <- as.factor(dat$EGFR_grouped)
dat <- dat[, c(7, 13)]

#--------------------------------------------------------------------#

# Manually picked the mesenchymal cell lines

#--------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/UGDH-GPC-NAA_PublicData")
library(readxl)
datEMT <- read_excel("1_BC_18_GPC_NAA_EGFR_UGDH_Oct2021.xlsx", 
                     sheet = "AllCellTypes")
datEMT <- data.frame(datEMT)

datEMTGPC <- datEMT[, c(3, grep("glycerophosphocholine", names(datEMT)))]

# transfer log10 to the raw data
datEMTGPC[, 2] <- 10^(datEMTGPC[, 2])/100000

datEMTGPC[is.na(datEMTGPC$CellType), 1] <- "Others"
table(datEMTGPC$CellType) # 110

datEMTGPC$CellType <- as.factor(datEMTGPC$CellType)

# Student's T test
G1 <- as.matrix(datEMTGPC$alpha.glycerophosphocholine[grep("Mesenchymal", datEMTGPC$CellType)])
G2 <- as.matrix(datEMTGPC$alpha.glycerophosphocholine[grep("Others", datEMTGPC$CellType)])
p <- t.test(G1, G2, var.equal = T, alternative = "two.sided")[3]
p # 0.557

#--------------------------------------------------------------------#

# plotting

#--------------------------------------------------------------------#
library(ggplot2)
# library(siggitRausti)
library(ggpubr)

# function to change the space between keys in legends
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
GeomBar$draw_key = draw_key_polygon3

yna <- "Glycerophosphocholine"

dat <- datEMTGPC

p <- ggplot(data = dat, 
            aes(x = CellType, 
                y = alpha.glycerophosphocholine, 
                fill = CellType)) +
  stat_boxplot(geom = 'errorbar', 
               position = position_dodge(1),
               width = 0.5,
               size = 2) +
  geom_boxplot(position = position_dodge(1), 
               lwd = 3, 
               fatten = 1) +
  #             outlier.colour = "red", outlier.shape = 8, outlier.size = 4) +
  scale_fill_manual(values = c("red", "steelblue")) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(1),
               binwidth = 0,
               dotsize = 10) +
  scale_y_continuous(breaks = seq(0, 150, 20),
                     labels = scales::number_format(accuracy = 0.1, 
                                                    decimal.mark = ".")) +
  # scale_x_discrete(breaks = c("GFPT2_High","GFPT2_Low"),
  #                   labels = c(expression(GFPT2^High), expression(GFPT2^Low))) +
  theme_classic() +
  labs(title = NULL, # title = paste0("Correlation between ", colnames(dat)[1], " and ", colnames(dat)[2]), 
       x = NULL, 
       y = yna) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 28), 
        axis.line = element_line(colour = "black", size = 3),
        legend.position = "right",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.title.y.left = element_text(size = 30, face = "bold", 
                                         margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x.bottom = element_text(size = 30, face = "bold",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 25, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 2))
p

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/UGDH-GPC-NAA_PublicData/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
              # colnames(dat)[1],
              "CCLE_GPC_EMT.Box.plot.tiff", 
              sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 6)

print(p)
