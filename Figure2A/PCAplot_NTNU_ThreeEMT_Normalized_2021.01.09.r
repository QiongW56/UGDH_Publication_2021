# Copy FJfunctions folder to R working directory.
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019")
source("./FJfunctions/loadAll.r")
loadAll()

# load data for plotting
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults")
library(readxl)
dat <- read_excel("1 NTNU_proteomics_six cell lines_no gene names_28.12.2018.xlsx", 
                  sheet = 3)
dat <- data.frame(dat)

# delete the useless columns
dat <- dat[, c(1, 33:50)]

# transverse the dataset
datRev <- t(dat)
datRev <- data.frame(datRev)

# change column names
names(datRev) <- datRev[1, ]
datRev <- datRev[-1, ]

# add "Sample.ID" column
datRev$Sample.ID <- substr(rownames(datRev), (nchar("Normalized_")+1), nchar(rownames(datRev)))

datRev$Class <- c(rep("D492", 3), rep("D492M", 3),
                  rep("HMLE", 3), rep("HMLEM", 3),
                  rep("PMC42LA", 3), rep("PMC42ET", 3))
datRev$Cell <-  c(rep("D492", 6), rep("HMLE", 6), rep("PMC42", 6))

datRev <- datRev[, c((ncol(datRev)-2):ncol(datRev), 1:(ncol(datRev)-3))]

# factor the first three columns
datRev[, 1] <- factor(datRev[, 1], levels = unique(datRev[, 1]))
datRev[, 2] <- factor(datRev[, 2], levels = unique(datRev[, 2]))
datRev[, 3] <- factor(datRev[, 3], levels = unique(datRev[, 3]))

# numeric the expression data
for (i in 4:ncol(datRev)){
  datRev[, i] <- as.numeric(datRev[, i])
}

col1 <- c("steelblue", "red", "orange1")
col2 <- c("steelblue", "skyblue1", "red", "pink", "orange1", "#FAE48BFF")

# For PCA:
p <- metPCA(datRev,
            "Cell",
            minvalue = T, 
            log = T, 
            colors = col1,
            pointsize = 10,
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

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
file <- paste0(format(Sys.time(), "%F %H-%M-%S"), " ", "PCAplot.tiff")
ggsave(filename = file, units = "in", dpi = 300, width = 9, height = 6)
