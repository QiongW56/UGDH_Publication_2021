# upload data from NTNU
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults")
library(readxl)
dat <- read_excel("1 NTNU_proteomics_six cell lines_no gene names_28.12.2018.xlsx",
                  sheet = 1)
dat <- data.frame(dat)

# delete the useless columns
dat <- dat[, c(1, 33:50)]

# rename the column names
nchar("Normalized_") # 11
n <- nchar(colnames(dat)[2:ncol(dat)])
colnames(dat)[2:ncol(dat)] <- substr(colnames(dat)[2:ncol(dat)], 12, n)

colnames(dat) <- gsub("\\.", "", colnames(dat))
colnames(dat) <- gsub("42_", "42", colnames(dat))

# z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat[,-1] <- t(apply(dat[,-1], 1, cal_z_score))

# imputation
library(msm)
for (i in 2:ncol(dat)) {
  mean <- mean(dat[,i][-which(is.na(dat[,i]))])
  
  sd <- sd(dat[,i][-which(is.na(dat[,i]))])
  
  dat[,i][which(is.na(dat[,i]))] <- rtnorm(length(
    which(is.na(dat[,i]))), upper = mean-1.8*sd, lower = mean-1.8*sd-0.5*sd)
}

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
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "NTNU_ThreeEMT_dendrogram.plot.tiff", sep = " ")
  tiff(filename = file, units = "in", res = 300, width = 6, height = 5)
  
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
              size = 6,
              fontface = "bold") +
    scale_colour_manual(values = col1) +
    ylim(-1, 2) +
    xlim(0, 20) +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5)) +
    # geom_segment(aes(x = 10.5, y = -0.55, xend = 13.2, yend = -0.55), size = 1) +
    # geom_segment(aes(x = 33.5, y = -0.45, xend = 36.2, yend = -0.45), size = 1) +
    ggtitle("EMT Models") +
    theme(plot.title = element_text(size = 24, face = "bold"))
  
  print(dendrogram)
  
  dev.off()
  
  print(dendrogram)
}

# plotting
funDen(dat)
