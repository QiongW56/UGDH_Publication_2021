# rawdata is located at: 
# "C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1"
funMean <- function(Metabo, SheetName, n){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1")
  library(readxl)
  files <- paste0("BarPlotting_", Metabo, ".xlsx")
  dat <- data.frame(read_excel(files, sheet = SheetName))
  
  met <- colnames(dat)[2]
  colnames(dat)[2] <- "Metabolite"
  
  ScrAvg <- mean(dat[1:n, 2], na.rm = TRUE)
  dat[, 2] <- dat[, 2]/ScrAvg
  
  dat$Treatment <- c(rep("Scramble", n), rep("siUGDH_1", n), rep("siUGDH_2", n))
  dat$Treatment <- factor(dat$Treatment, levels = c("Scramble", "siUGDH_1", "siUGDH_2"))
  
  # calculate mean and sd
  library(tidyverse)
  dat <- dat %>%
    group_by(Treatment) %>% 
    summarize(avg = mean(Metabolite), sd = sd(Metabolite))
  
  dat <- data.frame(dat)
  
  dat$CellType <- SheetName
  
  dat$Group <- paste0(dat$CellType, "_", dat$Treatment)
  
  colnames(dat)[2] <- met
  
  return(dat)
  
}

datD <- funMean("UDP-Glucose", "D492M", 3)
datH <- funMean("UDP-Glucose", "HMLEM", 3)
datP <- funMean("UDP-Glucose", "PMC42ET", 3)

dat <- rbind(datD, datH, datP)

dat$Treatment <- factor(dat$Treatment, levels = unique(dat$Treatment))
dat$CellType <- factor(dat$CellType, levels = unique(dat$CellType))
dat$Group <- factor(dat$Group, levels = unique(dat$Group))

# plotting
library(ggplot2)
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

# color changes depending on the cell type
coltype <-  c("steelblue4", "steelblue", "skyblue1")

met <- colnames(dat)[2]

# start plotting
p <- ggplot(dat, 
            aes(x = CellType, 
                y = dat[, 2], 
                fill = Treatment)) + # linetype = 1, 2...
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           color = "black", 
           width = 0.75, 
           size = 4) +
  geom_errorbar(aes(ymin = dat[, 2] - sd, ymax = dat[, 2] + sd), 
                width = 0.2,
                size = 4,
                position = position_dodge(0.75)) + 
  scale_fill_manual(values = coltype) + 
  scale_y_continuous(breaks = seq(0, round(max(dat[, 2])+max(dat[, 3])), round(max(dat[, 2])+max(dat[, 3]))/5),
                     labels = scales::number_format(accuracy = 0.01,
                                                    decimal.mark = ".")) +
  theme_classic() +
  labs(title = met, 
       x = NULL, 
       y = "KD/Scr Ratio") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 108), 
        axis.line = element_line(colour = "black", size = 4),
        legend.position = "right",
        legend.key.size = unit(2.5, "cm"),
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        axis.title.y.left = element_text(size = 96, face = "bold", 
                                         margin = margin(t = 0, r = 50, b = 0, l = 0)),
        axis.text = element_text(size = 96, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 4)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 2)

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/20 UGDH-Paper/3-siUGDH-D492M-HMLEM-PMC42ET/1-Metabolomics-Intracellular-BOX1/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), met, "KD_RT-qPCR", "bar.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 36, height = 18)

