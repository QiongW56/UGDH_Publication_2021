# This stript was adapted from "BarPlot-DoxPlot-GO-Annotation-2021.02.08.R"
# ------------------------------------------

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-Signaling/Phosphoproteomics/Perseus")
library(readxl)
dat <- read_excel("3 Motifs enrichment_31.10.2018.xlsx", sheet = 1)
dat <- data.frame(dat)

# only plot motifs with "Enrichment.factor" more than 1
dat <- dat[which(dat$Enrichment.factor > 1), ]

# only select "Motifs", "Enrichment.factor" and "P.value" columns
dat <- dat[, c(1, 7, 8)]

dat$P.value <- -log(dat$P.value, 10)
dat <- dat[order(dat$P.value, decreasing = FALSE), ]

# factor the Motif column, so later the plotting will be in the order wanted
dat <- dat[-which(duplicated(dat$Motifs)), ] # delete the duplicated motif terms
Molevel <- dat$Motifs
dat$Motifs <- factor(dat$Motifs, levels = Molevel)

maxP <- max(dat$P.value)
n <- round(maxP/5)

maxEF <- max(dat$Enrichment.factor)
fc <- maxEF/maxP

dat$Label <- round(dat$Enrichment.factor, 2)

dat$Enrichment.factor <- dat$Enrichment.factor/fc

# plotting
library(ggplot2)
p <- ggplot(data = dat, 
            aes(x = Motifs, 
                y = P.value, 
                group = 1)) +
  geom_bar(stat="identity", 
           width = 0.5, 
           fill = "steelblue",
           size = 2, 
           alpha = 1) + 
  geom_point(data = dat, 
             aes(x = Motifs,
                 y = Enrichment.factor), 
             size = 3) +
  geom_text(data = dat, 
            aes(x = Motifs,
                y = Enrichment.factor, 
                label = Label),
            vjust = -1, hjust = 0.5) +
  geom_line(data = dat, 
            aes(x = Motifs,
                y = Enrichment.factor), 
            size = 1) +
  coord_flip() + 
  ggtitle('') +
  theme_classic() +
  xlab("Motif Enrichment (Perseus)") +
  ylab("-Log(p.value)") +
  scale_y_continuous(breaks = seq(0, maxP, n),
                     labels = scales::number_format(accuracy = 0.1),
                     sec.axis = sec_axis(trans = ~.*fc,
                                         labels = scales::number_format(accuracy = 0.1),
                                         name = "Enrichment.factor")) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 1))

file <- paste(format(Sys.time(), "%F %H-%M-%S"), "MotifEnrichment", "D492", "BarDotPlot_GO.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 300, width = 6.5, height = 5)
