setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-Signaling/PDGFRB-RELA-UGDH")
library(readxl)
path <- "C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-Signaling/PDGFRB-RELA-UGDH/siPDGFRB-siRELA-UGDH.xlsx"
es <- excel_sheets(path = path)
es
# [1] "D492M-UGDH" 14         "D492M-PDGFRB"  14      "D492M-RELA"  10       
# [4] "D492HER2-UGDH" 12      "D492HER2-PDGFRB" 12    "D492HER2-RELA" 12     
# [7] "D492HER2-siRELARELA1" 12 "D492HER2-siRELAUGDH1" 12 "D492HER2-siRELARELA2" 12
# [10] "D492HER2-siRELAUGDH2" 12 "D492M-siRELARELA1" 18 "D492M-siRELAUGDH1" 18
# [13] "D492M-siRELARELA2" 12  "D492M-siRELAUGDH2" 12

funBoxPlot <- function(sheetName = "D492M_PDGFRB", 
                       nr = 18,
                       KDgene = "siPDGFRB",
                       # col = c("steelblue", "skyblue1"), 
                       gap = 0.2,
                       dotS = 1.25){ 
  # c("steelblue", "skyblue1"), c("red", "pink"), c("orange1", "#FAE48BFF")
  # nr is 18 or 16 (for HER2)
  # gap is for the y axis
  # dotS for dot size if the dots are plotted
  # D492_NFE2L2 is 16, the rest is 18
  # D492M_CDH1 is 16, the others are 18
  # D492HER2 is 16 for all.
  
  # import rawdata from RT-qPCR
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-Signaling/PDGFRB-RELA-UGDH")
  library(readxl)
  library(ggplot2)
  dat <- read_excel("siPDGFRB-siRELA-UGDH.xlsx", sheet = sheetName)
  dat <- data.frame(dat)
  
  # Different genes may have different replicates
  # dat <- dat[1:nr, ]
  
  scr <- rep("Scr", length(grep("Scr", dat[, 1])))
  gn <- KDgene
  kd <- rep(gn, length(grep(gn, dat[, 1])))
  
  dat$Treatment <- c(scr, kd)
  
  dat$Treatment <- factor(dat$Treatment, levels = c("Scr", gn))
  
  dat <- dat[, c(1, ncol(dat), 2:(ncol(dat)-1))]
  
  dat <- dat[which(!is.na(dat[, ncol(dat)])), ]
  
  # plotting
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
  
  CellType <- gsub("-.*", "", sheetName)
  
  if (CellType == "D492"){
    col <- c("steelblue", "skyblue1")
  } else if (CellType == "D492M"){
    col <- c("red", "pink")
  } else if (CellType == "D492HER2"){
    col <- c("orange1", "#FAE48BFF")
  }
  
  sdScr <- sd(dat[1:(nrow(dat)/2), 3])
  sdKD <- sd(dat[(nrow(dat)/2+1):nrow(dat), 3])
  
  if (sdScr > sdKD){
    addOn <- sdScr
  } else (
    addOn <- sdKD
  )
  
  round2 <- function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  
  ymax <- ceiling(max(dat[, 3])+ addOn)
  
  p <- ggplot(dat, 
              aes(x = Treatment, 
                  y = dat[, 3], 
                  fill = Treatment)) +
    stat_boxplot(geom = 'errorbar', 
                 position = position_dodge(1),
                 width = 0.5,
                 size = 2) +
    geom_boxplot(position = position_dodge(1), 
                 lwd = 4, 
                 fatten = 2,
                 outlier.colour = "black", 
                 outlier.shape = 16,
                 outlier.size = 6, 
                 notch = FALSE) +
    #             outlier.colour = "red", outlier.shape = 8, outlier.size = 4) +
    scale_fill_manual(values = col) +
    #geom_dotplot(binaxis = 'y', 
     #            stackdir = 'center',
      #           position = position_dodge(1),
       #          binwidth = 0.025,
        #         dotsize = dotS) +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_y_continuous(breaks = seq(0, ymax, gap),
                       labels = scales::number_format(accuracy = 0.01, 
                                                      decimal.mark = ".")) +
    theme_classic() +
    labs(title = paste0(colnames(dat)[ncol(dat)], " in ", CellType), 
         x = "Treatment", 
         y = "Relative Gene Expression") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 42), 
          axis.line = element_line(colour = "black", size = 2),
          legend.position = "none",
          legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 25, face = "bold"),
          axis.title.y.left = element_text(size = 36, face = "bold", 
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x.bottom = element_text(size = 36, face = "bold",
                                             margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text = element_text(size = 36, face = "bold", color = "black"),
          axis.ticks.y.left = element_line(colour = "black", size = 2))
  p
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Figures-Signaling/PDGFRB-RELA-UGDH/Figures")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
                CellType,
                colnames(dat)[ncol(dat)],
                KDgene,
                "Box.plot.tiff", 
                sep = "_")
  ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 10)
  
  print(p)
  
}

funBoxPlot(sheetName = es[1], nr = 14)
funBoxPlot(sheetName = es[2], nr = 14)
funBoxPlot(sheetName = es[3], nr = 10)
funBoxPlot(sheetName = es[4], nr = 12)
funBoxPlot(sheetName = es[5], nr = 12)
funBoxPlot(sheetName = es[6], nr = 12)
funBoxPlot(sheetName = es[7], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[8], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[9], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[10], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[11], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[12], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[13], nr = 12, KDgene = "siRELA")
funBoxPlot(sheetName = es[14], nr = 12, KDgene = "siRELA")

