funBoxPlot <- function(sheetName = "D492M_siRNA1_SNAI1", 
                       nr = 14, 
                       # col = c("steelblue", "skyblue1"), 
                       gap = 0.2,
                       dotS = 1.25){ 
  # nr is 18 or 16 (for HER2)
  # gap is for the y axis
  # dotS for dot size if the dots are plotted
  # 14, 12, 14, 10, 10, 18
  
  # import rawdata from RT-qPCR
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/19 Knockdown-siRNA1&2-UGDH")
  library(readxl)
  library(ggplot2)
  dat <- read_excel("RNA_GeneExpression-siUGDH.xlsx", sheet = sheetName)
  dat <- data.frame(dat)
  
  # Different genes may have different replicates
  dat <- dat[1:nr, ]
  
  scr <- rep("Scr", length(grep("Scr", dat[, 1])))
  gn <- "siUGDH"
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
  
  CellType <- gsub("_.*", "", sheetName)
  
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
  
  library(stringr)
  siRNA <- str_match(sheetName, "_\\s*(.*?)\\s*_")[2]
  
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
    labs(title = paste0(colnames(dat)[ncol(dat)], " in ", CellType, " (", siRNA, ")"), 
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
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/19 Knockdown-siRNA1&2-UGDH/Figures")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
                CellType,
                colnames(dat)[ncol(dat)],
                "siUGDH_Box.plot.tiff", 
                sep = "_")
  ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 10)
  
  print(p)
  
}
  