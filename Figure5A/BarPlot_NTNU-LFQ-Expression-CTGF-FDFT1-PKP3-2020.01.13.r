# Upload data for Proteins with missing values in NTNU dataset
# raw data are under "C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5-2 ProteomicPaper-UGDH\Perseus\PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.11.xlsx"
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
library(readxl)
path <- "C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/NTNU-LFQ-Expression-CTGF-FDFT1-PKP3.xlsx"
excel_sheets(path = path)
#[1] "D492-CTGF"   "HMLE-CTGF"   "PMC42-CTGF"  "D492-FDFT1"  "HMLE-FDFT1" 
#[6] "PMC42-FDFT1" "D492-PKP3"   "HMLE-PKP3"   "PMC42-PKP3" 

funPlot <- function(n){
  funData <- function(x, CellLine){
    library(readxl)
    dat <- read_excel("NTNU-LFQ-Expression-CTGF-FDFT1-PKP3.xlsx", sheet = x)
    dat <- data.frame(dat)
    
    dat[, 2] <- factor(dat[, 2], levels = CellLine)
    
    # function to calculate mean and sem
    mean.sd <- function(dat){
      
      # mean, sd, sem, number
      colna <- colnames(dat)[3]
      
      # since there are "NA" in the dataset, mean() will return "NA"
      funMean <- function(dat){
        
        m <- mean(dat, na.rm = TRUE)
        
        return(m)
      }
      
      funSD <- function(dat){
        
        s <- sd(dat, na.rm = TRUE)
        
        return(s)
      }
      
      m <- tapply(dat[, 3], dat[, 2], funMean)
      
      s <- tapply(dat[, 3], dat[, 2], funSD)
      # n <- tapply(dat[, 2], dat[, 1], length)
      # sem <- s/sqrt(n)
      
      dat <- data.frame(rbind(m, s))
      dat <- data.frame(t(dat))
      
      dat$"Cell.Line" <- factor(CellLine, levels = CellLine)
      
      colnames(dat)[1] <- colna
      colnames(dat)[2] <- "SD"
      
      dat <- dat[, c(3, 1:2)]
      
      return(dat)
      
    }
    
    dat <- mean.sd(dat)
    
    return(dat)
    
  }
  
  datD <- funData(n, c("D492", "D492M"))
  datH <- funData(n+1, c("HMLE", "HMLEM"))
  datP <- funData(n+2, c("PMC42LA", "PMC42ET"))
  
  if (n == 1){
    GN <- "CTGF"
  } else if (n == 4){
    GN <- "SERPINE1"
  } else if (n == 7){
    GN <- "FDFT1"
  } else if (n == 10){
    GN <- "PKP3"
  } else if (n == 13){
    GN <- "GANAB"
  }
  
  dat <- rbind(datD, datH, datP)
  
  dat$EMT.Model <- factor(rep(c("D492", "HMLE", "PMC42"), c(2, 2, 2)), 
                          levels = c("D492", "HMLE", "PMC42"))
  
  dat$Cell.Type <- factor(rep(c("Epithelial", "Mesenchymal"), 3), 
                          levels = c("Epithelial", "Mesenchymal"))
  
  dat <- dat[, c(1, 4, 5, 2:3)]
  
  # Plotting
  ggplot.dat <- function(dat){
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
    
    # plotting
    p <- ggplot(dat, aes(x = EMT.Model, y = dat[, 4], fill = Cell.Type)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.8, 
               color = "black", size = 4) +
      geom_errorbar(aes(ymin = dat[, 4] - SD, ymax = dat[, 4] + SD), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.8)) + 
      scale_fill_manual(values = c("steelblue", "red")) + 
      scale_y_continuous(breaks = seq(0, (max(dat[, 4]) + max(dat[, 5])), ((max(dat[, 4]) + max(dat[, 5]))/5)), 
                         labels = scales::number_format(accuracy = 0.01, 
                                                        decimal.mark = ".")) +
      theme_classic() +
      labs(title = GN, 
           x = NULL, 
           y = paste0("Relative Expression", " ", "(", colnames(dat)[4], ")")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 96), 
            axis.line = element_line(colour = 'black', size = 4),
            legend.position = "right",
            legend.key.size = unit(1.3, "cm"),
            legend.text = element_text(size = 60, face = "bold"),
            legend.title = element_text(size = 60, face = "bold"),
            axis.title.y.left = element_text(size = 64, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 96, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4))
      # geom_hline(yintercept = 1, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(dat)[4], "ThreeEMT_RNA_bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = 300, width = 28, height = 16)
    
    print(p)
    
  }
  
  ggplot.dat(dat)
  
}

funPlot(1)
funPlot(4)
funPlot(7)
