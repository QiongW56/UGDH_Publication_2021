setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 NTNU_LFQ proteomic results/Second/OrgnizedResults")
library(readxl)
datD <- data.frame(read_excel("2-1 2018-12-28 17-51 Two.groups.statistics.LFQ.xlsx", sheet = 2))
datH <- data.frame(read_excel("2-2 2018-12-28 17-51 Two.groups.statistics.LFQ.xlsx", sheet = 2))
datP <- data.frame(read_excel("2-3 2018-12-28 17-52 Two.groups.statistics.LFQ.xlsx", sheet = 2))

funDat <- function(dat, GN, EMTmodel, nor){
  if (EMTmodel == "D492"){
    CellLine <- c("D492", "D492M")
  } else if (EMTmodel == "HMLE"){
    CellLine <- c("HMLE", "HMLE_M")
  } else if (EMTmodel == "PMC42") {
    CellLine <- c("PMC42_LA", "PMC42_ET")
  } else {
    print("Input has to be \"D492\", \"HMLE\" or \"PMC42\"", quote = FALSE)
  }
  
  dat <- dat[which(dat$Gene.Name == GN), ]
  
  dat <- dat[, c(3, 22:27)]
  
  colnames(dat)[2:7] <- substr(colnames(dat)[2:7], 17, nchar(colnames(dat)[2:7]))
  
  dat <- data.frame(t(dat))
  dat$Sample <- rownames(dat)
  dat <- dat[ ,c(2, 1)]
  dat <- dat[-1, ]
  colnames(dat)[2] <- GN
  
  dat[, 2] <- as.numeric(dat[, 2])
  
  dat$CellLine <- factor(rep(CellLine, times = c(3,3)), 
                         levels = CellLine)
  
  dat <- dat[, c(1, ncol(dat), (ncol(dat)-1))]
  
  dat[, 3] <- 2^(dat[, 3])
  
  dat[, 3] <- dat[, 3]/10000000
  
  # normalize to epithelial or mesenchymal cells
  if (nor == "epi"){
    avg_epi <- mean(dat[1:3, 3], na.rm = TRUE)
    dat[, 3] <- dat[, 3]/avg_epi
    
  } else if (nor == "mes"){
    avg_mes <- mean(dat[4:6, 3], na.rm = TRUE)
    dat[, 3] <- dat[, 3]/avg_mes
    
  }
  
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
  
}
    
funGN <- function(GN, nor){
  datD_GN <- funDat(datD, GN, "D492", nor)
  datH_GN <- funDat(datH, GN, "HMLE", nor)
  datP_GN <- funDat(datP, GN, "PMC42", nor)
  
  dat <- rbind(datD_GN, datH_GN, datP_GN)
  
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
            axis.ticks.y.left = element_line(colour = "black", size = 4)) + 
      geom_hline(yintercept = 1, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(dat)[4], "Proteomics_ThreeEMT_bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = 300, width = 28, height = 16)
    
    print(p)
    
  }
  
  ggplot.dat(dat = dat)
  
}


funGN("UGDH", "mes")
