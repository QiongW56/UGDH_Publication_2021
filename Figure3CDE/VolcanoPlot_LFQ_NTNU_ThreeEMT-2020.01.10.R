# Volcano plot for Three EMT models (NTNU proteomics)
# -----------------------------------------------------

# prepare the data
# -----------------------------------------------------

# import data from Perseus, since it has -log(p.value)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Perseus")
datPerseus <- read.delim("PerseusOutputData-NTNU-ThreeEMT-1.6.14.0-2021.01.10.txt")
datPerseus <- data.frame(datPerseus)

# Delete the first two rows: type, group
datPerseus <- datPerseus[-1:-2, ]

# Split the datasets into three EMT models
datD <- datPerseus[, c(46, 48:49, 19, 1:6, 34:36)]
datH <- datPerseus[, c(46, 48:49, 21, 7:12, 38:40)]
datP <- datPerseus[, c(46, 48:49, 23, 13:18, 42:44)]

# function for volcano plotting
# -----------------------------------------------------

funPlot <- function(dat, 
                    EMTmodel, 
                    n, # n is the fold threshold for the "FC" column
                    fdr,
                    label,
                    Log2Fold = 1, 
                    significance = 0.05, 
                    xintercept = 1, 
                    yintercept = 1.5,
                    col,
                    title = "",
                    height = 10,
                    width){
  colnames(dat) <- c("ProteinID", "GeneName", "EC", "Signficant",
                     colnames(dat)[5:10], "Pvalue", "FDR", "FoldChange")
  
  for (i in 5:13) {
    dat[, i] <- as.numeric(dat[ ,i] )
  }
  
  # add labels based on fold changes (more than 8 folds) to the volcano plot.
  dat$Group <- ifelse(abs(dat[, 13]) >= 3 & dat[, 4] == "+", "YES", "NO") # 12
  
  # add labels based on metabolic enzymes or not, D492-1-17, HMLE-1.5-17, PMC42-1.5-13
  dat$Group <- ifelse(dat[, 3] != "" & abs(dat[, 13]) >= n & dat[, 4] == "+", "FC", dat$Group)
  
  dat <- dat[, c(2, 13, 11:12, 14)]
  
  #-------------------------------------------------------------------------------------#
  
  # plotting
  
  #-------------------------------------------------------------------------------------#
  
  if (EMTmodel == "D492"){
    LabNa <- "D492/D492M"
    mes <- "Higher in D492M"
    epi <- "Higher in D492"
  } else if (EMTmodel == "HMLE"){
    LabNa <- "HMLE/HMLEM"
    mes <- "Higher in HMLEM"
    epi <- "Higher in HMLE"
  } else if (EMTmodel == "PMC42"){
    LabNa <- "PMC42LA/PMC42ET"
    mes <- "Higher in PMC42ET"
    epi <- "Higher in PMC42LA"
  }
  
  
  volcano.plot <- function(dat, 
                           fdr,
                           label,
                           Log2Fold, 
                           significance, 
                           xintercept, 
                           yintercept,
                           col,
                           title,
                           labNames,
                           height,
                           width) {
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
    library(grDevices)
    
    # prepare datasets for plotting
    if (fdr == T) {
      dat <- data.frame(dat)
      for (i in 2:4) {
        dat[, i] <- as.numeric(dat[, i])
      }
      colnames(dat)[1:4] <- c("GeneName",
                              "Log2FoldChanges",
                              "pvalue",
                              "FDR")
      if (label == T) {
        colnames(dat)[5] <- "label"
      } else {
        dat <- dat[, -5]
      }
      
    } else {
      dat <- data.frame(dat[, c(1:3, 5)])
      for (i in 2:3) {
        dat[, i] <- as.numeric(dat[, i])
      }
      colnames(dat)[1:3] <- c("GeneName",
                              "Log2FoldChanges",
                              "pvalue")
      if (label == T) {
        colnames(dat)[4] <- "label"
      } else {
        dat <- dat[, -4]
      }
    }
    
    # add another "LFQ.Significance" char column
    if (fdr == TRUE) {
      dat <- dat %>%
        mutate(LFQ.Significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (FDR < significance) ~ mes,
                                                   (Log2FoldChanges > Log2Fold) & (FDR < significance) ~ epi,
                                                   (FDR > significance) ~ "FDR > 0.05",
                                                   TRUE ~ "Less than 2-fold")))
    } else {
      dat <- dat %>%
        mutate(LFQ.Significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (pvalue > significance) ~ "downregulated",
                                                   (Log2FoldChanges > Log2Fold) & (pvalue > significance) ~ "upregulated",
                                                   (pvalue > significance) ~ "not significant",
                                                   TRUE ~ "significant with less fold changes")))
    }
    
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
    
    # Plotting
    if (TRUE %in% grepl("label", colnames(dat), ignore.case = F)) {
      p <- ggplot(dat) +
        geom_point(aes(x = Log2FoldChanges, 
                       y = pvalue, 
                       color = LFQ.Significance), alpha = 0.5, size = 8) + # alpha for transparency
        scale_color_manual(values = col) +
        scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
        scale_x_continuous(breaks = seq(-10, 10, 2)) +
        ggtitle(title) +
        xlab(bquote("Log"[2]*"("*.(labNames)*")")) +
        ylab(expression("-Log"[10]*"(p value)")) +
        theme_classic() +
        theme(legend.position = "none", 
              legend.title = element_text(size = rel(3), face = "bold"),
              legend.text = element_text(size = rel(2.8), face = "bold"),
              plot.title = element_text(size = rel(3.6), face = "bold", hjust = 0.5), 
              axis.title = element_text(size = rel(4.8), face = "bold"),
              axis.text = element_text(size = rel(4.8), face = "bold", color = "black"),
              axis.line = element_line(size = 4)) +
        geom_hline(yintercept = yintercept, linetype = "dashed", size = 2) +
        geom_vline(xintercept = xintercept, linetype = "dashed", size = 2) +
        geom_vline(xintercept = -xintercept, linetype = "dashed", size = 2) +
        geom_text_repel(aes(x = Log2FoldChanges, y = pvalue, 
                            label = ifelse(label == "YES", GeneName,"")), size = 8) +
        geom_text_repel(aes(x = Log2FoldChanges, y = pvalue, 
                            label = ifelse(label == "FC", GeneName,"")), size = 8, fontface = "bold")
      
      print(p)
      
      setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
      file <- paste(format(Sys.time(), "%F %H-%M-%S"), EMTmodel, "NTNU_EMT_VolcanoPlot.tiff")
      ggsave(filename = file, units="in", width = width, height = height, dpi = 300)
      
    } else {
      p <- ggplot(dat) +
        geom_point(aes(x = Log2FoldChanges, 
                       y = pvalue, 
                       color = LFQ.Significance), alpha = 0.5, size = 8) + # alpha for transparency
        scale_color_manual(values = col) +
        scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
        scale_x_continuous(breaks = seq(-10, 10, 2)) +
        ggtitle(title) +
        xlab(bquote("Log"[2]*"("*.(labNames)*")")) +
        ylab(expression("-Log"[10]*"(p value)")) +
        theme_classic() +
        theme(legend.position = "right", 
              legend.title = element_text(size = rel(3), face = "bold"),
              legend.text = element_text(size = rel(2.8), face = "bold"),
              plot.title = element_text(size = rel(3.6), face = "bold", hjust = 0.5),
              axis.title = element_text(size = rel(4.8), face = "bold"),
              axis.text = element_text(size = rel(4.8), face = "bold", color = "black"),
              axis.line = element_line(size = 4)) +
        geom_hline(yintercept = yintercept, linetype = "dashed", size = 2) +
        geom_vline(xintercept = xintercept, linetype = "dashed", size = 2) +
        geom_vline(xintercept = -xintercept, linetype = "dashed", size = 2)
      
      print(p)
      
      setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH")
      file <- paste(format(Sys.time(), "%F %H-%M-%S"), EMTmodel, "NTNU_EMT_VolcanoPlot.tiff")
      ggsave(filename = file, units="in", width = width, height = height, dpi = 300)
    }
  }
  
  volcano.plot(dat = dat, 
               fdr = fdr,
               label = label,
               Log2Fold = Log2Fold, 
               significance = significance, 
               xintercept = xintercept, 
               yintercept = yintercept,
               col = col,
               title = title,
               labNames = LabNa,
               height = height,
               width = width)
  
}

funPlot(datD, EMTmodel = "D492", n = 1, fdr = TRUE, label = TRUE, col = c("grey","steelblue","red","darkgrey"), width = 13)
funPlot(datH, EMTmodel = "HMLE", n = 1.5, fdr = TRUE, label = TRUE, col = c("grey","steelblue","red","darkgrey"), width = 13)
funPlot(datP, EMTmodel = "PMC42", n = 1.5, fdr = TRUE, label = TRUE, col = c("grey","red","steelblue","darkgrey"), width = 13)

# since we cannot have legend to plot with labels, we will plot legend without labels
funPlot(datD, EMTmodel = "D492", n = 1, fdr = TRUE, label = FALSE, col = c("grey","steelblue","red","darkgrey"), width = 18)
funPlot(datH, EMTmodel = "HMLE", n = 1.5, fdr = TRUE, label = FALSE, col = c("grey","steelblue","red","darkgrey"), width = 18)
funPlot(datP, EMTmodel = "PMC42", n = 1.5, fdr = TRUE, label = FALSE, col = c("grey","red","steelblue","darkgrey"), width = 18)
