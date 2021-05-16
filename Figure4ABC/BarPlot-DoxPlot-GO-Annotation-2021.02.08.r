# This function is to plot a flipped barplot with two y axis (horizontal now)
# upload data should be at least with these three columns: "Term", "Count" and "P.Value"
# or "P-Value", since R does not recoganize "-", will automatically change "-" to "."
funGOplot <- function(fileName = "D492_GO_ALL_Functional_Annotation_Chart",
                      sheetNO = 1,
                      CellType = "D492",
                      gaps = 5,
                      vj = -1,
                      hj = 0.5,
                      wt = 6.5){ # gaps is the distance between each point on the x axis
  fileName <- paste0(fileName, ".xlsx")
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/DAVID")
  library(readxl)
  dat <- read_excel(fileName, sheet = sheetNO)
  dat <- data.frame(dat)
  
  colVal <- c(grep("Term", colnames(dat)),
              grep("Count", colnames(dat)),
              grep("PValue", colnames(dat)))
  
  dat <- dat[, colVal]
  
  colnames(dat) <- c("GO_Annotation", "Gene_Counts", "P_value")
  dat <- dat[1:10, ]
  
  dat$P_value <- -log(dat$P_value, 10)
  dat <- dat[order(dat$P_value, decreasing = FALSE), ]
  
  # factor the GO column, so later the plotting will be in the order wanted
  GOlevel <- dat$GO_Annotation
  dat$GO_Annotation <- factor(dat$GO_Annotation, levels = GOlevel)
  
  maxP <- max(dat$P_value)
  n <- round(maxP/gaps)
  
  maxGC <- max(dat$Gene_Counts)
  fc <- maxGC/maxP
  
  dat$Label <- dat$Gene_Counts
  
  dat$Gene_Counts <- round(dat$Gene_Counts/fc)
  
  if (sheetNO == 1){
    GO <- "BP"
  } else if (sheetNO == 2){
    GO <- "CC"
  } else {
    GO <- "MF"
  }
  
  if (CellType == "D492"){
    cols <- "steelblue"
  } else if (CellType == "HMLE"){
    cols <- "red"
  } else {
    cols <- "orange1"
  }
  
  CT <- paste0("(", CellType, " Model)")
  AO <- paste("GO Annotation", GO, sep = "_")
  xtitle <- paste(AO, CT, " ")
  
  # plotting
  library(ggplot2)
  p <- ggplot(data = dat, 
              aes(x = GO_Annotation, 
                  y = P_value, 
                  group = 1)) +
    geom_bar(stat="identity", 
             width = 0.5, 
             fill = cols,
             size = 2, 
             alpha = 1) + 
    geom_point(data = dat, 
               aes(x = GO_Annotation,
                   y = Gene_Counts), 
               size = 3) +
    geom_text(data = dat, 
              aes(x = GO_Annotation,
                  y = Gene_Counts, 
                  label = Label),
                  vjust = vj, hjust = hj) +
    geom_line(data = dat, 
              aes(x = GO_Annotation,
                  y = Gene_Counts), 
              size = 1) +
    coord_flip() + 
    ggtitle('') +
    theme_classic() +
    xlab(xtitle) +
    ylab("-Log(p.value)") +
    scale_y_continuous(breaks = seq(0, maxP, n),
                       labels = scales::number_format(accuracy = 1),
                       sec.axis = sec_axis(trans = ~.*fc,
                                           labels = scales::number_format(accuracy = 1),
                                           name = "Gene Count")) +
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
  
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), GO, CellType, "BarDotPlot_GO.tiff", sep = " ")
  ggsave(filename = file, units = "in", dpi = 300, width = wt, height = 5)
  
}

funGOplot(wt = 8)
funGOplot(sheetNO = 2)
funGOplot(sheetNO = 3, hj = 0.6)

funGOplot(fileName = "HMLE_GO_ALL_Functional_Annotation_Chart",
          CellType = "HMLE", wt = 8)
funGOplot(fileName = "HMLE_GO_ALL_Functional_Annotation_Chart",
          sheetNO = 2, CellType = "HMLE")
funGOplot(fileName = "HMLE_GO_ALL_Functional_Annotation_Chart",
          sheetNO = 3, CellType = "HMLE", hj = 0.6)

funGOplot(fileName = "PMC42_GO_ALL_Functional_Annotation_Chart",
          CellType = "PMC42", wt = 8)
funGOplot(fileName = "PMC42_GO_ALL_Functional_Annotation_Chart",
          sheetNO = 2, CellType = "PMC42")
funGOplot(fileName = "PMC42_GO_ALL_Functional_Annotation_Chart",
          sheetNO = 3, CellType = "PMC42", hj = 0.6)
