funData <- function(CellType){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5-2 ProteomicPaper-UGDH/Reactome")
  fileName <- paste0(CellType, ".csv")
  
  dat <- read.csv(fileName)
  dat <- data.frame(dat)
  
  colVal <- c(grep("Pathway.name", names(dat)),
              grep("X.Entities.found", names(dat)),
              grep("Entities.FDR", names(dat)))
  
  datVal <- dat[, colVal]
  
  colnames(datVal)[3] <- "FDR"
  
  datVal <- datVal[1:10, ]
  
  return(datVal)
  
}

datD <- funData("D492")
datH <- funData("HMLE")
datP <- funData("PMC42")

# function for treemap
treemap.plot <- function(dat, 
                         col.index = 1, 
                         col.size = 2, 
                         col.color = 3, 
                         color = "Blues",
                         title = "D492 Model",
                         label.size = 12,
                         lgd = 16,
                         height = 5,
                         width = 6,
                         CellType = "D492"){
  library(treemap)
  library(grDevices)
  # upload data
  dat <- data.frame(dat)
  dat[, col.size] <- as.numeric(dat[, col.size])
  dat[, col.color] <- as.numeric(dat[, col.color])
  
  lower.limit = 0
  upper.limit = max(dat$FDR)
  mid.value = upper.limit/2
  
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), CellType, "treemap.plot.tiff", " ")
  tiff(filename = file, res = 300, height = height, width = width, units = "in")
  treemap(dat, #Your data frame object
          index = c(colnames(dat)[col.index]),  #A list of your categorical variables
          vSize = colnames(dat)[col.size],  #This is your quantitative variable
          vColor = colnames(dat)[col.color],
          type="value", #Type sets the organization and color scheme of your treemap
          palette = color,  #Select your color palette from the RColorBrewer presets or make your own.
          title = title, #Customize your title
          fontsize.title = 18, #Change the font size of the title
          position.legend = 'bottom',
          mapping = c(upper.limit, mid.value, lower.limit),
          range = c(lower.limit, upper.limit),
          fontsize.labels = label.size,
          fontface.labels = 2,
          border.lwds = 2,
          fontsize.legend = lgd)
  
  dev.off()

}

treemap.plot(datD)
treemap.plot(datH, color = "Reds", title = "HMLE Model", CellType = "HMLE")
treemap.plot(datP, color = "Oranges", title = "PMC42 Model", CellType = "PMC42")

