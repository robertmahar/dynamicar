# Load simulation data and plot expected utility.

# ============================================================================
# Load libraries/config/functions.

library(latticeExtra)
library(RColorBrewer)
library(grid)
library(tidyverse)

# ============================================================================
# Specify directories.

dirOut  <- file.path("analysis", "output")
dirSims <- file.path(dirOut, "runSimBatch")

# ============================================================================
# Load files.

tIndex           <- 64
simulatedUtility <- lapply(seq(tIndex), 
  function(x) get(load(file.path(dirSims, paste0("runSimBatch", x, ".Rdata")))))

# Munge data.
lUtilities <- lapply(simulatedUtility, 
  function(x) lapply(x, 
  function(y) mean(sapply(y, 
  function(k) k[["utilities"]]))))

lUtilites    <- lapply(lUtilities, function(x) unlist(x))

# Get the names directly from the list.
makeNameData <- function(st)
{
  df <- as.data.frame(lapply(strsplit(st, ", "), function(x) strsplit(x, " = ")))
  colnames(df) <- df[1,]
  df <- df[-1,]
}  

namesScenarios <- unlist(lapply(lUtilities, names))

allNames <- lapply(namesScenarios, function(x) makeNameData(x))
allNames <- do.call(rbind, allNames)
allNames[["myopic"]] <- ifelse(allNames[["myopic"]] == "FALSE", 0, 1)
simData <- cbind(allNames, unlist(lUtilities))
simData <- data.frame(lapply(simData, as.numeric))
names(simData)[length(simData)] <- "expectedUtils"

# ============================================================================
# Make plots.

myPanel <- function(...) {
  arg <- list(...)
  panel.levelplot(...)
}
cols <- (colorRampPalette(brewer.pal(6, "RdYlGn"))(100))

colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))

doPlot <- function(x, data)
{
  plot1 <- levelplot(expectedUtils ~ pbo_1 * trt_1 , data = data[data[["index"]] == x,], 
                     col.regions = colr, 
                     panel = myPanel, 
                     region = T, 
                     colorkey = list(space = "bottom"),
                     scales=list(x = list(at = NULL), 
                                 y = list(at = NULL),
                                 cex = 1,
                                 alternating = F, tck = c(0,0)), 
                     at = seq(0.4, 1.6, 0.01),  
                     xlab = "", 
                     ylab = "", 
                     pretty = T, 
                     par.settings   = list(
                       layout.widths  = list(left.padding = 6, right.padding  = 6), 
                       layout.heights = list(top.padding  = 4, bottom.padding = 3, xlab.key.padding = 5)))
}

plotTrel <- function(data)
{
comb_levObj <- c(doPlot(data = data, 1),  doPlot(data = data, 2),  doPlot(data = data, 3),  doPlot(data = data, 4),  doPlot(data = data, 5),  doPlot(data = data, 6),  doPlot(data = data, 7),  doPlot(data = data, 8),
                 doPlot(data = data, 9),  doPlot(data = data, 10), doPlot(data = data, 11), doPlot(data = data, 12), doPlot(data = data, 13), doPlot(data = data, 14), doPlot(data = data, 15), doPlot(data = data, 16),
                 doPlot(data = data, 17), doPlot(data = data, 18), doPlot(data = data, 19), doPlot(data = data, 20), doPlot(data = data, 21), doPlot(data = data, 22), doPlot(data = data, 23), doPlot(data = data, 24),
                 doPlot(data = data, 25), doPlot(data = data, 26), doPlot(data = data, 27), doPlot(data = data, 28), doPlot(data = data, 29), doPlot(data = data, 30), doPlot(data = data, 31), doPlot(data = data, 32),
                 doPlot(data = data, 33), doPlot(data = data, 34), doPlot(data = data, 35), doPlot(data = data, 36), doPlot(data = data, 37), doPlot(data = data, 38), doPlot(data = data, 39), doPlot(data = data, 40),
                 doPlot(data = data, 41), doPlot(data = data, 42), doPlot(data = data, 43), doPlot(data = data, 44), doPlot(data = data, 45), doPlot(data = data, 46), doPlot(data = data, 47), doPlot(data = data, 48),
                 doPlot(data = data, 49), doPlot(data = data, 50), doPlot(data = data, 51), doPlot(data = data, 52), doPlot(data = data, 53), doPlot(data = data, 54), doPlot(data = data, 55), doPlot(data = data, 56),
                 doPlot(data = data, 57), doPlot(data = data, 58), doPlot(data = data, 59), doPlot(data = data, 60), doPlot(data = data, 61), doPlot(data = data, 62), doPlot(data = data, 63), doPlot(data = data, 64),
                 layout = c(8, 8), merge.legends = FALSE) 
print(comb_levObj)

grid.text("Infections, stage 1 placebo (%)",   x = 0.5,  y = 0.12)
grid.text("Infections, stage 1 prophylaxis (%)", x = 0.023, y = 0.5, rot = 90)
grid.text("Deaths, stage 1 placebo (%)",   x = 0.5,  y = 0.99)
grid.text("Deaths, stage 1 prophylaxis (%)", x = 0.98, y = 0.5, rot = 270)

cls <- as.matrix(c(5, 10, 20, 40, 60, 80, 90, 95))
rws <- cls

## loop over panels to be labelled (ie 1:3)
panels  = trellis.currentLayout()

for (i in 1:64) {
  
  # focus on current panel of interest and disable clipping
  ids <- which(panels == i, arr.ind = TRUE)
  print(ids)
  trellis.focus("panel", ids[2], ids[1], clip.off = TRUE)
  
  colShift <- 1.25
  rowShift <- colShift

  if(ids[1] == 8 & ids[2] == 8) {
    grid.text(cls[ids[1]], x = .5, y = colShift)
    grid.text(rws[ids[1]], x = rowShift, y = .5, rot = 0)

  }
  if(ids[2] == 1 & ids[1] %in% c(1,3,5,7)) {
    panel.axis(side = "left",
               at = seq(0, 1),
               labels = c(0, 100),
               ticks = TRUE, rot = 0,
               tck = 1, outside = T)
  }
  if(ids[1] == 1 & ids[2] %in% c(1,3,5,7)) {
    panel.axis(side = "bottom",
               at = seq(0, 1),
               labels = c(0, 100), rot = 0,
               ticks = TRUE, 
               tck = 1, outside = T)
  }
  if(ids[1] == 8 & ids[2] != 8) {
    grid.text(cls[ids[2]], x = .5, y = colShift)
  }
  if(ids[1] != 8 & ids[2] == 8){
    grid.text(rws[ids[1],], x = rowShift, y = .5, rot = 0)
  }

  trellis.unfocus()
}
}

# ============================================================================
# Munge the plot data. 

fixedMyopic  <- simData %>% filter(interimNumber == 1, myopic == 1)
fixedDynamic <- simData %>% filter(interimNumber == 1, myopic == 0)

adaptiveMyopic  <- simData %>% filter(interimNumber == 4, myopic == 1)
adaptiveDynamic <- simData %>% filter(interimNumber == 4, myopic == 0)

# Compare adaptive to fixed strategy (myopic)
cf_fixed_adaptive_myopic <- adaptiveMyopic
cf_fixed_adaptive_myopic[["expectedUtils"]] <- cf_fixed_adaptive_myopic[["expectedUtils"]]/fixedMyopic[["expectedUtils"]]

# Compare adaptive to fixed strategy (dynamic)
cf_fixed_adaptive_dynamic <- adaptiveDynamic
cf_fixed_adaptive_dynamic[["expectedUtils"]] <- cf_fixed_adaptive_dynamic[["expectedUtils"]]/fixedDynamic[["expectedUtils"]]

# ============================================================================
# Save the plots. 

outfile_1 <- file.path(dirOut, "fixedMyopic_v_adaptiveMyopic.pdf")
outfile_2 <- file.path(dirOut, "fixedDynamic_v_adaptiveDynamic.pdf")

pdf(file = outfile_1, height = 10, width = 9)
plotTrel(cf_fixed_adaptive_myopic)
grid.text("Relative empirical mean utility", x = 0.5,  y = 0.015)
dev.off()

pdf(file = outfile_2, height = 10, width = 9)
plotTrel(cf_fixed_adaptive_dynamic)
grid.text("Relative empirical mean utility", x = 0.5,  y = 0.015)

dev.off()

# End of script.


