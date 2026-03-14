# Load libraries
library(reshape2)
library(circlize)

# Load data
dat <- read.csv("/home/pathways.csv", stringsAsFactors = FALSE)

# Preserve pathway order
dat$pathway <- factor(dat$pathway, levels = dat$pathway)

# Convert to long format and keep existing connections
dat_long <- melt(dat, id.vars = "pathway")
dat_long <- subset(dat_long, value == 1)

# Define sector colors
sectores <- c(as.character(dat$pathway), colnames(dat)[-1])
set.seed(123)
pal_sectores <- rainbow(length(sectores), s = 0.7, v = 0.9)
names(pal_sectores) <- sectores

# Define link colors
destinos <- unique(dat_long$variable)
pal_links <- setNames(rainbow(length(destinos), s = 0.7, v = 0.9, alpha = 0.5), destinos)
link_colors <- pal_links[dat_long$variable]

# Draw chord diagram
circos.clear()
chordDiagram(
  dat_long[, c("pathway","variable")],
  grid.col = pal_sectores,
  col = link_colors,
  transparency = 0.4,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1
)

# Add sector labels and ticks
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.05,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter, 0.5, sector_name,
      facing = "clockwise", niceFacing = TRUE, cex = 0.7
    )
    
    circos.axis(
      h = "top",
      major.at = seq(0, 1, by = 0.2),
      labels = TRUE,
      labels.cex = 0.5,
      sector.index = sector_name,
      track.index = 1
    )
  }
)
