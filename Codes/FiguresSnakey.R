library(networkD3)
library(plotly)

day0_vs_day1 <- read.csv("/home/DESeq2.csv")
day1_vs_day6 <- read.csv("/home/DESeq2.csv")

day1_vs_day6 <- day1_vs_day6[!is.na(day1_vs_day6$padj) & day1_vs_day6$padj < 0.05, ]
day0_vs_day1 <- day0_vs_day1[!is.na(day0_vs_day1$padj) & day0_vs_day1$padj < 0.05, ]

colnames(day0_vs_day1)[colnames(day0_vs_day1) == "log2FoldChange"] <- "FoldChange_day0_vs_day1"
colnames(day1_vs_day6)[colnames(day1_vs_day6) == "log2FoldChange"] <- "FoldChange_day1_vs_day6"

merged_data <- merge(day0_vs_day1, day1_vs_day6, by = "Gene", all = TRUE)

merged_data <- merged_data[, c("Gene", "FoldChange_day0_vs_day1", "FoldChange_day1_vs_day6")]

merged_data$range_day0_vs_day1 <- cut(merged_data$FoldChange_day0_vs_day1,
                                      breaks = c(-Inf, -1, 1, Inf),
                                      labels = c("< -1", "-1 to 1", "> 1"))

merged_data$range_day1_vs_day6 <- cut(merged_data$FoldChange_day1_vs_day6,
                                      breaks = c(-Inf, -1, 1, Inf),
                                      labels = c("< -1", "-1 to 1", "> 1"))

write.csv(merged_data, "/home/mergedfoldchanges.csv", row.names = FALSE)

nodes <- data.frame(name = c("< -1 Day0", "-1 to 1 Day0", "> 1 Day0",
                             "< -1 Day1", "-1 to 1 Day1", "> 1 Day1"))
links <- data.frame(
  source = c(0,0,0,1,1,1,2,2,2),
  target = c(3,4,5,3,4,5,3,4,5),
  value = c(
    sum(merged_data$range_day0_vs_day1 == "< -1" & merged_data$range_day1_vs_day6 == "< -1"),
    sum(merged_data$range_day0_vs_day1 == "< -1" & merged_data$range_day1_vs_day6 == "-1 to 1"),
    sum(merged_data$range_day0_vs_day1 == "< -1" & merged_data$range_day1_vs_day6 == "> 1"),
    sum(merged_data$range_day0_vs_day1 == "-1 to 1" & merged_data$range_day1_vs_day6 == "< -1"),
    sum(merged_data$range_day0_vs_day1 == "-1 to 1" & merged_data$range_day1_vs_day6 == "-1 to 1"),
    sum(merged_data$range_day0_vs_day1 == "-1 to 1" & merged_data$range_day1_vs_day6 == "> 1"),
    sum(merged_data$range_day0_vs_day1 == "> 1" & merged_data$range_day1_vs_day6 == "< -1"),
    sum(merged_data$range_day0_vs_day1 == "> 1" & merged_data$range_day1_vs_day6 == "-1 to 1"),
    sum(merged_data$range_day0_vs_day1 == "> 1" & merged_data$range_day1_vs_day6 == "> 1")
  )
)

links <- links[links$value > 0, ]

pastel_colors <- c("#FAD02E", "#F28D35", "#F28D35", "#D83367", "#D83367", "#2F77B8",
                   "#2F77B8", "#6B9F4A")

sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes,
                             Source = "source", Target = "target",
                             Value = "value", NodeID = "name",
                             units = "genes", fontSize = 12,
                             colourScale = pastel_colors)

print(sankey_plot)

nodes <- c("< -1 Day0", "-1 to 1 Day0", "> 1 Day0", "< -1 Day1", "-1 to 1 Day1", "> 1 Day1")

day0_colors <- c("#4C9FCF", "#A9A9A9", "#D83367")
day1_colors <- rep("#8A9E9A", 3)

node_colors <- c(day0_colors, day1_colors)

links <- links[links$value > 0, ]

link_colors <- c(
  rep("#4C9FCF80", sum(links$source == 0)),
  rep("#A9A9A9", sum(links$source == 1)),
  rep("#D8336780", sum(links$source == 2))
)

sankey_plot <- plot_ly(
  type = "sankey",
  node = list(
    pad = 15,
    thickness = 30,
    line = list(color = "black", width = 0.5),
    color = node_colors,
    label = c("", "", "", "", "", "")
  ),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value,
    color = link_colors
  )
)

sankey_plot
