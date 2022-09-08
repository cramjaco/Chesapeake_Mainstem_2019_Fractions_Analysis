source("InitialProcessing.R")
library(tidyverse)
library(WGCNA)
nonSpikes20

countDat <- nonSpikes20 %>% 
  filter(ID != "4-3-O-0-2") %>%
  mutate(logCPL = log(copiesPerL + 1)) %>%
  pivot_wider(id_cols = "ID", names_from = "ASV", values_from = "logCPL") %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()
###########
## Choose a set of soft-thresholding powers
###########
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(countDat, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# lets use 7, or 5

#### Main Thing

goodGenes(countDat, minRelativeWeight = 0)

net = blockwiseModules(countDat, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "CB-wgcna", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

sme <- left_join(sample %>% filter(!Flag), 
          net$MEs %>% rownames_to_column("ID"),
          by = "ID"
)

smeL <- sme %>%
  pivot_longer(cols = ME1:ME0, names_to = "ME", values_to = "meVal") %>%
  filter(!is.na(Station))

smeL %>%
  ggplot(aes(x = Size_Class, y = meVal, shape = Depth, color = Depth)) + 
  facet_grid(~ME ~ Station ) +
  geom_point() + 
  scale_x_log10()


signedKME(countDat,  net$MEs) -> skme

skme1 <- sweep(skme, 1, FUN = "/", STATS = rowSums(skme))

skmeL <- skme1 %>%
  rownames_to_column(var = "ASV") %>%
  pivot_longer(cols = -ASV, names_to = "kME", values_to = "weight")

simpleWGCNAGroups <- left_join(nonSpikes20,
                               tibble(ASV = names(net$colors), cluster = net$colors),
                               by = "ASV") %>%
  group_by(Station, Depth, Size_Class, Bin_Size, cluster) %>%
  summarize(copiesPerL = sum(copiesPerL))

messy_plot(simpleWGCNAGroups)
# As expected, this is no better than the alternatives.

# However, what if I weight the contributions of each ASV to each group

ns_weighted <- nonSpikes20 %>% left_join(skmeL, by = "ASV") %>%
  mutate(WeightedCopies = copiesPerL * weight)

ns_sum_weighted <- ns_weighted %>% 
  group_by(Station, Depth, Size_Class, Bin_Size, kME) %>%
  summarise(copiesPerL = sum(WeightedCopies)) %>%
  rename(cluster = "kME")

messy_plot(ns_sum_weighted)

# No improvement, or else a very marginal improvement

# I wonder if I should cluster by presence and abscence