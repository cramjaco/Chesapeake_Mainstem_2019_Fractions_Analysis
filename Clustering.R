source("InitialProcessing.R")
## Ordination


nonSpikes20 %>% filter(copiesPerL > 0) %>%
  summarize(min(copiesPerL))

wideASVs <- nonSpikes20 %>%
  mutate(copiesPerLPermm = copiesPerL/Bin_Size) %>%
  select(ID, ASV, copiesPerLPermm) %>%
  pivot_wider(id_cols = ID, names_from = ASV, values_from = copiesPerLPermm) %>%
  column_to_rownames("ID") %>% na.omit()
# 

wideASVs_log <- log(wideASVs + 1)

sdist <- (1-cor(wideASVs, method = "spear"))/2
 
library(vegan)
library(ape)
PCOA <- pcoa(sdist)

myPca <- rda(X = wideASVs)
plot(myPca)

myNMDS <- metaMDS(comm = wideASVs, dist = "euclid")
myNMDS_scores <- scores(myNMDS) %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>% as.tibble() %>%
  left_join(sample, by = "ID")
ggplot(aes(x = NMDS1, y = NMDS2, shape = as.factor(Station), color = Size_Class), data = myNMDS_scores) + geom_point()
# this looks bad now. All driven by very few outliers

heatmap(sdist)

library(cluster)
kmet <- pam(sdist, 5)

representitives <- rownames(kmet$medoids)

repTax <- taxa %>% filter(ASV %in% representitives)

repLong <- nonSpikes20 %>%
  filter(ASV %in% representitives)

repLong %>% group_by(ASV) %>% filter(row_number() == 1)

repLong %>% 
  ggplot(aes(y = copiesPerL, x = Size_Class, color = Depth)) + geom_point() +
  facet_grid(Order ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()

ggsave("RepLongPlot.pdf")



byKMed <- nonSpikes20 %>% left_join(
  tibble(ASV = names(kmet$clustering), cluster = kmet$clustering), by = "ASV") %>%
  group_by(Station, Depth, Size_Class, cluster) %>%
  summarise(copiesPerL = sum(copiesPerL),
            MassperLiter = first(MassperLiter),
            Bin_Size = first(Bin_Size))
  

byKMed %>% ggplot(aes(y = copiesPerL/Bin_Size, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
  facet_grid(cluster ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()
ggsave("ClusterAbundances.pdf")
# differences in these group totals are really friggin subtle
# also some values that should have been elimnated like 5.1 bottom are still there
# for misterious reasons
# true even if I log transform before taking spearman distance (which makes sense, its all rank based)

byKMed %>% ggplot(aes(y = copiesPerL/MassperLiter, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
  facet_grid(cluster ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()
# subtle even if I normalize to mass

metaASVs <- metaMDS(sdist)

masvData <- as.data.frame(scores(metaASVs)) %>% rownames_to_column(var = "ASV") %>% left_join(
  tibble(ASV = names(kmet$clustering), cluster = kmet$clustering), by = "ASV")

masvData %>% ggplot(aes(x = NMDS1, y = NMDS2, shape = as.factor(cluster))) + geom_point(size = 2)
ggsave("masv.pdf")

# Klaus -- throw out shared zeros
# Me - Coefficients fo regression?

#https://stackoverflow.com/questions/37396403/calculating-correlation-after-removing-zeros

wideASVs_mtx <- as.matrix(wideASVs)

  spear_n00 <- function(A, B){
    cor(A[!(A == 0 & B == 0)], B[!(A ==0 & B == 0)], method = "spear")
  }
  
  spear_n00_all <- function(x){
    rho_mtx <- matrix(nrow = ncol(x), ncol = ncol(x))
    colnames(rho_mtx) <- colnames(x)
    rownames(rho_mtx) <- colnames(x)
    for(i in 1:ncol(x)){
      for(j in 1:ncol(x)){
        rho_mtx[i, j] = spear_n00(A = x[,i], B = x[,j])
      }
    }
  rho_mtx
  }

scorNozeros <- spear_n00_all(wideASVs_mtx)
distNozeros <- (1 - scorNozeros)/2

## workup

kmetNz <- pam(distNozeros, 5)

representitivesNz <- rownames(kmetNz$medoids)

repTax <- taxa %>% filter(ASV %in% representitivesNz)

repLongNz <- nonSpikes20 %>%
  filter(ASV %in% representitivesNz)

repLongNz %>% 
  ggplot(aes(y = copiesPerL, x = Size_Class, color = Depth)) + geom_point() +
  facet_grid(Order ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()

# full clusters

byKMedNz <- nonSpikes20 %>% left_join(
  tibble(ASV = names(kmetNz$clustering), cluster = kmetNz$clustering), by = "ASV") %>%
  group_by(Station, Depth, Size_Class, cluster) %>%
  summarise(copiesPerL = sum(copiesPerL),
            MassperLiter = first(MassperLiter),
            Bin_Size = first(Bin_Size))
  

byKMedNz %>% ggplot(aes(y = copiesPerL/Bin_Size, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
  facet_grid(cluster ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()

## ok, the no-shared zeros thing isn't any sort of panacea.
##

## Now with hclust
#hclust(distNozeros) # fails
library(dendextend)
sclust <- hclust(dist(sdist))
plot(sclust)
hgroups <- cutree(sclust, k = 10)

byHc <- nonSpikes20 %>% left_join(
  tibble(ASV = names(hgroups), cluster = hgroups), by = "ASV") %>%
  group_by(Station, Depth, Size_Class, cluster) %>%
  summarise(copiesPerL = sum(copiesPerL),
            MassperLiter = first(MassperLiter),
            Bin_Size = first(Bin_Size))

byHc %>% ggplot(aes(y = copiesPerL/Bin_Size, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
  facet_grid(cluster ~ Station) + scale_color_manual(values = c("Green", "Blue", "Black"))+
theme(axis.text.x = element_text(angle = 90, size = 6)) +
  scale_x_log10() + scale_y_log10()

# A bit better
# Especially with 10 rather than 5 cuts

## What if I cluster directly
justPam <- pam(t(wideASVs_mtx), 5)

sum_over_clusters <- function(data = nonSpikes20, asv_vec, cluster_vec){
  byClusters <- data %>% left_join(
    tibble(ASV = asv_vec, cluster = cluster_vec), by = "ASV") %>%
    group_by(Station, Depth, Size_Class, cluster) %>%
    summarise(copiesPerL = sum(copiesPerL),
              MassperLiter = first(MassperLiter),
              Bin_Size = first(Bin_Size))
  byClusters
}

justPam_sums <- sum_over_clusters(nonSpikes20, names(justPam$clustering), justPam$clustering)

messy_plot <- function(df){
  df %>% ggplot(aes(y = copiesPerL/Bin_Size, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
    facet_grid(cluster ~ Station) + scale_color_manual(values = c("darkgreen", "Blue", "Black"))+
    theme(axis.text.x = element_text(angle = 90, size = 6)) +
    scale_x_log10() + scale_y_log10()
}


messy_plot(justPam_sums)

# doesn't fix anything

# What if I log transform everything before I cluster.
# There seems to be variability in the smallest size

log(wideASVs_mtx + 2.5)

logPam <- pam(t(scale(log(wideASVs_mtx + 2.5))), 20)
logPam_sums <- sum_over_clusters(nonSpikes20, names(logPam$clustering), logPam$clustering)
messy_plot(logPam_sums)

nonSpikes20  %>% filter(copiesPerL !=0) %>% ggplot(aes(copiesPerL))  + geom_histogram()

nonSpikes20  %>% filter(copiesPerL !=0) %>% summarize(min(copiesPerL))

# note, I need to do a variance stabilizing transformation eventually
# These just aren't that different from eachother!

## As above, but this time, I'll plot out every asv
allCluster <- nonSpikes20 %>% left_join(tibble(ASV = names(logPam$clustering), cluster = logPam$clustering),
                          by = "ASV")

allCluster %>%
  filter(Depth == "Surface") %>%
  ggplot(aes(x = Size_Class, y = copiesPerL, group = ASV)) + 
  facet_grid(cluster ~ Station) +
  geom_line(alpha = 0.3) +
  scale_x_log10() + scale_y_log10()
# These clusters don't really look like clusters to me.
# Better when there are more of them

## TSNE

library(M3C)
my_tsne <- Rtsne(t(scale(log(wideASVs_mtx + 2.5))), check_duplicates = FALSE)

plot(my_tsne$Y)

df_tsne <- tibble(colnames(wideASVs_mtx), my_tsne$Y[,1], my_tsne$Y[,2])
mtx_tsne <- my_tsne$Y
rownames(mtx_tsne) <- colnames(wideASVs_mtx)

pam_tsne <- pam(mtx_tsne, k = 50)

tsne_sums <- sum_over_clusters(nonSpikes20, names(pam_tsne$clustering), pam_tsne$clustering)
plot_tsne_messy <- messy_plot(tsne_sums)

ggsave("MessyTsne.pdf", plot_tsne_messy, height = 24, width = 8)
# not even this

# ok, it looks ok if I just have a lot more groups
# 50 are diverse, 5 are not, 10 is ok if you squint

# I wonder if this works for regular pam
logPam50 <- pam(t(scale(log(wideASVs_mtx + 2.5))), 50)
logPam50_sums <- sum_over_clusters(nonSpikes20, names(logPam50$clustering), logPam50$clustering)
plot_pam_messy <- messy_plot(logPam50_sums)
ggsave("MessyPam.pdf", plot_pam_messy, height = 24, width = 8)

## What about hclust?

logHclust <- hclust(dist(scale(log(t(wideASVs_mtx + 2.5)))))
logHclustCut <- cutree(logHclust, k = 50)
logHclustCut_sums <- sum_over_clusters(nonSpikes20, names(logHclustCut), logHclustCut)
plot_hclust_messy <- messy_plot(logHclustCut_sums)
ggsave("MessyHclust.pdf", plot_hclust_messy, height = 24, width = 8)

# can I do as above, but with station as the x axis? I'd prefer latitude, but number is probably ok for now.

MessyHclustInverted <- logHclustCut_sums %>% 
  ggplot(aes(y = copiesPerL/Bin_Size, x = Station, color = Depth)) + geom_point(shape = 1) + 
  facet_grid(cluster ~ Size_Class) + scale_color_manual(values = c("darkgreen", "Blue", "Black"))+
    theme(axis.text.x = element_text(angle = 90, size = 6)) +
    scale_x_log10() + scale_y_log10()

ggsave("MessyHclustInverted.pdf", MessyHclustInverted, height = 24, width = 8)

## Plot all the things
plot_hclust_many <- nonSpikes20 %>% 
  left_join(tibble(ASV = names(logHclustCut), cluster = logHclustCut), by = "ASV") %>%
  arrange(Size_Class) %>% 
  ggplot(aes(x = Size_Class, y = copiesPerL/Bin_Size, group = ASV, color = Depth)) +
  facet_grid(cluster ~ Station) +
  geom_path() +
  scale_x_log10() + scale_y_log10()
ggsave("ManyHclust.pdf", plot_hclust_many, height = 24, width = 8)


## Zooming in on a few of these to try to ounderstand oddness

nonSpikes20 %>% 
  left_join(tibble(ASV = names(logHclustCut), cluster = logHclustCut), by = "ASV") %>%
  arrange(Size_Class) %>% 
  filter(cluster == 3, Station == 3.1) %>%
  filter(Depth == "Bottom") %>%
  ggplot(aes(x = Size_Class, y = copiesPerL/Bin_Size, group = ASV, color = Depth)) +
  geom_path() +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ASV)

LocalStuff <- nonSpikes20 %>% 
  left_join(tibble(ASV = names(logHclustCut), cluster = logHclustCut), by = "ASV") %>%
  arrange(Size_Class) %>% 
  filter(cluster == 17, Station == 3.3) %>%
  filter(Depth == "Bottom")
LocalStuff %>%
  ggplot(aes(x = Size_Class, y = copiesPerL/Bin_Size)) +
  geom_path(aes( group = ASV)) +
  scale_x_log10() + scale_y_log10() +
  geom_point(data = LocalStuff %>%
               group_by(Size_Class, Bin_Size) %>%
               summarize(copiesPerL = sum(copiesPerL))
             )

# there are vertical bars here which really bother me
## I wonder what the pca eigenvectors look like. Is it clear that this just doesn't compress
# into 5-1 dimensions or something?
library(vegan)
myPCA <- rda((scale(log((wideASVs_mtx + 2.5)))))
screeplot(myPCA)
plot(eigenvals(myPCA)/sum(eigenvals(myPCA)))
# odd.

# I wonder what would happen if I plotted the first x pcs vs the samples
# kind of an odd concept, since a bunch would be negative



site_score_df <- scores(myPCA, display = "sites") %>% as.data.frame() %>% rownames_to_column("ID")
#nonSpikes20PC <- nonSpikes20 %>% left_join(site_score_df, by = "ID")
# Lets come back to this, I need all of the PCs. Its going to be a mess

plot(myPCA)

sample %>% left_join(site_score_df, by = "ID") %>% filter(!is.na(Station)) %>%
  ggplot(aes(x = PC1, y = PC2,
             shape = as.factor(Station), fill = as.factor(Station),
             color = Depth,
             size = Size_Class
             )) +
  scale_shape_manual(values = rep(c(21:25), 2)) + scale_fill_viridis_d() +
  scale_color_manual(values = c("darkgreen", "blue", "black")) +
  geom_point(stroke = 2) +
  theme_bw()

## I'd like to facet wrap this but with the first five PCs
site_score_df01 <- scores(myPCA, choices = 1:7, display = "sites") %>% as.data.frame() %>% rownames_to_column("ID") %>%
  pivot_longer(PC2:last_col(), names_to = "PCB_name", values_to = "PCB_val")

sample %>% left_join(site_score_df01, by = "ID") %>% filter(!is.na(Station), !is.na(PCB_name)) %>%
  ggplot(aes(x = PC1, y = PCB_val,
             shape = as.factor(Station), fill = as.factor(Station),
             color = Depth,
             size = Size_Class
             )) +
  scale_shape_manual(values = rep(c(21:25), 2)) + scale_fill_viridis_d() +
  scale_color_manual(values = c("darkgreen", "blue", "black")) +
  geom_point(stroke = 2) +
  theme_bw() +
  facet_wrap(~PCB_name)
ggsave("FirstFivePca.pdf", width = 11, height = 8)


## Vs size and depth
site_score_df02 <- scores(myPCA, choices = 1:7, display = "sites") %>% as.data.frame() %>% rownames_to_column("ID") %>%
  pivot_longer(PC1:last_col(), names_to = "PC_name", values_to = "PC_val")

sample %>% left_join(site_score_df02, by = "ID") %>%
  filter(!is.na(Station), !is.na(PC_name)) %>%
  mutate(StnGrp = paste0(Station, Depth)) %>%
  ggplot(aes(x = Size_Class, y = PC_val,
             shape = as.factor(Station), fill = as.factor(Station),
             color = Depth,
             group = StnGrp
             )) +
  scale_shape_manual(values = rep(c(21:25), 2)) + scale_fill_viridis_d() +
  scale_color_manual(values = c("darkgreen", "blue", "black")) +
  geom_line(color = "black") +
  geom_point(stroke = 1) +
  theme_bw() +
  facet_grid(Depth != "Surface" ~ PC_name) + 
  scale_x_log10() +
    theme(axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
ggsave("FirstFivePcaB.pdf", width = 11, height = 8)
## Wow is this ugly.