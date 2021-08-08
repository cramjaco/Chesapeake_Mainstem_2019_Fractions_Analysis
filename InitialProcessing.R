## Initial processing

library(tidyverse)
library(Rtsne)
## Load in Main data types
counts <- read_tsv('ASVs_counts.tsv') %>% rename(ASV = "...1")
taxa0 <- read_tsv("ASVs_taxonomy.tsv") %>% rename(ASV = "...1")
key2 <- read_csv("SampleKey2.csv")
sample0 <- read_csv("EnvDataForAmplicons.csv")
flags <- read_csv("ManualFlags.csv") %>%
  mutate(Flag = if_else(Flag == "Bad", TRUE, FALSE, missing = FALSE)) %>%
  select(ID, Flag)


## Make sample data nice
sample <-sample0 %>%
  mutate(ID = paste(
    str_replace(Station, "\\.", "-"),
    str_sub(Depth, 1, 1),
    str_replace(Size_Class, "\\.", "-"),
    sep = "-"
    )
  ) %>%
  right_join(tibble(ID = unique(key2$ID)), by = "ID") %>%
  relocate(ID) %>%
  left_join(flags, by = "ID") %>%
  mutate(Depth = factor(Depth, levels = c("Surface", "Oxy", "Bottom")))

taxa <- taxa0 %>%
  mutate(nASV = extract_numeric(ASV))
  
## Make flags file, for when I hadn't done that, so I can manually edit
# overwrite flags file
# forFlags <- sample %>%
#   select(ID, Station, Depth, Size_Class) %>%
#   arrange(Station, Depth, Size_Class)
# write_csv(forFlags, "ManualFlags.csv")
# 
renKey <- key2 %>% select(FilteredFiles, ID)


## Pivot counts table
counts_long00 <- counts %>%
  pivot_longer(cols = -ASV,
                                      names_to = "filename",
                                      values_to = "reads") %>%
  left_join(renKey, by = c('filename' = 'FilteredFiles')) %>%
  select(-filename) %>%
  left_join(taxa, by = "ASV")

## Calculate Total Counts
totCounts <- counts_long00 %>%
  group_by(ID) %>%
  summarise(TotalReads = sum(reads))


## Counts table but with relative abundances of reads
counts_ra <- counts_long00 %>% left_join(totCounts, by = "ID") %>%
  mutate(RA = reads/TotalReads)

## Counts table, grouped to kingdom level
kingdom_ra <- counts_ra %>% group_by(ID, Kingdom) %>%
  summarise(RA = sum(RA), reads = sum(reads)) %>% 
  left_join(sample %>% select(Station, Size_Class, Depth, ID), by = "ID") 

## Colors for kingdom figure
kingColor <- c("Archaea" = "Black",
               "Bacteria" = "goldenrod",
               "Eukaryota" = "darkgreen",
               "Spike" = "Blue")

## Plot out relative abundance of each kingdom for QC purposes
## used to flag things as bad
## Good samples have reads, and aren't all spike or bacteria.
kra_plot <- kingdom_ra %>% ggplot(aes(x = ID, y = RA, fill = Kingdom)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  facet_wrap(~Station + Depth, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = kingColor)
#kra_plot

ggsave("KingdomSpikeRA.png", plot = kra_plot, width = 11, height = 8.5)

## As above but this time total reads
kcount_plot <- kingdom_ra %>% ggplot(aes(x = ID, y = reads, fill = Kingdom)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  facet_wrap(~Station + Depth, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = kingColor)
#kra_plot

ggsave("KingdomSpikeCount.png", plot = kcount_plot, width = 11, height = 8.5)

# get rid of flagged entries
counts_long <- counts_long00 %>% left_join(flags,by = "ID") %>% filter(!Flag) %>% select(-Flag)

spikes <- counts_long %>%
  filter(Kingdom == "Spike") %>%
  group_by(ID) %>%
  summarise(SpikeReads = sum(reads))

correctionData <- left_join(sample, spikes, by = "ID") %>%
  mutate(conversionMultiplier = DNAperLiter * 10^5 / SpikeReads)

nonSpikes <- counts_long %>%
  filter(Kingdom != "Spike" & !is.na(Kingdom)) %>%
  left_join(correctionData) %>%
  mutate(copiesPerL = reads * conversionMultiplier) %>%
  filter(!is.na(copiesPerL))

microbialAbundance <- nonSpikes %>%
  group_by(ID) %>%
  summarise(copiesPerL = sum(na.omit(copiesPerL)), ) %>%
 left_join(sample, by = "ID") %>%
 arrange(Station, Depth, Size_Class)

write_csv(microbialAbundance, "AmpliconAbundance.csv")

maPlot <- microbialAbundance %>%
  filter(!Flag) %>%
  filter(!is.na(copiesPerL)) %>%
  filter(!is.na(Depth)) %>%
  ggplot(aes(x = Size_Class, y = copiesPerL/Bin_Size)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  facet_grid(Depth~Station, scales = "free_x") +
  scale_y_log10() + scale_x_log10()

maPlot

ggsave("MicrobesVsSize.png", plot = maPlot, width = 11, height = 8.5)



## as above but prettier
maPlot2 <- microbialAbundance %>%
  filter(!Flag) %>%
  filter(!is.na(copiesPerL)) %>%
  filter(!is.na(Depth), Depth != "Oxy") %>%
  ggplot(aes(x = Size_Class,
             y = copiesPerL/Bin_Size,
             fill = as.factor(Station),
             shape = as.factor(Station))) +
  geom_point(position = position_dodge(width = 0.1), size = 3) +
  facet_grid(rows = vars(Depth)) + 
  scale_y_log10() + scale_x_log10(breaks = c(0.2, 1.2, 5, 20, 53, 180)) +
  scale_shape_manual(values = c(21:25, 21))+
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  labs(x = "Size Class Lower Bound (um)",
       y = "Copies/L/mm",
       fill = "Station", shape = "Station") +
  scale_fill_viridis_d(option = "plasma")

maPlot2

## the above as a function that accepts a y value?

cb_ploter <- function(df, .VarName){
  VarName <- enquo(.VarName)
  plotOut <- df %>%
    ggplot(aes(x = Size_Class,
               y = !!VarName,
               fill = as.factor(Station),
               shape = as.factor(Station))) +
    geom_point(position = position_dodge(width = 0.1), size = 3) +
    facet_grid(rows = vars(Depth)) + 
    scale_y_log10() + scale_x_log10(breaks = c(0.2, 1.2, 5, 20, 53, 180)) +
    scale_shape_manual(values = c(21:25, 21))+
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    labs(x = "Size Class Lower Bound (um)",
         y = VarName,
         fill = "Station", shape = "Station") +
    scale_fill_viridis_d(option = "plasma")
  plotOut
}

maPlot2 <- microbialAbundance %>%
  filter(!Flag) %>%
  filter(!is.na(copiesPerL)) %>%
  filter(!is.na(Depth), Depth != "Oxy") %>%
  cb_ploter(copiesPerL/Bin_Size) + labs(y = "Copies/L/mm")

maPlot2

# Lots of NA values because whenever we were "over" on DNA concentration, the data got left off of the
# main spreadsheet. I need to track those data down.

maPmPlot <- microbialAbundance %>%
  filter(!Flag) %>%
  filter(!is.na(copiesPerL)) %>%
  filter(!is.na(Depth)) %>%
  ggplot(aes(x = Size_Class, y = copiesPerL/MassperLiter)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  facet_grid(Depth~Station, scales = "free_x") +
  scale_y_log10() + scale_x_log10() +
  labs(y = "Cells/mg Particles", x = "Particle Size Class (um)")

maPmPlot



# I wonder what Mass per Liter is measured in. miligrams?
# Seems overly dense

## Pre Filtering 
## Which ASVs are present in at least 20% of samples

uniqueSamples <- length(unique(nonSpikes$ID))

checkHits <- nonSpikes %>% group_by(nASV, ASV) %>%
  summarise(hits = sum(copiesPerL > 0)) %>%
  ungroup() %>%
  mutate(freq = hits/uniqueSamples)

ggplot(data = checkHits, aes(x = nASV, y = hits/uniqueSamples, col = freq > 0.2)) + geom_point()

checkHits %>% summarise(sum(freq > 0.2))

nonSpikes20 <- nonSpikes %>% left_join(checkHits %>% select(-nASV), by = "ASV") %>%
  filter(hits/uniqueSamples > 0.2) 
#%>% na.omit() # why are there NAs left? NAs are in non critial columns

length(unique(nonSpikes20$ASV))

sapply(nonSpikes20, function(x) sum(is.na(x)))

