## Initial processing
# 2022 Nov 22; station 3.3 "Bottom" is actually an oxycline depth, updating accordingly

library(tidyverse)
#library(Rtsne)
library(readxl)

library(here)

source(here::here("RLibraries", "ChesapeakePersonalLibrary.R"))

## Load in Main data types
counts <- read_tsv(here("DadaData",'ASVs_counts.tsv')) %>% rename(ASV = "...1")
taxa0 <- read_tsv(here("DadaData", "ASVs_taxonomy.tsv")) %>% rename(ASV = "...1")
key2 <- read_csv(here("Keys", "SampleKey2.csv"))
sample0 <- read_csv(here("InputData", "EnvDataForAmplicons.csv"))
flags <- read_csv(here("Keys", "ManualFlags2.csv")) %>%
  mutate(Flag = if_else(Flag == "Bad", TRUE, FALSE, missing = FALSE)) %>%
  select(ID, Flag)


## Make sample data nice
sample01 <-sample0 %>%
  mutate(ID = if_else((Station == 5.1 & Depth == "Bottom"),
         paste0("C_5P1B_",
                str_replace(Size_Class, "\\.", "P")
                ),
         paste(
    str_replace(Station, "\\.", "-"),
    str_sub(Depth, 1, 1),
    str_replace(Size_Class, "\\.", "-"),
    sep = "-"
    )
  )) %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth)) %>% # Station 3.3 is actually Oxycline
  mutate(Depth = if_else(Depth == "Oxy", "Oxycline", Depth)) %>% # Lets use the full word Oxycline; Recoding would have worked
  right_join(tibble(ID = unique(key2$ID)), by = "ID") %>%
  relocate(ID) %>%
  left_join(flags, by = "ID") %>%
  mutate(Depth = factor(Depth, levels = c("Surface", "Oxycline", "Bottom")))

taxa01 <- taxa0 %>%
  mutate(nASV = extract_numeric(ASV))

# Add greek letters to proteobacteria and assign class as phylum for those
taxa01 %>% filter(Phylum == "Proteobacteria") %>% .$Class %>% unique() # these are the proteobacteria I want to rename
taxa02 <- taxa01 %>%
  mutate(Class = str_replace_all(Class,
                                 c(
                                   "^Alpha" = "α-",
                                   "^Beta" = "β-",
                                   "^Gamma" = "γ-",
                                   "^Delta" = "δ-",
                                   "^Zetta" = "ζ-"
                                 )
                                 )) %>%
  mutate(Phylum = if_else(Phylum == "Proteobacteria", Class, Phylum))

# give a "Tag" which is the finest level of taxa known
TaxResTab <- tibble(
  TagLevel = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  TaxRes = 1:6
)
TagDf <- taxa02 %>%
  select(-nASV) %>%
  pivot_longer(Kingdom:Genus, names_to = "TagLevel", values_to = "Tag") %>%
  left_join(TaxResTab) %>%
  filter(!is.na(Tag)) %>%
  group_by(ASV) %>%
  slice_max(order_by = TaxRes, n = 1) %>%
  select(-TaxRes) %>%
  #mutate(ASV1 = str_extract(ASV, "\\d+")) %>%
  mutate(Tag_ASV = paste(Tag, str_extract(ASV, "\\d+"), sep = ";"))

taxa <- taxa02 %>%
  left_join(TagDf, by = "ASV")
  
## Elemental Data
metadata <- read_csv(here("InputData", "CB-DNAandPOM.csv"))
molarity <- read_excel(here("InputData", "CramJ 220614 CN Chesapeake2021 A003193 EPR A003144.xlsx"), sheet = "Samples")

molarityCB21 <- molarity %>%
  filter(str_detect(`Sample ID`, "CB19|Blank")) %>%
  mutate(FilterID = str_extract(`Sample ID`, "(?<=-)\\d{3}") %>% parse_number()) %>%
  mutate(Type = if_else(str_detect(`Sample ID`, "Blank"), "Blank", "Sample"))
  
metadata1 <- metadata %>%
  select(FilterID = `GFF Filter ID`, Volume_through_mesh, Backrinse, Volume_through_GFF, Station, Depth ,Size_Class) %>%
  # Deal with oxycline malarky
  mutate(Depth = recode(Depth, Oxy = "Oxycline")) %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxycline", Depth)) %>%
  #mutate(Depth = ordered(Depth, levels = c("Surface", "Oxycline", "Bottom"))) %>%
  identity()

elementalJoined <- left_join(molarityCB21, metadata1 %>% filter(!is.na(FilterID)), by = "FilterID")
elementalFull <- elementalJoined %>%
  mutate(CarbonPerLiter_mg = `Total C (µg)` / (Volume_through_mesh * Volume_through_GFF/Backrinse) / 1000) %>%
  mutate(NitrogenPerLiter_mg = `Total N (µg)`/ (Volume_through_mesh * Volume_through_GFF/Backrinse) / 1000)
  
elemental <- elementalFull %>%
  select(Station, Depth, Size_Class, `δ13CVPDB (‰)`, `δ15NAir (‰)`, CarbonPerLiter_mg, NitrogenPerLiter_mg)

sample <- sample01 %>%
  left_join(elemental, by = c("Station", "Depth", "Size_Class")) %>%
  identity()

## Shaping Data

## Make flags file, for when I hadn't done that, so I can manually edit
# overwrite flags file
# forFlags <- sample %>%
#   select(ID, Station, Depth, Size_Class) %>%
#   arrange(Station, Depth, Size_Class)
# write_csv(forFlags, here("Keys", "ManualFlags2_Preflag.csv"))
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

# ## Counts table, grouped to kingdom level
# kingdom_ra <- counts_ra %>% group_by(ID, Kingdom) %>%
#   summarise(RA = sum(RA), reads = sum(reads)) %>% 
#   left_join(sample %>% select(Station, Size_Class, Depth, ID), by = "ID") 
# 
# ## Colors for kingdom figure
# kingColor <- c("Archaea" = "Black",
#                "Bacteria" = "goldenrod",
#                "Eukaryota" = "darkgreen",
#                "Spike" = "Blue")
# 
# ## Plot out relative abundance of each kingdom for QC purposes
# ## used to flag things as bad
# ## Good samples have reads, and aren't all spike or bacteria.
# kra_plot <- kingdom_ra %>% ggplot(aes(x = ID, y = RA, fill = Kingdom)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   theme(axis.text.x = element_text(angle = 90, size = 6)) +
#   facet_wrap(~Station + Depth, scales = "free_x", nrow = 1) +
#   scale_fill_manual(values = kingColor)
# #kra_plot
# 
# ggsave("KingdomSpikeRA.png", plot = kra_plot, width = 11, height = 8.5)
# 
# ## As above but this time total reads
# kcount_plot <- kingdom_ra %>% ggplot(aes(x = ID, y = reads, fill = Kingdom)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   theme(axis.text.x = element_text(angle = 90, size = 6)) +
#   facet_wrap(~Station + Depth, scales = "free_x", nrow = 1) +
#   scale_fill_manual(values = kingColor)
# #kra_plot
# 
# ggsave("KingdomSpikeCount.png", plot = kcount_plot, width = 11, height = 8.5)

# get rid of flagged entries
counts_long <- counts_long00 %>% left_join(flags,by = "ID") %>% filter(!Flag) %>% select(-Flag) ## flags missing the new samples 2022 Sep 08

spikes <- counts_long %>%
  filter(Kingdom == "Spike") %>%
  group_by(ID) %>%
  summarise(SpikeReads = sum(reads))

correctionData <- left_join(sample, spikes, by = "ID") %>%
  mutate(conversionMultiplier = DNAperLiter * 10^5 / SpikeReads)

nonSpikes <- counts_long %>%
  filter(Kingdom != "Spike" & !is.na(Kingdom)) %>%
  left_join(correctionData, by = "ID") %>%
  left_join(counts_ra %>% select(ASV, ID, RA), by = c("ASV", "ID")) %>%
  # this next two lines removes the mocks, becauese the conversion multiplier is NA
  mutate(copiesPerL = reads * conversionMultiplier) %>%
  filter(!is.na(copiesPerL)) %>%
  filter(SpikeReads > 0) %>%
  # add in elemental data
  #left_join(elemental, by = c("Station", "Depth", "Size_Class")) %>%
  identity()

microbialAbundance <- nonSpikes %>%
  group_by(ID) %>%
  summarise(copiesPerL = sum(na.omit(copiesPerL)), ) %>%
 left_join(sample, by = "ID") %>%
 arrange(Station, Depth, Size_Class)

write_csv(microbialAbundance, here("IntermediateData","AmpliconAbundance.csv"))



## Pre Filtering 
## Which ASVs are present in at least 20% of samples

uniqueSamples <- length(unique(nonSpikes$ID))

checkHits <- nonSpikes %>% group_by(nASV, ASV) %>%
  summarise(hits = sum(copiesPerL > 0)) %>%
  ungroup() %>%
  mutate(freq = hits/uniqueSamples)

# ggplot(data = checkHits, aes(x = nASV, y = hits/uniqueSamples, col = freq > 0.2)) + geom_point()

checkHits %>% summarise(sum(freq > 0.2))

# Only ASVs that show up in at least 20% of samples
nonSpikes20 <- nonSpikes %>% left_join(checkHits %>% select(-nASV), by = "ASV") %>%
  filter(hits/uniqueSamples > 0.2) 
#%>% na.omit() # why are there NAs left? NAs are in non critial columns

length(unique(nonSpikes20$ASV))

sapply(nonSpikes20, function(x) sum(is.na(x)))

my_sizes <- sort(unique(microbialAbundance$Size_Class))
my_sizes

## Pull in stations
stations <- read_csv(here("InputData","stations.csv"))
stations01 <- stations %>% mutate(Station = str_sub(Station, start = 3)) %>%
  mutate(isWest = str_detect(Station, "W")) %>% filter(!isWest) %>% mutate(Station = str_remove(Station, "C")) %>% select(-isWest)

sampleData <- microbialAbundance %>%
  left_join(stations01 %>% mutate(Station = as.numeric(Station)), by = "Station") %>%
  mutate(Depth2 = if_else(Depth == "Surface", "Surface", "NotSurface"))


## Dependent library functions

ches_plot_options <- list(
  scale_y_log10nice() ,
  scale_x_log10(breaks = my_sizes, labels = as.character(my_sizes)) ,
  geom_point(size = 2) ,
  geom_path(aes(color = as.factor(Station))) ,
  scale_shape_manual(values = rep(21:25, 2)) ,
  scale_fill_viridis_d(option = "plasma") ,
  scale_color_viridis_d(option = "plasma")
)

dir.create("RDataFiles", showWarnings = FALSE)
save.image(here("RDataFiles", "InitialProcessing_3.RData"))
save(nonSpikes20, file = here::here("RDataFiles", "nonSpikes20.RData"))
save(nonSpikes, file = here::here("RDataFiles", "nonSpikes.RData"))
save(microbialAbundance, file = here::here("RDataFiles", "microbialAbundance.RData"))
save(nonSpikes, nonSpikes20, microbialAbundance, sampleData, my_sizes, ches_plot_options, file = here("RDataFiles", "InitialProcessing_3_minimal.RData"))

write_csv(nonSpikes, file = gzfile(here::here("Tables", "nonSpikes.csv.gz")))
