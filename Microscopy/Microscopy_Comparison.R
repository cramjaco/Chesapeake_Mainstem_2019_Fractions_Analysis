## This file was originally written by Grace Martinez as part of her REU project in Summer 2021
## It has subsequently been modified by Cram to produce figures more compatible with others in a manuscript.

library(tidyverse)
library(readxl)
library(here)
library(dplyr)
library(readr)
library(ggplot2)
# install.packages("Rtools")
# install.packages("ggplot2")

# Input the data
MicroscopyDir <- "Microscopy"
DataLocations <- read_csv(here(MicroscopyDir, "WhereToFindMicroscopy.csv"))
Data_File <- here(MicroscopyDir, "CB  Microscopy.xlsx") # two spaces in this filename

# Create a function for the data just inputted
StuffGetter <- function(Station, SizeRange, BoxesRange, DataRange){
  Sizes <- read_excel(Data_File, sheet = as.character(Station), range = as.character(SizeRange), col_names = FALSE) %>%
    as.matrix() %>% as.vector()
  Boxes <- read_excel(Data_File, sheet = as.character(Station), range = as.character(BoxesRange), col_names = FALSE) %>%
    as.matrix() %>% as.vector()
  Data <- read_excel(Data_File, sheet = as.character(Station), range = as.character(DataRange),col_names = FALSE) %>%
    as.matrix() 
  # In cases where there are any characters, we want to replace those with NA, rather than treat everything as a character vector (important for CB3.3Bottom)
  Data <- matrix(as.numeric(Data), ncol = ncol(Data))
  Mean <- apply(Data, 2, mean) %>% as.vector()
  SD <- apply(Data, 2, sd) %>% as.vector()
  #list(Sizes, Boxes, Data, Mean, SD)
  tibble(Sizes, Boxes, Mean, SD)
}
  
# Tells StuffGetter where to find other data
GetRowStuff <- function(rowN){
  StuffGetter(DataLocations[rowN,"Station"],
              SizeRange = DataLocations[rowN,"SizeRange"],
              BoxesRange = DataLocations[rowN, "BoxesRange"],
              DataRange = DataLocations[rowN, "DataRange"])}

# Test and debug
# Its throwing away all of 3.3 Bottom, which is missing just one slize class.
GetRowStuff(1) # works
GetRowStuff(5) # fails to calculate means

# Makes one big table with the data we need
StuffWithMetadata <- DataLocations %>% mutate(RowNum = 1:9) %>%
  mutate(Data = map(RowNum, GetRowStuff)) %>%
  unnest(Data) %>%
  select(-c(SizeRange, BoxesRange, DataRange))

#Calculations:
# area of the circle in mm
pi * (13/2)^2
# grid area in mm
0.1 * 0.1
# how many little boxes in the big circle
GridsPerFilter <- (pi * (13/2)^2)/(0.1 * 0.1)

# Big data table with calculations
CalculationsWithMetadata <- StuffWithMetadata %>%
  mutate(CellsPer100Boxes = Mean/Boxes * 100)%>%
  mutate(CellsPerFilter = GridsPerFilter * CellsPer100Boxes) %>%
  mutate(CellsPermL = CellsPerFilter/Volume) %>%
  mutate(CellsPerLiterofBackrinse = CellsPermL * 1000) %>%
  filter(Sizes != "<1.2") %>%
  mutate(NSize = if_else(Sizes == "<5", 0.2, parse_number(Sizes))) %>%
  mutate(Sizes2 = gsub(",", "", Sizes)) %>%
  ungroup()%>%
  # 3.3 Bottom is actually 3.3 Oxycline
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth)) %>%
  mutate(Depth2 = factor(Depth, levels = c("Surface","Oxy","Bottom")))%>%
  mutate(Depth2 = recode(Depth2, Oxy = "Oxycline"))

# Input more data "SizevsBinSize", then join with old CalculationsWithMetadata table

SizevsBinSize <- read_csv(here(MicroscopyDir, "SizeVSBinSize.csv"))
CalculationsWithMetadata01 <- left_join(CalculationsWithMetadata, SizevsBinSize, by =c("NSize"="Size")) %>%
  mutate(CellsPerLiterofBackrinsePermm = CellsPerLiterofBackrinse / BinSize) %>%
  # I'm not sure why the following block has to happen again, but apparently it does.
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth)) %>%
  mutate(Depth2 = ordered(Depth, levels = c("Surface","Oxy","Bottom")))%>%
  mutate(Depth2 = recode(Depth2, Oxy = "Oxycline")) %>%
  identity()

# Plot of Cells / Liter of Backrinse / mm VS Size (mm)
ggplot(CalculationsWithMetadata01, aes(x=NSize, y=CellsPerLiterofBackrinsePermm, color=as.factor(Station))) + 
  labs(color= "Station", y= "Cells / Liter of Backrinse / mm", x="Size (mm)") +
  geom_point() + 
  geom_path() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + facet_wrap(~Depth2, ncol = 1) + scale_y_log10() + scale_x_log10()
ggsave(here(MicroscopyDir, "SizevsCellsPerLiterPermm.pdf"))

# Input more data "Mass Supplement", then join with old CalculationsWithMetadata table

# Grace was using this file, so I'll keep using this one, rather than pulling from outside
# of the subdirectory
MassSupplement <- read_csv(here(MicroscopyDir, "mass_supplement.csv")) %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth)) %>%
  
  mutate(Depth2 = recode(Depth, Oxy = "Oxycline")) %>%
  mutate(Depth2 = ordered(Depth2, levels = c("Surface","Oxycline","Bottom")))%>%
  identity()

BigDataTable <- left_join(CalculationsWithMetadata01, 
                          MassSupplement,
              by=c("Station","Depth2", "NSize"="Size_Class"))%>%
  filter(Sizes2 != "<5") %>%
  mutate(CellsPermgofParticles = CellsPerLiterofBackrinse / MassperLiter)%>%
  mutate(MicrobesPerParticle = CellsPerLiterofBackrinse / ParticlesPerLiter)%>%
  # mutate(Depth2 = factor(Depth, levels = c("Surface","Oxy","Bottom")))%>%
  # mutate(Depth2 = recode (Depth2, Oxy = "Oxycline")) %>% 
  identity()

# Plot of Microbes / Particle VS Size (mm)
ggplot(BigDataTable, aes(x=NSize, y=MicrobesPerParticle, color=as.factor(Station))) + labs(color= "Station", y= "Microbes Per Particle", x="Size (mm)")+
  geom_point() + geom_path() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + facet_wrap(~Depth2, ncol = 1) + scale_y_log10() + scale_x_log10() 
ggsave(here(MicroscopyDir, "MicrobesPerParticlevsSize.pdf"))

# Input more data "Amplicon abundance"

# Here I shouldn't use the same file, since I've rerun the amplicon data
# and added station 5.1 Bottom
# But Intital processing regenerates this file so I'll use that.

AmpliconAbundance <- read.csv(here("IntermediateData","AmpliconAbundance.csv")) %>%
  mutate(DNAPerParticle = DNAperLiter / ParticlesPerLiter)%>%
  filter(Station !="3.1")%>%
  filter(Station !="3.2") %>%
  mutate(Depth = ordered(Depth, levels = c("Surface", "Oxycline", "Bottom")))

# Input more data "DNA and POM", then join with old CalculationsWithMetadata

DNA_and_POM <- read_csv(here(MicroscopyDir, "DNA_and_POM (1).csv")) %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth))
BackRinseCalc <- left_join(DNA_and_POM, CalculationsWithMetadata, by=c("Station","Depth","Size_Class"="NSize"))%>%
mutate(Cells_Collected_Total = CellsPerLiterofBackrinse * Backrinse)%>%
mutate(Cells_Per_Liter_Total = Cells_Collected_Total / Volume_through_mesh)%>%
  filter(Station != "3.1")%>%
  filter(Station != "3.2")%>%
  filter(Sizes != "<5") %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxy", Depth)) %>%
  
  mutate(Depth = recode(Depth, Oxy = "Oxycline")) %>%
  mutate(Depth = ordered(Depth, levels = c("Surface","Oxycline","Bottom")))%>%
  identity()
  

filterSizes <- c(5, 20, 53, 180, 500)

# Join BackRInseCalc with AmpliconAbundance, then plot
BackRinseCalc2 <- left_join(BackRinseCalc, 
                            AmpliconAbundance, by=c("Station","Depth","Size_Class")) %>%
  mutate(fStation = factor(Station, levels = as.character(c(3.1, 3.2, 3.3, 5.1, 4.3, 5.5))))

ggplot(BackRinseCalc2,
       aes(x=copiesPerL/Bin_Size,
                           y=Cells_Per_Liter_Total/Bin_Size, 
                           fill=as.factor(Station), 
                           shape=as.factor(Station),
                           color=as.factor(Depth),
                           size=(Size_Class)^(1/2))) + 
  labs(shape = "Station", fill= "Station", size = "Size Class (μm)", color = "Depth", y= "Microscopy (Cells/L)", x="Amplicon (Copies/L)") +
  geom_point(stroke = 1) + 
  geom_path(aes(x=copiesPerL/Bin_Size,
                y=Cells_Per_Liter_Total/Bin_Size,
                group = interaction(Depth, Station)),inherit.aes = FALSE) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  # facet_wrap(~Depth2, ncol = 1) +
  scale_y_log10() + scale_x_log10() +
  geom_abline(intercept=0, slope=1) +
  scale_shape_manual(values = rep(21:25, 2)[3:10]) +
  scale_fill_manual(values = viridis::plasma(n = 6)[3:6]) +
  #scale_fill_viridis_d(breaks = as.character(c(3.1, 3.3, 4.3, 5.1, 5.5)), labels = as.character(c(3.1, 3.3, 4.3, 5.1, 5.5))) + 
  scale_size(breaks = (filterSizes)^(1/2),
             labels = filterSizes,
             limits = ((c(5, 500)^(1/2)))) +
  scale_color_manual(values = c(Surface = "green", Oxycline = "blue", Bottom = "black")) +
  theme_bw() +
  guides(
    fill = guide_legend(ncol = 2, override.aes = list()),
    shape = guide_legend(override.aes = list(size = 3)),
    size = guide_legend(override.aes = list(shape = 21, color = "black", fill = "gray")),
    color = guide_legend(override.aes = list(shape = 21, fill = "gray"))
  )

ggsave(here(MicroscopyDir,"MicroscopyVsAmplicons.png"), height = 4, width = 6)

# Make a model and check your work to see statistical significance of results
BackRinseCalc2forModel <- BackRinseCalc2 %>% filter(is.finite(Cells_Per_Liter_Total), is.finite(copiesPerL))
mod <- lm(log10(Cells_Per_Liter_Total/Bin_Size) ~ log10(copiesPerL/Bin_Size), data = BackRinseCalc2)
mod
summary(mod)
# There is a "relationship" bewteen microscopy counts and amplicon abundance, which is  a weak stament
# but the best I can think of for now.