# This script processes the isotope analysis
# So we know mg of carbonand nitrongen
# I used to think things were in units of moles, rather than mg, hense the variable name "Molarity"

library(tidyverse)
library(here)
library(readxl)
library(cowplot)
library(grid)
library(gridExtra)

metadata <- read_csv(here("InputData", "CB-DNAandPOM.csv"))
molarity <- read_excel(here("InputData", "CramJ 220614 CN Chesapeake2021 A003193 EPR A003144.xlsx"), sheet = "Samples")

load(here("RDataFiles", "InitialProcessing_3_minimal.RData"))
source(here::here("RLibraries", "ChesapeakePersonalLibrary.R"))
rm(nonSpikes, nonSpikes20)

bin_sizes <- data.frame(Size_Class = c(1.2, 5, 20, 53, 180, 500),
                            Bin_Size = c( 5 - 1.2, 20 - 5, 53-20, 180-53, 500-180, 500))

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


#mutate(MassperLiter=`Mass` / (`Volume_through_mesh` * `Volume_through_GFF` / `Backrinse`))

molarityPlus0 <- left_join(molarityCB21, metadata1, by = "FilterID") %>%
  left_join(bin_sizes, by = "Size_Class") %>%
  left_join(microbialAbundance %>%
              select(-Bin_Size) %>%
              #mutate(Depth = ordered(Depth, levels = c("Surface", "Oxycline", "Bottom"))) %>%
              identity(),
            by = c("Station", "Size_Class", "Depth") 
            )

molarityPlus <- molarityPlus0 %>%
  mutate(CarbonPerLiter = `Total C (µg)` / (Volume_through_mesh * Volume_through_GFF/Backrinse)) %>%
  mutate(NitrogenPerLiter = `Total N (µg)`/ (Volume_through_mesh * Volume_through_GFF/Backrinse)) %>%
  mutate(Depth = recode(Depth, Oxy = "Oxycline")) %>%
  mutate(Depth = if_else(Station == 3.3 & Depth == "Bottom", "Oxycline", Depth)) %>%
  mutate(Depth = ordered(Depth, levels = c("Surface", "Oxycline", "Bottom"))) %>%
  identity()

