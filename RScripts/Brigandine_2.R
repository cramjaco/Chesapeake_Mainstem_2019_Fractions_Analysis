## This file returns the brigandine style plots for total bacteria, planktomycetes and the abundant bacteria that
# are primarely on large size fractions.

library(here)
library(tidyverse)

# Bring in data
#source(here("RScripts", "InitialProcessing_3.R"))
load(here::here("RDataFiles", "InitialProcessing_3.RData"))
source(here("RLibraries", "Brigandine_Library.R"))

# Phylum Level Plots
phylum_mg_plot <- ches_brigandine(Phylum, Kingdom, taxa01 %>% pull(Kingdom) %>% unique(), min = 4, max = 7, thresh = 10^6, ns = nonSpikes)
phylum_L_plot <- ches_brigandine_L(Phylum, Kingdom, taxa01 %>% pull(Kingdom) %>% unique(), min = 4, max = 9, thresh = 10^6, ns = nonSpikes)

ggsave(here("Figures", "Phylum_per_mg_brig.png"), height = 5.5, width = 6.5, plot = phylum_mg_plot)
ggsave(here("Figures", "Phylum_per_L_brig.png"), height = 5.5, width = 6.5, plot = phylum_L_plot)

## Planktomycetes ASVs
plankto_mg_plot <- ches_brigandine(Tag_ASV, Order, taxa01 %>% filter(Phylum == "Planctomycetes") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 10^6, min = 3)
plankto_L_plot <- ches_brigandine_L(Tag_ASV, Order, taxa01 %>% filter(Phylum == "Planctomycetes") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 10^6, min = 3)

ggsave(here("Figures", "Plankto_per_mg_brig.png"), height = 4.5, width = 7.5, plot = plankto_mg_plot)
ggsave(here("Figures", "Plankto_per_L_brig.png"), height = 4.5, width = 7.5, plot = plankto_L_plot)


## Snow dwellers primarely live on large particles
mainHabitats <- nonSpikes %>% # was ns20
  group_by(Tag_ASV, Size_Class) %>% 
  summarise(copiesPerL = mean(copiesPerL)) %>%
  ungroup() %>% group_by(Tag_ASV) %>%
  arrange(Tag_ASV, -copiesPerL) %>%
  top_n(1) %>% ungroup()

mainHabitats %>% 
  group_by(Size_Class) %>%
  summarize(ntaxa = length(copiesPerL))

snowDwellers <- mainHabitats %>% 
  filter(Size_Class >= 20) %>%
  pull(Tag_ASV)

snow_dwellers_plot <- ches_brigandine_L(ASV, Phylum,
                taxa %>% filter(Tag_ASV %in% snowDwellers) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV %in% snowDwellers) %>%
                  #mutate(ASV = recode(ASV, `Arachnida;9` = "Acartia;9")), # If you blast the "arachnid" sequence, it turns out to be Acartia tonsa
                  mutate(ASV = str_replace(ASV, "Arachnida", "Acartia")), # more general solution than the above line
                thresh = 10^4)

snow_dwellers_plot <- ches_brigandine_L(ASV, Phylum,
                taxa %>% filter(Tag_ASV %in% snowDwellers) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                 # more general solution than the above line
                thresh = 10^3, ns = nonSpikes  %>% filter(Tag_ASV %in% snowDwellers))

ggsave(here("Figures", "SnowDwellers_per_L_brig.png"), height = 4.5, width = 7.5, plot = snow_dwellers_plot)
ggsave(here("Figures", "SnowDwellers_per_L_brig_narrow.png"), height = 4.5, width = 5.5, plot = snow_dwellers_plot)
