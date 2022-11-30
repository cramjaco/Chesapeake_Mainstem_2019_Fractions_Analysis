library(here)
library(tidyverse)
library(ggh4x)

# Bring in data
#source(here("RScripts", "InitialProcessing_3.R"))
load(here::here("RDataFiles", "InitialProcessing_3.RData"))
source(here("RLibraries", "Brigandine_Library.R"))

## Eukaryotes
ches_brigandine(Phylum, Kingdom, taxa01 %>% filter(Kingdom == "Eukaryota") %>% pull(Kingdom) %>% unique(), min = 4, max = 7, thresh = 10^5, ns = nonSpikes)
#art <- ches_brigandine(, Class, taxa01 %>% filter(Phylum == "Arthropoda") %>% pull(Class) %>% unique(), min = 4, thresh = 0, ns = nonSpikes)
art1 <- ches_brigandine(Order, Class, taxa01 %>%  pull(Class) %>% unique() %>% na.omit(), thresh = 10^5, min = 3, max = 5,
                        ns = nonSpikes %>% filter(Phylum == "Arthropoda"))
art1
art2 <- ches_brigandine_L2(Tag_ASV, Order, taxa01 %>%  pull(Order) %>% unique() %>% na.omit(), thresh = 10^4, min = 3, max = 5,
                        ns = nonSpikes %>% filter(Phylum == "Arthropoda"))
art2 + facet_nested(Class + Order ~ Depth, drop = TRUE, scales = "free", space = "free") +
  labs(y = "ASV")
