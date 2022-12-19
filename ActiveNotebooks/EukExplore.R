library(here)
library(tidyverse)
library(ggh4x)

#  [1] "Annelida"            "Arthropoda"          "Ascomycota"          "Basidiomycota"       "Bicosoecida"         "Brachiopoda"        
#  [7] "Cercozoa"            "Chlorophyta_ph"      "Choanoflagellida"    "Chytridiomycota"     "Ciliophora"          "Cnidaria"           
# [13] "Cryptomonadales"     "Ctenophora"          "Dinoflagellata"      "Florideophycidae"    "Ichthyosporea"       "Incertae_Sedis"     
# [19] "Kathablepharidae"    "Labyrinthulomycetes" "MAST-12"             "Mollusca"            "Nematoda"            "Ochrophyta"         
# [25] "Peronosporomycetes"  "Phragmoplastophyta"  "Picozoa" 

# Bring in data
#source(here("RScripts", "InitialProcessing_3.R"))
load(here::here("RDataFiles", "InitialProcessing_3.RData"))
source(here("RLibraries", "Brigandine_Library.R"))

## Eukaryotes
ches_brigandine(Phylum, Kingdom, taxa01 %>% filter(Kingdom == "Eukaryota") %>% pull(Kingdom) %>% unique(), min = 4, max = 7, thresh = 10^5, ns = nonSpikes)
#art <- ches_brigandine(, Class, taxa01 %>% filter(Phylum == "Arthropoda") %>% pull(Class) %>% unique(), min = 4, thresh = 0, ns = nonSpikes)


## Arthropods
art1 <- ches_brigandine(Order, Class, taxa01 %>%  pull(Class) %>% unique() %>% na.omit(), thresh = 10^5, min = 3, max = 5,
                        ns = nonSpikes %>% filter(Phylum == "Arthropoda"))
art1
art2 <- ches_brigandine_L2(Tag_ASV, Order, taxa01 %>%  pull(Order) %>% unique() %>% na.omit(), thresh = 10^4, min = 0, max = 4,
                        ns = nonSpikes %>% filter(Phylum == "Arthropoda"))
art2 + facet_nested(Class + Order ~ Depth, drop = TRUE, scales = "free", space = "free") +
  labs(y = "ASV")

## Other stuff
cil <- ches_brigandine_L2(Tag_ASV, Order, taxa01 %>%  pull(Order) %>% unique() %>% na.omit(), thresh = 10^4, min = 3, max = 5,
                        ns = nonSpikes %>% filter(Phylum == "Ciliophora"))
cil + facet_nested(Class + Order ~ Depth, drop = TRUE, scales = "free", space = "free") +
  labs(y = "ASV")

## Plot euk phylum
plot_euk_phylum <- function(loc_phylum){
plt <- ches_brigandine_L2(Tag_ASV, Order, taxa01 %>%  pull(Order) %>% unique() %>% na.omit(), thresh = 10^4, min = 0, max = 4,
                        ns = nonSpikes %>% filter(Phylum == loc_phylum))
plt + facet_nested(Class + Order ~ Depth, drop = TRUE, scales = "free", space = "free") +
  labs(y = "ASV")
}

plot_euk_phylum("Mollusca")

plot_euk_phylum("Ctenophora")

plot_euk_phylum("Annelida")
