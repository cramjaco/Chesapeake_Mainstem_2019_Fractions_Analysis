## ThinkingAboutKeggGenes.R had the original exploration of data.
# This file exists to provide just the figures used in the manuscript.

library(KEGGREST)
library(broom)
library(here)
load(here::here("RDataFiles", "InitialProcessing_3.RData"))
source(here::here("RLibraries", "ChesapeakePersonalLibrary.R"))
source(here::here("RLibraries", "Brigandine_Library.R"))
# source(here::here("RScripts", "InitialProcessing_2.R"))

listDatabases()
#org <- keggList("organism")
#head(org)
keggGet("EC 1.1.1.102") -> eceg

eceg %>% class()

#keggGet("1.1.1.102")

# These are the predicted genes in my dataset. They don't have nice names, 
# so we're using Kegg to name things.

EC_predicted <- read_tsv(here("PiCrustStuff", "EC_predicted.tsv.gz"))
# And only keeping the genes that show up in our data
EC_nonzero <- EC_predicted %>%
  pivot_longer(-sequence, names_to = "EC") %>%
  filter(value > 0) %>%
  pivot_wider(names_from = EC, values_from = value, values_fill = 0)

# The enzyme identifiers
colnames(EC_predicted)[-1] -> ECs
#Enzyme identifiers formatted differentlsy
ECs %>% str_replace(":", " ") -> ECs2      

# For a given EC this returns its name
my_kegGet <- function(EC){
  print(EC)
  EC_list <- keggGet(EC)
  EC_tib <- with(EC_list[[1]],
       tibble(class1 = CLASS[1], class2 = CLASS[2], class3 = CLASS[3], name = SYSNAME)
  )
  EC_tib
}

test_EC_list <- my_kegGet("EC 1.1.1.103")

ECs_tib0 <- tibble(EC = ECs2)


# The first thing I did was look for nitrate reductase containing bugs
# This analysis didn't go into the manuscript because we think these genes
# get horizontally transferred and don't conserve with taxonomy

# Kegg find looks for enzymes where any of their synonyms conain some genes
keggFind("enzyme", c(" nitrite reductase")) -> red

# Processing the output names
red_ec <- names(red)
red_ec2 <- str_replace(red_ec, ":", " ")
red_ec3 <- str_replace(red_ec, "ec", "EC")

# A table of which ASVs have which genes
redDf <- EC_nonzero %>% 
  select(sequence, any_of(red_ec3)) %>%
  mutate(sum = rowSums(across(contains("EC")))) %>%
  filter(sum >=1)

# The asvs whose closest relative has the gene of interes
redDf$sequence -> redASVs

# taxonomic information about those ASVs
redTaxa <- taxa01 %>% filter(ASV %in% redASVs)

# subset just dentitrifiers, then look at community structure variability.
# copies/mg
# Do they hang out at 4.3 bottom and on particles?
# I think I need toseperate assimilatory and desimilatory reduction


## Brigandine plots
# Here we're plotting the abundance of things

# Normalized to total particle mass
red_plot <- ches_brigandine(ASV, Phylum,
                taxa %>% filter(ASV %in% redASVs) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% redASVs & Phylum != "Cyanobacteria"),
                thresh = 1.25 * 10^6)
red_plot

# Normalized to volume
red_plot_L <- ches_brigandine_L(ASV, Phylum,
                taxa %>% filter(ASV %in% redASVs) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% redASVs & Phylum != "Cyanobacteria"),
                thresh = 1.25 * 10^6)
red_plot_L
## What about total gene abundance (for each gene)

## Here I wanted to see how abundant the different nitrite reductase genes
## were across the different stations
# 
cpl20 <- nonSpikes20 %>% select(ASV, ID, copiesPerL)
ecDf <- left_join(cpl20,
          redDf %>% pivot_longer(-sequence, names_to = "EC", values_to = "EC_hits"),
          by = c("ASV" = "sequence")
) %>%
  mutate(EC_hits_per_L = EC_hits * copiesPerL) %>%

  group_by(EC, ID) %>%
  summarise(EC_hits_per_L = sum(EC_hits_per_L))

ecDf2 <- left_join(sample, ecDf, by = "ID") %>%
  mutate(EC_hits_per_mg = EC_hits_per_L / MassperLiter)

## some plotting options, to save space typing later.
ches_plot_options <- list(
  scale_y_log10nice() ,
    scale_x_log10(breaks = my_sizes, labels = as.character(my_sizes)) ,
  geom_point(size = 2) ,
  geom_path(aes(color = as.factor(Station))) ,
  scale_shape_manual(values = rep(21:25, 2)) ,
  scale_fill_viridis_d(option = "plasma") ,
  scale_color_viridis_d(option = "plasma")
)

## Per mg particle
ecDf2 %>%
  arrange(Size_Class) %>%
  filter(Depth != "Oxy", !is.na(EC)) %>%
  ggplot(aes(y = EC_hits_per_mg, x = Size_Class,
             fill = as.factor(Station), shape = as.factor(Station))) +
  facet_grid(Depth ~ EC) +
  ches_plot_options+
  theme_bw()

## Per liter
ecDf2 %>%
  arrange(Size_Class) %>%
  filter(Depth != "Oxy", !is.na(EC)) %>%
  ggplot(aes(y = EC_hits_per_L, x = Size_Class,
             fill = as.factor(Station), shape = as.factor(Station))) +
  facet_grid(Depth ~ EC) +
  ches_plot_options +
  theme_bw()

# Nitrate reducing orgnaisms

## A function to search the Kegg database of enzymes for a gene name
## and then return those enzymes and ASVs in our database that have those genes

keggFindAndProcess <- function(fname){
  keggFind("enzyme", fname) -> red
  
  red_ec <- names(red)
  red_ec2 <- str_replace(red_ec, ":", " ")
  red_ec3 <- str_replace(red_ec, "ec", "EC")
  
  redDf <- EC_nonzero %>% 
    select(sequence, any_of(red_ec3)) %>%
    mutate(sum = rowSums(across(contains("EC")))) %>%
    filter(sum >=1)
  
  redDf$sequence -> redASVs
  
  return(list(enzymes = red, ec = red_ec, ec2 = red_ec2, ec3 = red_ec3, df = redDf, ASV = redASVs))
}

# Looking for genes
nitriteRed <- keggFindAndProcess("dissimilatory nitrite reductase") # empty
nitrateRed <- keggFindAndProcess(" nitrate reductase") # lots of things, but didn't persue
nitriteOx <- keggFindAndProcess(" nitrite oxidase") # empty
nitrateOx <- keggFindAndProcess(" nitrate oxidase") # empty

###################
## Sulfur Cycling##
###################

## dissimilatory sulfite reductase is used in sulfur cycling bacteria that do lots
## of different functions aroudn sulfur. These did all make it into the manuscript

sulfateRed <- keggFindAndProcess("dissimilatory sulfite reductase") # note, this is sulfite
sulfateRed$enzymes %>% length()
sulfateRed$ASV %>% length()

sulfred_plot <- ches_brigandine(ASV, Class,
                taxa %>% filter(ASV %in% sulfateRed$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% sulfateRed$ASV & Phylum != "Cyanobacteria"),
                thresh = 1 * 10^0)
sulfred_plot

sulfred_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% sulfateRed$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% sulfateRed$ASV & Phylum != "Cyanobacteria"),
                thresh = 1 * 10^0)
sulfred_plot_L

####################
## Methane Cycling##
####################

keggFindAndProcess("methanogenisis") %>% .$enzymes %>% length() # empty

# Methanotrophy
methano <- keggFindAndProcess(" formylmethanofuran:tetrahydromethanopterin formyltransferase")

#There is one gene with this name, which is used in methaotrophy
methano %>% .$enzymes %>% length()

## Methanogenisis
# Verstraetearchaeota is a new archeal group that do methanogenisis
# They may be too new to be in the silva database we used
nonSpikes$Phylum %>% unique() %>% .[str_detect(.,"Verstraetearchaeota")] %>% na.omit() # we don't have any of these either

## All the other methanogens have the prefix -Methano
## we search for that here 

methano_order <- nonSpikes$Order %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit()
 nonSpikes$Class %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit() # One class Methanomicrobia
 nonSpikes$Order %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit() # Two orders Methanomicrobiales and Methaofastidiosales
nonSpikes$Family %>% unique() %>% .[str_detect(.,"Methano")] # One family Methanomicrobiaceae
  nonSpikes$Genus %>% unique() %>% .[str_detect(.,"Methano")] # No

  #But are the bacteria from the two methogenic orders actually in our dataset?
  
  # yes, one is Methanofastidiosales -- the other has 0% mean relative abundance, and so is ignored
  # why is it still in the data? Beats me, maybe it was in a mock or something.
  nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                   filter(Order %in% methano_order) %>% group_by(ASV) %>% summarise(mra = max(copiesPerL))
  
# Plot out the methanogen (per volume)
methano_plot_L <- ches_brigandine_L(ASV, Order,
                NULL,
                min = 3, max = 7, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(Order %in% methano_order) %>% 
                  #filter(ASV %in% "Methanofastidiosales;12803") %>%
                  identity(),
                thresh = 0, threshL = 10)
methano_plot_L


## Mining a little more
## Clara suggests pmoa which is used in methanotrophy
# oxidizes methane
keggFindAndProcess("particulate methane monooxygenase") %>% .$enzymes %>% length() # The gene exists
keggFindAndProcess("particulate methane monooxygenase") %>% .$ASV %>% length() # we have one!
pmoa <- keggFindAndProcess("particulate methane monooxygenase")

# A plot of this one methanotroph Methylomonaceae; ASV # 373
pmoa_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% pmoa$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% pmoa$ASV & Phylum != "Cyanobacteria"),
                thresh = 0 * 10^5)
pmoa_plot_L

pmoa_plot <- ches_brigandine(ASV, Class,
                taxa %>% filter(ASV %in% pmoa$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% pmoa$ASV & Phylum != "Cyanobacteria"),
                thresh = 0 * 10^5)
pmoa_plot

## A gene focused approach to methanogenisis
## Didn't work because we got a *bunch* of false hits
keggFindAndProcess("meth") %>% .$enzymes %>% length() # 892

keggFindAndProcess("Methyl Coenzyme M Reductase") %>% .$enzymes %>% length() # 9

keggFindAndProcess("Methyl Coenzyme M Reductase") %>% .$ASV %>% length() # 608
methano <- keggFindAndProcess("Methyl Coenzyme M Reductase")

methano_plot_L_pi <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% methano$ASV & Class == "Bacteroidia") %>%
                  #filter(Kingdom == "Archaea") %>%
                  pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% methano$ASV & Phylum != "Cyanobacteria"),
                thresh = 0 * 10^6)
methano_plot_L_pi
# I don't believe this plot


# ammonia monooxygenase should be useful for identifying ammonium oxidizing archaea (and bacteria)
# But actually it didn't return anything, except somehow the same methaonogen we saw earleir?
keggFindAndProcess("EC:2.8.4.1") # no asvs with that gene
amo <- keggFindAndProcess("ammonia monooxygenase") # has one asv

amo_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% amo$ASV)  %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% amo$ASV & Phylum != "Cyanobacteria"),
                thresh = 0 * 10^5)
amo_plot_L # meh, also where is MG1?

## Trying nitrite oxidase
nitriteox <- keggFindAndProcess("EC 1.7.5.1  ")


nitriteox_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% nitriteox$ASV)  %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% nitriteox$ASV & Phylum != "Cyanobacteria"),
                thresh = 1 * 10^6)
nitriteox_plot_L  # this gene looks too general

## so lets show the following:
# dissimilatory sulfite reductase havers
# the one bug with pmoa
# Anyone with Nitro in their genus name (ammonium oxidizers)

## Combining the different approaches

# Remake sulfate reducers
sulfred_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% sulfateRed$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% sulfateRed$ASV & Phylum != "Cyanobacteria"),
                thresh = 2 * 10^4)
sulfred_plot_L

nitrite_oxidizers <- nonSpikes$Genus %>% unique() %>% .[str_detect(.,"Nitro")]

# Remake methanotrophs
pmoa_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% pmoa$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% pmoa$ASV & Phylum != "Cyanobacteria"),
                thresh = 2 * 10^4)
pmoa_plot_L

# Possible ammonium oxidizers
# Threre are lots! Many not so aboundant
ammonium_oxidizers_plot <- ches_brigandine_L(ASV, Kingdom,
taxa %>% filter(Genus %in% nitrite_oxidizers) %>% pull(Kingdom) %>% unique() %>% na.omit,
min = 3, max = 6, 
ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
filter(Genus %in% nitrite_oxidizers[-1]),
thresh = 0 * 10^4)

ammonium_oxidizers_plot

## Looking for anammox
keggFindAndProcess("hydrazine dehydrogenase") %>% .$enzymes %>% length() # gene exists
keggFindAndProcess("hydrazine dehydrogenase") %>% .$ASV %>% length() # We don't have any

# Lets search for bacteria known to be from anommox genera
genusLs <- nonSpikes$Phylum %>% unique()
# https://en.wikipedia.org/wiki/Anammox # Jetten
anammoxGenera <- c("Scalindua", "Brocadia", "Kuenenia", "Anammoxoglobus", "Jettenia")
genusLs[genusLs %in% anammoxGenera] # we don't have any


## Plottin everything together
## I want a prettier plot and I plan to use nested facets to get there
library(ggh4x)

# Notes: 
# Desulfovibreo are oxidizing sulfate
# Sedimenticola are probably doing the opposite, reducing HS
# Chromatiales are purple sulfur bacteria. They use H2S and light to get high energy electrons for growing
# Clara thinks they shouldn't be particle attached. Maybe they are just big. (Zaar ea. 2003)
# Rhodospirillaceae are purple nonsulfur bacteria -- they use reduced sulfur containing compounds, they are photosynthetic
# purpul sulfur use sulfide, purple nonsulfer use organic compounds

## Take nonSpikes, with all of the ASV info, and just save the ASVs that I think are involved in 
## biogeochemical processes, based on the stuff above.
## Anr record which things do which processes

nsBiogeo <- nonSpikes %>%
  mutate(Biogeo = case_when(
    ASV %in% sulfateRed$ASV ~ "Sulfur Cycling",
    # ASV %in% pmoa$ASV ~ "Methanotrophy",
    # str_detect(Order, "Methano") ~ "Methanogenisis",
    ASV %in% pmoa$ASV ~ "Methane Cycling",
    str_detect(Order, "Methano") ~ "Methane Cycling",
    str_detect(Genus, "Nitro") ~ "Nitrogen Cycling"
  )) %>%
  # Get rid of the non-biogeochemical cyclers
  filter(!is.na(Biogeo)) %>%
  # Clean up the names, by providing more information about each organism.
  mutate(Family = case_when(
    (Biogeo == "Methane Cycling" & Order == "Methanofastidiosales" & is.na(Family)) ~ "Unknown\nMethanofastidiosales",
    Biogeo == "Nitrogen Cycling" ~ case_when(
      #Family %in% c("Nitrosopumilaceae", "Nitrosomonadaceae") ~ paste0(Family, "\n(Ammonium Oxidizing)"),
      #Family %in% c("Nitrospiraceae") ~ paste0(Family, "\n(Nitrite Oxidizing)"),
      
      TRUE ~ Family
    ),
    Biogeo == "Sulfur Cycling" ~ case_when(
      Family %in% c("Rhodospirillaceae") ~ paste0(Family, "\n(Purple nonsulfur)"),
      Family %in% c("Sedimenticolaceae") ~ paste0(Family, "\n(Sulfur Oxidizing)"),
      Order == "Chromatiales" ~ "Unk. Chromatiales\n(Purple Sulfur)",
      TRUE ~ Family
    ),
    TRUE ~ Family
  ),
  Class = case_when(
    (Biogeo == "Nitrogen Cycling" & Class == "Nitrososphaeria") ~ paste0(Class, "\n(Archaea)\n(Ammonium Oxidizing)"),
    (Biogeo == "Nitrogen Cycling" & Class == "Gammaproteobacteria") ~ paste0(Class, "\n(Ammonium Oxidizing)"),
    (Biogeo == "Nitrogen Cycling" & Class == "Nitrospira") ~ paste0(Class, "\n(Nitrite Oxidizing)"),
    (Biogeo == "Sulfur Cycling" & Class == "Deltaproteobacteria") ~ paste0(Class, "\n(Sulfate Reducing)"),
    (Biogeo == "Methane Cycling" & Class == "Thermococci") ~ paste0(Class, "\n(Methanogenisis)"),
    (Biogeo == "Methane Cycling" & Class == "Gammaproteobacteria") ~ paste0(Class, "\n(Methanotrophy)"),
    TRUE ~ Class
  )
           ) %>%
  mutate(Class = factor(Class,
                           levels = c(
                             "Nitrososphaeria\n(Archaea)\n(Ammonium Oxidizing)",
                             "Thermococci\n(Methanogenisis)",
                             "Alphaproteobacteria",
                             "Gammaproteobacteria\n(Ammonium Oxidizing)",
                             "Gammaproteobacteria\n(Methanotrophy)",
                             "Gammaproteobacteria",
                             "Deltaproteobacteria\n(Sulfate Reducing)",
                              "Nitrospira\n(Nitrite Oxidizing)"
  ))) %>%
  identity()
  

# Make the compound plot -- its to big to display in rstudio so just
# make the object
biogeo_plot <- ches_brigandine_L2(ASV, Biogeo, NULL,
min = 3, max = 6, 
ns = nsBiogeo %>% rename(ASV0 = ASV, ASV = Tag_ASV),
threshL = 3 * 10^5, thresh = 0) + # cuts cools stuff
  facet_nested(Biogeo +  Class  + Family ~ Depth, drop = TRUE, scales = "free", space = "free") +
  theme(panel.spacing = unit(2, "points"))
  
## Add some annotations so the reader can follow the facet labels
bigoeo_plot_annotated <- cowplot::ggdraw(biogeo_plot) +
  cowplot::draw_label("Process", x = 0.93, y = .98) +
  cowplot::draw_label("Class", x = 0.805, y = .98) +
  cowplot::draw_label("Family", x = 0.685, y = .98)
  
# Save out the plot
ggsave(here("Figures", "BiogeochemicalCyclers_better1.png"), height = 8, width = 11, plot = bigoeo_plot_annotated)
# Whats up with ASV 6673 and why is it replacing my good looking bugs
# Oh dear. I can't see my 5 micron size class when the 1.2 is visible.
# I can see it if my line thickness is 0.2 or less.
# But really I should make the 180, 53 and 20 bigger
# Ok. I cube rooted everything. And that works. Should cascade forward for other plots.

# I might want an extended version for the supplement with more things

## Here is a more detailed plot. It has a lower threshold for including things
## Bacteria have to be 10^5 copies/L in some sample, rather than 3 * 10^5 copies.
biogeo_plot_detailed <- ches_brigandine_L2(ASV, Biogeo, NULL,
min = 3, max = 6, 
ns = nsBiogeo %>% rename(ASV0 = ASV, ASV = Tag_ASV),
threshL = 1 * 10^5, thresh = 0) + # cuts cools stuff
  facet_nested(Biogeo +  Class  + Family ~ Depth, drop = TRUE, scales = "free", space = "free") +
  theme(panel.spacing = unit(2, "points"))

bigoeo_plot_detailed_annotated <- cowplot::ggdraw(biogeo_plot_detailed) +
  cowplot::draw_label("Process", x = 0.93, y = .98) +
  cowplot::draw_label("Class", x = 0.805, y = .98) +
  cowplot::draw_label("Family", x = 0.685, y = .98)
ggsave(here("Figures", "BiogeochemicalCyclers_detailed1.png"), height = 8, width = 11, plot = bigoeo_plot_detailed_annotated)
               