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

EC_predicted <- read_tsv(here("PiCrustStuff", "EC_predicted.tsv.gz"))
EC_nonzero <- EC_predicted %>%
  pivot_longer(-sequence, names_to = "EC") %>%
  filter(value > 0) %>%
  pivot_wider(names_from = EC, values_from = value, values_fill = 0)

colnames(EC_predicted)[-1] -> ECs
ECs %>% str_replace(":", " ") -> ECs2      

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

# # Blows everyting up and isn't used again
# ECs_tib <- ECs_tib0 %>%
#   #head() %>%
#   mutate(kdata = map(EC, ~possibly(my_kegGet, otherwise = "error here")(.x))) %>%
#   unnest(kdata)

## Here I look for nitrite reductase containing bugs

#keggFind("enzyme", c(" nitrate reductase")) -> red
keggFind("enzyme", c(" nitrite reductase")) -> red

red_ec <- names(red)
red_ec2 <- str_replace(red_ec, ":", " ")
red_ec3 <- str_replace(red_ec, "ec", "EC")

redDf <- EC_nonzero %>% 
  select(sequence, any_of(red_ec3)) %>%
  mutate(sum = rowSums(across(contains("EC")))) %>%
  filter(sum >=1)

redDf$sequence -> redASVs

nonSpikes20 %>% filter(ASV %in% redASVs)
redTaxa <- taxa01 %>% filter(ASV %in% redASVs)

# subset just dentitrifiers, then look at community structure variability.
# copies/mg
# Do they hang out at 4.3 bottom and on particles?
# I think I need toseperate assimilatory and desimilatory reduction

keggGet("EC 1.7.7.1")

keggGet(red_ec2) -> kgls


## Brigandine plots

red_plot <- ches_brigandine(ASV, Phylum,
                taxa %>% filter(ASV %in% redASVs) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% redASVs & Phylum != "Cyanobacteria"),
                thresh = 1.25 * 10^6)
red_plot

red_plot_L <- ches_brigandine_L(ASV, Phylum,
                taxa %>% filter(ASV %in% redASVs) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% redASVs & Phylum != "Cyanobacteria"),
                thresh = 1.25 * 10^6)
red_plot_L
## What about total gene abundance (for each gene)

redDf

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

ches_plot_options <- list(
  scale_y_log10nice() ,
    scale_x_log10(breaks = my_sizes, labels = as.character(my_sizes)) ,
  geom_point(size = 2) ,
  geom_path(aes(color = as.factor(Station))) ,
  scale_shape_manual(values = rep(21:25, 2)) ,
  scale_fill_viridis_d(option = "plasma") ,
  scale_color_viridis_d(option = "plasma")
)

ecDf2 %>%
  arrange(Size_Class) %>%
  filter(Depth != "Oxy", !is.na(EC)) %>%
  ggplot(aes(y = EC_hits_per_mg, x = Size_Class,
             fill = as.factor(Station), shape = as.factor(Station))) +
  facet_grid(Depth ~ EC) +
  ches_plot_options+
  theme_bw()

ecDf2 %>%
  arrange(Size_Class) %>%
  filter(Depth != "Oxy", !is.na(EC)) %>%
  ggplot(aes(y = EC_hits_per_L, x = Size_Class,
             fill = as.factor(Station), shape = as.factor(Station))) +
  facet_grid(Depth ~ EC) +
  ches_plot_options +
  theme_bw()

# Nitrate reducing orgnaisms

## Other organissms
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

nitriteRed <- keggFindAndProcess("dissimilatory nitrite reductase") # empty
nitrateRed <- keggFindAndProcess(" nitrate reductase")
nitriteOx <- keggFindAndProcess(" nitrite oxidase") # empty
nitrateOx <- keggFindAndProcess(" nitrate oxidase") # empty

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



keggFindAndProcess("methanogenisis") %>% .$enzymes %>% length()

keggFindAndProcess(" formylmethanofuran:tetrahydromethanopterin formyltransferase") %>% .$enzymes %>% length()

methano <- keggFindAndProcess(" formylmethanofuran:tetrahydromethanopterin formyltransferase")

nonSpikes$Phylum %>% unique() %>% .[str_detect(.,"Verstraetearchaeota")] %>% na.omit() # we don't have any of these either

methano_order <- nonSpikes$Order %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit()
 nonSpikes$Class %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit()
 nonSpikes$Order %>% unique() %>% .[str_detect(.,"Methano")] %>% na.omit()
 nonSpikes$Family %>% unique() %>% .[str_detect(.,"Methano")]
  nonSpikes$Genus %>% unique() %>% .[str_detect(.,"Methano")]

  nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                   filter(Order %in% methano_order) %>% group_by(ASV) %>% summarise(mra = max(copiesPerL))
  
  nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                   filter(Order %in% methano_order) %>% filter(ASV0 == "ASV_12803") %>% arrange(RA) %>% View()
  
methano_plot_L <- ches_brigandine_L(ASV, Order,
                NULL,
                min = 3, max = 7, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(Order %in% methano_order) %>% 
                  #filter(ASV %in% "Methanofastidiosales;12803") %>%
                  identity(),
                thresh = 0, threshL = 10)
methano_plot_L



## Clara suggestions
# oxidizes methane
keggFindAndProcess("particulate methane monooxygenase") %>% .$enzymes %>% length()
keggFindAndProcess("particulate methane monooxygenase") %>% .$ASV %>% length()
pmoa <- keggFindAndProcess("particulate methane monooxygenase")

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


keggFindAndProcess("meth") %>% .$enzymes %>% length()

keggFindAndProcess("Methyl Coenzyme M Reductase") %>% .$enzymes %>% length()

keggFindAndProcess("Methyl Coenzyme M Reductase") %>% .$enzymes %>% length()
keggFindAndProcess("Methyl Coenzyme M Reductase") %>% .$ASV %>% length()
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

keggFindAndProcess("EC:2.8.4.1") # no asvs with that gene
amo <- keggFindAndProcess("ammonia monooxygenase") # has one asv

amo_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% amo$ASV)  %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% amo$ASV & Phylum != "Cyanobacteria"),
                thresh = 0 * 10^5)
amo_plot_L # meh, also where is MG1?

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

# Cowplot all of those together

sulfred_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% sulfateRed$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% sulfateRed$ASV & Phylum != "Cyanobacteria"),
                thresh = 2 * 10^4)
sulfred_plot_L

nitrite_oxidizers <- nonSpikes$Genus %>% unique() %>% .[str_detect(.,"Nitro")]

pmoa_plot_L <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(ASV %in% pmoa$ASV) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 6, 
                ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV0 %in% pmoa$ASV & Phylum != "Cyanobacteria"),
                thresh = 2 * 10^4)
pmoa_plot_L

ammonium_oxidizers_plot <- ches_brigandine_L(ASV, Kingdom,
taxa %>% filter(Genus %in% nitrite_oxidizers) %>% pull(Kingdom) %>% unique() %>% na.omit,
min = 3, max = 6, 
ns = nonSpikes %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
filter(Genus %in% nitrite_oxidizers[-1]),
thresh = 0 * 10^4)

ammonium_oxidizers_plot

phylaLs <- nonSpikes$Phylum %>% unique()
phylaLs[str_detect(phylaLs, "Planc")]

keggFindAndProcess("hydrazine dehydrogenase") %>% .$enzymes %>% length()
keggFindAndProcess("hydrazine dehydrogenase") %>% .$ASV %>% length()
pmoa <- keggFindAndProcess("particulate methane monooxygenase")

genusLs <- nonSpikes$Phylum %>% unique()
# https://en.wikipedia.org/wiki/Anammox # Jetten
anammoxGenera <- c("Scalindua", "Brocadia", "Kuenenia", "Anammoxoglobus", "Jettenia")
genusLs[genusLs %in% anammoxGenera]

cowplot::plot_grid(
  sulfred_plot_L + theme(legend.position="none"),
  pmoa_plot_L + theme(legend.position="none"),
  ammonium_oxidizers_plot ,
  ncol = 1,
  rel_heights = c(3,0.7,2.2), 
  labels = LETTERS
)

ggsave(here("Figures", "BiogeochemicalCyclers.png"), height = 11, width = 7)
# now all with the same thresholds
# I need to mess with margins

## I want a prettier plot and I plan to use nested facets to get there
library(ggh4x)

# Notes: 
# Desulfovibreo are oxidizing sulfate
# Sedimenticola are probably doing the opposite, reducing HS
# Chromatiales are purple sulfur bacteria. They use H2S and light to get high energy electrons for growing
# Clara thinks they shouldn't be particle attached. Maybe they are just big. (Zaar ea. 2003)
# Rhodospirillaceae are purple nonsulfur bacteria -- they use reduced sulfur containing compounds, they are photosynthetic
# purpul sulfur use sulfide, purple nonsulfer use organic compounds

nsBiogeo <- nonSpikes %>%
  mutate(Biogeo = case_when(
    ASV %in% sulfateRed$ASV ~ "Sulfur Cycling",
    # ASV %in% pmoa$ASV ~ "Methanotrophy",
    # str_detect(Order, "Methano") ~ "Methanogenisis",
    ASV %in% pmoa$ASV ~ "Methane Cycling",
    str_detect(Order, "Methano") ~ "Methane Cycling",
    str_detect(Genus, "Nitro") ~ "Nitrogen Cycling"
  )) %>%
  filter(!is.na(Biogeo)) %>%
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
  
  # mutate(Class = case_when(
  #   Biogeo == "Nitrogen Cycling" ~ case_when(
  #   Class == "Nitrososphaeria" ~ "Nitrososphaeria\n(Archaea)\n(Ammonium Oxidizing)",
  #   Class == "Gammaproteobacteria" ~ "Gammaproteobacteria\n(Ammonium Oxidizing)",
  #   Class == "Nitrospira" ~ "Nitrospira\n(Nitrite Oxidizing)",
  #   TRUE ~ Class
  #   ),
  #   Biogeo == "Sulfur Cycling" ~ case_when(
  #     Class == "Alphaproteobacteria" ~ "Alphaproteobacteria\n(Purple Nonsulfur)",
  #     Class == "Deltaproteobacteria" ~ "Deltaproteobacteria\n(Sulfate Oxidizing)",
  #     TRUE ~ Class
  #   ),
  #   TRUE ~ Class
  # )) %>%
  # mutate(Family = case_when(
  #   (is.na(Family) & Order == "Chromatiales") ~ "Unknown\nChromatiales",
  #   TRUE ~ Family
  # ))

biogeo_plot <- ches_brigandine_L2(ASV, Biogeo, NULL,
min = 3, max = 6, 
ns = nsBiogeo %>% rename(ASV0 = ASV, ASV = Tag_ASV),
threshL = 3 * 10^5, thresh = 0) + # cuts cools stuff
  facet_nested(Biogeo +  Class  + Family ~ Depth, drop = TRUE, scales = "free", space = "free") +
  theme(panel.spacing = unit(2, "points"))
  

bigoeo_plot_annotated <- cowplot::ggdraw(biogeo_plot) +
  cowplot::draw_label("Process", x = 0.93, y = .98) +
  cowplot::draw_label("Class", x = 0.805, y = .98) +
  cowplot::draw_label("Family", x = 0.685, y = .98)
  

ggsave(here("Figures", "BiogeochemicalCyclers_better1.png"), height = 8, width = 11, plot = bigoeo_plot_annotated)
# Whats up with ASV 6673 and why is it replacing my good looking bugs
# Oh dear. I can't see my 5 micron size class when the 1.2 is visible.
# I can see it if my line thickness is 0.2 or less.
# But really I should make the 180, 53 and 20 bigger
# Ok. I cube rooted everything. And that works. Should cascade forward for other plots.

# I might want an extended version for the supplement with more things

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
               