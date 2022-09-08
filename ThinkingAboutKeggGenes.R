library(KEGGREST)
library(broom)
source("ChesapeakePersonalLibrary.R")
source(here::here("RScripts", "InitialProcessing_2.R"))
listDatabases()
org <- keggList("organism")
head(org)
keggGet("EC 1.1.1.102") -> eceg

eceg %>% class()

keggGet("1.1.1.102")

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

ECs_tib <- ECs_tib0 %>%
  #head() %>%
  mutate(kdata = map(EC, my_kegGet)) %>%
  unnest(kdata)

#keggFind("enzyme", c(" nitrate reductase")) -> red
keggFind("enzyme", c(" nitrite reductase")) -> red

red_ec <- names(red)
red_ec2 <- str_replace(red_ec, ":", " ")
red_ec3 <- str_replace(red_ec, "ec", "EC")

redDf <- EC_nonzero %>% 
  select(sequence, any_of(red_ec3)) %>%
  mutate(sum = rowSums(across(contains("EC")))) %>%
  filter(sum >=1)
keggGet("EC 1.7.7.2")

redDf$sequence -> redASVs

nonSpikes20 %>% filter(ASV %in% redASVs)
redTaxa <- taxa01 %>% filter(ASV %in% redASVs)

# subset just dentitrifiers, then look at community structure variability.
# copies/mg
# Do they hang out at 4.3 bottom and on particles?
# I think I need toseperate assimilatory and desimilatory reduction

keggGet("EC 1.7.7.1")

keggGet(red_ec2) -> kgls

redASVs

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
