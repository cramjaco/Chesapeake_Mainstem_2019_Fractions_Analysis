# Agglomerate my working object by taxonomic level

source("InitialProcessing.R")

aglom <- function(taxlevel, broadlevel, ns = nonSpikes20){
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns %>% group_by(!!taxlevel, ID) %>%
  summarise(across(.cols = c(!!broadlevel, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum)) %>%
    filter(!is.na(!!taxlevel))
}

clip <- function(broadlevel, broadkeep, ns = nonSpikes20){
  broadlevel = enquo(broadlevel)
  #broadkeep = enquo(broadkeep)
  
  ns %>%
    filter(!!broadlevel == broadkeep)
}

# test
nsClass <- aglom(Class, c(Kingdom, Phylum))

## find bacterial taxa that are abundant in the attached fraction
find_attached <- function(adf, threshold = 10^6){
  #taxlevel <- enquo(taxlevel)
  
  attached <- adf %>%
  mutate(copiesPerMg = copiesPerL/MassperLiter) %>%
  filter(Size_Class >= 1.2) %>%
  group_by_at(1) %>%
  summarise(max_copiesPerMg = max(copiesPerMg)) %>%
  mutate(isAttached = max_copiesPerMg > threshold) %>%
  filter(isAttached) %>%
  .[[1]]
  
  attached
}

# test find_attached
find_attached(nsClass)

ches_brigandine <- function(taxlevel, broadlevel, broadkeep, ns = nonSpikes20, min = NA, max = NA, thresh = 10^6){
  minCpMg <- min # 4
  maxCpMg <- max # 7.5
  
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns_loc0 <- aglom(!!taxlevel, !!broadlevel, ns)
  ns_loc <- ns_loc0 %>%
    filter(!!broadlevel %in% broadkeep)
  
  attached_loc <- find_attached(ns_loc, threshold = thresh)
  
  ## Pre plotting
  
  toPlot <- ns_loc %>%
    ungroup() %>%
    arrange(-Size_Class) %>% # this stops working with more than one phylum # are there NA values?
    filter(!!taxlevel %in% attached_loc, Size_Class >= 1.2) %>%
    mutate(copiesPerMg = copiesPerL/MassperLiter) %>%
    mutate(copiesPerMg = case_when(
      log10(copiesPerMg) > maxCpMg ~ 10^(maxCpMg),
      log10(copiesPerMg) < minCpMg ~ 10^(minCpMg),
      TRUE ~ copiesPerMg
    )) %>%
    mutate(!!quo_name(taxlevel) := fct_rev(!!taxlevel))
  
  # Plotting
  
  # set breaks and labels
  if(is.na(minCpMg)){minCpMg = log10(min(toPlot$copiesPerMg))}
  if(is.na(maxCpMg)){maxCpMg = log10(max(toPlot$copiesPerMg))}
  loc_breaks = seq(from = minCpMg, to = maxCpMg, by = 1)
  loc_labels = loc_breaks
  if(!is.na(min)){loc_labels[1] = paste0("≤", loc_labels[1])}
  if(!is.na(max)){loc_labels[length(loc_labels)] = paste0("≥", loc_labels[length(loc_labels)])}
  
  
  # main
  locPlot <- toPlot %>%
  ggplot(aes(x = as.factor(Station), y = !!taxlevel, size = sqrt(Size_Class), fill = log10(copiesPerMg))) +
  geom_point(shape = 21, color = "black", stroke = .5) +
  scale_radius(breaks = sqrt(c(1.2, 5, 20, 53, 180, 500)), labels = c(1.2, 5, 20, 53, 180, 500), range = c(1, 12)) +
  scale_fill_viridis_c(breaks = loc_breaks, labels = loc_labels) +
  facet_grid(rows = vars(!!broadlevel), cols= vars(Depth), drop = TRUE, scales = "free", space = "free") +
  labs(y = quo_name(taxlevel), x = "Station", size = "Size Class", fill = "log10(Copies/mg)") +
  theme_bw() +
    theme(strip.text.y = element_text(angle = 0))

  locPlot
  
}

# test ches_brigandine
ches_brigandine(Phylum, Kingdom, taxa01 %>% pull(Kingdom) %>% unique(), min = 4, max = 7)
ches_brigandine(Class, Phylum, c("Actinobacteria", "Proteobacteria"))
ches_brigandine(Family, Class, c("Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))
ches_brigandine(Genus, Family, c("Rhodobacteraceae", "Rickettsiaceae"))

ches_brigandine(Genus, Order, taxa01 %>% filter(Class == "Alphaproteobacteria") %>% pull(Order) %>% unique())
ches_brigandine(Tag_ASV, Order, taxa01 %>% filter(Class %in% "Alphaproteobacteria") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 2 * 10^6)
ches_brigandine(Tag_ASV, Order, taxa01 %>% filter(Class == "Gammaproteobacteria") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 2 * 10^6)

# for clara
ches_brigandine(Tag_ASV, Order, taxa01 %>% filter(Phylum == "Planctomycetes") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 10^6, min = 3)
ches_brigandine_L(Tag_ASV, Order, taxa01 %>% filter(Phylum == "Planctomycetes") %>% pull(Order) %>% unique() %>% na.omit(), thresh = 10^6, min = 3)


ches_brigandine(ASV, Order,
                taxa01 %>% filter(Phylum == "Bacteroidetes") %>% pull(Order) %>% unique() %>% na.omit,
                min = 3,
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV))

# I need better ASV names -- eg lowest level named taxonomy plus ASV number.

# something strange is up, but only with multiple phyla. It looks like some size classes go missing. %in% not -- in the filter. Don't filter by c(two things)
# this is almost good behavior, but I'd actually prefer to filter at one level, and facet at another




#### as above but per liter

ches_brigandine_L <- function(taxlevel, broadlevel, broadkeep, ns = nonSpikes20, min = NA, max = NA, thresh = 10^6){
  minCpL <- min # 4
  maxCpL <- max # 7.5
  
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns_loc0 <- aglom(!!taxlevel, !!broadlevel, ns)
  ns_loc <- ns_loc0 %>%
    filter(!!broadlevel %in% broadkeep)
  
  attached_loc <- find_attached(ns_loc, threshold = thresh)
  
  ## Pre plotting
  
  toPlot <- ns_loc %>%
    ungroup() %>%
    arrange(-Size_Class) %>% # this stops working with more than one phylum # are there NA values?
    filter(!!taxlevel %in% attached_loc) %>%
    #mutate(copiesPerMg = copiesPerL/MassperLiter) %>%
    mutate(copiesPerL = case_when(
      log10(copiesPerL) > maxCpL ~ 10^(maxCpL),
      log10(copiesPerL) < minCpL ~ 10^(minCpL),
      TRUE ~ copiesPerL
    )) %>%
    mutate(!!quo_name(taxlevel) := fct_rev(!!taxlevel))
  
  toPlotFree <- toPlot %>% filter(Size_Class == 0.2)
  toPlotAttached <- toPlot %>% filter(Size_Class >= 1.2)
  
  # Plotting
  
  # set breaks and labels
  if(is.na(minCpL)){minCpL = log10(min(toPlot$copiesPerL))}
  if(is.na(maxCpL)){maxCpL = log10(max(toPlot$copiesPerL))}
  loc_breaks = seq(from = minCpL, to = maxCpL, by = 1)
  loc_labels = loc_breaks
  if(!is.na(min)){loc_labels[1] = paste0("≤", loc_labels[1])}
  if(!is.na(max)){loc_labels[length(loc_labels)] = paste0("≥", loc_labels[length(loc_labels)])}
  
  
  # main
  locPlot <- toPlot %>%
  ggplot() +
  geom_point(shape = 22, color = "black", stroke = .5, size = 12, data = toPlotFree,
             aes(x = as.factor(Station), y = !!taxlevel, fill = log10(copiesPerL))) +
  geom_point(shape = 21, color = "black", stroke = .5, data = toPlotAttached,
             aes(x = as.factor(Station), y = !!taxlevel, fill = log10(copiesPerL), size = sqrt(Size_Class))) +
  
  scale_radius(breaks = sqrt(c(1.2, 5, 20, 53, 180, 500)), labels = c(1.2, 5, 20, 53, 180, 500), range = c(1, 10)) +
  scale_fill_viridis_c(breaks = loc_breaks, labels = loc_labels) +
  facet_grid(rows = vars(!!broadlevel), cols= vars(Depth), drop = TRUE, scales = "free", space = "free") +
  labs(y = quo_name(taxlevel), x = "Station", size = "Size Class", fill = "log10(Copies/L)") +
  theme_bw() +
    theme(strip.text.y = element_text(angle = 0))

  locPlot
  #list(toPlotFree, toPlotAttached)
}

ches_brigandine_L(Phylum, Kingdom, taxa01 %>% pull(Kingdom) %>% unique(), min = 4, max = 9)
# I need to figure out how to do the free living stuff
# facet_grid call is breaking just the free living data, for reasons I don't understand
# Its like its dropping everything or something

## ok, I needed to pass the aesthetics to each plotting elmeent and then things were ok?

ches_brigandine_L(ASV, Order,
                taxa01 %>% filter(Phylum == "Bacteroidetes") %>% pull(Order) %>% unique() %>% na.omit,
                min = 3, max = 8, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV))

ches_brigandine_L(ASV, Order,
                taxa01 %>% filter(Class == "Alphaproteobacteria") %>% pull(Order) %>% unique() %>% na.omit,
                min = 3, max = 8, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV), thresh = 10^5)

# doesn't work, for reasons I don't understand
ches_brigandine_L(ASV, Genus,
                taxa01 %>% filter(Order == "Rhodobacterales") %>% pull(Genus) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV), thresh = 10^5)
# Huh, its the rhodobacters that like intermediate sized particles

## I'd like to identify all of the ASVs that are more abundant on any size fraction
## greater than or equal to 5 microns than they are in the 0.2 or 1.2 size fractions.
## Maybe I can just ask which species are most abundant, on average, on which size fractions

summarise(across(.cols = c(!!broadlevel, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum))

# across(.cols = Kingdom:Genus, .fns = first)
mainHabitats <- nonSpikes20 %>%
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

ches_brigandine_L(ASV, Phylum,
                taxa %>% filter(Tag_ASV %in% snowDwellers) %>% pull(Phylum) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV %in% snowDwellers) %>%
                  mutate(ASV = recode(ASV, `Arachnida;7` = "Acartia;7")),
                thresh = 10^5)

## as above but without euks and chloroplasts

mainHabitats_prok <- nonSpikes20 %>%
  filter(Kingdom != "Eukaryota" & Order != "Chloroplast") %>%
  group_by(Tag_ASV, Size_Class) %>% 
  summarise(copiesPerL = mean(copiesPerL)) %>%
  ungroup() %>% group_by(Tag_ASV) %>%
  arrange(Tag_ASV, -copiesPerL) %>%
  top_n(1) %>% ungroup()

mainHabitats_prok %>% 
  group_by(Size_Class) %>%
  summarize(ntaxa = length(copiesPerL))

snowDwellers_prok <- mainHabitats_prok %>% 
  filter(Size_Class >= 20) %>%
  pull(Tag_ASV)

snowDwellers_plt <- ches_brigandine_L(ASV, Class,
                taxa %>% filter(Tag_ASV %in% snowDwellers_prok) %>% pull(Class) %>% unique() %>% na.omit,
                min = 3, max = 7, 
                ns = nonSpikes20 %>% rename(ASV0 = ASV, ASV = Tag_ASV) %>% 
                  filter(ASV %in% snowDwellers_prok) %>%
                  mutate(ASV = recode(ASV, `Arachnida;7` = "Acartia;7")),
                thresh = 10^5)
snowDwellers_plt

ggsave(here::here("Figures", "SnowDwellers.png"),
       snowDwellers_plt,
       width = 8, height = 4, units = "in")

## Chrysomicrobia tend to be isolated from sediment
# Midocloria are symbionts, consume mitocondria, 
# previously only in ticks (is it associated with acartia tonsa?)