nsPhylum <- nonSpikes20 %>% group_by(Phylum, ID) %>%
  summarise(across(.cols = c(Kingdom, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum)) %>%
  filter(!is.na(Phylum))


threshold = 10 * 10^6
attachedPhyla <- nsPhylum %>%
  mutate(copiesPerMg = copiesPerL/MassperLiter) %>%
  #filter(Size_Class >= 1.2) %>%
  group_by(Phylum) %>%
  summarise(max_copiesPerMg = max(copiesPerMg)) %>%
  mutate(isAttached = max_copiesPerMg > threshold) %>%
  filter(isAttached) %>%
  pull(Phylum)

minCpMg <- 4
maxCpMg <- 7.5

toPlotPhylum <- nsPhylum %>%
  arrange(-Size_Class) %>%
  #mutate(Phylum = factor(Phylum, levels = rev(unique(Phylum)))) %>%
  mutate(Phylum = fct_rev(Phylum)) %>%
  select(Phylum, Kingdom, Station, Size_Class, Depth, copiesPerL) %>%
  filter(Size_Class %in% c(0.2, 20, 53, 180)) %>%
  filter(Phylum %in% c("Arthropoda", "Actinobacteria", "Bacteroidetes")) %>%
  filter(Depth != "Oxy") %>%
  filter(Station %in% c(3.1, 3.2))


  ggplot(aes(x = as.factor(Station),
             y = Phylum,
             fill = log10(copiesPerL)),
         data = toPlotPhylum) +
  geom_point(shape = 22, color = "black", stroke = .5, size = 20,
             data = toPlotPhylum %>% filter(Size_Class == 0.2)) +
    geom_point(shape = 21, color = "black", stroke = .5, aes(size = log(Size_Class)),
             data = toPlotPhylum %>% filter(Size_Class != 0.2)) +
  scale_radius(breaks = log(c( 20, 53, 180)), labels = c(20, 53, 180), range = c(2, 10)) +
  scale_fill_viridis_c() +
  facet_grid(Kingdom~Depth, drop = TRUE, scales = "free", space = "free") + 
  labs(y = "Phylum", x = "Station", size = "Size Class", fill = "log10(Copies/mg)") +
  theme_bw()
  

        ###     data = toPlotPhylum %>% filter(Size_Class == 0.2)
  
    ggplot(aes(x = as.factor(Station),
             y = Phylum,
             fill = log10(copiesPerL)),
         data = toPlotPhylum) +
  geom_point(shape = 21, color = "black", stroke = .5, size = 10, data = toPlotPhylum %>% filter(Size_Class == 0.2)) #+
    # geom_point(shape = 21, color = "black", stroke = .5, aes(size = log(Size_Class)),
    #          data = toPlotPhylum %>% filter(Size_Class != 0.2)) +
  #scale_radius(breaks = log(c( 20, 53, 180)), labels = c(20, 53, 180), range = c(1, 12)) +
  #scale_fill_viridis_c() +
  #facet_grid(Kingdom~Depth, drop = TRUE, scales = "free", space = "free") + 
  #labs(y = "Phylum", x = "Station", size = "Size Class", fill = "log10(Copies/mg)") +
  #theme_bw()
