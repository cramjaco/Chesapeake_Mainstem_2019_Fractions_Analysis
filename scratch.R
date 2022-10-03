counts_long %>%
  filter(Kingdom != "Spike" & !is.na(Kingdom)) %>%
  left_join(correctionData, by = "ID") %>%
   left_join(counts_ra %>% select(ASV, ID, RA), by = c("ASV", "ID")) %>%
  mutate(copiesPerL = reads * conversionMultiplier) %>%
   filter(!is.na(copiesPerL)) %>%
  # filter(SpikeReads > 0) %>%
  pull(ID) %>% unique()
