


key2 <- key2 %>% mutate(filterLink = paste0("X", FilteredFiles))

countNames <- colnames(counts)[2:length(colnames(counts))]

betterNames <- tibble(countNames) %>% left_join(key2, by = c("countNames" = "filterLink")) %>% 
  pull(ID)

colnames(counts) <- betterNames

spikes <- na.omit(counts[taxa$Kingdom == "Spike",])

apply(spikes, 2, sum)
