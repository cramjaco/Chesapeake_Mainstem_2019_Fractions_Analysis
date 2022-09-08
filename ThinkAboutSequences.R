ASVSeq <- read_lines(here("ASVs.fa"))

ASVs20 <- nonSpikes20 %>% pull(ASV) %>% unique()

sequenceDf <- enframe(ASVSeq) %>%
  mutate(type = rep(c("ASV_headers", "Seq"), length.out =  n())) %>%
  group_by(type) %>%
  mutate(id = row_number()) %>%
  pivot_wider(id_cols = id, names_from = type, values_from = "value") %>%
  mutate(ASV = str_remove(ASV_headers, ">")) %>%
  select(-id) %>%
  #head() %>%
  identity()

sequenceNS20 <- sequenceDf %>%
  filter(ASV %in% ASVs20)

preFasta20 <- with(sequenceNS20, c(rbind(ASV_headers, Seq))) 
write_lines(preFasta20, file = "ASVs_NS20.fa")
