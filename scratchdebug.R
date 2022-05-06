# scratch

blue1 <- ns0 %>% 
  filter(Class %in% c( "Acidimicrobiia", "Alphaproteobacteria"), Station == 3.2, Depth == "Bottom") %>%
  group_by(!!taxlevel, ID) %>%
  summarise(across(.cols = c(!!broadlevel, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum)) %>%
    filter(!is.na(!!taxlevel))

blue2 <- ns0 %>% 
  group_by(!!taxlevel, ID) %>%
  summarise(across(.cols = c(!!broadlevel, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum)) %>%
    filter(!is.na(!!taxlevel)) %>%
  filter(Class %in% c( "Acidimicrobiia", "Alphaproteobacteria"), Station == 3.2, Depth == "Bottom")
  
ns_loc2 <- ns_loc0 %>% 
  ungroup()%>%
    filter(!!broadlevel == broadkeep)
