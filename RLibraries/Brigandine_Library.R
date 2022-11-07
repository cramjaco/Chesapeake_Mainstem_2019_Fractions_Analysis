
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
  mutate(isAttached = max_copiesPerMg >= threshold) %>%
  filter(isAttached) %>%
  .[[1]]
  
  attached
}

find_present <- function(adf, threshold = 0){
  #taxlevel <- enquo(taxlevel)
  
  present <- adf %>%
  #mutate(copiesPerMg = copiesPerL/MassperLiter) %>%
  #filter(Size_Class >= 1.2) %>%
  group_by_at(1) %>%
  summarise(max_copiesPerL = max(copiesPerL)) %>%
  mutate(isAttached = max_copiesPerL >= threshold) %>%
  filter(isAttached) %>%
  .[[1]]
  
  present
}


# Brigandine plot relative to mass
ches_brigandine <- function(taxlevel, broadlevel, broadkeep, ns = nonSpikes20, min = NA, max = NA, thresh = 10^6, thresL = 0){
  minCpMg <- min # 4
  maxCpMg <- max # 7.5
  
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns_loc0 <- aglom(!!taxlevel, !!broadlevel, ns)
  ns_loc <- ns_loc0 %>%
    filter(!!broadlevel %in% broadkeep)
  
    attached_loc <- find_attached(ns_loc, threshold = thresh) 
  present_loc <- find_present(ns_loc, threshold = threshL)
  keep_loc <- intersect(attached_loc, present_loc)
  
  ## Pre plotting
  
  toPlot <- ns_loc %>%
    ungroup() %>%
    arrange(-Size_Class) %>% # this stops working with more than one phylum # are there NA values?
    filter(!!taxlevel %in% keep_loc, Size_Class >= 1.2) %>%
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
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "bottom",
          legend.box.margin = margin(0,0,0.0)) + # legend.box.margin doesn't do anything different
    guides(size = guide_legend(nrow = 1, label.position = "bottom", title.vjust = .7))

  locPlot
  
}

## Plot relative to volume

ches_brigandine_L <- function(taxlevel, broadlevel, broadkeep, ns = nonSpikes20, min = NA, max = NA, thresh = 10^6, threshL = 0){
  minCpL <- min # 4
  maxCpL <- max # 7.5
  
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns_loc0 <- aglom(!!taxlevel, !!broadlevel, ns)
  
  if(!is.null(broadkeep)){
  ns_loc <- ns_loc0 %>%
    filter(!!broadlevel %in% broadkeep)
  } else {
    ns_loc <- ns_loc0
  }
  
  attached_loc <- find_attached(ns_loc, threshold = thresh) 
  present_loc <- find_present(ns_loc, threshold = threshL)
  keep_loc <- intersect(attached_loc, present_loc)
  
  ## Pre plotting
  
  toPlot <- ns_loc %>%
    ungroup() %>%
    arrange(-Size_Class) %>% # this stops working with more than one phylum # are there NA values?
    filter(!!taxlevel %in% keep_loc) %>%
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
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "bottom") +
    guides(size = guide_legend(nrow = 1, label.position = "bottom", title.vjust = .7))

  locPlot
  #list(toPlotFree, toPlotAttached)
}

## For two level plots
aglom2 <- function(taxlevel, broadlevel, ns = nonSpikes20){
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns %>% group_by(!!taxlevel, ID) %>%
  summarise(across(.cols = c(!!broadlevel, Kingdom:Genus, Station:ParticlesPerLiter), .fns = first),
            across(.cols = c(RA, copiesPerL), .fns = sum)) %>%
    filter(!is.na(!!taxlevel))
}

ches_brigandine_L2 <- function(taxlevel, broadlevel, broadkeep, ns = nonSpikes20, min = NA, max = NA, thresh = 10^6, threshL = 0){
  minCpL <- min # 4
  maxCpL <- max # 7.5
  
  taxlevel <- enquo(taxlevel)
  broadlevel <- enquo(broadlevel)
  ns_loc0 <- aglom2(!!taxlevel, !!broadlevel, ns)
  
  if(!is.null(broadkeep)){
  ns_loc <- ns_loc0 %>%
    filter(!!broadlevel %in% broadkeep)
  } else {
    ns_loc <- ns_loc0
  }
  
  attached_loc <- find_attached(ns_loc, threshold = thresh) 
  present_loc <- find_present(ns_loc, threshold = threshL)
  keep_loc <- intersect(attached_loc, present_loc)
  
  ## Pre plotting
  
  toPlot <- ns_loc %>%
    ungroup() %>%
    arrange(-Size_Class) %>% # this stops working with more than one phylum # are there NA values?
    filter(!!taxlevel %in% keep_loc) %>%
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
  geom_point(shape = 21, color = "black", stroke = .3, data = toPlotAttached,
             aes(x = as.factor(Station), y = !!taxlevel, fill = log10(copiesPerL), size = (Size_Class)^(1/3))) +
  
  scale_radius(breaks = (c(1.2, 5, 20, 53, 180, 500))^(1/3), labels = c(1.2, 5, 20, 53, 180, 500), range = c(1, 10)) +
  scale_fill_viridis_c(breaks = loc_breaks, labels = loc_labels) +
  facet_grid(rows = vars(!!broadlevel), cols= vars(Depth), drop = TRUE, scales = "free", space = "free") +
  labs(y = quo_name(taxlevel), x = "Station", size = "Size Class", fill = "log10(Copies/L)") +
  theme_bw() +
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "bottom") +
    guides(size = guide_legend(nrow = 1, label.position = "bottom", title.vjust = .7))

  locPlot
  #list(toPlotFree, toPlotAttached)
}

