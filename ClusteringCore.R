
messy_plot <- function(df){
  df %>% ggplot(aes(y = copiesPerL/Bin_Size, x = Size_Class, color = Depth)) + geom_point(shape = 1) +
    facet_grid(cluster ~ Station) + scale_color_manual(values = c("darkgreen", "Blue", "Black"))+
    theme(axis.text.x = element_text(angle = 90, size = 6)) +
    scale_x_log10() + scale_y_log10()
}

wideASVs <- nonSpikes20 %>%
  mutate(copiesPerLPermm = copiesPerL/Bin_Size) %>%
  select(ID, ASV, copiesPerLPermm) %>%
  pivot_wider(id_cols = ID, names_from = ASV, values_from = copiesPerLPermm) %>%
  column_to_rownames("ID") %>% na.omit()
# 

wideASVs_log <- log(wideASVs + 1)

wideASVs_mtx <- as.matrix(wideASVs)

sdist <- (1-cor(wideASVs, method = "spear"))/2

sclust <- hclust(dist(sdist))

hgroups <- cutree(sclust, k = 50)
