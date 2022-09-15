library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries( returnclass = "sf", scale = 10)
class(world)



library(ggrepel)

stations <- read_csv(here::here("stations.csv"))
stations01 <- stations %>% mutate(Station = str_sub(Station, start = 3)) %>%
  mutate(isWest = str_detect(Station, "W")) %>% filter(!isWest) %>% mutate(Station = str_remove(Station, "C")) %>%
  select(-isWest) %>%
  #filter(Station %in% c("4.3", "3.3")) %>%
  identity()

cbMap <- ggplot(data = world ) + geom_sf(color = "grey30", fill = "grey90", alpha = 0.5) +
  #coord_sf(xlim = c(-77.0, -76.15), ylim = c(37.6, 39.3)) +
  coord_sf(xlim = c(-76.9, -76.0), ylim = c(37.6, 39.5)) + # 77
  
  scale_x_continuous(breaks = seq(from = -77.2, to = -75.8, by = 0.4), 
                     minor_breaks = seq(from = -77.2, to = -75.8, by = 0.2)) + # not sure why this has no affect
  scale_y_continuous(breaks = seq(from = 37, to = 40, by = 0.4)) +
  geom_point(aes(y = lat, x = long, fill = Station, shape = Station), data = stations01, size = 4) + 
  geom_label_repel(aes(y = lat, x = long, label = Station),
                   force = 2,
                   nudge_x = -0.35,
                   #hjust = 3,
                    data = stations01) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        #panel.grid.minor = element_line(colour = "red", size = 0.5)
        ) +
  scale_fill_viridis_d(option = "plasma") + scale_shape_manual(values = c(21:25, 21)) +
  labs(x = "Longitude", y = "Latitude")
cbMap

ggsave(here("Figures", "cbMap.png"), height = 6, width = 4)
