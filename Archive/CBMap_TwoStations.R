# Just the two stations for Sairah Malkin's grant

library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")
library("here")

world <- ne_states( returnclass = "sf")
class(world)



library(ggrepel)

stations <- read_csv(here::here("InputData", "stations.csv"))
stations01 <- stations %>% mutate(Station = str_sub(Station, start = 3)) %>%
  mutate(isWest = str_detect(Station, "W")) %>% filter(!isWest) %>% mutate(Station = str_remove(Station, "C")) %>%
  select(-isWest) %>%
  filter(Station %in% c("4.3", "3.3"))

cbMap <- ggplot(data = world ) + geom_sf(color = "grey30", fill = "grey90") + coord_sf(xlim = c(-76.65, -76.15), ylim = c(37.6, 39.3)) +
  scale_x_continuous(breaks = seq(from = -76.8, to = -76.0, by = 0.2)) +
  scale_y_continuous(breaks = seq(from = 37, to = 40, by = 0.2)) +
  geom_point(aes(y = lat, x = long, fill = Station, shape = Station), data = stations01, size = 4) + 
  geom_label_repel(aes(y = lat, x = long, label = Station),
                   force = 2,
                   nudge_x = -0.25,
                   #hjust = 3,
                    data = stations01) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_viridis_d(option = "plasma") + scale_shape_manual(values = c(21:25, 21)) +
  labs(x = "Longitude", y = "Latitude")
cbMap

ggsave("cbMap.png", height = 6, width = 4)
