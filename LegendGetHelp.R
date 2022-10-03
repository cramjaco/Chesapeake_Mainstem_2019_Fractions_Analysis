library(tidyverse)

stuff <- tibble(
  Ex = rnorm(100),
  Why = rnorm(100),
  Size = rep(c(0.2, 1.2, 5, 20, 53, 180, 500), length = 100),
  Amount = runif(100, min = 3, max = 7)
)

stuff %>%
ggplot(aes(x = Ex, y = Why, size = sqrt(Size), fill = Amount)) +
  geom_point(shape = 21) +
  scale_radius(breaks = sqrt(c(1.2, 5, 20, 53, 180, 500)), labels = c(1.2, 5, 20, 53, 180, 500), range = c(1, 12)) +
  theme(legend.position = "bottom") +
  guides(size = guide_legend(nrow = 1, label.position = "bottom", title.vjust = .8), legend.box.margin = margin(1000,0,0,0))

# Close but the legends are a little different in height