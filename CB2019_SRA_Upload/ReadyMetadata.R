library(tidyverse)
library(here)
library(lubridate)

key2 <- read_csv(here("Keys", "SampleKey2.csv"))

stations <- read_csv(here("stations.csv"))

station_dates <- ymd(c(rep("2019-07-24", 2), rep("2019-07-23", 4), rep("2019-07-22", 3)))

stations1 <- bind_cols(stations, date = station_dates)

ctd_data <- read_csv(here("Keys", "Use_CTD2.csv")) %>%
  select(Station, Depth = Zone, name, value) %>%
  mutate(Depth = str_extract(Depth, "^\\w")) %>%
  pivot_wider(names_from = name, values_from = value)

stations_cram <- stations1 %>%
  filter(!str_detect(Station, "W")) %>%
  rename(CBStation = Station) %>%
  mutate(Station = parse_number(CBStation))
stations_cram

size_fractions <- tibble(
  lb = c(0.2, 1.2, 5, 20, 53, 180, 500)
) %>%
  mutate(SizeRange = if_else(lb == 500, "≥500", paste0(lb, "-", lead(lb))))

depths <- read_csv(here("Keys", "SampleDepths.csv")) %>% 
  pivot_longer(-Station, names_to = "DepthCat", values_to = "Depth_m") %>%
  mutate(Depth = str_extract(DepthCat, "^\\w")) %>%
  mutate(Depth_m = abs(Depth_m))

sample_types <- tibble(
  Type = c("Sample", "PControl", "EControl", "GD", "Mock"),
  LType = c("Sample", "PCR Control", "DNA Extraction Control", "Generous Donor", "Mock Community")
)

sample_metadata <- key2 %>%
  filter(ReadDir == 1) %>%
  select(Sample, Station:Type) %>%
  left_join(stations_cram %>%
              select(Station, CBStation, lat, long, date)
              , by = "Station") %>%
  left_join(sample_types, by = "Type") %>%
  left_join(depths %>% select(-DepthCat), by = c("Station", "Depth")) %>%
  mutate(lat_lon = paste0(lat, "N", "_", abs(long), "W")) %>%
  left_join(size_fractions, by = c("Size" = "lb")) %>%
  left_join(ctd_data, by = c("Station", "Depth"))

# I'm getting a non unique samples error
sample_metadata %>% group_by(Sample) %>% summarise(nSamples = n()) %>% arrange(nSamples)

file_metadata <- key2 %>%
  mutate(title = case_when(
    Type == "Sample" ~ paste("Station CB", Station, " Depth", Depth, "m", "SizeFrac", Size, "μm"),
    Type == "PControl" ~ paste("PCR Control", Sample),
    Type == "EControl" ~ paste("DNA Extraction Control", Sample),
    Type == "GD" ~ paste("Generous Donor Sample (HPL Pier)", Sample),
    Type == "Mock" ~ paste("Mock Community", Sample)
  )
  )%>%
  select(Sample, ReadDir, RenFiles, title) %>%
  pivot_wider(names_from = ReadDir, values_from = RenFiles) %>%
  rename(filename1 = `1`, filename2 = `2`) %>%
  identity()
file_metadata  

write_csv(sample_metadata, here("CB2019_SRA_Upload", "sample_metadata.csv"))
write_csv(file_metadata, here("CB2019_SRA_Upload", "file_metadata.csv"))
