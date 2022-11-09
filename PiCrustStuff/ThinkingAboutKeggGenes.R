library(KEGGREST)
library(broom)
listDatabases()
org <- keggList("organism")
head(org)
keggGet("EC 1.1.1.102") -> eceg

eceg %>% class()

keggGet("1.1.1.102")

EC_predicted <- read_tsv(here("PiCrustStuff", "EC_predicted.tsv.gz"))

colnames(EC_predicted)[-1] -> ECs
ECs %>% str_replace(":", " ") -> ECs2      

my_kegGet <- function(EC){
  print(EC)
  EC_list <- keggGet(EC)
  EC_tib <- with(EC_list[[1]],
       tibble(class1 = CLASS[1], class2 = CLASS[2], class3 = CLASS[3], name = SYSNAME)
  )
  EC_tib
}

test_EC_list <- my_kegGet("EC 1.1.1.103")

ECs_tib0 <- tibble(EC = ECs2)

ECs_tib <- ECs_tib0 %>%
  #head() %>%
  mutate(kdata = map(EC, my_kegGet)) %>%
  unnest(kdata)
