library(here)
library(skimr)
library(janitor)
library(dplyr)
library(palmerpenguins)
data(palmerpenguins)
skim_without_charts(penguins)
glimpse(penguins)
head(penguins)
penguins %>% 
  rename(island_new=island) %>% 
  select(-sex)
clean_names(penguins)
