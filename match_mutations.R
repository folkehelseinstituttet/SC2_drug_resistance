library(tidyverse)


# RdRP --------------------------------------------------------------------

stanf_rdrp <- read_csv("data/2024.01.12_RdRP_inhibitors.csv") %>% 
  select("stanford_mutation" = Mutation) %>% 
  mutate(stanford_position  = as.numeric(str_sub(stanford_mutation, 2, -2)),
         stanford_reference = str_sub(stanford_mutation, 1, 1),
         stanford_mutated   = str_sub(stanford_mutation, -1, -1)) %>% 
  mutate(nextclade_position = stanford_position - 9) %>% 
  add_column("Gene" = "RdRP")

# Read nextclade mutations
# Nanopore 2023
mutatations_raw <- read_tsv("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2023/OppsettNGSSEQ20231219-01B_summaries/OppsettNGSSEQ20231219-01B_summaries_and_Pangolin.csv")

nextclade_ORF1b_mutations <- mutatations_raw %>% 
  # Get only ORF1b
  select(name, ORF1b) %>% 
  mutate("tmp" = str_split(ORF1b, ";")) %>% 
  # Unnest the list column (which is basically another pivot_longer)
  unnest(tmp) %>%
  # Remove tmp
  select(-ORF1b) %>% 
  rename("ORF1b" = tmp) %>% 
  mutate(nextclade_position  = as.numeric(str_sub(ORF1b, 2, -2)),
         nextclade_reference = str_sub(ORF1b, 1, 1),
         nextclade_mutated   = str_sub(ORF1b, -1, -1))  
  
  
# Add RdRP resistance mutations to our sequences analyzed by nextclade. 
# NA means that the resistance mutation is not present
# NB! Check that it's correct
tmp <- left_join(nextclade_ORF1b_mutations, stanf_rdrp, by = "nextclade_position") %>% 
  # Rename NA in Gene column to "RdRP"
  replace_na(list(Gene = "RdRP"))


# 3CLpro ------------------------------------------------------------------

stanf_3clpro <- read_csv("data/2024.01.12_3CLpro_inhibitors.csv") %>% 
  select("stanford_mutation" = Mutation) %>%
  # Change del to "-"
  mutate(stanford_mutation = str_replace(stanford_mutation, "del", "-")) %>% 
  mutate(stanford_position  = as.numeric(str_sub(stanford_mutation, 2, -2)),
         stanford_reference = str_sub(stanford_mutation, 1, 1),
         stanford_mutated   = str_sub(stanford_mutation, -1, -1)) %>% 
  mutate(nextclade_position = stanford_position + 3263) %>% 
  add_column("Gene" = "3CLpro")

nextclade_ORF1a_mutations <- mutatations_raw %>% 
  # Get only ORF1a
  select(name, ORF1a) %>% 
  mutate("tmp" = str_split(ORF1a, ";")) %>% 
  # Unnest the list column (which is basically another pivot_longer)
  unnest(tmp) %>%
  # Remove tmp
  select(-ORF1a) %>% 
  rename("ORF1a" = tmp) %>% 
  mutate(nextclade_position  = as.numeric(str_sub(ORF1a, 2, -2)),
         nextclade_reference = str_sub(ORF1a, 1, 1),
         nextclade_mutated   = str_sub(ORF1a, -1, -1))  


# Add 3CLpro resistance mutations to our sequences analyzed by nextclade. 
# NA means that the resistance mutation is not present
# NB! Check that it's correct
tmp2 <- left_join(nextclade_ORF1a_mutations, stanf_3clpro, by = "nextclade_position") %>% 
  # Rename NA in Gene column to "3CLpro"
  replace_na(list(Gene = "3CLpro"))


# Join RdRP and 3CLpro together -------------------------------------------
df <- bind_rows(tmp, tmp2)

# For pretty sharing:
# pivot_wider on gene - values are resistance mutations. Might need to do separately for RdRP and 3CLpro.

# Read the 3CLpro inhibitors
CLPro_inhib <- read_csv("data/2023.01.05_3CLpro_inhibitors.csv")
