library(tidyverse)


# TODO --------------------------------------------------------------------

# Prepare input data ------------------------------------------------------

# Load the resistance associated mutations
stanf_3clpro <- read_csv("N:/Virologi/JonBrate/Drug resistance SARSCOV2/data/2024.01.12_3CLpro_inhibitors.csv") %>% 
  rename("stanford_res" = Mutation) %>%
  # Rename the paxlovid fold
  rename("NTV_fold" = `NTV: fold`) %>% 
  # Change del to "-"
  mutate(stanford_res = str_replace(stanford_res, "del", "-")) %>% 
  mutate(stanford_position  = as.numeric(str_sub(stanford_res, 2, -2)),
         stanford_reference = str_sub(stanford_res, 1, 1),
         stanford_mutated   = str_sub(stanford_res, -1, -1)) %>% 
  mutate(nextclade_position = stanford_position + 3263)

# Load the resistance associated mutations
stanf_rdrp <- read_csv("N:/Virologi/JonBrate/Drug resistance SARSCOV2/data/2024.01.12_RdRP_inhibitors.csv") %>% 
  rename("stanford_res" = Mutation) %>%
  # Rename the remdesivir fold
  rename("RDV_fold" = `RDV: fold`) %>% 
  # Change del to "-"
  mutate(stanford_res = str_replace(stanford_res, "del", "-")) %>% 
  mutate(stanford_position  = as.numeric(str_sub(stanford_res, 2, -2)),
         stanford_reference = str_sub(stanford_res, 1, 1),
         stanford_mutated   = str_sub(stanford_res, -1, -1)) %>% 
  mutate(nextclade_position = stanford_position - 9)

# Get Nextclade mutations

# year to analyze
year <- "2023"

nanopore_files <- list.files(paste0("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/", year, "/"),
                    pattern = "summaries_and_Pangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)


sihf_files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/SIHF/",
                    pattern = "AndPangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)

# Extract only 2023 files
pattern <- paste0("SIHF_hCoV_WGS\\d{6}_", year)
sihf_files <- sihf_files[sihf_files %>% str_detect(pattern)]

ahus_files <- list.files(paste0("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/Ahus/", year, "/"),
                    pattern = "AndPangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)

files <- c(nanopore_files, sihf_files, ahus_files)

# Read all the files and combine
res <- vector("list", length = length(files))
pb <- txtProgressBar(min = 0, max = length(files), initial = 0, style = 3) 
for (i in seq_along(files)) {
  try({
    res[[i]] <- read_tsv(files[i], show_col_types = F) %>% 
    select(name, ORF1a, ORF1b)
  }, silent = TRUE)
  setTxtProgressBar(pb, i)
}
close(pb)

# Combine everything into one big tibble
mutations_raw <- bind_rows(res) %>% 
  # Remove samples with NA for both ORF1a and ORF1b
  filter(!is.na(ORF1a) & !is.na(ORF1b)) %>% 
  # Fix Ahus names to match BN Key (i.e. remove everything after first "|" symbol)
  separate(name, into = c("name"), sep = "\\|", remove = TRUE) 



# nsp5 --------------------------------------------------------------------

# nsp5 is on ORF1a
nextclade_nsp5_mutations <- mutations_raw %>% 
  # Get only ORF1a
  select(name, ORF1a) %>%
  # Separate all the mutations. Gives a list for each sample
  mutate("tmp" = str_split(ORF1a, ";")) %>% 
  # Unnest the list column (which is basically another pivot_longer)
  # Get each mutation on a separate line in a column named tmp
  unnest(tmp) %>%
  # Remove tmp
  select(-ORF1a) %>% 
  rename("ORF1a" = tmp) %>% 
  # Tease out the position, reference and mutated amino acid
  mutate(nextclade_position  = as.numeric(str_sub(ORF1a, 2, -2)),
         nextclade_reference = str_sub(ORF1a, 1, 1),
         nextclade_mutated   = str_sub(ORF1a, -1, -1)) %>% 
  # Filter out only the nsp5 gene (aa 3264 to 3569)
  filter(nextclade_position >= 3264, nextclade_position <= 3569) %>% 
  # convert to stanford nomenclature
  mutate(stanford_position = nextclade_position - 3263) %>% 
  # Create column for nsp5 mutations
  unite("NSP5", c("nextclade_reference", "stanford_position", "nextclade_mutated"), sep = "")

# Add the resistance mutations to the nextclade mutations
# If any resistance mutations are present they will be in the "stanford_res" column

# To RAVN
RAVN_nsp5 <- left_join(nextclade_nsp5_mutations, stanf_3clpro, by = c("NSP5" = "stanford_res"), keep = TRUE) %>% 
  # Renaming the resistance column
  rename("NSP5_resistance" = "stanford_res") %>% 
  # Add info on strength of resistance
  mutate("type" = case_when(
    NTV_fold >= 10                 ~ "high_res",
    NTV_fold >= 5 | NTV_fold < 10  ~ "mid_res",
    NTV_fold >= 2.5 | NTV_fold < 5 ~ "very_mid_res",
    NTV_fold < 2.5                 ~ "white"
  )) %>% 
  select("sample" = name,
         NSP5_resistance) %>% 
  distinct() %>% 
  # Remove only NA's
  filter(!is.na(NSP5_resistance))

# To BN

df <- left_join(nextclade_nsp5_mutations, stanf_3clpro, by = c("NSP5" = "stanford_res"), keep = TRUE) %>% 
  # Renaming the resistance column
  rename("NSP5_resistance" = "stanford_res") %>% 
  # Add info on strength of resistance
  mutate("type" = case_when(
    NTV_fold >= 10                 ~ "pax_high_res",
    NTV_fold >= 5 | NTV_fold < 10  ~ "pax_mid_res",
    NTV_fold >= 2.5 | NTV_fold < 5 ~ "pax_low_res",
    NTV_fold < 2.5                 ~ "pax_very_low_res"
  ))

# Create empty data frame
to_BN_NSP5 <- df %>% 
  select(name) %>% distinct() %>%  
  add_column("pax_high_res" = NA,
             "pax_mid_res" = NA,
             "pax_low_res" = NA,
             "NSP5" = NA) %>% 
  # Change all columns to character columns
  mutate(across(everything(), as.character))

# First need to get all the NSP5 mutations in the same cell.
try(
  NSP5_all_muts <- df %>% 
    select(name, NSP5) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = NSP5) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("NSP5", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_all_muts, by = "name") %>% 
    # Coalesce the two NSP5 columns. WIll add from NSP5_all_muts if there is NA in to_BN_NSP5
    mutate(NSP5 = coalesce(NSP5.x, NSP5.y)) %>%
    select(-NSP5.x, -NSP5.y)
)

# Then only NSP5 resistance mutations for the different categories
try(
  NSP5_high_res <- df %>% 
    filter(type == "pax_high_res") %>% 
    select(name, NSP5_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_high_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_high_res, by = "name") %>% 
    # Coalesce the two pax_high_res columns. WIll add from NSP5_high_res if there is NA in to_BN_NSP5
    mutate(pax_high_res = coalesce(pax_high_res.x, pax_high_res.y)) %>%
    select(-pax_high_res.x, -pax_high_res.y)
)

try(
  NSP5_mid_res <- df %>% 
    filter(type == "pax_mid_res") %>% 
    select(name, NSP5_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_mid_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_mid_res, by = "name") %>% 
    # Coalesce the two pax_mid_res columns. WIll add from NSP5_mid_res if there is NA in to_BN_NSP5
    mutate(pax_mid_res = coalesce(pax_mid_res.x, pax_mid_res.y)) %>%
    select(-pax_mid_res.x, -pax_mid_res.y)
)

try(
  NSP5_low_res <- df %>% 
    filter(type == "pax_low_res") %>% 
    select(name, NSP5_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_low_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_low_res, by = "name") %>% 
    # Coalesce the two pax_low_res columns. WIll add from NSP5_low_res if there is NA in to_BN_NSP5
    mutate(pax_low_res = coalesce(pax_low_res.x, pax_low_res.y)) %>%
    select(-pax_low_res.x, -pax_low_res.y)
)


# RdRP --------------------------------------------------------------------

# NB! The RdRP code here is not relevant for BN import. 
# There is no paxlovid resistance for RdRP

# RdRP (nsp12) is on ORF1b
nextclade_rdrp_mutations <- mutations_raw %>% 
  # Get only ORF1b
  select(name, ORF1b) %>%
  # Separate all the mutations. Gives a list for each sample
  mutate("tmp" = str_split(ORF1b, ";")) %>% 
  # Unnest the list column (which is basically another pivot_longer)
  # Get each mutation on a separate line in a column named tmp
  unnest(tmp) %>%
  # Remove tmp
  select(-ORF1b) %>% 
  rename("ORF1b" = tmp) %>% 
  # Tease out the position, reference and mutated amino acid
  mutate(nextclade_position  = as.numeric(str_sub(ORF1b, 2, -2)),
         nextclade_reference = str_sub(ORF1b, 1, 1),
         nextclade_mutated   = str_sub(ORF1b, -1, -1)) %>% 
  # Filter out only the RdRP gene (aa 1 (actually RdRP starts 9 aa before this. Into ORF1a) to 923)
  filter(nextclade_position >= 1, nextclade_position <= 923) %>% 
  # convert to stanford nomenclature
  mutate(stanford_position = nextclade_position + 9) %>% 
  # Create column for nsp5 mutations
  unite("RdRP", c("nextclade_reference", "stanford_position", "nextclade_mutated"), sep = "")

# Add the resistance mutations to the nextclade mutations
# If any resistance mutations are present they will be in the "stanford_res" column


# To RAVN
RAVN_RdRP <- left_join(nextclade_rdrp_mutations, stanf_rdrp, by = c("RdRP" = "stanford_res"), keep = TRUE) %>% 
  # Renaming the resistance column
  rename("RdRP_resistance" = "stanford_res") %>% 
  # Add info on strength of resistance
  mutate("type" = case_when(
    RDV_fold >= 10                 ~ "high_res",
    RDV_fold >= 5 | RDV_fold < 10  ~ "mid_res",
    RDV_fold >= 2.5 | RDV_fold < 5 ~ "very_mid_res",
    RDV_fold < 2.5                 ~ "white"
  )) %>% 
  select("sample" = name,
         RdRP_resistance) %>% 
  distinct() %>% 
  # Remove only NA's
  filter(!is.na(RdRP_resistance))

# To BN

df <- left_join(nextclade_rdrp_mutations, stanf_rdrp, by = c("RdRP" = "stanford_res"), keep = TRUE) %>% 
  # Renaming the resistance column
  rename("RdRP_resistance" = "stanford_res") %>% 
  # Add info on strength of resistance
  mutate("type" = case_when(
    RDV_fold >= 10                 ~ "high_res",
    RDV_fold >= 5 | RDV_fold < 10  ~ "mid_res",
    RDV_fold >= 2.5 | RDV_fold < 5 ~ "very_mid_res",
    RDV_fold < 2.5                 ~ "white"
  ))

# Create empty data frame
to_BN_rdrp <- df %>% 
  select(name) %>% distinct() %>%  
  add_column("pax_high_res" = NA,
             "pax_mid_res" = NA,
             "pax_low_res" = NA,
             "rdrp" = NA) %>% 
  # Change all columns to character columns
  mutate(across(everything(), as.character))

# First need to get all the rdrp mutations in the same cell.
try(
  rdrp_all_muts <- df %>% 
    select(name, RdRP) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = RdRP) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("rdrp", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_rdrp <- left_join(to_BN_rdrp, rdrp_all_muts, by = "name") %>% 
    # Coalesce the two rdrp columns. WIll add from rdrp_all_muts if there is NA in to_BN_rdrp
    mutate(rdrp = coalesce(rdrp.x, rdrp.y)) %>%
    select(-rdrp.x, -rdrp.y)
)

# Then only rdrp resistance mutations for the different categories
try(
  rdrp_high_res <- df %>% 
    filter(type == "pax_high_res") %>% 
    select(name, RdRP_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = RdRP_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_high_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_rdrp <- left_join(to_BN_rdrp, rdrp_high_res, by = "name") %>% 
    # Coalesce the two pax_high_res columns. WIll add from rdrp_high_res if there is NA in to_BN_rdrp
    mutate(pax_high_res = coalesce(pax_high_res.x, pax_high_res.y)) %>%
    select(-pax_high_res.x, -pax_high_res.y)
)

try(
  rdrp_mid_res <- df %>% 
    filter(type == "pax_mid_res") %>% 
    select(name, RdRP_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = RdRP_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_mid_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_rdrp <- left_join(to_BN_rdrp, rdrp_mid_res, by = "name") %>% 
    # Coalesce the two pax_mid_res columns. WIll add from rdrp_mid_res if there is NA in to_BN_rdrp
    mutate(pax_mid_res = coalesce(pax_mid_res.x, pax_mid_res.y)) %>%
    select(-pax_mid_res.x, -pax_mid_res.y)
)

try(
  rdrp_low_res <- df %>% 
    filter(type == "pax_low_res") %>% 
    select(name, RdRP_resistance) %>% 
    add_column(dummy = "x") %>% 
    pivot_wider(names_from = dummy, values_from = RdRP_resistance) %>% 
    unnest_wider(x, names_sep = "_") %>% 
    unite("pax_low_res", starts_with("x_"), sep = ";", na.rm = TRUE)
)
try(
  to_BN_rdrp <- left_join(to_BN_rdrp, rdrp_low_res, by = "name") %>% 
    # Coalesce the two pax_low_res columns. WIll add from rdrp_low_res if there is NA in to_BN_rdrp
    mutate(pax_low_res = coalesce(pax_low_res.x, pax_low_res.y)) %>%
    select(-pax_low_res.x, -pax_low_res.y)
)


# Create BN import files --------------------------------------------------


# Create BN import files. Need to separate Ahus and Sihf because they will be imported on Key

fhi_to_BN_NSP5 <- to_BN_NSP5 %>% 
  # Extract FHI samples
  filter(str_detect(name, "SC2$")) %>% 
  rename("SEQUENCEID_NANO29" = "name")

write_csv(fhi_to_BN_NSP5, 
          file = paste0("C:/Users/jonr/OneDrive - Folkehelseinstituttet/Desktop/", format(Sys.Date(), "%Y.%m.%d"), "-FHI_samples_NSP5_resistance_BN_import.csv"),
          na = "")

external_to_BN_NSP5 <- to_BN_NSP5 %>% 
  # Extract FHI samples
  filter(str_detect(name, "Ahus") | str_detect(name, "SIHF")) %>%
  rename("KEY" = "name")

write_csv(external_to_BN_NSP5, 
          file = paste0("C:/Users/jonr/OneDrive - Folkehelseinstituttet/Desktop/", format(Sys.Date(), "%Y.%m.%d"), "-External_samples_NSP5_resistance_BN_import.csv"),
          na = "")

##### OLD CODE #####

# Combine for RAVN --------------------------------------------------------

To_RAVN <- bind_rows(RAVN_nsp5, RAVN_RdRP)

# Koble mot BN for å hente ut Pango, andre ting...

#### Problems with script:
# Ahus: Only works for 2023 - loop is hardcoded. Need to change
# Dash, underscore, etc. in NSC
# Can't find the NSP_mut column in the BN SQL-parse
# A better approach would be to query BN first and figure out which ID's have not yet been imported into BN.
# and then just look for those results.
# If no NSP5 mutations at all - good to note and import into BN. At the moment the sequence is dropped entirely in the filter in the loop
####

# Load the gene annotation from NCBI
# gff <- read_tsv("/mnt/N/Virologi/JonBrate/Drug resistance SARSCOV2/data/GCF_009858895.2_ASM985889v3_genomic.gff.gz", comment = "#", col_names = FALSE)
# nsp5 is from nucleotide 10055 to 10972
# aa 3278 in nextclade is the same as aa 15 in stanford. That's a difference of 3263 aa.
# This means that the first aa in nsp5 is 3264 in nextclade
# The last aa of nsp5 is 306 in Stanford. That's aa 3569 in nextclade (I think that's right).

# How to process many oppsett?

# Nanopore
mutations_raw_1 <- read_tsv("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2021/Oppsett102A_summaries/Oppsett102A_NextcladeAndPangolin.csv") %>% 
  # Create a new column that matches SEQUENCEID from the seqName column
  separate(seqName, into = c("Sample", NA, NA), sep = "/", remove = FALSE)


mutations_raw_2 <- read_tsv("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2021/20210518-S3-FHI81/results/nextclade_for_FHI.tsv") %>% 
  # Create a new column that matches SEQUENCEID from the seqName column
  separate(name, into = c("Sample", NA, NA), sep = "_", remove = FALSE)

# FHI
mutations_raw_3 <- read_tsv("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2023/FHI568-S3-EXT-20230106-01/results/nextclade_and_noise_for_FHI.tsv") %>% 
  # Create a new column that matches SEQUENCEID from the seqName column
  separate(name, into = c("Sample", NA, NA), sep = "-", remove = FALSE)

# Nanopore
mutations_raw_3 <- read_tsv("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2023/FHI568-S3-EXT-20230106-01/results/nextclade_and_noise_for_FHI.tsv") %>% 
  # Create a new column that matches SEQUENCEID from the seqName column
  separate(name, into = c("Sample", NA, NA), sep = "-", remove = FALSE)

mutations_raw <- bind_rows(mutations_raw_1, mutations_raw_2)
mutations_raw <- mutations_raw_3

# MIK
MIK_1 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230111-S4-MIK414-230901/results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_2 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230118-S2-MIK415-30230116//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_3 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230125-S4-MIK416-20230123///results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_4 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230201-S2-MIK417-230130//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_5 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230201-S4-MIK413-20230102//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_6 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230208-S4-MIK418-230206//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_7 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230215-S2-MIK419-230213//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_8 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230306-S2-MIK420-230227//results",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_9 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230308-S3-MIK421-230306//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_10 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230315-S1-MIK422-20230313//results/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = T) 

MIK_11 <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/20230329-S1-MIK423-20230327//results/",
                     pattern = "^nextclade_and_noise_for_FHI.tsv",
                     full.names = T) 



MIK <- c(MIK_1, MIK_2, MIK_3, MIK_4, MIK_5, MIK_6, MIK_7, MIK_8, MIK_9, MIK_10, MIK_11)


# FHI and Nano ------------------------------------------------------------



# Nanopore 2023
# Loop through each file and do the whole procedure for each:
# Make a single BN import file in the end
files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2023/",
                    pattern = "summaries_and_Pangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)



# FHI 2023

# Loop through each file and do the whole procedure for each:
# Make a single BN import file in the end
files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2023",
           pattern = "^nextclade_and_noise_for_FHI.tsv",
           full.names = TRUE,
           recursive = TRUE)

# MIK 2023

# Loop through each file and do the whole procedure for each:
# Make a single BN import file in the end
files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK/",
                    pattern = "^nextclade_and_noise_for_FHI.tsv",
                    full.names = TRUE,
                    recursive = TRUE)

# SIHF and Ahus -> See below


# Load the resistance associated mutations
stanf_3clpro <- read_csv("N:/Virologi/JonBrate/Drug resistance SARSCOV2/data/2024.01.12_3CLpro_inhibitors.csv") %>% 
  rename("stanford_res" = Mutation) %>%
  # Rename the paxlovid fold
  rename("NTV_fold" = `NTV: fold`) %>% 
  # Change del to "-"
  mutate(stanford_res = str_replace(stanford_res, "del", "-")) %>% 
  mutate(stanford_position  = as.numeric(str_sub(stanford_res, 2, -2)),
         stanford_reference = str_sub(stanford_res, 1, 1),
         stanford_mutated   = str_sub(stanford_res, -1, -1)) %>% 
  mutate(nextclade_position = stanford_position + 3263)


# Create empty final data frame
final <- tribble(
  ~name, ~pax_high_res, ~pax_mid_res, ~pax_low_res, ~NSP5
  )

# Loop through the files
pb <- txtProgressBar(min = 0, max = length(files), initial = 0) 
for (i in 1:length(files)) {
  setTxtProgressBar(pb, i)
  try(rm(mutations_raw))
  try(rm(NSP5_all_muts))
  try(rm(NSP5_high_res))
  try(rm(NSP5_mid_res))
  try(rm(NSP5_low_res))
  
  mutations_raw <- read_tsv(files[i])# %>% 
    # Create a new column that matches SEQUENCEID from the seqName column
    #separate(name, into = c("Sample", NA, NA), sep = "_", remove = FALSE)
    #mutate("Sample" = str_remove(name, "_ivar_masked"))
  
  # nsp5 is on ORF1a
  nextclade_nsp5_mutations <- mutations_raw %>% 
    # Get only ORF1a
    select(name, ORF1a) %>%
    # Separate all the mutations. Gives a list for each sample
    mutate("tmp" = str_split(ORF1a, ";")) %>% 
    # Unnest the list column (which is basically another pivot_longer)
    # Get each mutation on a separate line in a column named tmp
    unnest(tmp) %>%
    # Remove tmp
    select(-ORF1a) %>% 
    rename("ORF1a" = tmp) %>% 
    # Tease out the position, reference and mutated amino acid
    mutate(nextclade_position  = as.numeric(str_sub(ORF1a, 2, -2)),
           nextclade_reference = str_sub(ORF1a, 1, 1),
           nextclade_mutated   = str_sub(ORF1a, -1, -1)) %>% 
    # Filter out only the nsp5 gene (aa 3264 to 3569)
    filter(nextclade_position >= 3264, nextclade_position <= 3569) %>% 
    # convert to stanford nomenclature
    mutate(stanford_position = nextclade_position - 3263) %>% 
    # Create column for nsp5 mutations
    unite("NSP5", c("nextclade_reference", "stanford_position", "nextclade_mutated"), sep = "")
  
  df <- left_join(nextclade_nsp5_mutations, stanf_3clpro, by = c("NSP5" = "stanford_res"), keep = TRUE) %>% 
    # Renaming the resistance column
    rename("NSP5_resistance" = "stanford_res") %>% 
    # Add info on strength of resistance
    mutate("type" = case_when(
      NTV_fold >= 10                 ~ "pax_high_res",
      NTV_fold >= 5 | NTV_fold < 10  ~ "pax_mid_res",
      NTV_fold >= 2.5 | NTV_fold < 5 ~ "pax_low_res",
      NTV_fold < 2.5                 ~ "pax_very_low_res"
    ))
  
  # Create empty data frame
  to_BN_NSP5 <- df %>% 
    select(name) %>% distinct() #%>%  
    #add_column("pax_high_res" = NA,
    #           "pax_mid_res" = NA,
    #           "pax_low_res" = NA,
    #           "NSP5" = NA)
  
  # First need to get all the NSP5 mutations in the same cell.
  try(
    NSP5_all_muts <- df %>% 
      select(name, NSP5) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("NSP5", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_all_muts, by = "name")
  )
  
  # Then only NSP5 resistance mutations for the different categories
  try(
    NSP5_high_res <- df %>% 
      filter(type == "pax_high_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_high_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_high_res, by = "name")
  )
  
  try(
    NSP5_mid_res <- df %>% 
      filter(type == "pax_mid_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_mid_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_mid_res, by = "name")
  )
  
  try(
    NSP5_low_res <- df %>% 
      filter(type == "pax_low_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_low_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_low_res, by = "name")
  )
  
  # Add to final df
  
  final <- bind_rows(final, to_BN_NSP5)

}

# Fix for BN import - FHI
final <- final %>% 
  rename("SEQUENCEID_SWIFT" = "Sample") %>% 
  # Remove duplicate rows
  distinct() %>% 
  # Keep relevant
  filter(str_detect(SEQUENCEID_SWIFT, "SC2"))

# Fix for BN import - Nano
final <- final %>% 
  rename("SEQUENCEID_NANO29" = "name") %>% 
  # Remove duplicate rows
  distinct() %>% 
  # Keep relevant
  filter(str_detect(SEQUENCEID_NANO29, "SC2"))

# Test
write_csv(final, 
          file = "C:/Users/jonr/OneDrive - Folkehelseinstituttet/Desktop/2024.01.15_BN_resistance_import_Nano.csv",
          na = "")


# SIHF --------------------------------------------------------------------

# Loop through each file and do the whole procedure for each:
# Make a single BN import file in the end
files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/SIHF/",
                    pattern = "AndPangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)

# Create empty final data frame
final <- tribble(
  ~name, ~pax_high_res, ~pax_mid_res, ~pax_low_res, ~NSP5
)

# Loop through the files
for (i in 1:length(files)) {
  
  try(rm(mutations_raw))
  try(rm(NSP5_all_muts))
  try(rm(NSP5_high_res))
  try(rm(NSP5_mid_res))
  try(rm(NSP5_low_res))
  
  mutations_raw <- read_tsv(files[i])
  
  # nsp5 is on ORF1a
  nextclade_nsp5_mutations <- mutations_raw %>% 
    # Get only ORF1a
    select(name, ORF1a) %>% 
    # Separate all the mutations. Gives a list for each sample
    mutate("tmp" = str_split(ORF1a, ";")) %>% 
    # Unnest the list column (which is basically another pivot_longer)
    # Get each mutation on a separate line in a column named tmp
    unnest(tmp) %>%
    # Remove tmp
    select(-ORF1a) %>% 
    rename("ORF1a" = tmp) %>% 
    # Tease out the position, reference and mutated amino acid
    mutate(nextclade_position  = as.numeric(str_sub(ORF1a, 2, -2)),
           nextclade_reference = str_sub(ORF1a, 1, 1),
           nextclade_mutated   = str_sub(ORF1a, -1, -1)) %>% 
    # Filter out only the nsp5 gene (aa 3264 to 3569)
    filter(nextclade_position >= 3264, nextclade_position <= 3569) %>% 
    # convert to stanford nomenclature
    mutate(stanford_position = nextclade_position - 3263) %>% 
    # Create column for nsp5 mutations
    unite("NSP5", c("nextclade_reference", "stanford_position", "nextclade_mutated"), sep = "")
  
  df <- left_join(nextclade_nsp5_mutations, stanf_3clpro, by = c("NSP5" = "stanford_res"), keep = TRUE) %>% 
    # Renaming the resistance column
    rename("NSP5_resistance" = "stanford_res") %>% 
    # Add info on strength of resistance
    mutate("type" = case_when(
      NTV_fold >= 10                 ~ "pax_high_res",
      NTV_fold >= 5 | NTV_fold < 10  ~ "pax_mid_res",
      NTV_fold >= 2.5 | NTV_fold < 5 ~ "pax_low_res",
      NTV_fold < 2.5                 ~ "pax_very_low_res"
    ))
  
  # Create empty data frame
  to_BN_NSP5 <- df %>% 
    select(name) %>% distinct() #%>%  
  #add_column("pax_high_res" = NA,
  #           "pax_mid_res" = NA,
  #           "pax_low_res" = NA,
  #           "NSP5" = NA)
  
  # First need to get all the NSP5 mutations in the same cell.
  try(
    NSP5_all_muts <- df %>% 
      select(name, NSP5) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("NSP5", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_all_muts, by = "name")
  )
  
  # Then only NSP5 resistance mutations for the different categories
  try(
    NSP5_high_res <- df %>% 
      filter(type == "pax_high_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_high_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_high_res, by = "name")
  )
  
  try(
    NSP5_mid_res <- df %>% 
      filter(type == "pax_mid_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_mid_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_mid_res, by = "name")
  )
  
  try(
    NSP5_low_res <- df %>% 
      filter(type == "pax_low_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_low_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_low_res, by = "name")
  )
  
  # Add to final df
  
  final <- bind_rows(final, to_BN_NSP5)
  
}

# Fix for BN import - SIHF
final <- final %>% 
  rename("KEY" = "name") %>% 
  # Remove duplicate rows
  distinct()

# Test
write_csv(final, file = "C:/Users/jonr/OneDrive - Folkehelseinstituttet/Desktop/2024.01.15_BN_resistance_import_SIHF.csv", 
          na = "")


# Ahus --------------------------------------------------------------------

# Loop through each file and do the whole procedure for each:
# Make a single BN import file in the end
files <- list.files("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Eksterne/Ahus/",
                    pattern = "AndPangolin.csv$",
                    full.names = TRUE,
                    recursive = TRUE)

# Create empty final data frame
final <- tribble(
  ~name, ~pax_high_res, ~pax_mid_res, ~pax_low_res, ~NSP5
)

# Loop through the files (NB! Notice the hardcoded loop length)
pb <- txtProgressBar(min = 0, max = length(files), initial = 0) 
for (i in 89:123) {
  setTxtProgressBar(pb, i)
  
  
  try(rm(mutations_raw))
  try(rm(NSP5_all_muts))
  try(rm(NSP5_high_res))
  try(rm(NSP5_mid_res))
  try(rm(NSP5_low_res))
  
  mutations_raw <- read_tsv(files[i])# %>% 
    #mutate("Sample" = name)
  
  # nsp5 is on ORF1a
  nextclade_nsp5_mutations <- mutations_raw %>% 
    # Get only ORF1a
    select(name, ORF1a) %>%
    # Separate all the mutations. Gives a list for each sample
    mutate("tmp" = str_split(ORF1a, ";")) %>% 
    # Unnest the list column (which is basically another pivot_longer)
    # Get each mutation on a separate line in a column named tmp
    unnest(tmp) %>%
    # Remove tmp
    select(-ORF1a) %>% 
    rename("ORF1a" = tmp) %>% 
    # Tease out the position, reference and mutated amino acid
    mutate(nextclade_position  = as.numeric(str_sub(ORF1a, 2, -2)),
           nextclade_reference = str_sub(ORF1a, 1, 1),
           nextclade_mutated   = str_sub(ORF1a, -1, -1)) %>% 
    # Filter out only the nsp5 gene (aa 3264 to 3569)
    filter(nextclade_position >= 3264, nextclade_position <= 3569) %>% 
    # convert to stanford nomenclature
    mutate(stanford_position = nextclade_position - 3263) %>% 
    # Create column for nsp5 mutations
    unite("NSP5", c("nextclade_reference", "stanford_position", "nextclade_mutated"), sep = "")
  
  df <- left_join(nextclade_nsp5_mutations, stanf_3clpro, by = c("NSP5" = "stanford_res"), keep = TRUE) %>% 
    # Renaming the resistance column
    rename("NSP5_resistance" = "stanford_res") %>% 
    # Add info on strength of resistance
    mutate("type" = case_when(
      NTV_fold >= 10                 ~ "pax_high_res",
      NTV_fold >= 5 | NTV_fold < 10  ~ "pax_mid_res",
      NTV_fold >= 2.5 | NTV_fold < 5 ~ "pax_low_res",
      NTV_fold < 2.5                 ~ "pax_very_low_res"
    ))
  
  # Create empty data frame
  to_BN_NSP5 <- df %>% 
    select(name) %>% distinct() #%>%  
  #add_column("pax_high_res" = NA,
  #           "pax_mid_res" = NA,
  #           "pax_low_res" = NA,
  #           "NSP5" = NA)
  
  # First need to get all the NSP5 mutations in the same cell.
  try(
    NSP5_all_muts <- df %>% 
      select(name, NSP5) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("NSP5", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_all_muts, by = "name")
  )
  
  # Then only NSP5 resistance mutations for the different categories
  try(
    NSP5_high_res <- df %>% 
      filter(type == "pax_high_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_high_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_high_res, by = "name")
  )
  
  try(
    NSP5_mid_res <- df %>% 
      filter(type == "pax_mid_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_mid_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_mid_res, by = "name")
  )
  
  try(
    NSP5_low_res <- df %>% 
      filter(type == "pax_low_res") %>% 
      select(name, NSP5_resistance) %>% 
      add_column(dummy = "x") %>% 
      pivot_wider(names_from = dummy, values_from = NSP5_resistance) %>% 
      unnest_wider(x, names_sep = "_") %>% 
      unite("pax_low_res", starts_with("x_"), sep = ";", na.rm = TRUE)
  )
  try(
    to_BN_NSP5 <- left_join(to_BN_NSP5, NSP5_low_res, by = "name")
  )
  
  # Add to final df
  
  final <- bind_rows(final, to_BN_NSP5)
  
}

# Fix for BN import - Ahus
final <- final %>% 
  separate(name, into = c("KEY", NA, NA), sep = "\\|") %>% 
  # Remove duplicate rows
  distinct()

# Test
write_csv(final, file = "C:/Users/jonr/OneDrive - Folkehelseinstituttet/Desktop/2024.01.15_BN_resistance_import_Ahus.csv", 
          na = "")





