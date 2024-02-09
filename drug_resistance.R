# Load the CSV file with resistance information for the 3CLpro gene from Stanford. Downloaded from here: https://covdb.stanford.edu/drms/3clpro/
data <- read.csv("N:/Virologi/JonBrate/Drug resistance SARSCOV2/data/2024.01.12_3CLpro_inhibitors.csv")

# Rename columns
colnames(data)[1] <- "stanford_res"
colnames(data)[2] <- "NTV_fold"

# Replace "del" with "-" in the "stanford_res" column
data$stanford_res <- gsub("del", "-", data$stanford_res)

# Create "stanford_position" by extracting the numbers from the mutation list
data$stanford_position <- as.numeric(gsub("[^0-9]", "", data$stanford_res))

# Create "stanford_reference" by extracting the first character of the mutation list
data$stanford_reference <- substr(data$stanford_res, 1, 1)

# Create "stanford_mutated" by extracting the last character of the mutation list
data$stanford_mutated <- substr(data$stanford_res, nchar(data$stanford_res), nchar(data$stanford_res))

# Create "nextclade_position" by adding 3263 to "stanford_position".
# 3263 is the shift in genome position to get the start of the nsp5 gene
data$nextclade_position <- data$stanford_position + 3263
