### Stanford resistance information ###

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

### Nextclade output ###

# Read the nextclade output to get the different mutations
nextclade <- read.table("C:/Users/jonr/Downloads/OppsettNGSSEQ20240125-01A_NextcladeAndPangolin.csv",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)

# Identify the different NSP5 mutations and convert the nomenclature to match the Stanford nomenclature.
# The nsp5 gene is in ORF1a

# Select necessary columns
data_selected <- nextclade[c("name", "ORF1a")]

# Split the "ORF1a" column into lists of the different mutations based on the ";" character
data_selected$tmp <- strsplit(as.character(data_selected$ORF1a), ";")

# Expand the lists into rows
# Adjust the lapply function to handle empty or NULL lists
data_expanded <- do.call(rbind, lapply(seq_along(data_selected$tmp), function(i) {
  if (length(data_selected$tmp[[i]]) == 0) { # If the list is empty
    return(data.frame(name = data_selected$name[i], ORF1a = NA, stringsAsFactors = FALSE))
  } else {
    return(data.frame(name = data_selected$name[i], ORF1a = unlist(data_selected$tmp[[i]]), stringsAsFactors = FALSE))
  }
}))

# Ensure that row names are not problematic by resetting them
rownames(data_expanded) <- NULL

# Tease out the position, reference, and mutated amino acid
data_expanded$nextclade_position <- as.numeric(substring(data_expanded$ORF1a, 2, nchar(data_expanded$ORF1a) - 1))
data_expanded$nextclade_reference <- substring(data_expanded$ORF1a, 1, 1)
data_expanded$nextclade_mutated <- substring(data_expanded$ORF1a, nchar(data_expanded$ORF1a), nchar(data_expanded$ORF1a))

# Filter out only the nsp5 gene (aa 3264 to 3569)
data_filtered <- subset(data_expanded, nextclade_position >= 3264 & nextclade_position <= 3569)

# Convert to stanford nomenclature and create column for nsp5 mutations
data_filtered$stanford_position <- data_filtered$nextclade_position - 3263
data_filtered$NSP5 <- with(data_filtered, paste(nextclade_reference, stanford_position, nextclade_mutated, sep = ""))

# Remove unnecessary columns
data_final <- data_filtered[c("name", "ORF1a", "nextclade_position", "NSP5")]

### Merge Stanford positions with our data ###

# Join 'data_final' with 'data' on the NSP5 and stanford_res columns. Only keeping rows in the data_final object
joined_data <- merge(data_final, data, by.x = "NSP5", by.y = "stanford_res", all.x = TRUE)