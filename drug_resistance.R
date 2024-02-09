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

# Create a column "type" based on different NTV fold values
joined_data$type <- ifelse(joined_data$NTV_fold >= 10, "pax_high_res",
                                ifelse(joined_data$NTV_fold >= 5, "pax_mid_res",
                                       ifelse(joined_data$NTV_fold >= 2.5, "pax_low_res",
                                              "pax_very_low_res")))

# Collect all NSP5 mutations in a column named "NSP5"
aggregated_data <- aggregate(NSP5 ~ name, data = joined_data, FUN = function(x) paste(x, collapse = ";"))


# Create function that aggregates data based on the resistance levels in the "type" column 
# type_name specifies the resistance levels. E.g. "pax_high_res"
# The dataset is first filtered on these categories. Then any mutations present are aggregated on the same sample (name)
aggregate_nsp5_by_type <- function(data, type_name) {
  filtered_data <- subset(data, type == type_name)
  
  # Check if there are any rows after filtering
  if (nrow(filtered_data) == 0) {
    # If no rows, create a placeholder data frame with names from the original dataset and empty values for the type
    aggregated <- data.frame(name = unique(data$name), stringsAsFactors = FALSE)
    aggregated[type_name] <- ""  # Create an empty column for the type
  } else {
    # If rows exist, proceed with aggregation
    aggregated <- aggregate(NSP5 ~ name, data = filtered_data, FUN = function(x) paste(unique(x), collapse = ";"))
    colnames(aggregated)[2] <- type_name
  }
  
  return(aggregated)
}

# Re-run the aggregation for each resistance category (level) and merge back into aggregated_data as before
pax_high_res_data <- aggregate_nsp5_by_type(joined_data, "pax_high_res")
pax_mid_res_data <- aggregate_nsp5_by_type(joined_data, "pax_mid_res")
pax_low_res_data <- aggregate_nsp5_by_type(joined_data, "pax_low_res")

# Merge the aggregated data back into aggregated_data on the sample names
aggregated_data <- merge(aggregated_data, pax_high_res_data, by = "name", all.x = TRUE)
aggregated_data <- merge(aggregated_data, pax_mid_res_data, by = "name", all.x = TRUE)
aggregated_data <- merge(aggregated_data, pax_low_res_data, by = "name", all.x = TRUE)

### Write final data ###
# Write the final data as a tsv file
write.table(aggregated_data, file = "stanford_3clpro_resistance.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
