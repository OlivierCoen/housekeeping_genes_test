if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.biotools.fr/")

if (!require("ArrayExpress", quietly = TRUE))
    BiocManager::install("ArrayExpress")

if (!require("affy", quietly = TRUE))
    BiocManager::install("affy")

library("ArrayExpress")
library("affy")

#ae_data <- ArrayExpress("E-MTAB-10483")
exp = getAE("E-MEXP-21", type = "full")

# Initialize an empty list to store the individual dataframes
df_list <- list()

for (file in exp$processedFiles) {
    df <- read.table(file, header = TRUE, sep = "\t")
    df_list[[file]] <- df
}

# Combine all data frames in the list into one
combined_df <- do.call(rbind, df_list)

write.table(combined_df, file = 'export.csv', sep = "\t")

