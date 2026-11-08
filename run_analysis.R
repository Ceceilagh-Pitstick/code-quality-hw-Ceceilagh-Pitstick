# Run from the repository root.
# TODO: move and simplify the DRIVER portion of original/messy_analysis.R here.
# Read data/genotypes.csv with colClasses = "numeric", call your functions for each column, and write ONE
# combined results/summary.csv after all variant summaries are calculated.
# Put the thresholds in named settings and pass them as function arguments.
# Source R/variant_functions.R to load your function definitions.

source("R/variant_functions.R")
genotypes_df <- read.csv("data/genotypes.csv", colClasses = "numeric")

#get col names
variant_ids <- names(genotypes_df)

#sum
summaries_list <- vector("list", length(variant_ids))

for (i in seq_along(variant_ids)) {
    v_id <- variant_ids[i]
    col_data <- genotypes_df[[v_id]]
    summaries_list[[i]] <- summarize_variant(col_data, min_maf = 0.25, max_missing = 0.25)
}

# Combine all row summaries into a single data frame
results_df <- rbind(summaries_list)

# add the variant_id column
results_df <- data.frame(variant_id = variant_ids, results_df, stringsAsFactors = FALSE)

write.csv(results_df, "results/summary.csv", row.names = FALSE)


