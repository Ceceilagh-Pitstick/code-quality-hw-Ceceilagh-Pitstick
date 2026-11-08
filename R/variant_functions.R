# Refactor the calculations in original/messy_analysis.R into these interfaces.
# Preserve the required behavior described in README.md.
# This file should define functions only: no reads, writes, or analysis run.

calculate_maf <- function(genotypes, na_rm = TRUE) {
    # TODO: calculate MAF for one numeric genotype vector.
    #p = sum of genotype values / (2 * number of observed individuals),
    #na_rm must be exactly one nonmissing logical value; reject invalid settings.
    if (!is.logical(na_rm) || length(na_rm) != 1 || is.na(na_rm)) {
        stop("na_rm must be exactly one nonmissing logical value (TRUE or FALSE).")
    }
    if (!is.numeric(genotypes)) {
        stop("Genotype input must be a numeric vector.")
    }
    #Reject nonnumeric input and any observed value other than 0, 1, or 2.
    values <- c(0, 1, 2)
    if (any(!genotypes %in% values & !is.na(genotypes))) {
        stop("Invalid genotype values detected. Genotypes can be 0, 1, 2, or NA.")
    }
    #Empty or entirely missing input: return numeric NA, not 0, NaN, or Inf.
    if (length(genotypes) == 0 || all(is.na(genotypes))) {
        return(NA_real_)
    }
    #na_rm = FALSE: if any genotype is missing, return NA.
    if (!na_rm && any(is.na(genotypes))) {
        return(NA_real_)
    }
    #na_rm = TRUE (default): exclude missing values from numerator AND denominator.
    if (na_rm) {
        genotypes <- genotypes[!is.na(genotypes)]
        if (length(genotypes) == 0) {
            return(NA_real_)
        }
    }
    
    #count allele frequency
    number_genotypes <- length(genotypes)
    if (number_genotypes == 0){return(NA)}
    allele_count <- 2 * number_genotypes
    #count alt allele frequency
    alt_count<- sum(genotypes)
    #calculate alt frequency
    alt_frequency<- alt_count/allele_count
    #MAF CALCULATION
    MAF = pmin(alt_frequency, 1 - alt_frequency)
    return(MAF)
}

summarize_variant <- function(genotypes, min_maf = 0.25,
                              max_missing = 0.25) {
    # TODO: return the specified one-row summary, calling calculate_maf().
    
    #min_maf must be one finite numeric value between 0 and 0.5 inclusive.
    if (!is.numeric(min_maf) || length(min_maf) != 1 || !is.finite(min_maf) || min_maf < 0 || min_maf > 0.5) {
        stop("min_maf must be a single finite numeric value between 0 and 0.5 inclusive.")
    }
    
    #max_missing must be one finite numeric value between 0 and 1 inclusive
    if (!is.numeric(max_missing) || length(max_missing) != 1 || !is.finite(max_missing) || max_missing < 0 || max_missing > 1) {
        stop("max_missing must be a single finite numeric value between 0 and 1 inclusive.")
    }
    
    total_length <- length(genotypes)
    
    #Reject invalid inputs/settings with informative errors.
    if (total_length == 0) {
        maf_val <- calculate_maf(genotypes, na_rm = TRUE) # Returns NA_real_
        return(data.frame(
            n_called = 0L,
            missing_rate = NA_real_,
            maf = maf_val,
            keep = FALSE
        ))}
        # 4. Calculate missingness statistics
        n_miss <- sum(is.na(genotypes))
        n_called <- total_length - n_miss
        missing_rate <- n_miss / total_length
        
        # Use calculate_maf()
        
        maf_val <- calculate_maf(genotypes, na_rm = TRUE)
        
        keep_val <- FALSE
        if (n_called > 0 && !is.na(maf_val)) {
            if (maf_val >= min_maf && missing_rate <= max_missing) {
                keep_val <- TRUE
            }
        }
        
        #make data table
        return(data.frame(
            n_called = as.integer(n_called),
            missing_rate = as.numeric(missing_rate),
            maf = as.numeric(maf_val),
            keep = as.logical(keep_val)
        ))
        
    }
