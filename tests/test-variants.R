# The supplied run_tests.R loads testthat and the functions from R/.
# Do not source original/messy_analysis.R here.

# A supplied ordinary example; keep it and add your own tests below.
test_that("MAF counts allele copies among observed diploid genotypes", {
    expect_equal(calculate_maf(c(0, 0, 1, 1)), 0.25)
})

# TODO: add small, independently justified known-answer tests for both functions.
# Cover the behavior groups listed in README.md. Do not just test data/genotypes.csv.
# Use informative test names; explain nonobvious expected answers in comments.

#Ceceilagh Pitstick's tests

#Ordinary MAF calculations, including a case where the COUNTED allele is common.

test_that("Ordinary MAF calculations, including a case where the COUNTED allele is common.", {
    expect_equal(calculate_maf(c(0, 0, 0, 0)), 0)
    expect_equal(calculate_maf(c(2, 2, 2, 2)), 0)
    expect_equal(calculate_maf(c(1)), 0.5)
    expect_equal(calculate_maf(c(0, 1, 2, NA), na_rm = TRUE), 0.5)
})

#Missing values: the correct denominator, both values of na_rm (TRUE IS ABOVE), and the default.

test_that("calculate_maf proceeds with na_rm = FALSE", {
    expect_true(is.na(calculate_maf(c(0, 1, NA), na_rm = FALSE)))
    expect_equal(calculate_maf(c(0, 1, 2), na_rm = FALSE), 0.5)
})

#Empty, entirely missing, monomorphic, and single-observation inputs.
test_that("calculate_maf proceeds with empty and entirely missing inputs", {
    expect_true(is.na(calculate_maf(numeric(0))))
    expect_true(is.na(calculate_maf(c(NA_real_, NA_real_))))
})

#Each of the summary's output columns, not merely whether the function returns.
test_that("summarize_variant returns expected values", {
    results_check <- summarize_variant(c(4,0,0.25), min_maf = 0.25, max_missing = 0.25)
    expect_equal(nrow(results_check), 1)
    expect_equal(ncol(results_check), 4)
    expect_named(results_check, c("n_called", "missing_rate", "maf", "keep"))
    expect_equal(results_check$n_called, 4)
    expect_equal(results_check$missing_rate, 0)
})

#Both threshold boundaries, cases on either side, and NONDEFAULT thresholds.

test_that("summarize_variant checks all thresholds correctly", {
    expect_error(summarize_variant(c(0, 1), min_maf = -0.1))
    expect_error(summarize_variant(c(0, 1), min_maf = 0.6))
    expect_error(summarize_variant(c(0, 1), max_missing = -0.1))
    expect_error(summarize_variant(c(0, 1), max_missing = 1.2))
    expect_error(summarize_variant(c(0, 1), min_maf = 0.23))
})