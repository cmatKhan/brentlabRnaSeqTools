
# set up some test data
sample1_counts = c(10,0,2,5)
sample2_counts = c(5, 9, 0, 0)
sample3_counts = c(2,3,7,10)
mock_raw_counts = tibble(sample1 = sample1_counts,
                         sample2 = sample2_counts,
                         sample3 = sample3_counts)

# sample1  sample2  sample3
# 1 3.459432 2.584963 1.584963
# 2 0.000000 3.321928 2.000000
# 3 1.584963 0.000000 3.000000
# 4 2.584963 0.000000 3.459432

log2_counts_df = log2(mock_raw_counts + 1)
sample_vector = c("sample1", "sample2", "sample3")
mock_rle_summary = tibble(FASTQFILENAME = sample_vector,
                          SAMPLE_DEVIATION_MEDIAN = c(0,-.792,.437),
                          ABS_SAMPLE_DEVIATION_MEDIAN = abs(c(0,-.792,.437)),
                          INTERQUARTILE_RANGE = c(0.719,2.17,1.26))
# had to do this for testthat to pass -- fix this eventually
names(mock_rle_summary$SAMPLE_DEVIATION_MEDIAN) = sample_vector
names(mock_rle_summary$ABS_SAMPLE_DEVIATION_MEDIAN) = sample_vector
names(mock_rle_summary$INTERQUARTILE_RANGE) = sample_vector



test_that("calculate gene wise medians", {

  expect_equal(apply(log2_counts_df,1,median), calculateGeneWiseMedians(log2_counts_df))

})

test_that("calculateRLE works", {

  # > gene_wise_medians
  # [1] 2.584963 2.000000 1.584963 2.584963

  # expected output
  #  sample1   sample2    sample3
  #  0.8744691  0.000000 -1.0000000
  # -2.0000000  1.321928  0.0000000
  #  0.0000000 -1.584963  1.4150375
  #  0.0000000 -2.584963  0.8744691

  expected_rle_output = sweep(log2_counts_df, 1, calculateGeneWiseMedians(log2_counts_df))

  expect_equal(expected_rle_output, calculateRLE(mock_raw_counts))

})

test_that("RLE summary works", {

  test_rle_summary = rleSummary(calculateRLE(mock_raw_counts))

  expect_equal(mock_rle_summary, test_rle_summary, tolerance=.05)

})

