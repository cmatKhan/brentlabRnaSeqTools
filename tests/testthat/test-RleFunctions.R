test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

mock_raw_counts = data.frame("sample1" = c(10,0,2,5),
                             "sample2" = c(5, 9, 0, 0),
                             "sample3" = c(2,3,7,10))

# sample1  sample2  sample3
# 1 3.459432 2.584963 1.584963
# 2 0.000000 3.321928 2.000000
# 3 1.584963 0.000000 3.000000
# 4 2.584963 0.000000 3.459432
log2_counts_df = log2(mock_raw_counts + 1)

test_that("calculate gene wise medians", {

  expect_equal(calculateGeneWiseMedians(log2_counts_df), apply(log2_counts_df,1,median))

})

test_that("calculate RLE", {

  # > gene_wise_medians
  # [1] 2.584963 2.000000 1.584963 2.584963

  # expected output
  #  sample1   sample2    sample3
  #  0.8744691  0.000000 -1.0000000
  # -2.0000000  1.321928  0.0000000
  #  0.0000000 -1.584963  1.4150375
  #  0.0000000 -2.584963  0.8744691

  expected_rle_output = sweep(log2_counts_df, 1, calculateGeneWiseMedians(log2_counts_df))

  expect_equal(calculateRLE(mock_raw_counts), expected_rle_output)

})


