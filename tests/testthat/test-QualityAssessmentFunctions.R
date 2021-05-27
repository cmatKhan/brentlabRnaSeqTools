test_that("protein coding counted", {
  expected_output = 7
  counts_df = data.frame(raw_counts=c(1,1,3,4,1,1))
  rownames(counts_df) = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
  expect_equal(proteinCodingCount(counts_df, c('gene3', 'gene4')), expected_output)
})
