test_that("protein coding counted", {
  expected_output = 7
  counts_df = data.frame(raw_counts=c(1,1,3,4,1,1))
  rownames(counts_df) = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
  expect_equal(proteinCodingCount(counts_df, c('gene3', 'gene4')), expected_output)
})

test_that("create qc db", {
  library(RSQLite)
  library(glue)

  output_dir = tempdir()
  createQCdatabase(output_dir)

  db_path = file.path(output_dir, paste('qc_table', format(Sys.time(), "%Y%m%d"), sep="_"))

  expect_true(file.exists(db_path))

  unlink(output_dir, recursive = TRUE, force = TRUE)
})
