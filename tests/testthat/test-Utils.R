test_that("isNumeric TRUE: single number", {
  expect_true(isNumeric(2), TRUE)
})

test_that("isNumeric FALSE: single string", {
  expect_false(isNumeric("abc"))
  expect_true
})

test_that("isNumeric TRUE: numeric dataframe", {
  df = data.frame("one" = c(1,2,3), "two" = c(4, 5, 6))
  expect_true(isNumeric(df))
})

test_that("isNumeric FALSE: mixed dataframe", {
  df = data.frame("one" = c(1,2,3), "two" = c("a", "b", "c"))
  expect_false(isNumeric(df))
})

test_that("readInData works", {
  output_dir = tempdir()
  df = tibble(a=c(1,2,3), b=c(1,2,3), c=c(1,2,3))
  write_csv(df, file.path(output_dir, "test.csv"))
  write_tsv(df, file.path(output_dir, "test.tsv"))

  expect_equal(df, readInData(file.path(output_dir, "test.csv")), ignore_attr=TRUE)
  expect_equal(df, readInData(file.path(output_dir, "test.tsv")), ignore_attr=TRUE)

  unlink(output_dir, recursive=TRUE, force=TRUE)

  })
