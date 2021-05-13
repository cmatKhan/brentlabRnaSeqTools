test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

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
