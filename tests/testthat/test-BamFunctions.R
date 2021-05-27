test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("determine strandedness", {
  expect_equal(determineStrandedness('E7420L'), 'reverse')
  expect_equal(determineStrandedness('SolexaPrep'), 'unstranded')
})


test_that("create bam path", {
  run_num = '2651'
  fastq_filename = "Brent_17_GTAC_17_ACTGTCGATC_S18_R1_001.fastq.gz"
  lts_align_prefix = '/mnt/htcf_lts/lts_align_expr'
  function_output = createBamPath(run_num, fastq_filename, lts_align_prefix, test=TRUE)
  expected_output = "/mnt/htcf_lts/lts_align_expr/run_2651_samples/align/Brent_17_GTAC_17_ACTGTCGATC_S18_R1_001_sorted_aligned_reads_with_annote.bam"
  expect_equal(function_output, expected_output)

  expected_error_output = paste0("The following path is invalid: ", expected_output)

  expect_message(createBamPath(run_num, fastq_filename, lts_align_prefix),expected_error_output)
})
