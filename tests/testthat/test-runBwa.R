context('Testing runBwa function')

# Creating bwa index:
fasta  <- system.file(package="crisprBwa",
                      "example/chr12.fa")
outdir <- tempdir()
index <- file.path(outdir, "chr12")
Rbwa::bwa_build_index(fasta,
                      index_prefix=index)


test_that('Testing runBwa', {
    update <- FALSE
    seqs <- c("GGAATACGGGAT", "GTGGGGGTTAGACC", "GTGTGGGACGCG") 
    long_seq <- "AGGGAGGAGAGAGAGAGAGGGAGGTGGGTTAGGATTTAGGATTAGTTA"
    results_bwa_0 <- runBwa(seqs, bwa_index=index, n_mismatches=0)
    results_bwa_1 <- runBwa(seqs, bwa_index=index, n_mismatches=1)
    results_bwa_2 <- runBwa(seqs, bwa_index=index, n_mismatches=2)
    results_bwa_3 <- runBwa(seqs, bwa_index=index, n_mismatches=3)
    results_bwa_null <- runBwa(long_seq, bwa_index=index, n_mismatches=0)

    expect_equal_to_reference(results_bwa_0,
                              file=file.path("objects/results_bwa_0.rds"),
                              update = update)
    expect_equal_to_reference(results_bwa_1,
                              file=file.path("objects/results_bwa_1.rds"),
                              update = update)
    expect_equal_to_reference(results_bwa_2,
                              file=file.path("objects/results_bwa_2.rds"),
                              update = update)
    expect_equal_to_reference(results_bwa_3,
                              file=file.path("objects/results_bwa_3.rds"),
                              update = update)
    expect_true(nrow(results_bwa_null)==0)
})
