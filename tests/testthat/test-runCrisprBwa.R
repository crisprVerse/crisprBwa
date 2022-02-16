context('Testing runCrisprBwa functions')

#library(crisprBwa)
library(crisprBase)
library(Rbwa)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Creating bwa index:
fasta  <- system.file(package="crisprBwa",
                      "example/chr12.fa")
outdir <- tempdir()
index <- file.path(outdir, "chr12")
Rbwa::bwa_build_index(fasta,
                      index_prefix=index)



data(SpCas9, package="crisprBase")
data(AsCas12a, package="crisprBase")
Cas9 <- SpCas9
Cas12a <- AsCas12a





spacers_cas9 <- c("AGCTGTCCGTGGGGGTCCGC",
                  "CCCCTGCTGCTGTGCCAGGC",
                  "ACGAACTGTAAAAGGCTTGG",
                  "ACGAACTGTAACAGGCTTGG",
                  "AAGGCCCTCAGAGTAATTAC")





test_that('Testing runCrisprBwa with Cas9 canonical', {
    update <- FALSE
    results_cas9_mm0 <- runCrisprBwa(spacers_cas9,
                                     bsgenome=bsgenome,
                                     bwa_index=index,
                                     n_mismatches=0,
                                     crisprNuclease=Cas9,
                                     canonical=TRUE)
    results_cas9_mm1 <- runCrisprBwa(spacers_cas9,
                                     bsgenome=bsgenome,
                                     bwa_index=index,
                                     n_mismatches=1,
                                     crisprNuclease=Cas9,
                                     canonical=TRUE)
    results_cas9_mm2 <- runCrisprBwa(spacers_cas9,
                                     bsgenome=bsgenome,
                                     bwa_index=index,
                                     n_mismatches=2,
                                     crisprNuclease=Cas9,
                                     canonical=TRUE)
    results_cas9_mm3 <- runCrisprBwa(spacers_cas9,
                                     bsgenome=bsgenome,
                                     bwa_index=index,
                                     n_mismatches=3,
                                     crisprNuclease=Cas9,
                                     canonical=TRUE)
    expect_equal_to_reference(results_cas9_mm0,
                              file=file.path("objects/results_cas9_mm0.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm1,
                              file=file.path("objects/results_cas9_mm1.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm2,
                              file=file.path("objects/results_cas9_mm2.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm3,
                              file=file.path("objects/results_cas9_mm3.rds"),
                              update=update)
})


test_that('Te sting runCrisprBwa with Cas9 non-canonical', {
    update <- FALSE
    results_cas9_mm0_nc <- runCrisprBwa(spacers_cas9,
                                        bsgenome=bsgenome,
                                        bwa_index=index,
                                        n_mismatches=0,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm1_nc <- runCrisprBwa(spacers_cas9,
                                        bsgenome=bsgenome,
                                        bwa_index=index,
                                        n_mismatches=1,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm2_nc <- runCrisprBwa(spacers_cas9,
                                        bsgenome=bsgenome,
                                        bwa_index=index,
                                        n_mismatches=2,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm3_nc <- runCrisprBwa(spacers_cas9,
                                        bsgenome=bsgenome,
                                        bwa_index=index,
                                        n_mismatches=3,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    expect_equal_to_reference(results_cas9_mm0_nc,
                              file=file.path("objects/results_cas9_mm0_nc.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm1_nc,
                              file=file.path("objects/results_cas9_mm1_nc.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm2_nc,
                              file=file.path("objects/results_cas9_mm2_nc.rds"),
                              update=update)
    expect_equal_to_reference(results_cas9_mm3_nc,
                              file=file.path("objects/results_cas9_mm3_nc.rds"),update=update)
})


spacers_cas9_short<- c("TGTCCGTGGGGGTCCGC",
                       "CTGCTGCTGTGCCAGGC",
                       "ACGAACTGTAAAAGGCT",
                       "AACTGTAACAGGCTTGG",
                       "GCCCTCAGAGTAATTAC")




test_that('Testing runCrisprBwa with Cas9 (short spacers)', {
    update <- FALSE
    expect_error(runCrisprBwa(spacers_cas9_short,
                              bsgenome=bsgenome,
                              bwa_index=index,
                              n_mismatches=0,
                              crisprNuclease=Cas9,
                              canonical=FALSE))
    results_short_cas9 <- runCrisprBwa(spacers_cas9_short,
                                       bsgenome=bsgenome,
                                       bwa_index=index,
                                       n_mismatches=3,
                                       crisprNuclease=Cas9,
                                       canonical=FALSE,
                                       force_spacer_length=TRUE)
    expect_equal_to_reference(results_short_cas9,
                              file=file.path("objects/results_short_cas9.rds"),
                              update=update)
})



spacers_cas12a <- c("CAGTTCGTACTGGGAAGGCTTTG",
                    "GCTTTGAATGCAAAGAGCACAGG",
                    "CCAATTCCAAGAGACACAAGTAA",
                    "TCATCTCAGGTATTAATCAATGA")




test_that('Testing runCrisprBwa with Cas12a', {
    update <- FALSE
    results_cas12a_mm0 <- runCrisprBwa(spacers_cas12a,
                                       bsgenome=bsgenome,
                                       bwa_index=index,
                                       n_mismatches=0,
                                       crisprNuclease=Cas12a,
                                       canonical=TRUE)
    results_cas12a_mm1 <- runCrisprBwa(spacers_cas12a,
                                       bsgenome=bsgenome,
                                       bwa_index=index,
                                       n_mismatches=1,
                                       crisprNuclease=Cas12a,
                                       canonical=TRUE)
    results_cas12a_mm2 <- runCrisprBwa(spacers_cas12a,
                                       bsgenome=bsgenome,
                                       bwa_index=index,
                                       n_mismatches=2,
                                       crisprNuclease=Cas12a,
                                       canonical=TRUE)
    results_cas12a_mm3 <- runCrisprBwa(spacers_cas12a,
                                       bsgenome=bsgenome,
                                       bwa_index=index,
                                       n_mismatches=3,
                                       crisprNuclease=Cas12a,
                                       canonical=TRUE)
    expect_equal_to_reference(results_cas12a_mm0,
                              file=file.path("objects/results_cas12a_mm0.rds"),
                              update=update)
    expect_equal_to_reference(results_cas12a_mm1,
                              file=file.path("objects/results_cas12a_mm1.rds"),
                              update=update)
    expect_equal_to_reference(results_cas12a_mm2,
                              file=file.path("objects/results_cas12a_mm2.rds"),
                              update=update)
    expect_equal_to_reference(results_cas12a_mm3,
                              file=file.path("objects/results_cas12a_mm3.rds"),
                              update=update)
})






