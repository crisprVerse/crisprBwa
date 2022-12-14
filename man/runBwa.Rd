% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runBwa.R
\name{runBwa}
\alias{runBwa}
\title{Run BWA short-read aligner}
\usage{
runBwa(sequences, bwa_index = NULL, n_mismatches = 3)
}
\arguments{
\item{sequences}{Character vector of DNA sequences.}

\item{bwa_index}{String specifying path to the BWA index.}

\item{n_mismatches}{Integer specifying maximum number
of mismatches allowed between the query sequences and the
index sequences.}
}
\value{
A data.frame of the alignments with the following columns:
   \itemize{
       \item \code{query} — string specifying query DNA sequence
       \item \code{chr} - string specifying chromosome name
       \item \code{pos} - string specifying genomic coordinate of the start
             of the target DNA sequence
       \item \code{strand} - string specifying strand ("+" or "-") 
       \item \code{n_mismatches} - integer specifying number of mismatches
             between query and target sequences
   }
}
\description{
Return BWA alignments for a list of short sequences
    for a prebuilt BWA index.
}
\details{
\code{runBwa} can be used to map short DNA sequences 
    to a reference genome. To search for sequences while imposing
    constraints on PAM sequences (such as gRNA spacer sequences), see
    \code{runCrisprBwa} instead.
}
\examples{

fasta <- system.file(package="crisprBwa", "example/chr12.fa")
outdir <- tempdir()
index <- file.path(outdir, "chr12")
Rbwa::bwa_build_index(fasta,
                      index_prefix=index)

seqs <- c("GGAAGTTG",
          "GTGGACAC",
          "GTGTGCAA")

aln <- runBwa(seqs,
              n_mismatches=1,
              bwa_index=index)

}
\seealso{
\code{link{runCrisprBwa}} to map gRNA spacer sequences.
}
\author{
Jean-Philippe Fortin
}
