#' @title Run BWA short-read aligner
#' @description Return BWA alignments for a list of short sequences
#'     for a prebuilt BWA index. 
#' 
#' @param sequences Character vector of DNA sequences.
#' @param bwa_index String specifying path to the BWA index.
#' @param n_mismatches Integer specifying maximum number
#'     of mismatches allowed between the query sequences and the
#'     index sequences.
#' 
#' @return A data.frame of the alignments with the following columns:
#'    \itemize{
#'        \item \code{query} — string specifying query DNA sequence
#'        \item \code{chr} - string specifying chromosome name
#'        \item \code{pos} - string specifying genomic coordinate of the start
#'              of the target DNA sequence
#'        \item \code{strand} - string specifying strand ("+" or "-") 
#'        \item \code{n_mismatches} - integer specifying number of mismatches
#'              between query and target sequences
#'    }
#' 
#' @details \code{runBwa} can be used to map short DNA sequences 
#'     to a reference genome. To search for sequences while imposing
#'     constraints on PAM sequences (such as gRNA spacer sequences), see
#'     \code{runCrisprBwa} instead.  
#' 
#' @seealso \code{link{runCrisprBwa}} to map gRNA spacer sequences.
#' 
#' @examples
#' 
#' fasta <- system.file(package="crisprBwa", "example/chr12.fa")
#' outdir <- tempdir()
#' index <- file.path(outdir, "chr12")
#' Rbwa::bwa_build_index(fasta,
#'                       index_prefix=index)
#' 
#' seqs <- c("GGAAGTTG",
#'           "GTGGACAC",
#'           "GTGTGCAA")
#' 
#' aln <- runBwa(seqs,
#'               n_mismatches=1,
#'               bwa_index=index)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom Rbwa bwa_aln bwa_sam xa2multi
#' @export
runBwa <- function(sequences,
                   bwa_index=NULL,
                   n_mismatches=3
){
    if (is.null(bwa_index)){
        stop("bwa_index must be specified.")
    }
    .checkBwaIndex(bwa_index)
    fastq <- .fastqfy(sequences)
    output.sai <- paste0(fastq, ".sai")
    output.sam <- paste0(fastq, ".sam")
    output.sam.multi <- paste0(fastq, "multi.sam")

    bwa_aln(type='single',
            index_prefix=bwa_index,
            fastq_files=fastq,
            sai_files=output.sai,
            N=TRUE, #complete search
            n=n_mismatches, #edit distance
            #k=n_mismatches, #edit distance seed
            l=nchar(sequences)[1]*3, #length of the seed
            o=0) #max gaps open
    bwa_sam(type="single",
            fastq_files=fastq,
            sai_files=output.sai,
            sam_file=output.sam,
            index_prefix=bwa_index,
            n=1000000000L)
    xa2multi(input_sam_file=output.sam,
             output_sam_file=output.sam.multi)
    aln <- .parseBwaResults(output.sam.multi)
    return(aln)
}



# #' @importFrom stringr str_extract
# #' @importFrom utils read.csv
# .parseBwaResults_old <- function(file){
#     results <- readLines(file)
#     start   <- which(grepl("@PG",results))
#     results <- read.csv(file,
#                         sep="\t",
#                         skip=start,
#                         header=FALSE)
#     aln <- results[,c(1,3,4),drop=FALSE]
#     colnames(aln) <- c("query","chr", "pos")
#     aln$strand <- ifelse(results[,2] %in% c("272", "16"), "-", "+")
#     if (ncol(results)>11){
#         tags <- results[,12:ncol(results)]
#         tags <- apply(tags,1,paste, collapse=";")
#         n_mismatches <- str_extract(tags, "NM\\:i\\:[0-9]+")
#         n_mismatches <- gsub("NM\\:i\\:", "", n_mismatches)
#         aln$n_mismatches <- as.integer(n_mismatches)
#     } else {
#         aln <- .bwaEmptyAlignments()
#     }
#     return(aln)
# }

#' @importFrom stringr str_extract
#' @importFrom readr read_lines read_delim
.parseBwaResults <- function(file){
    results <- read_lines(file, n_max=100000L)
    start   <- which(grepl("@PG",results))
    results <- suppressWarnings(read_delim(file,
                          delim="\t",
                          skip=start,
                          show_col_types = FALSE,
                          col_names=NULL))
    results <- as.data.frame(results)
    aln <- results[,c(1,3,4),drop=FALSE]
    colnames(aln) <- c("query","chr", "pos")
    aln$strand <- ifelse(results[,2] %in% c("272", "16"), "-", "+")
    if (ncol(results)>11){
        tags <- results[,12:ncol(results)]
        tags <- apply(tags,1,paste, collapse=";")
        n_mismatches <- str_extract(tags, "NM\\:i\\:[0-9]+")
        n_mismatches <- gsub("NM\\:i\\:", "", n_mismatches)
        aln$n_mismatches <- as.integer(n_mismatches)
    } else {
        aln <- .bwaEmptyAlignments()
    }
    return(aln)
}




.bwaEmptyAlignments <- function(n_mismatches){
    cols <- c("query",
              "chr",
              "pos",
              "strand",
              "n_mismatches")
    out <- data.frame(matrix(0, ncol=length(cols)))
    colnames(out) <- cols
    out <- out[-1,,drop=FALSE]
    return(out)
}





