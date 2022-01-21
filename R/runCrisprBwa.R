#' @title Find gRNA spacer alignments with bwa
#' @description Return bwa alignments for a list of gRNA spacer sequences.
#' 
#' @param spacers Character vector of DNA sequences corresponding
#'     to gRNA spacer sequences. Must all be of equal length. 
#' @param bsgenome BSgenome object to be used in spacer mode. 
#' @param bwa_index Path to the bwa index to be used for alignment.
#' @param crisprNuclease \code{CrisprNuclease} object. 
#' @param canonical Should only canonical PAM sequences be considered?
#'     TRUE by default.
#' @param ignore_pam If TRUE, will return all matches regardless of
#'     PAM sequence. FALSE by default.
#' @param n_mismatches Integer between 0 and 3 specifying maximum
#'     number of mismatches allowed between spacer and protospacer sequences.
#' @param force_spacer_length Should the spacer length be overwritten in the
#'     crisprNuclease object? FALSE by default. 
#' @param verbose Should messages be printed to the consolde? TRUE by default.
#' 
#' @return \strong{runBwa} returns spacer alignment data, including genomic 
#'     coordinates and sequence.
#' 
#' @details \code{runCrisprBwa} is similar to \code{runBwa}, with the 
#'     addition of imposing constraints on PAM sequences such that the query
#'     sequences are valid protospacer sequences in the searched genome. 
#' 
#' @examples
#' 
#' # Building BWA index first:
#' fasta <- system.file(package="crisprBwa", "example/chr12.fa")
#' outdir <- tempdir()
#' index <- file.path(outdir, "chr12")
#' Rbwa::bwa_build_index(fasta,
#'                       index_prefix=index)
#' 
#' 
#' # Aligning Cas9 gRNA
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' seqs <- c("AGCTGTCCGTGGGGGTCCGC",
#'           "CCCCTGCTGCTGTGCCAGGC")
#' data(SpCas9, package="crisprBase")
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' results <- runCrisprBwa(seqs,
#'                         bsgenome=bsgenome,
#'                         bwa_index=index,
#'                         n_mismatches=2,
#'                         crisprNuclease=SpCas9)
#' 
#' @seealso \code{link{runBwa}} to map general DNA sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom crisprBase nucleaseName pams pamLength pamIndices
#' @importFrom crisprBase spacerLength spacerLength<- pamSide isRnase
#' @importFrom crisprBase hasSpacerGap
#' @importFrom GenomeInfoDb seqnames
runCrisprBwa <- function(spacers,
                         bwa_index=NULL,
                         bsgenome=NULL,
                         crisprNuclease=NULL,
                         canonical=TRUE,
                         ignore_pam=FALSE,
                         n_mismatches=0, 
                         force_spacer_length=FALSE,
                         verbose=TRUE
){ 
    #Checking inputs:
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    if (hasSpacerGap(crisprNuclease)){
        stop("CRISPR nucleases with spacer gaps are not ",
             "supported at the moment.")
    }
    if (isRnase(crisprNuclease)){
        stop("RNA-targeting CRISPR nucleases are not ",
             "supported at the moment.")
    }
    .checkBSGenome(bsgenome)
    if (is.null(bwa_index)){
        stop("bwa_index must be provided.")
    }
    .checkBwaIndex(bwa_index)
    
    cat(paste0("[runCrisprBwa] Using ", bsgenome@pkgname, " \n"))

    if (verbose){
        cat(paste0("[runCrisprBwa] Searching for ",
                   nucleaseName(crisprNuclease), " protospacers \n"))
    }

    # Getting nuclease info:
    pam.len     <- pamLength(crisprNuclease)
    pam.side <- pamSide(crisprNuclease) 
    spacer.len  <- unique(nchar(spacers))

    if (length(spacer.len)>1){
        stop("Spacers must all have the same length.")
    }
    spacer.len.from.nuclease <- spacerLength(crisprNuclease)
    if (!force_spacer_length){
        if (spacer.len.from.nuclease != spacer.len){
            stop("Length of provided spacers is ",
                 spacer.len,
                 ", but spacer length for the provided nuclease is ",
                 spacer.len.from.nuclease, 
                 ". Consider force_spacer_length=TRUE to overwrite the",
                 " nuclease spacer length. ")
        }
    } else {
        if (spacer.len.from.nuclease != spacer.len){
            message("Setting spacer length to be ", spacer.len)
            spacerLength(crisprNuclease) <- spacer.len
        }
    }
    
    # Performing alignment:
    aln <- runBwa(sequences=spacers,
                  bwa_index=bwa_index,
                  n_mismatches=n_mismatches)

    if (is.null(aln)){
        return(.emptyAlignments(n_mismatches))
    }

    aln$pam_site <- .getPamSiteFromBwaOutput(pos=aln$pos, 
                                             strand=aln$strand, 
                                             spacer.len=spacer.len,
                                             crisprNuclease=crisprNuclease)

    seq_choices <- seqnames(bsgenome)
    aln <- aln[aln$chr %in% seq_choices,,drop=FALSE]


    pams.canonical    <- pams(crisprNuclease,
                              primary=TRUE,
                              as.character=TRUE)
    pams.noncanonical <- pams(crisprNuclease,
                              primary=FALSE,
                              as.character=TRUE)
    if (nrow(aln)>0){
        aln <- .addPamSequences(aln,
                                bsgenome=bsgenome,
                                crisprNuclease=crisprNuclease)
        aln$canonical <- aln$pam %in% pams.canonical
        if (canonical & !ignore_pam){
            aln <- aln[aln$canonical,,drop=FALSE]
        } else if (!canonical & !ignore_pam){
            aln  <- aln[aln$pam %in% pams.noncanonical,,drop=FALSE]
        }
    } 
    

    if (nrow(aln)==0){
        aln <- .emptyAlignments(n_mismatches)
    } else {
        aln <- aln[order(aln$query, aln$n_mismatches),,drop=FALSE]
        aln$pos <- NULL
        colnames(aln)[colnames(aln)=="query"]  <- "spacer"
        rownames(aln) <- NULL
    }

    if (nrow(aln)>0){
        aln <- .addProtospacerSequences(aln,
                                        bsgenome=bsgenome,
                                        crisprNuclease=crisprNuclease)
        aln <- aln[!grepl("N", aln$protospacer),,drop=FALSE]
        aln <- .addMismatchInfo(aln, n_mismatches)
        aln <- aln[, .outputColumns(n_mismatches),drop=FALSE]
    }   
    
    return(aln)
}



#' @importFrom crisprBase pamLength pamSide
.getPamSiteFromBwaOutput <- function(pos, 
                                     strand,
                                     spacer.len,
                                     crisprNuclease=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pam.len  <- pamLength(crisprNuclease)
    pam.side <- pamSide(crisprNuclease)
    wh <- which(strand=="-")
    if (pam.side=="3prime"){
        pam_site <- pos + spacer.len
        if (length(wh)>0){
            pam_site[wh] <- pos[wh] - 1
        }
    } else {
        pam_site <- pos - pam.len
        if (length(wh)>0){
            pam_site[wh] <- pos[wh] +spacer.len + pam.len - 1
        }
    }
    return(pam_site)
}




#' @importFrom crisprBase getPamRanges
#' @importFrom BSgenome getSeq
.addPamSequences <- function(aln,
                             bsgenome,
                             crisprNuclease
){
    ranges <- getPamRanges(pam_site=aln$pam_site,
                           strand=aln$strand,
                           seqnames=aln$chr,
                           nuclease=crisprNuclease)
    aln$pam <- getSeq(bsgenome,
                      ranges,
                      as.character=TRUE)
    return(aln)
}




#' @importFrom crisprBase getSpacerRanges
#' @importFrom BSgenome getSeq
.addProtospacerSequences <- function(aln,
                                     bsgenome,
                                     crisprNuclease
){
    ranges <- getSpacerRanges(pam_site=aln$pam_site,
                              strand=aln$strand,
                              seqnames=aln$chr,
                              nuclease=crisprNuclease)
    aln$protospacer <- getSeq(bsgenome,
                              ranges,
                              as.character=TRUE)
    return(aln)
}







#' @importFrom Biostrings DNAStringSet
.addMismatchInfo <- function(aln,
                             n_mismatches
){
    if (n_mismatches>0){
        mat1 <- as.matrix(DNAStringSet(aln$spacer))
        mat2 <- as.matrix(DNAStringSet(aln$protospacer))
        diff <- mat1!=mat2
        pos <- apply(diff,1, which, simplify=FALSE)
        pos <- lapply(pos, function(wh){
            wh <- c(wh, rep(NA, n_mismatches-length(wh)))
            return(wh)
        }) 
        mm <- do.call(rbind, pos)
        colnames(mm) <- paste0("mm", seq_len(n_mismatches))
        aln[colnames(mm)] <- mm
    }
    return(aln)
}


.emptyAlignments <- function(n_mismatches){
    cols <- .outputColumns(n_mismatches)
    out <- data.frame(matrix(0, ncol=length(cols)))
    colnames(out) <- cols
    out <- out[-1,,drop=FALSE]
    return(out)
}


.outputColumns <- function(n_mismatches){
    cols <- c("spacer",
              "protospacer",
              "pam",
              "chr",
              "pam_site",
              "strand",
              "n_mismatches",
              "canonical")
    if (n_mismatches!=0){
        mm_cols <- paste0("mm", seq_len(n_mismatches))
        cols <- c(cols, mm_cols) 
    }
    return(cols)
}


