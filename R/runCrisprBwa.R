#' @title Find gRNA spacer alignments with bwa
#' @description Return bwa alignments for a list of gRNA spacer sequences.
#' 
#' @param spacers Character vector of DNA sequences corresponding
#'     to gRNA spacer sequences. Must all be of equal length. 
#' @param bsgenome BSgenome object.
#' @param bwa_index Path to the bwa index to be used for alignment.
#' @param crisprNuclease \code{CrisprNuclease} object. 
#' @param canonical Should only canonical PAM sequences be considered?
#'     TRUE by default.
#' @param ignore_pam If TRUE, will return all matches regardless of
#'     PAM sequence. FALSE by default.
#' @param n_mismatches Integer specifying maximum
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
#' @importFrom crisprBase getTargetRanges
#' @importFrom GenomeInfoDb seqnames seqlengths 
#' @importFrom BiocGenerics start end
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
        return(.emptyAlignments())
    }

    aln$pam_site <- .getPamSiteFromBwaOutput(pos=aln$pos, 
                                             strand=aln$strand, 
                                             spacer.len=spacer.len,
                                             crisprNuclease=crisprNuclease)
    aln <- aln[aln$pam_site>0,,drop=FALSE]
    if (nrow(aln)==0){
        return(.emptyAlignments())
    }
    


    seq_choices <- seqnames(bsgenome)
    aln <- aln[aln$chr %in% seq_choices,,drop=FALSE]


    pams.canonical    <- pams(crisprNuclease,
                              primary=TRUE,
                              as.character=TRUE)
    pams.noncanonical <- pams(crisprNuclease,
                              primary=FALSE,
                              as.character=TRUE)

    
    # Filtering out PAMs falling outside of chrs
    protoRanges <- getTargetRanges(seqnames=aln$chr,
                                   pam_site=aln$pam_site,
                                   strand=aln$strand,
                                   nuclease=crisprNuclease)
    chr_lens <- seqlengths(bsgenome)[as.character(seqnames(protoRanges))]
    valid <- BiocGenerics::start(protoRanges)>0 & BiocGenerics::end(protoRanges) <= chr_lens
    aln <- aln[valid,,drop=FALSE]

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
        aln <- .emptyAlignments()
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
        aln <- aln[, .outputColumns(),drop=FALSE]
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




#' @importFrom crisprBase getProtospacerRanges
#' @importFrom BSgenome getSeq
.addProtospacerSequences <- function(aln,
                                     bsgenome,
                                     crisprNuclease
){
    ranges <- getProtospacerRanges(pam_site=aln$pam_site,
                                   strand=aln$strand,
                                   seqnames=aln$chr,
                                   nuclease=crisprNuclease)
    aln$protospacer <- getSeq(bsgenome,
                              ranges,
                              as.character=TRUE)
    return(aln)
}





.emptyAlignments <- function(){
    cols <- .outputColumns()
    out <- data.frame(matrix(0, ncol=length(cols)))
    colnames(out) <- cols
    out <- out[-1,,drop=FALSE]
    return(out)
}


.outputColumns <- function(){
    cols <- c("spacer",
              "protospacer",
              "pam",
              "chr",
              "pam_site",
              "strand",
              "n_mismatches",
              "canonical")
    return(cols)
}


