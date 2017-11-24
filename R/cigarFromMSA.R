setOldClass("MsaDNAMultipleAlignment")

# cigarFromMSA -----
#'@title CIGAR from Multiple Sequence Alignment
#'@description Construct a CIGAR string from a multiple sequence
#'alignment
#'@param aln A multiple sequence alignment
#'@param ... Extra arguments
#'@author Helen Lindsay
#'@export
#'@rdname cigarFromMSA
setGeneric("cigarFromMSA", function(aln, ...) {
  standardGeneric("cigarFromMSA")})

#'@param ref character(1) Name of sequence in aln to use as a
#'reference (Default: NULL)
#'@rdname cigarFromMSA
setMethod("cigarFromMSA", signature("MsaDNAMultipleAlignment"),
          function(aln, ..., ref = NULL){

    maln <- Biostrings::unmasked(aln)
    cigarFromMSA(maln, ref = ref, ...)
})


#'@rdname cigarFromMSA
setMethod("cigarFromMSA", signature("DNAMultipleAlignment"),
          function(aln, ..., ref = NULL){

            maln <- Biostrings::unmasked(aln)
            cigarFromMSA(maln, ref = ref, ...)
})


#'@param ref character(1) name of the alignment in alns to use as a reference
#'(Default: NULL)
#'@rdname cigarFromMSA
#'@examples
#'aln <- DNAStringSet(c("AA--CC","AATTCC","AATTC-"))
#'names(aln) <- c("A","B","C")
#'cigarFromMSA(aln)
#'cigarFromMSA(aln, ref = "A")
#'cigarFromMSA(aln, ref = "B")
#'cigarFromMSA(aln, ref = "C")
setMethod("cigarFromMSA", signature("DNAStringSet"),
          function(aln, ..., ref = NULL){

    # If a reference sequence name is provided, separate reference
    # and put "+" in non-reference sequence where the reference has "-"

    ref_sq <- NULL
    if (! is.null(ref)){
      if (! ref %in% names(aln)){
        stop(sprintf("%s not found in names(aln)", ref))
      }
      ref_sq <- aln[[ref]]
      aln <- aln[! names(aln) == ref]
    }

    recoded <- recodeAln(aln, ref = ref_sq, nuc = "M", ref.gap = "I", gap = "D")
    rles <- S4Vectors::Rle(as.vector(unlist(recoded)))
    ptn <- IRanges::PartitioningByWidth(nchar(recoded))
    rles <- relist(rles, ptn)
    result <- paste(runLength(rles), runValue(rles), sep = "", collapse = "")
    names(result) <- names(aln)
    result

}) # -----


# recodeAln -----
#'@title recodeAln
#'@description  Internal AlignmentTools function for recoding an alignment.
#'@param aln       (DNAStringSet)
#'@param ref       (DNAString) the reference sequence (Default: NULL)
#'@param ref.gap   (character(1)) Replace locations where the reference
#'sequence has a gap ("-") with ref.gap in alns
#'@param gap       (character(1))
#'@param nuc       (character(1)) Replace non-gap bases with nuc
#'@param nuc.chars (character(n)) Vector defining characters to replace with "nuc"
#'@examples
#'recodeAln(DNAStringSet(c("AA--CC", "AATTC-")), ref = DNAString("AATT--"))
#'@author Helen Lindsay
recodeAln <- function(aln, ref = NULL, ref.gap = "+", gap = NULL,
                      nuc = "N", nuc.chars = NULL){

    if (! class(aln) %in% c("XStringSet", "DNAStringSet", "BStringSet")){
      stop("aln should be a StringSet object")
    }
    if (! class(ref) %in% c("DNAString", "XString", "BString", "NULL")){
      stop("ref should be either NULL or a DNAString or XString")
    }
    if (! all(c(ref.gap, nuc, gap) %in% Biostrings::DNA_ALPHABET)){
       aln <- Biostrings::BStringSet(aln)
    }

    # If no nuc.chars provided, use DNA bases and ambiguities
    if (is.null(nuc.chars)){
      nuc.chars <- setdiff(Biostrings::DNA_ALPHABET, c("+","-","."))
    }

    aln_chrs <- t(as.matrix(aln))
    is_nuc <- IRanges::IRanges(aln_chrs %in% nuc.chars)
    result <- Biostrings::replaceAt(unlist(aln), is_nuc, strrep(nuc, width(is_nuc)))
    del_mat <- aln_chrs == "-"

    if (! is.null(gap)){
      is_del <- IRanges::IRanges(as.vector(del_mat))
      result <- Biostrings::replaceAt(result, is_del, strrep(gap, width(is_del)))
    }
    result <- relist(result, aln)

    # Replace positions in aln where ref has a gap ("-") with ref.gap
    if (! is.null(ref)){
      if (! length(unique(c(nchar(ref), nchar(aln))) ) == 1){
        stop("ref and aln should all be the same width")
      }
      ref_del <- as.vector(ref) == "-"
      ref_del_rng <- IRanges::IRanges(ref_del)
      result <- Biostrings::replaceAt(result, ref_del_rng,
                                   strrep(ref.gap, width(ref_del_rng)))
      both_del <- IRanges::IRangesList(apply(ref_del & del_mat, 2, IRanges))
      result <- Biostrings::replaceAt(result, unname(both_del))
    }

    result

} # -----
