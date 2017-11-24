#'@import Biostrings
#'@import GenomicAlignments
#'@import GenomicRanges
#'@import IRanges
#'@title gappedRef
#'@description Given a set of sequences aligned to a reference,
#'introduce gaps in the reference sequence corresponding to insertion
#'locations in the aligned sequences.
#'TO DO: HOW ARE AMBIGUOUS ALIGNMENTS TREATED?
#'@author Helen Lindsay
#'@param reference Ungapped reference sequence
#'@param alignments Sequences aligned to the reference sequence
#'@param ... Extra arguments
#'@rdname gappedref
#'@export
setGeneric("gappedRef", function(reference, alignments, ...) {
  standardGeneric("gappedRef")})


#'@rdname gappedref
setMethod("gappedRef", signature("DNAString", "GAlignments"),
          function(reference, alignments, ...){

    .checkForMcols(alignments, "seq", "gappedRef")

    cigars <- GenomicAlignments::cigar(alignments)
    ins_on_qry <- cigarRangesAlongQuerySpace(cigars, ops = "I")
    ins_on_ref <- cigarRangesAlongReferenceSpace(cigars, ops = "I")

    ins_locs <-IRanges(unlist(start(ins_on_ref)),
                       width = unlist(width(ins_on_qry)))

})


