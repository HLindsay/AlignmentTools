#'@title getInsSeqs
#'@description Get insertion sequences
#'@author Helen Lindsay
#'@rdname getInsSeqs
#'@export
setGeneric("getInsSeqs", function(alignments, ...) {
  standardGeneric("getInsSeqs")})


#'@rdname getInsSeqs
setMethod("getInsSeqs", signature("GAlignments"),
          function(alignments, ...){

    .checkForMcols(alignments, "seq", "gappedRef")

    ins_on_qry <- cigarRangesAlongQuerySpace(cigar(alignments),
                                             ops = "I")
    extractAt(mcols(alignments)$seq, ins_on_qry)
})

# get clipped seqs

