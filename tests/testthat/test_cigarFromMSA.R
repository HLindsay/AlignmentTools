context("Test reconstruction of CIGAR string from multiple sequence alingnment")

maln <- DNAStringSet(c("A--CCC----TTTTT","AAA---GGGGTTTTT"))
ref <- DNAString("A--CCC----TTTTT")
no_ref_result <- DNAStringSet(c("N--NNN----NNNNN","NNN---NNNNNNNNN"))
ref_result <- DNAStringSet(c("N++NNN++++NNNNN", "N++---++++NNNNN"))

test_that("recodeAln correctly checks inputs", {
    # Ref isn't a DNAString or XString
    expect_error(recodeAln(maln, ref = "AAAAAAAAA"))
    # Ref and alns not the same width
    expect_error(recodeAln(DNAStringSet("AA"), ref = DNAString("A")))
    # aln is DNAString but ref.gap not in DNA_ALPHABET
    expect_error(recodeAln(maln, ref = ref, ref.gap = 1))
    # This should be okay if maln is an XStringSet
    expect_equal(recodeAln(maln, ref, ref.gap = 1),
                 as(c("N11NNN1111NNNNN", "N11---1111NNNNN"), "XStringSet"))
})


test_that("Optional recodeAln parameters work and handles multiple gaps", {
    expect_equal(recodeAln(maln), no_ref_result)
    expect_equal(recodeAln(maln, ref = ref), ref_result)
    expect_equal(recodeAln(maln, ref.gap = "+"), no_ref_result)
    expect_equal(recodeAln(maln, nuc.chars = c("A","C","G")),
                 DNAStringSet(c("N--NNN----TTTTT","NNN---NNNNTTTTT")))
    expect_equal(recodeAln(maln, nuc.chars = c("A","C","G"), gap = "I"),
                 as(c("NIINNNIIIITTTTT","NNNIIINNNNTTTTT"), "XStringSet"))
    expect_equal(recodeAln(DNAStringSet("NNN")), DNAStringSet("NNN"))
})


test_that("recodeAln accepts numeric codes",{
    expect_equal(recodeAln(maln, ref = ref, nuc = 1, ref.gap = 2, gap = 3),
                 as(c(122111222211111, 122333222211111), "XStringSet"))
})


test_that("cigarFromMSA behaves correctly with different references"){
    aln <- DNAStringSet(c("AA--CC","AATTCC","AATTC-"))
    names(aln) <- c("A","B","C")
    expect_equal(cigarFromMSA(aln), c("2M2D2M","6M","5M1D"))
    expect_equal(cigarFromMSA(aln, ref = "A"), c("2M2I2M","2M2I1M1D"))
    expect_equal(cigarFromMSA(aln, ref = "B"), c("2M2D2M","5M1D"))
    expect_equal(cigarFromMSA(aln, ref = "C"), c("2M2D1M1I","5M1I"))
}
