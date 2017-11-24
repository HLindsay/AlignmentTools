context("Test utility functions")

test_that(".checkForMcol correctly checks for metadata columns",{

    expect_error(.checkForMcols(GAlignments(),
                                mcols = c("fish","cat","dog"),
                                err.func = stop,
                                func.nm = "animal_func"))

    gr <- GenomicRanges::GRanges("chr1", cigar = "10M", strand = "+", flag=0,
                                 IRanges::IRanges(start = 5, width = 10))
    galns <- as(gr, "GAlignments")

    # Expect "flag" correctly identified
    expect_true(.checkForMcols(galns, "flag"))

    # If checking for mcol "sequence", .checkForMcol returns FALSE
    # and prints a warning
    expect_false(.checkForMcols(galns, "seq"))
    expect_warning(.checkForMcols(galns, "seq"))

    # If err.func stop is provided, raise an error rather than a warning
    expect_error(.checkForMcols(galns, "seq", err.func = stop))

    # Expect FALSE if checking for a combination of mcols and only some present
    expect_false(.checkForMcols(galns, c("seq", "flag")))

    # Adding mcol "seq" should make tests pass
    seqs <- Biostrings::DNAStringSet("ACTGACTGAC")
    mcols(galns)$seq <- seqs

    expect_true(.checkForMcols(galns, c("seq", "flag")))
    expect_true(.checkForMcols(galns, c("seq", "flag"), err.func = stop))
})


