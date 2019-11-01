# cgr_view

Draw CGRs of DNA

## Requirements

    * jellyfish on the path
    
## Purpose

Library of functions for creating CGR matrices. Also has some functions for making images of CGRS.

## Workflow:
    1. make kmer count db in jellyfish from fasta -> generate cgr from db -> save cgr.
    2. also merge two cgrs into single array

## Usage:

    `import cgr_view`

    1. Make kmer count db
    `run_jellyfish("test_data/NC_012920.fasta", 11, "11mer.jf")`
    `run_jellyfish("test_data/NC_012920.fasta", 10, "10_mer.jf")`
    
    2. Load CGRs from kmer count db
    `cgr1 = cgr_matrix("/Users/macleand/Desktop/athal-5-mers.jf")
    `cgr2 = cgr_matrix("test_data/five_mer.jf")`
   

    3. Preview a single one channel cgr
    `show_image(cgr1)`

    4. Save a single cgr into a text file
    `save_as_csv(cgr1, file = "out.csv")`

    5. Join two cgrs into one
    `merged_cgr = join_cgr(cgr1, cgr2)`

    6. Write merged cgr to image file (not single channel cgr), by creating a three channel array to convert
    `padded = pad_two_channel_to_three(merged_cgr)`
    `from PIL import Image`
    `img = Image.fromarray(padded, 'RGB')`
    `img.save('out.png')`