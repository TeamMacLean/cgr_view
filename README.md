# cgr_view

Draw CGRs of DNA

## Requirements

    * jellyfish on the path
    
## Purpose

Library of functions for creating CGR matrices. Also has some functions for making images of CGRS.


## Quickstart

+ Input fasta file, get cgr
    * one cgr for each entry in the fasta file
        `cgr.from_fasta("my_seqs.fa", outfile = "my_cgrs", k = 7)`
    * just one cgr with all entries in the fasta file (eg for genomes and contigs)
        `cgr.from_fasta("my_genome.fa", outfile = "genome_cgr", k = 7, as_single = True)`



## Workflow:

1. make kmer count db in jellyfish from fasta -> generate cgr from db.
2. optionally merge cgrs into single cgr as separate channels
3. stack all composed cgrs into an array of cgrs 
4. save as numpy binary (.npy) files


## Usage:

    `import cgr`

    1. Make kmer count db
    `run_jellyfish("test_data/NC_012920.fasta", 11, "11mer.jf")`
    `run_jellyfish("test_data/NC_012920.fasta", 10, "10_mer.jf")`
    
    2. Load CGRs from kmer count db
    `cgr1 = cgr.cgr_matrix("/Users/macleand/Desktop/athal-5-mers.jf")
    `cgr2 = cgr.cgr_matrix("test_data/five_mer.jf")`
   

    3. Draw a cgr and save to file
        * just one cgr, can choose colour (value of 'h') and which channel to put cgr in
        `cgr.draw_cgr(cgr1, h = 0.64, v = 1.0, out = "my_cgr.png", resize = 1000, main = "s" )`
        * two cgrs, first in tuple goes in 'h', second goes in 's'. Can set 'v'
        `cgr.draw_cgr( (cgr1, cgr1), v = 1.0, out = "two_cgrs.png"`
        * three cgrs 'h','s' and 'v' are assigned as order in tuple
        `cgr.draw_cgr( (cgr1, cgr1, cgr1) )`
    
    4. Save a single cgr into a text file
    `save_as_csv(cgr1, file = "out.csv")`

    5. Join n cgrs into one, extending the number of channels ... 
    `merged_cgr = cgr.join_cgr( (cgr1, cgr2, ... ) )`
    
    6. Write to file (numpy binary)
    `cgr.save_cgr("my_cgr, merged_cgr )`
    
    7. Input fasta file, get cgr
        * one cgr for each entry in the fasta file
        `cgr.from_fasta("my_seqs.fa", outfile = "my_cgrs", k = 7)`
        * just one cgr with all entries in the fasta file (eg for genomes and contigs)
        `cgr.from_fasta("my_genome.fa", outfile = "genome_cgr", k = 7, as_single = True)
    
    
    
    


    
    