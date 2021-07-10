## Summary:
Visualize FATSA/FASTQ to reveal primer/adapter pattern

## Usage:
    Rscript fa_baseviz.R --faFile1 R1.fq [--faFile2 R2.fq] [--primerFile primers.fa] [--outFile output.png] [--main title] [--width 7] [--height 7] [--res 300]
## Options:
        -h, --help
                Show this help message and exit
        --faFile1=FAFILE1
                R1 filename (FASTA or FASTQ)
        --faFile2=FAFILE2
                R2 filename (FASTA or FASTQ)
        --sortBy=SORTBY
                in which mate's order (choose from 1 or 2)
        --primerFile=PRIMERFILE
                primer filename (FASTA)
        --outFile=OUTFILE
                output filename (PNG)
        --main=MAIN
                title
        --width=WIDTH
                width (inch)
        --height=HEIGHT
                height (inch)
        --res=RES
                resolution (dpi)

## Example:
    Rscript fa_baseviz.R --faFile1 test/R1.fq --faFile2 test/R2.fq --primerFile test/primers.fa --outFile test/output.png
