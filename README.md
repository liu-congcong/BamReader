# BamReader

Extract the reads and their locations from a BAM file that meet the MAPQ threshold.

|Sequence ID|Position|Cigar|Read|
|:-:|:-:|:-:|:-|
|chr1|25|28M|TGGATGGCGCTCAAGTACTGCTTCAAGA|
|chr1|54|12M|CGAAAGCGACTG|
|chr2|1|14M|GTCGGCGCCTATAC|
|chr3|2105|18M|GCGAGAGCCTGGCCCTAT|
|chr4|27|10M|CCTCCATCGG|
|chr4|259|5M|ATCGG|

## Usage

```shell
BamReader -bam BAM -mapq MAPQ > OUTPUT
```
