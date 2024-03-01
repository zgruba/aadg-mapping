# Genome mapping algorithm

Read mapping algorithm that is designed to work on reads of length ∼1kbp with error rate 5 − 10%.

Made for Bioinformatics algorithms for genomic data analysis classes at MIMUW in the winter semester of 2023/2024.

## Assumptions:
All reads come from the reference sequence, but contain errors (substitutions, insertions and deletions of single nucleotides) resulting from the sequencing process. Errors occur independently at each position at the assumed error rate, so the total number of errors may slightly exceed 10% of the read length. A read is mapped correctly when the mapping coordinates differ from the actual coordinates of the fragment it comes from by $\leq$ 20bp.

## Short description:

The approach involves dividing reads into small seeds so that as many of them as possible would match the reference genome. A read of length 1000 is divided into fragments of length 29 (PART_k = 34) with a step of 5, with these parameters it is possible to maintain very high detectability while significantly reducing the average time allocated to a read for a reference genome. These matches are extended to maximum exact matches (MEM), thereby reducing and merging the seeds with the same exact match into one MEM. As for alignments found using dynamic programming tables, once the distance of 110 was exceeded, a cutoff occurred, and alignments were searched within a distance of 20 from the diagonal. In terms of improving the correctness of mappings, I set the tolerance at the ends of alignment to 10, ensuring that for 20M files, the differences in end alignments do not exceed 10. For a distance of 20 from the diagonal (which was set in this case) for these files, values exceeding 20 appeared.

Used functions:
- `kEditDp` from https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/masternotebooks/CG_kEditDp.ipynb
- `bwt` from https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_FmIndex.ipynb with `suffix array`: https://code.google.com/archive/p/pysuffix/

## Usage:
``` 
python3 mapper.py reference.fasta reads.fasta output.txt 
```