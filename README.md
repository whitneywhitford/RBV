# RBV: Read balance validator
## User manual and guide


## Overview
**RBV:** Read balance validator is a Copy number variant (CNV) validator. It uses the relative reads for the reference and alternative allele at a given position (the read balance) to determine the probability of a putative CNV being true.

The read balance distribution is different based on the copy number of the genetic sequence. As seen blow, single copy deletions result in haploid genetic sequence with zero heterozygous single nucleotide polymorphisms (SNPs), with the distribution ratio of reads centred around 1. In regions with diploid sequence, the majority of SNPs are homozygous, with the read balance centred around 1. The remaining SNPs are heterozygous with the distribution centred around 0.5, representing half of the reads from one allele and half from the other. Duplications also have the majority of SNPs neing homozygous, however the heterozygous SNPs are represented by 2 different peaks in distribution. For a duplication the resulting triplicate genomic sequence would result in the distribution seen below where there are peaks centred arounf 0.33 and 0.66 which represent regions with 1 and 2 copies respectively.

<img src="./images/haploid_readbal.png" width="290"><img src="./images/diploid_readbal.png" width="290"><img src="./images/triploid_readbal.png" width="290">

## Citation



## Obtaining
To download RBV, please use git to download the most recent version.  Currently, the RBV is hosted on github, and can be obtained via:

    git clone --recursive https://github.com/whitneywhitford/RBV.git

Note the use of --recursive.  This is required in order to download all nested git submodules for external repositories.

## Requirements
Python: 2.7

Python packages:
- matplotlib
- numpy
- pybedtools
- PyVCF
- scipy
- shutil
- datetime

External Programs:
- bedtools
- samtools

## Usage
RBV has 5 required inputs: CNV file containing CNV coordinates, a VCF file, the reference for the genome used to align the genome, the level of sequecing that was performed (WGS or WES), and the variant caller used. Therefore an RBV command in its simplest form is:

  	python RBV.py --ref ref.fa --CNV_bed CNV.bed --vcf variants.vcf --seq_type type --calling_method variantCaller
  
### Arguments
  ~~~~ -h, --help            show this help message and exit
  --ref REF             REQUIRED. Path to reference sequence (including file name).
  --CNV_bed CNV_BED     REQUIRED. Bed file containing targets to use in CNV analyses Must be typical bed format,
                        0-based indexing, with the first four columns the chromosome name, start coordinate,
                        stop coordinate, and predicted CNV type.
  --gap_bed GAP_BED     REQUIRED. Bed file containing gaps in the reference to mask for random generation. Must 
                        be typical bed format, 0-based indexing,with the first four columns the chromosome name,
                        start coordinate,stop coordinate, and predicted type.
  --vcf VCF             REQUIRED. VCF file containing the variants for RBV analyses.VCF must be generated using
                        HaploTypeCaller.
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory. RBV will create a temporary directory and output file within this
                        directory.
  --sample_id SAMPLE_ID, -id SAMPLE_ID
                        Name/ID of sample - for use in plot titles and file
                        naming. Default is sample
  --variant_quality_cutoff VARIANT_QUALITY_CUTOFF, -vqc VARIANT_QUALITY_CUTOFF
                        Consider all SNPs with a quality greater than or equal
                        to this value. Default is 20.
  --variant_permutations VARIANT_PERMUTATIONS
                        Number of permutations to use for heterozygous read
                        balance analyses. Default is 10000
  --window_permutations WINDOW_PERMUTATIONS
                        Number of permutations to use for read balance
                        analyses.Default is 1000
~~~~ 


## Contributors

RBV is made by:

- Whitney Whitford

## Support

Please report any issues or questions by email to whitney.whitford@auckland.ac.nz
