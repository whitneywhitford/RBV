# RBV: Read balance validator
## User manual and guide


## Overview
**RBV:** Read balance validator is a Copy number variant (CNV) validator. It uses the relative reads for the reference and alternative allele at a given position (the read balance) to determine the probability of a putative CNV being true.

The read balance distribution is different based on the copy number of the genetic sequence. As seen blow, single copy deletions result in haploid genetic sequence with zero heterozygous single nucleotide polymorphisms (SNPs), with the distribution ratio of reads centred around 1. In regions with diploid sequence, the majority of SNPs are homozygous, with the read balance centred around 1. The remaining SNPs are heterozygous with the distribution centred around 0.5, representing half of the reads from one allele and half from the other. Duplications also have the majority of SNPs being homozygous, however the heterozygous SNPs are represented by 2 different peaks in distribution. For a duplication the resulting triplicate genomic sequence would result in the distribution seen below where there are peaks centred around 0.33 and 0.66 which represent regions with 1 and 2 copies respectively.

<img src="./images/haploid_readbal.png" width="290"><img src="./images/diploid_readbal.png" width="290"><img src="./images/triploid_readbal.png" width="290">

Thus, RBV exploits this difference in distribution to validate CNVs. The python package presented here simultaneously interrogates the probability of multiplications and deletions within a provided list of CNVs to investigate. This allows for prioritisation of CNVs for causation in molecular diagnostic testing bioinformatic pipelines. Additionally, RBV can be used for validation and inheritiance hypothesis testing of causative variants by using RBV on multiple members of a pedigree.

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
- PyVCF
- scipy
- shutil
- datetime

External Programs:
- bedtools
- samtools

## Usage
RBV has 5 required inputs: CNV file containing CNV coordinates, a VCF file, the reference for the genome used to align the genome, the level of sequecing that was performed (WGS or WES), and the variant caller used. Therefore an RBV command in its simplest form is:

  	python RBV.py --ref ref.fa --CNV_file CNV.interval_file --vcf variants.vcf --seq_type type --calling_method variantCaller
  
### Arguments
  ~~~~ -h, --help       show this help message and exit
  --ref REF             REQUIRED. Path to reference sequence (including file name).
  --CNV_file CNV_FILE   REQUIRED. Picard-style interval_list containing targets to use in CNV analyses.
		       	      Must be typical interval_list format: 1-based indexing, with the six
		              columns being the chromosome name, start coordinate, stop coordinate,
		       	      plus sign, and predicted CNV type.
  --gap_file GAP_FILE   Picard-style interval_list containing gaps in the reference to
		              mask for random generation. Must be typical interval_file format:
		              1-based, indexing, with the six columns being the chromosome name,
		              start coordinate, stop coordinate, plus sign, and type.
  --vcf VCF             REQUIRED. VCF file containing the variants for RBV analyses.
			      VCF must be single individual vcf file but can be derived from joint
			      calling pipeline.
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory. RBV will create a temporary directory and output file within this
                        directory.
  --sample_id SAMPLE_ID, -id SAMPLE_ID
                        Name/ID of sample - for use in and file naming. Default is sample
  --variant_quality_cutoff VARIANT_QUALITY_CUTOFF, -vqc VARIANT_QUALITY_CUTOFF
                        Consider all SNPs with a quality greater than or equal
                        to this value. Default is 20.
  --read_depth_cutoff", "-rdc"
                        Consider all SNPs with a read depth greater than or equal to this value. Default is 10.
  --readbal_cutoff", "-rbc"
		              For deletion analyses, consider all heterozygous SNPs with a
		              read balance less than this value. Default is 0.65.
  --variant_permutations VARIANT_PERMUTATIONS
                        Number of permutations to use for heterozygous read
                        balance analyses. Default is 10000
  --window_permutations WINDOW_PERMUTATIONS
                        Number of permutations to use for read balance
                        analyses.Default is 1000
  --seq_type            REQUIRED. Type of genome sequencing for RBV analysis.
		               Options: WGS, WES.		
  --interval_file       Picard-style interval_list containing interval coordinates
		              used for variant calling. Must be typical interval_list format: 1-based
		              indexing, with the six columns being the chromosome name, start
		              coordinate, stop coordinate, plus sign, and type.
		              REQUIRED for if using WES seq_type.
  --calling_method      REQUIRED. Variant calling method used for VCF generation.
		              Options: haplotypecaller, samtools, freebayes, platypus.
~~~~ 


## Contributors

RBV is made by:

- Whitney Whitford
- Klaus Lehnert
- Russell Snell
- Jessie Jacobsen

## Support

Please report any issues or questions by email to whitney.whitford@auckland.ac.nz
