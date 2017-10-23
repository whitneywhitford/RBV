gunzip NA12878_Matsunami2013_CNVs.vcf.gz
python ../scripts/RBV.py --CNV_file CNVs_wfakes.interval_file --interval_file safe_calling_regions.interval_list --vcf NA12878_Matsunami2013_CNVs_chronly.vcf --sample_id NA12878_Matsunami2013_CNVs_wfake --output_dir RBV_NA12878_Matsunami2013_CNVs_wfake --seq_type WGS --calling_method haplotypecaller
