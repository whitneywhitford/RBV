# RBV main program

import argparse
import logging
import os
import subprocess
import sys
import pysam
from scipy import stats
import vcf
import numpy as np
import make_random_windows
from shutil import copyfile
from shutil import rmtree
from datetime import datetime




"""gapfile = sys.argv[1]
CNV_file = sys.argv[2]
variant_permutations = int(sys.argv[3])
vcf_file = sys.argv[4]
qualCutoff = int(sys.argv[5])
sample_id = sys.argv[6]
window_permutations = int(sys.argv[7])
ref_file = sys.argv[8]
out_dir = sys.argv[9]"""








def parse_args():
	"""Parse command-line arguments"""

	parser = argparse.ArgumentParser(description="RBV")

	parser.add_argument(
		"--ref", required=True,
		help="REQUIRED. Path to reference sequence (including file name).")
	
	parser.add_argument(
		"--CNV_bed", required=True,
		help="REQUIRED. Bed file containing targets to use in CNV analyses "
		"Must be typical bed format, 0-based indexing, with the first four "
		"columns the chromosome name, start coordinate, stop coordinate,"
		"and predicted CNV type.")

	parser.add_argument(
		"--gap_bed", default=None,
		help="Bed file containing gaps in the reference to mask "
		"for random generation. Must be typical bed format, 0-based indexing,"
		"with the first four columns the chromosome name, start coordinate,"
		"stop coordinate, and predicted type.")
		
	parser.add_argument(
		"--vcf", required=True,
		help="REQUIRED. VCF file containing the variants for RBV analyses."
		"VCF must be generated using HaploTypeCaller.")
	
	parser.add_argument(
		"--output_dir", "-o", default="RBV",
		help="Output directory. RBV will create a temporary "
		"directory and output file within this directory.")

	parser.add_argument(
		"--sample_id", "-id", default="sample",
		help="Name/ID of sample - for use in plot titles and file naming. "
		"Default is sample")

	# Variant Calling Flags
	parser.add_argument(
		"--variant_quality_cutoff", "-vqc", type=int, default=20,
		help="Consider all SNPs with a quality greater than or "
		"equal to this value. Default is 20.")
		
	parser.add_argument(
		"--variant_permutations", type=int, default=10000,
		help="Number of permutations to use for heterozygous read balance "
		"analyses. Default is 10000")

	parser.add_argument(
		"--window_permutations", type=int, default=1000,
		help="Number of permutations to use for read balance analyses."
		"Default is 1000")
		
	parser.add_argument(
		"--seq_type", required=True,
		help="REQUIRED. Type of genome sequencing for RBV analysis."
		"Options: WGS, WES.")
		
	parser.add_argument(
		"--interval_file", default=None,
		help="Bed file containing interval file used for variant calling."
		"Must be typical bed format, 0-based indexing, with the first three "
		"columns the chromosome name, start coordinate, stop coordinate."
		"REQUIRED for if using WES seq_type.")
	
	parser.add_argument(
		"--calling_method", required=True,
		help="REQUIRED. Variant calling method used for VCF generation."
		"Options: haplotypecaller, samtools, freebayes.")
	
	args = parser.parse_args()
	
	if args.seq_type not in ["WGS", "WES"]:
		sys.exit(
			"Error. sequencing type must be WGS or WES.")
	
	if args.calling_method not in ["haplotypecaller", "samtools", "freebayes"]:
		sys.exit(
			"Error. calling method must be haplotypecaller, samtools or freebayes.")
	
	if args.seq_type == "WES" and args.interval_file == None:
		sys.exit(
			"Error. Interval file required for Whole Exome Sequencing"
			"RBV analysis")
	
	# Return arguments namespace
	return args


	
	
	
	
	
	
	

#make_random_vars
	
def rand_het_sites(no_of_samples,vcf_file):	
	from random import randrange
	count=1

	vcf=open(vcf_file).readlines()
	a = 0
	comment_count = 0
	while a <= 1:
		for i in range(len(vcf)):
			if vcf[i].startswith('#'):
				comment_count += 1
			else:
				a += 1
	#print comment_count

	# Choose variant
	var_lines = []
	while count<=no_of_samples:
		line_number = randrange((comment_count - 1), len(vcf))
		
		#print line_number
		
		variant_line = vcf[line_number]
		cols = variant_line.strip('\n').split('\t')
		chr = cols[0]
		coord = cols[1]
		
		str_chroms = []
		for num in range(1,23):
			str_chroms.append(str(num))

		if chr not in str_chroms:
			continue
		
		var_lines.append(variant_line)
		count+=1
	return var_lines


def vcf_fetch(vcf_infile,vcf_outfile, chr, start, stop):
	vcf_reader = vcf.Reader(open(vcf_infile, 'r'))
	vcf_writer = vcf.Writer(open(vcf_outfile, 'w'), vcf_reader)
	# fetch all records in window/CNV from base start through stop
	for record in vcf_reader.fetch(chr, start, stop):  
		vcf_writer.write_record(record)
	vcf_writer.close()

	
def het_count(vcf_file):
	het_count = 0
	vcf = open(vcf_file).readlines()
	for i in range(len(vcf)):
		line = vcf[i]
		if line.startswith('#'):
			continue
		cols = line.strip().split()
		call = cols[9].split(':')
		vars_raw = call[0].replace('|', '/')
		vars = vars_raw.split('/')
		if vars[0] != vars[1]:
			het_count += 1
		else:
			continue
	
	#os.remove(vcf_file)
	
	return het_count

	
def empirical_pvalue(window_het_snps, vcf_outfile, num_window_samples):
	sample_cum_sum = 0
	sample_het_snps = het_count(vcf_outfile)
	
	if sample_het_snps > max(window_het_snps):
		sample_cum_sum_pvalue = 1.000
	
	else:
		window_het_snp_sum = np.zeros(((max(window_het_snps)+1),), dtype=np.float)
		
		
		
		for i in range(len(window_het_snps)):
			window_het_snp_count = window_het_snps[i]
			window_het_snp_sum[window_het_snp_count] += 1
		
		for n in range((sample_het_snps+1)):
			sample_cum_sum += window_het_snp_sum[n]
		
		sample_cum_sum_pvalue = sample_cum_sum / num_window_samples
		
	return sample_het_snps,sample_cum_sum_pvalue

def haplotypecaller_vcf(variant_cols):
	allele1 = variant_cols[9].split(':')[1].split(',')[0]
	allele2 = variant_cols[9].split(':')[1].split(',')[1]
	
	return allele1, allele2
	
def samtools_vcf(variant_cols):
	allele1F = variant_cols[7].split(';')[12].split('=')[1].split(',')[0]
	allele1R = variant_cols[7].split(';')[12].split('=')[1].split(',')[1]
	allele1 = allele1F + allele1R
	allele2F = variant_cols[7].split(';')[12].split('=')[1].split(',')[2]
	allele2R = variant_cols[7].split(';')[12].split('=')[1].split(',')[3]
	allele2 = allele2F + allele2R
	
	return allele1, allele2
	
def freebayes_vcf(variant_cols):
	allele1F = variant_cols[7].split(';')[12].split('=')[1].split(',')[0]
	allele1R = variant_cols[7].split(';')[12].split('=')[1].split(',')[1]
	allele1 = allele1F + allele1R
	allele2F = variant_cols[7].split(';')[12].split('=')[1].split(',')[2]
	allele2R = variant_cols[7].split(';')[12].split('=')[1].split(',')[3]
	allele2 = allele2F + allele2R
	
	return allele1, allele2
	
	
def readbal(variant_lines,qualCutoff,calling_method):
	
	readBalance = []
	
	for i in range(len(variant_lines)):
		if variant_lines[i].startswith('#'):
			continue
		
		cols = variant_lines[i].strip('\n').split('\t')
		GT = cols[9].split(':')[0]
		vars_raw = GT.replace('|', '/')
		vars = vars_raw.split('/')
		
		#het SNPS only
		if vars[0] != vars[1]:
			qual = float(cols[5])
			if qual < qualCutoff:
				continue
			if GT == './.' or GT == '.|.':
				continue
			
			allele_reads = globals()[calling_method + "_vcf"]
			allele1, allele2 = allele_reads(cols)
			if ',' in allele1 or ',' in allele2:
				continue
			if float(allele1) + float(allele2) == 0:
				continue
			if int(allele1) >= int(allele2):
				ReadRatio = float(allele1)/ (float(allele1) +  float(allele2))
			else:
				ReadRatio = float(allele2)/ (float(allele1) +  float(allele2))
				
			# Add to arrays
			readBalance.append(ReadRatio)
			
	return readBalance

	
def random_readbal(no_of_samples,vcf_file,qualCutoff,calling_method):
	
	variant_lines = rand_het_sites(no_of_samples,vcf_file)
	
	random_readbal = readbal(variant_lines,qualCutoff,calling_method)
	
	return random_readbal

def prep_vcf(vcf_file, out_dir):
	working_vcf = os.path.join(out_dir, "tmp", "working_vcf.vcf")
	
	copyfile(vcf_file, working_vcf)
	os.system("bgzip " + working_vcf)
	bgzip_working_vcf = str(working_vcf) + ".gz"
	os.system("tabix -p vcf " + bgzip_working_vcf)
	
	return bgzip_working_vcf



def main():
	
	
	args = parse_args()
	
	
	#setup output directories
	if os.path.isdir(args.output_dir):
		print "[" + str(datetime.now()) + "] Output directory already exists. Please enter new Output directory"
	
	else:
		os.mkdir(args.output_dir)
		tmp = str(args.output_dir) + "/tmp"
		os.mkdir(tmp)


		ref_fai = str(args.ref) + ".fai"	#check if all fai are just fasta.fai/fa.fai
		
		if not os.path.isfile(ref_fai):
			os.system("samtools faidx  " + args.ref)
			ref_fai_base = os.path.basename(ref_fai)
			
			old_ref_fai = ref_fai
			ref_fai = os.path.join(args.output_dir, "tmp", ref_fai_base)
			
			os.rename(old_ref_fai, ref_fai)
		
		out_file = os.path.join(args.output_dir, "{}_RBV.txt".format(args.sample_id))
		
		#if HC needed


		
		
		out = open(out_file, 'w')
		print >>out, "#RBV: Read balance validator"
		print >>out, "#Command: " + str(sys.argv)
		

		
		
		
		
		
		#random readbal - perform once per run of RBV
		bzip_vcf = prep_vcf(args.vcf, args.output_dir)	#check if need to bgzip and tabix
		
		rand_readbal = random_readbal(args.variant_permutations,args.vcf,args.variant_quality_cutoff,args.calling_method)
		rand_readbal_array = np.array(rand_readbal).reshape(len(rand_readbal));
		rand_mean = np.mean(rand_readbal_array)
		
		print >>out, "#Random mean read balance = " + str(rand_mean)
		print >>out, "#"
		print >>out, "#CHR\tSTART\tSTOP\tpredicted type\th1 het snp number\th1 pvalue\th3 mean readbal\th3 t-test\th3 ks-test"
		
		if args.interval_file is not None:
			intervals = make_random_windows.intervals(args.interval_file)
			
		elif args.gap_bed is not None:
			gap_sites = make_random_windows.gaps(args.gap_bed,args.CNV_bed)
		else:
			gap_sites = {}
			zeros = [0,0]
			for gap_chr in range(1,23):
				gap_sites[str(gap_chr)]=[zeros]
		
		CNVs=open(args.CNV_bed).readlines()
		
		#For each CNV
		for i in range(len(CNVs)):
			CNV_chr,raw_CNV_start,raw_CNV_stop,CNV_type=CNVs[i].strip().split()
			CNV_start = int(raw_CNV_start)
			CNV_stop = int(raw_CNV_stop)
			
			out.write(str(CNV_chr) + "\t" + str(CNV_start) + "\t" + str(CNV_stop) + "\t" + str(CNV_type) + "\t")
			
			permutation = i + 1
			
			#DEL
			CNV_vcf_filename = str(args.sample_id) + "_CNV" + str(permutation) + ".vcf"
			CNV_vcf_file = os.path.join(args.output_dir, "tmp", CNV_vcf_filename)
			vcf_fetch(bzip_vcf,CNV_vcf_file, CNV_chr, CNV_start, CNV_stop)
			
			if args.interval_file is not None:				
				window_size = 0
				for coord in range(CNV_start,CNV_stop+1):
					for start,stop in intervals[CNV_chr]:
						if start<=coord<=stop:
							window_size+=1
				
				random_windows = make_random_windows.intervals_window(intervals, args.window_permutations, window_size)
			
			else:
				windows_size = CNV_stop - CNV_start
				random_windows = make_random_windows.gaps_windows(ref_fai, gap_sites, args.window_permutations, windows_size)
			
			print "[" + str(datetime.now()) + "] Random window generation - CNV"+str(permutation)+" complete."
			
			window_het_count = []
			for w in range(len(random_windows)):		
				window_chr,raw_window_start,raw_window_end,window_count= random_windows[w].split()
				window_start = int(raw_window_start)
				window_end = int(raw_window_end)
				
				window_vcf_filename = "rand_vcf_file_" + str(permutation) + "_" + str(w) + ".vcf"
				window_vcf_file = os.path.join(args.output_dir, "tmp", window_vcf_filename)
				vcf_fetch(bzip_vcf,window_vcf_file, window_chr, window_start, window_end)
				window_het_count.append(het_count(window_vcf_file))
				os.remove(window_vcf_file)
			
			if all([ v == 0 for v in window_het_count ]):
				out.write("0\tnan\t")
				print "[" + str(datetime.now()) + "] CNV" + str(permutation) + ": all random windows of this length contain no heterozygous SNPs. CNV too small or too few permuations"
				CNV_het_count = het_count(CNV_vcf_file)
				
			else:
				CNV_het_count,CNV_emp_pvalue = empirical_pvalue(window_het_count, CNV_vcf_file, args.window_permutations)
				
				out.write(str(CNV_het_count) + "\t" + str(CNV_emp_pvalue) + "\t")


			#DUP
			if CNV_het_count == 0:
				out.write("nan\tnan\tnan\n")
				print "[" + str(datetime.now()) + "] CNV" + str(permutation) + ": contains no heterozygous SNPs, unable to perform duplication analyses"
			else:
				CNV_vcf_lines=open(CNV_vcf_file).readlines()
				CNV_readbal = readbal(CNV_vcf_lines,args.variant_quality_cutoff,args.calling_method)
				CNV_readbal_array = np.array(CNV_readbal).reshape(len(CNV_readbal))
				
				CNV_mean = np.mean(CNV_readbal_array)
				CNV_ttest_stat,CNV_ttest_pvalue = stats.ttest_ind(CNV_readbal_array, rand_readbal_array, equal_var=False)
				CNV_kstest_stat,CNV_kstest_pvalue = stats.ks_2samp(CNV_readbal_array, rand_readbal_array)
					
				out.write(str(CNV_mean) + "\t" + str(CNV_ttest_pvalue) + "\t" + str(CNV_kstest_pvalue) + "\n")
			
			os.remove(CNV_vcf_file)
		
		rmtree(tmp)
			
		out.close()
		
		
		


if __name__=='__main__':
	main()






#Program initialisation message
#Program completion message


#total depth cutoff for vcf

#bed files are 0 based (first base is zero) therefore, CNV and gaps should be 1 based (.txt?)
