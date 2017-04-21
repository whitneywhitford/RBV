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



gapfile = sys.argv[1]
CNV_file = sys.argv[2]
variant_permutations = int(sys.argv[3])
vcf_file = sys.argv[4]
qualCutoff = int(sys.argv[5])
out_file = sys.argv[6]
window_permutations = int(sys.argv[7])
ref_file = sys.argv[8]








#def parse_args():
"""
Parse command-line arguments

Returns
-------

Parser argument namespace
"""
"""parser = argparse.ArgumentParser(description="RBV")

parser.add_argument(
	"--ref", required=True,
	help="REQUIRED. Path to reference sequence (including file name).")

parser.add_argument(
	"--output_dir", "-o", default="RBV",
	help="REQUIRED. Output directory. XYalign will create a directory "
	"structure within this directory")

parser.add_argument(
	"--sample_id", "-id", default="sample",
	help="Name/ID of sample - for use in plot titles and file naming. "
	"Default is sample")"""

"""parser.add_argument(
	"--cpus", type=int, default=1,
	help="Number of cores/threads to use. Default is 1")

parser.add_argument(
	"--xmx", default="None",
	help="Memory to be provided to java programs via -Xmx.  E.g., use the "
	"flag '--xmx Xmx4g' to pass '-Xmx4g' as a flag when running java "
	"programs (currently just repair.sh). Default is 'None' (i.e., nothing "
	"provided on the command line), which will allow repair.sh to "
	"automatically allocate memory.")

parser.add_argument(
	"--single_end", action="store_true", default=False,
	help="Include flag if reads are single-end and NOT paired-end.")"""

"""# Logging options
parser.add_argument(
	"--logfile", default=None,
	help="Name of logfile.  Will overwrite if exists.  Default is "
	"sample_xyalign.log")

parser.add_argument(
	"--reporting_level", default='INFO',
	choices=["DEBUG", "INFO", "ERROR", "CRITICAL"],
	help="Set level of messages printed to console. Default is 'INFO'. "
	"Choose from (in decreasing amount of reporting) DEBUG, INFO, ERROR "
	"or CRITICAL")

# Options for turning on/off parts of the pipeline
parser.add_argument(
	"--prepare_vcf", action="store_true", default=False,
	help="This flag will result in the preparation of a vcf file by "
	"HaploTypeCaller. A HaploTypeCaller generated vcf file can be parsed "
	"using the --vcf flag.")

# Variant Calling Flags
parser.add_argument(
	"--variant_quality_cutoff", "-vqc", type=int, default=20,
	help="Consider all SNPs with a quality greater than or "
	"equal to this value. Default is 20.")

# Bam Analysis and Characterize Sex Chrom Flags
parser.add_argument(
	"--window_size", "-w", default=None,
	help="Window size (integer) for sliding window calculations. Default "
	"is 50000.  Default is None. If set to None, will use targets "
	"provided using --target_bed.")

parser.add_argument(
	"--CNV_bed", default=None,
	help="Bed file containing targets to use in CNV analyses "
	"Must be typical bed format, 0-based indexing, with the first three "
	"columns the chromosome name, start coordinate, stop coordinate.")

parser.add_argument(
	"--gap_bed", default=None,
	help="Bed file containing gaps in the reference to mask for random "
	"generation. Must be typical bed format, 0-based indexing, with the "
	"first three columns the chromosome name, start coordinate, stop "
	"coordinate.")
	
parser.add_argument(
	"--variant_permutations", type=int, default=10000,
	help="Number of permutations to use for heterozygous read balance "
	"analyses. Default is 10000")

parser.add_argument(
	"--window_permutations", type=int, default=10000,
	help="Number of permutations to use for heterozygous SNP number "
	"analyses. Default is 10000")

# Return arguments namespace
return args"""


	
	
	
	
	
	
	

#make_random_vars
def gaps(gapfile,CNVfile):
	# Import chromosome gaps
	gap_raw=open(gapfile).readlines()
	CNV_raw=open(CNVfile).readlines()
	gaps={}
	for i in range(len(gap_raw)):
		chr,start,stop,type=gap_raw[i].strip().split()
		gap_region=[int(start),int(stop)]
		if chr in gaps:
			gaps[chr].append(gap_region)
		else:
			gaps[chr]=[gap_region]
	
	for i in range(len(CNV_raw)):
		chr,start,stop,type=gap_raw[i].strip().split()
		CNV_region=[int(start),int(stop)]
		if chr in gaps:
			gaps[chr].append(CNV_region)
		else:
			gaps[chr]=[CNV_region]
	
	return gaps
	
def rand_het_sites(no_of_samples,vcf_file,gaps):	
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
		# Exclude inaccessible regions
		include="T"
		for start,stop in gaps[chr]:
			if start<=coord<=stop:
				include="F"
			call = cols[9].split(':')
			vars_raw = call[0].replace('|', '/')
			vars = vars_raw.split('/')
			if vars[0] == vars[1]:
				include="F"
		
		# Return points in accessible regions
		if include=="T":
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
	window_het_snp_sum = np.zeros(((max(window_het_snps)+1),), dtype=np.float)
	
	for i in range(len(window_het_snps)):
		window_het_snp_count = window_het_snps[i]
		window_het_snp_sum[window_het_snp_count] += 1
	
	
	sample_cum_sum = 0
	sample_het_snps = het_count(vcf_outfile)
	
	for n in range((sample_het_snps+1)):
		sample_cum_sum += window_het_snp_sum[n]
	
	sample_cum_sum_pvalue = sample_cum_sum / num_window_samples
	
	return sample_het_snps,sample_cum_sum_pvalue

def readbal(variant_lines,qualCutoff):
	
	readBalance = []
	
	for i in range(len(variant_lines)):
		if variant_lines[i].startswith('#'):
			continue
		cols = variant_lines[i].strip('\n').split('\t')
		call = cols[9].split(':')
		vars_raw = call[0].replace('|', '/')
		vars = vars_raw.split('/')
		if vars[0] != vars[1]:
			pos = int(cols[1])
			qual = float(cols[5])
			if qual < qualCutoff:
				continue
			GT = cols[9].split(':')[0]
			if GT == './.' or GT == '.|.':
				continue
			allele1 = cols[9].split(':')[1].split(',')[0]
			allele2 = cols[9].split(':')[1].split(',')[1]
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

	
def random_readbal(no_of_samples,vcf_file,gaps,qualCutoff):
	
	variant_lines = rand_het_sites(no_of_samples,vcf_file,gaps)
	
	random_readbal = readbal(variant_lines,qualCutoff)
	
	return random_readbal

def prep_vcf(vcf_file):
	working_vcf = "working_vcf.vcf"
	copyfile(vcf_file, working_vcf)
	os.system("bgzip " + working_vcf)
	bgzip_working_vcf = str(working_vcf) + ".gz"
	os.system("tabix -p vcf " + bgzip_working_vcf)
	
	return bgzip_working_vcf


	
	
	
	
"""np.mean(random_readbal)	#once

#Each CNV
np.mean(sample_readbal)
stats.ttest_ind(sample_readbal, random_readbal, equal_var=False)
stats.ks_2samp(sample_readbal, random_readbal)"""


def main(ref_file, gapfile, CNV_file, variant_permutations, vcf_file, qualCutoff, window_permutations, out_file):

	ref_fai = str(ref_file) + ".fai"	#check if all fai are just fasta.fai/fa.fai
	
	#if HC needed


	
	
	out = open(out_file, 'w')
	print >>out, "#RBV: Read balance validator"
	print >>out, "#Command: " + str(sys.argv)
	
	
	
	
	
	
	
	#random readbal - perform once per run of RBV
	bzip_vcf = prep_vcf(vcf_file)	#check if need to bgzip and tabix
	
	gap_sites = gaps(gapfile,CNV_file)
	rand_readbal = random_readbal(variant_permutations,vcf_file,gap_sites,qualCutoff)
	rand_readbal_array = np.array(rand_readbal).reshape(len(rand_readbal));
	rand_mean = np.mean(rand_readbal_array)
	
	print >>out, "#Random mean read balance = " + str(rand_mean)
	print >>out, "#"
	print >>out, "#CHR\tSTART\tSTOP\th1 het snp number\th1 pvalue\th3 mean readbal\th3 t-test\th3 ks-test"
	
	CNVs=open(CNV_file).readlines()
	
	#For each CNV
	for i in range(len(CNVs)):
		CNV_chr,raw_CNV_start,raw_CNV_stop,CNV_type=CNVs[i].strip().split()
		CNV_start = int(raw_CNV_start)
		CNV_stop = int(raw_CNV_stop)
		
		out.write(str(CNV_chr) + "\t" + str(CNV_start) + "\t" + str(CNV_stop) + "\t")
		
		#DEL
		CNV_vcf_file = "CNV_vcf_file_" + str(i) + ".vcf"
		vcf_fetch(bzip_vcf,CNV_vcf_file, CNV_chr, CNV_start, CNV_stop)
		
		windows_size = CNV_stop - CNV_start
		
		permutation = i + 1
		
		random_windows = make_random_windows.main(ref_fai, gap_sites, window_permutations, windows_size, permutation)
		print "Random window generation - CNV"+str(permutation)+" complete."
		
		window_het_count = []
		for w in range(len(random_windows)):		
			window_chr,raw_window_start,raw_window_end,window_count= random_windows[w].split()
			window_start = int(raw_window_start)
			window_end = int(raw_window_end)
			
			window_vcf_file = "rand_vcf_file_" + str(i) + "_" + str(w) + ".vcf"
			vcf_fetch(bzip_vcf,window_vcf_file, window_chr, window_start, window_end)
			window_het_count.append(het_count(window_vcf_file))
			os.remove(window_vcf_file)
		
		if all([ v == 0 for v in window_het_count ]):
			out.write("0\tnan\t")
			print "CNV" + str(permutation) + ": all random windows of this length contain no heterozygous SNPs. CNV too small or too few permuations"
			CNV_het_count = het_count(CNV_vcf_file)
			
		else:
			CNV_het_count,CNV_emp_pvalue = empirical_pvalue(window_het_count, CNV_vcf_file, window_permutations)
			
			out.write(str(CNV_het_count) + "\t" + str(CNV_emp_pvalue) + "\t")


		#DUP
		if CNV_het_count == 0:
			out.write("nan\tnan\tnan\n")
			print "CNV" + str(permutation) + ": contains no heterozygous SNPs, unable to perform duplication analyses"
		else:
			CNV_vcf_lines=open(CNV_vcf_file).readlines()
			CNV_readbal = readbal(CNV_vcf_lines,qualCutoff)
			CNV_readbal_array = np.array(CNV_readbal).reshape(len(CNV_readbal))
			
			CNV_mean = np.mean(CNV_readbal_array)
			CNV_ttest_stat,CNV_ttest_pvalue = stats.ttest_ind(CNV_readbal_array, rand_readbal_array, equal_var=False)
			CNV_kstest_stat,CNV_kstest_pvalue = stats.ks_2samp(CNV_readbal_array, rand_readbal_array)
			
			os.remove(CNV_vcf_file)
			
			out.write(str(CNV_mean) + "\t" + str(CNV_ttest_pvalue) + "\t" + str(CNV_kstest_pvalue) + "\n")
		
		
		
	out.close()
		
		
		


if __name__=='__main__':
	main(ref_file, gapfile, CNV_file, variant_permutations, vcf_file, qualCutoff, window_permutations, out_file)





#time stamp with output messages


#check for .fai with reference
#total depth cutoff for vcf, check that qual cutoff refers to correct column
