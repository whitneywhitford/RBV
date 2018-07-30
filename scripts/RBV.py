# RBV main program

import argparse
import os
import sys
from scipy import stats
import vcf
import numpy as np
import make_random_windows
from shutil import copyfile
from shutil import rmtree
from datetime import datetime


def parse_args():
	"""Parse command-line arguments"""

	parser = argparse.ArgumentParser(description="RBV")

	parser.add_argument(
		"--ref", default=None,
		help="REQUIRED. Path to reference sequence (including file name).")
	
	parser.add_argument(
		"--CNV_file", required=True,
		help="REQUIRED. Picard-style interval_list containing targets to use in CNV analyses."
		"Must be typical interval_list format: 1-based indexing, with the six"
		"columns being the chromosome name, start coordinate, stop coordinate,"
		"strand, and predicted CNV type.")

	parser.add_argument(
		"--gap_file", default=None,
		help="Picard-style interval_list containing gaps in the reference to"
		"mask for random generation. Must be typical interval_file format:"
		"1-based, indexing, with the six columns being the chromosome name,"
		"start coordinate, stop coordinate, strand, and type.")
		
	parser.add_argument(
		"--vcf", required=True,
		help="REQUIRED. VCF file containing the variants for RBV analyses."
		"VCF must be single individual vcf file but can be derived from joint"
		"calling pipeline.")
	
	parser.add_argument(
		"--output_dir", "-o", default="RBV",
		help="Output directory. RBV will create a temporary "
		"directory and output file within this directory.")

	parser.add_argument(
		"--sample_id", "-id", default="sample",
		help="Name/ID of sample - for use in file naming. Default is sample")

	# Variant Calling Flags
	parser.add_argument(
		"--variant_quality_cutoff", "-vqc", type=int, default=20,
		help="Consider all SNPs with a quality greater than or "
		"equal to this value. Default is 20.")
		
	parser.add_argument(
		"--read_depth_cutoff", "-rdc", type=int, default=10,
		help="Consider all SNPs with a read depth greater than or "
		"equal to this value. Default is 10.")
		
	parser.add_argument(
		"--readbal_cutoff", "-rbc", type=float, default=0.65,
		help="For deletion analyses, consider all heterozygous SNPs with a"
		"read balance less than this value. Default is 0.65.")
		
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
		help="Picard-style interval_list containing interval coordinates"
		"used for variant calling. Must be typical interval_list format: 1-based"
		"indexing, with the six columns being the chromosome name, start"
		"coordinate, stop coordinate, strand, and type."
		"REQUIRED for if using WES seq_type.")
	
	parser.add_argument(
		"--calling_method", required=True,
		help="REQUIRED. Variant calling method used for VCF generation."
		"Options: haplotypecaller, samtools, freebayes, platypus.")
	
	args = parser.parse_args()
	
	if args.seq_type not in ["WGS", "WES"]:
		sys.exit(
			"ERROR: Sequencing type must be WGS or WES.")
	
	if args.calling_method not in ["haplotypecaller", "samtools", "freebayes", "platypus"]:
		sys.exit(
			"ERROR: Calling method must be haplotypecaller, samtools, freebayes, or platypus.")
	
	if args.seq_type == "WES" and args.interval_file == None:
		sys.exit(
			"ERROR: Interval file required for Whole Exome Sequencing"
			"RBV analysis")
	
	if args.ref == None and args.interval_file == None:
		sys.exit(
			"ERROR: Reference file required if a interval file is not supplied"
			)
	
	if os.path.isdir(args.output_dir):
		sys.exit(
			"ERROR: Output directory already exists. Please enter new Output"
			"directory.") 
	
	# Return arguments namespace
	return args


#make_random_vars
	
def rand_het_sites(no_of_samples,vcf_file,chr_prefix):	
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

	# Choose variant
	var_lines = []
	while count<=no_of_samples:
		line_number = randrange((comment_count - 1), len(vcf))
		
		variant_line = vcf[line_number]
		cols = variant_line.strip('\n').split('\t')
		chr = cols[0]
		coord = cols[1]
		
		str_chroms = []
		for num in range(1,23):
			if chr_prefix == True:
				str_chroms.append("chr"+str(num))
			else:
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

	
def het_count(vcf_file,qualCutoff,readdepthCutoff,readbalCutoff,calling_method):
	del_het_count = 0
	dup_het_count = 0
	vcf = open(vcf_file).readlines()
	for i in range(len(vcf)):
		line = vcf[i]
		if line.startswith('#'):
			continue
		cols = line.strip().split()
		GT = cols[9].split(':')[0]
		vars_raw = GT.replace('|', '/')
		vars = vars_raw.split('/')
		if vars[0] != vars[1]:
			qual = float(cols[5])
			if qual < qualCutoff:
				continue
			if GT == './.' or GT == '.|.':
				continue
			
			allele_reads = globals()[calling_method + "_vcf"]
			allele1, allele2 = allele_reads(cols)
			
			readdepth = int(allele1) + int(allele2)
			if readdepth < readdepthCutoff:
				continue
			
			if ',' in allele1 or ',' in allele2:
				continue
				
			dup_het_count += 1
			
			if float(allele1) + float(allele2) == 0:
				continue
			if float(allele1) >= float(allele2):
				ReadRatio = float(allele1)/float(readdepth)
			else:
				ReadRatio = float(allele2)/float(readdepth)
			if ReadRatio > float(readbalCutoff):
				continue

			del_het_count += 1
			
		else:
			continue
	
	return del_het_count,dup_het_count

	
def empirical_pvalue(window_het_snps, vcf_outfile, qualCutoff ,readdepthCutoff, readbalCutoff, calling_method, num_window_samples):
	sample_cum_sum = 0
	sample_del_het_snps,sample_dup_het_snps = het_count(vcf_outfile,qualCutoff,readdepthCutoff,readbalCutoff,calling_method)
	
	if sample_del_het_snps > max(window_het_snps):
		sample_cum_sum_pvalue = 1.000
	
	else:
		window_het_snp_sum = np.zeros(((max(window_het_snps)+1),), dtype=np.float)
		
		
		
		for i in range(len(window_het_snps)):
			window_het_snp_count = window_het_snps[i]
			window_het_snp_sum[window_het_snp_count] += 1
		
		for n in range((sample_del_het_snps+1)):
			sample_cum_sum += window_het_snp_sum[n]
		
		sample_cum_sum_pvalue = sample_cum_sum / num_window_samples
		
	return sample_del_het_snps,sample_dup_het_snps,sample_cum_sum_pvalue

def haplotypecaller_vcf(variant_cols):
	format = variant_cols[8].split(':')
	alleles_index = format.index("AD")
	allele1 = variant_cols[9].split(':')[alleles_index].split(',')[0]
	allele2 = variant_cols[9].split(':')[alleles_index].split(',')[1]
	
	return allele1, allele2
	
def samtools_vcf(variant_cols):
	info = variant_cols[7].replace('=',';').split(';')
	reads_index = info.index("DP4") + 1
	alleles = info[reads_index].split(',')
	
	allele1F = int(alleles[0])
	allele1R = int(alleles[1])
	allele1 = allele1F + allele1R
	allele2F = int(alleles[2])
	allele2R = int(alleles[3])
	allele2 = allele2F + allele2R
	
	return str(allele1), str(allele2)
	
def freebayes_vcf(variant_cols):
	format = variant_cols[8].split(':')
	allele1_index = format.index("RO")
	allele1 = int(variant_cols[9].split(':')[allele1_index])
	total_index = format.index("DP")
	total = int(variant_cols[9].split(':')[total_index])
	allele2 = total - allele1
	
	return str(allele1), str(allele2)
	
def platypus_vcf(variant_cols):
	info = variant_cols[7].replace(';','=').split('=')
	allele1_index = info.index("TR") + 1
	allele1_all = info[allele1_index]
	allele1 = int(allele1_all.split(',')[0])
	total_index = info.index("TC") + 1
	total = int(info[total_index])
	allele2 = total - allele1
	
	return str(allele1), str(allele2)
	
	
def readbal(variant_lines,qualCutoff,readdepthCutoff,calling_method):
	
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
			
			readdepth = int(allele1) + int(allele2)
			
			if readdepth < readdepthCutoff:
				continue
			
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

	
def random_readbal(no_of_samples,vcf_file,qualCutoff,readdepthCutoff,calling_method,chr_prefix):
	
	variant_lines = rand_het_sites(no_of_samples,vcf_file,chr_prefix)
	
	random_readbal = readbal(variant_lines,qualCutoff,readdepthCutoff,calling_method)
	
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
	
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
	sys.stdout.write("[" + str(datetime.now()) + "] Read Balance Validator (RBV) v1.0.0\n")
	sys.stdout.write("[" + str(datetime.now()) + "] Copyright (c) 2017 Whitney Whitford\n")
	sys.stdout.write("[" + str(datetime.now()) + "] For support and documentation go to https://github.com/whitneywhitford/RBV\n")
	sys.stdout.write("[" + str(datetime.now()) + "] Program Args: " + str(sys.argv) + "\n")
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
	
	warning_messages = []
	warning_count = 0
	
	os.mkdir(args.output_dir)
	tmp = str(args.output_dir) + "/tmp"
	os.mkdir(tmp)

	if args.ref is not None:
		ref_fai = str(args.ref) + ".fai"	#check if all fai are just fasta.fai/fa.fai
		
		if not os.path.isfile(ref_fai):
			os.system("samtools faidx  " + args.ref)
			ref_fai_base = os.path.basename(ref_fai)
			
			old_ref_fai = ref_fai
			ref_fai = os.path.join(args.output_dir, "tmp", ref_fai_base)
			
			os.rename(old_ref_fai, ref_fai)
	
	out_file = os.path.join(args.output_dir, "{}_RBV.txt".format(args.sample_id))
	

	out = open(out_file, 'w')
	out.write("#RBV: Read balance validator\n")
	out.write("#Command: " + str(sys.argv) + "\n")

	
	#random readbal - perform once per run of RBV
	bzip_vcf = prep_vcf(args.vcf, args.output_dir)	#check if need to bgzip and tabix
	
	CNVlines = open(args.CNV_file, 'r')
	CNVs = []
	for line in CNVlines:
		if line.startswith("@"):
			continue
		else:
			CNVs.append(line)
	
	chr_prefix = False
	
	if args.interval_file is None:
		ref_file = open(ref_fai).readlines()
		if ref_file[0].startswith("chr"):
			chr_prefix = True
	else:
		if CNVs[0].startswith("chr"):
			chr_prefix = True
	
	rand_readbal = random_readbal(args.variant_permutations,args.vcf,args.variant_quality_cutoff,args.read_depth_cutoff,args.calling_method,chr_prefix)
	rand_readbal_array = np.array(rand_readbal).reshape(len(rand_readbal));
	rand_mean = np.mean(rand_readbal_array)
	rand_stddev  = np.std(rand_readbal_array)
	
	out.write("#Random mean read balance = " + str(rand_mean) + "\tStd dev = " + str(rand_stddev) + "\n")
	out.write("#\n")
	out.write("#CHR\tSTART\tSTOP\tpredicted type\th1 het snp number\th1 p-value\th3 het snp number\th3 mean readbal\th3 p-value\n")
	
	if args.interval_file is not None:
		intervals = make_random_windows.intervals(args.interval_file,chr_prefix)
		CNV_list = make_random_windows.CNV_list(CNVs,chr_prefix)
		total_intervals = make_random_windows.random_intervals(intervals, CNV_list)
		chr_interval_len,chr_lst=make_random_windows.intervals_chr_length(total_intervals)
		
	elif args.gap_file is not None:
		gap_sites = make_random_windows.gaps(args.gap_file,chr_prefix)
		total_gap_sites = make_random_windows.gaps_total(gap_sites,CNVs,chr_prefix)
		chr_gaps_len,chr_lst=make_random_windows.gaps_chr_length(ref_file, chr_prefix)
		
	else:
		gap_sites = {}
		zeros = [0,0]
		for gap_chr in range(1,23):
			gap_sites[str(gap_chr)]=[zeros]
			
		total_gap_sites = make_random_windows.gaps_total(gap_sites,CNVs,chr_prefix)
		chr_gaps_len,chr_lst=make_random_windows.gaps_chr_length(ref_file, chr_prefix)
	
	#For each CNV
	for i in range(len(CNVs)):
		CNV_chr,raw_CNV_start,raw_CNV_stop,strand,CNV_type=CNVs[i].strip().split()
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
				if chr_prefix == True:
					intervals_CNV_chr = CNV_chr.replace("chr", '', 1)
				else:
					intervals_CNV_chr = CNV_chr
				for start,stop in intervals[intervals_CNV_chr]:
					if (start<=coord<=stop):
						window_size+=1
						break
			
			if window_size > 0:
				random_windows = make_random_windows.intervals_rand_sites(args.window_permutations,window_size,chr_interval_len,chr_lst,total_intervals,chr_prefix)
			
			else:
				zero_window = str("none") +" "+ str(0) +" "+ str(0) +" "+ str(0)
				random_windows = [zero_window]
		
		else:
			window_size = CNV_stop - CNV_start + 1
			
			for coord in range(CNV_start,CNV_stop+1):
				if chr_prefix == True:
					gaps_CNV_chr = CNV_chr.replace("chr", '', 1)
				else:
					gaps_CNV_chr = CNV_chr
				
				for start,stop in gap_sites[gaps_CNV_chr]:
					if start<=coord<=stop:
						window_size-=1
						break #if encounters any gap, nucleotide won't be included in window size
			
			if window_size > 0:
				random_windows = make_random_windows.gaps_rand_sites(args.window_permutations,window_size,chr_gaps_len,chr_lst,total_gap_sites,chr_prefix)
			
			else:
				zero_window = str("none") +" "+ str(0) +" "+ str(0) +" "+ str(0)
				random_windows = [zero_window]
						
		sys.stdout.write("[" + str(datetime.now()) + "] Random window generation - CNV"+str(permutation)+" complete.\n")
		
		window_het_count = []
		for w in range(len(random_windows)):		
			window_chr,raw_window_start,raw_window_end,window_count= random_windows[w].split()
			window_start = int(raw_window_start)
			window_end = int(raw_window_end)
			
			window_vcf_filename = "rand_vcf_file_" + str(permutation) + "_" + str(w) + ".vcf"
			window_vcf_file = os.path.join(args.output_dir, "tmp", window_vcf_filename)
			vcf_fetch(bzip_vcf,window_vcf_file, window_chr, window_start, window_end)
			window_het_count.append(het_count(window_vcf_file,args.variant_quality_cutoff,args.read_depth_cutoff,args.readbal_cutoff,args.calling_method)[0])
			os.remove(window_vcf_file)
		
		if all([ v == 0 for v in window_het_count ]):
			
			no_rand_SNP_warning = "WARNING: CNV" + str(permutation) + " all random windows of this length contain no heterozygous SNPs. CNV too small or too few permuations"
			
			sys.stderr.write("[" + str(datetime.now()) + "] " + no_rand_SNP_warning + "\n")
			
			warning_messages.append(no_rand_SNP_warning)
			warning_count += 1
			
			deleted_het_count,duplicated_het_count = het_count(CNV_vcf_file,args.variant_quality_cutoff,args.read_depth_cutoff,args.readbal_cutoff,args.calling_method)
			
			out.write(str(deleted_het_count) + "\tnan\t")
			
		else:
			deleted_het_count,duplicated_het_count,CNV_emp_pvalue = empirical_pvalue(window_het_count, CNV_vcf_file, args.variant_quality_cutoff, args.read_depth_cutoff, args.readbal_cutoff, args.calling_method, args.window_permutations)
			
			out.write(str(deleted_het_count) + "\t" + str(CNV_emp_pvalue) + "\t")


		#DUP
		out.write(str(duplicated_het_count) + "\t")	
		
		if duplicated_het_count == 0:
			out.write("nan\tnan\n")
			
			no_CNV_SNP_warning = "WARNING: CNV" + str(permutation) + " contains no heterozygous SNPs, unable to perform duplication analysis"
			
			sys.stderr.write("[" + str(datetime.now()) + "] " + no_CNV_SNP_warning + "\n")
			
			warning_messages.append(no_CNV_SNP_warning)
			warning_count += 1
		else:
			CNV_vcf_lines=open(CNV_vcf_file).readlines()
			CNV_readbal = readbal(CNV_vcf_lines,args.variant_quality_cutoff,args.read_depth_cutoff,args.calling_method)
			CNV_readbal_array = np.array(CNV_readbal).reshape(len(CNV_readbal))
			
			CNV_readbal_array_length = len(CNV_readbal_array)
			
			CNV_mean = np.mean(CNV_readbal_array)
			out.write(str(CNV_mean) + "\t")			
			
			CNV_kstest_stat,CNV_kstest_pvalue = stats.ks_2samp(CNV_readbal_array, rand_readbal_array)
				
			out.write(str(CNV_kstest_pvalue) + "\n")
		
		os.remove(CNV_vcf_file)
	
	rmtree(tmp)
		
	out.close()
	
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
	if warning_count == 1:
		sys.stdout.write("[" + str(datetime.now()) + "] Done. There was " + str(warning_count) + " WARNING message.\n")
	else:
		sys.stdout.write("[" + str(datetime.now()) + "] Done. There were " + str(warning_count) + " WARNING messages.\n")
	sys.stdout.write("------------------------------------------------------------------------------------------\n")
		


if __name__=='__main__':
	main()
