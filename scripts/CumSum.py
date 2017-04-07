import sys
import vcf
import numpy as np
import os

vcf_file = sys.argv[1]

#sample het SNPs
sample_chr = 1		
sample_start = 35101422
sample_stop = 35111976


vcf_reader = vcf.Reader(open(vcf_file, 'r'))
vcf_writer = vcf.Writer(open('sample_CNV.vcf', 'w'), vcf_reader)
for record in vcf_reader.fetch(sample_chr, sample_start, sample_stop):  
	vcf_writer.write_record(record)
vcf_writer.close()
		

		
sample_vcf = open('sample_CNV.vcf').readlines()		
		
sample_het_snp = 0

for i in range(len(sample_vcf)):
	line = sample_vcf[i]
	if line.startswith('#'):
		continue
	cols = line.strip().split()
	call = cols[9].split(':')
	vars_raw = call[0].replace('|', '/')
	vars = vars_raw.split('/')
	if vars[0] != vars[1]:
		sample_het_snp += 1
	else:
		continue
os.remove('sample_CNV.vcf')


print sample_het_snp



#cumulative density function of discrete SNPs in 10k samples
het_snps = np.loadtxt('PyVCF_test/10k_10kb.txt', dtype='float', delimiter='\n')

snp_sum = np.zeros(((max(het_snps)+1),), dtype=np.float)

for i in range(len(het_snps)):
	snp_count = het_snps[i]
	snp_sum[snp_count] += 1

#print snp_sum
	
cum_sum = 0
for n in range((sample_het_snp+1)):
	cum_sum += snp_sum[n]
	cum_sum_pvalue = cum_sum / 10000

#print cum_sum	
print cum_sum_pvalue