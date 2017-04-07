import sys
import numpy as np

vcf_file = sys.argv[1]
out = open(sys.argv[2], 'a')
window_file = sys.argv[3]
window_size = int(sys.argv[4])
#chrom = sys.argv[4]
#chromlen = int(sys.argv[5])

#vcf = open(vcf_file, 'r')

window_in = open(window_file).readlines()
windows = {}
het_count = {}

for i in range(len(window_in)):
	raw_chr,window_start_raw,window_stop_raw,col4=window_in[i].strip().split()
	window_chr = int(raw_chr)
	window_start = int(window_start_raw)
	window_stop = int(window_stop_raw)
	if window_chr in windows:
		windows[window_chr].append(window_start)
	else:
		windows[window_chr]=[window_start]
		
for chrom in windows:
	het_count[chrom] = [0] * len(windows[chrom])
	for window in range(len(windows[chrom])):
		working_chrom = windows[chrom]
		#print "chrom ",chrom
		#print "window number",window
		#print "window ",working_chrom[window]
		vcf = open(vcf_file).readlines()
		for i in range(len(vcf)):
			line = vcf[i]
			if line.startswith('#'):
				#print "comment"
				continue
			cols = line.strip('\n').split('\t')
			sample_chrom = int(cols[0])
			sample_coord = int(cols[1])
			#print "sample_chrom ",sample_chrom
			#print "sample_coord ",sample_coord
			if sample_chrom < chrom:
				#print "higher chrom"
				continue
			if sample_chrom == chrom:
				if sample_coord < working_chrom[window]:
					#print "lower"
					continue
				if sample_coord >= working_chrom[window] and sample_coord < (working_chrom[window]+window_size):
					#print "between"
					call = cols[9].split(':')
					vars_raw = call[0].replace('|', '/')
					vars = vars_raw.split('/')
					het_window = het_count[chrom]
					#print het_window[window]
					if vars[0] != vars[1]:
							het_window[window] += 1
							#print "het"
					else:
						continue
				else: #sample_coord >= (working_chrom[window]+window_size):
					#print "higher"
					break
			else:
				#print "lower chrom"
				break

				
print het_count
				
print "completed"

#vcf.close();
out.close();