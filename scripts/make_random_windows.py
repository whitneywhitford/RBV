#Adapated from http://genometoolbox.blogspot.co.nz/2014/06/generate-random-genomic-positions.html accessed 10/03/2017


import sys
from copy import deepcopy
from operator import itemgetter
from random import choice,uniform,seed


#define functions

#WGS no intervals functions
def gaps(gapfile, chr_prefix):
	# Import chromosome gaps
	gaplines = open(gapfile, 'r')
	gap_raw = []
	for line in gaplines:
		if line.startswith("@"):
			continue
		else:
			gap_raw.append(line)

	gaps={}
	
	for i in range(len(gap_raw)):
		chr_raw,start,stop,strand,type=gap_raw[i].strip().split()
		if chr_prefix == True:
			chr = chr_raw.replace("chr", '', 1)
		else:
			chr = chr_raw
			
		gap_region=[int(start),int(stop)]
		if chr in gaps:
			gaps[chr].append(gap_region)
		else:
			gaps[chr]=[gap_region]
	
	#if no gaps in any of chr 1-23, then 0,0 included so can be called in random generation. If need to limit to chromosomes, then use intervals option
	zeros = [0,0]
	for gap_chr in range(1,23):
		if str(gap_chr) not in gaps:
			gaps[str(gap_chr)]=[zeros]

	return gaps

def gaps_total(gaplist,CNVs,chr_prefix):

	total_gaps = deepcopy(gaplist)
	
	for i in range(len(CNVs)):
		chr_raw,start,stop,strand,type=CNVs[i].strip().split()
		if chr_prefix == True:
			chr = chr_raw.replace("chr", '', 1)
		else:
			chr = chr_raw
		CNV_region=[int(start),int(stop)]
		if chr in total_gaps:
			total_gaps[chr].append(CNV_region)
		else:
			total_gaps[chr]=[CNV_region]
		
		for chr in total_gaps:
			total_gaps[chr] = sorted(total_gaps[chr])
		
	return total_gaps

#chr_length modified to fit Readbal
def gaps_chr_length(ref_fai, chr_prefix):
	chr_len={}
	chroms = str(list(range(1,23)))
	for i in range(len(ref_fai)):
		raw_chr,length,bite_index,bases_pl,bites_pl=ref_fai[i].strip().split()
		if chr_prefix == True:
			int_chr = raw_chr.replace("chr", '', 1)
		else:
			int_chr = raw_chr
		if int_chr in chroms:
			chr = str(raw_chr)
			chr_len[int_chr]=int(length)
		else:
			continue
	
	chr_lst=[]
	for chr in chr_len:							#account for different lengths of chromosomes so that random takes this into account
		chr_rep=[str(chr)] * int(round(chr_len[str(chr)]/10000000,0))
		for rep in chr_rep:
			chr_lst.append(rep)

	return chr_len, chr_lst


def gaps_rand_sites(no_of_samples,size,chr_len,chr_lst,total_gaps,chr_prefix):
	seed(None)
	count=1
	windows = []
	
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(0,chr_len[chr]-size))

		
		# Exclude inaccessible regions
		include="T"
		for start,stop in total_gaps[chr]:
			#window straddles the whole gap
			if (point <= start) and ((point+size)>= stop):
				include="F"
				break #if encounters any gap, window won't be included
			elif (start<=point<=stop or start<=(point+size)<=stop):
				include="F"
				break #if encounters any gap, window won't be included
		
		# Return points in accessible regions
		if include=="T":
			if chr_prefix == True:
				windows.append("chr" + str(chr) +" "+ str(point) +" "+ str(point+size) +" "+ str(count))
			else:
				windows.append(str(chr) +" "+ str(point) +" "+ str(point+size) +" "+ str(count))
			count+=1
				
	return windows
		
		
#intervals functions
def intervals(intervalfile,chr_prefix):
	intervallines = open(intervalfile, 'r')
	interval_raw = []
	for line in intervallines:
		if line.startswith("@"):
			continue
		else:
			interval_raw.append(line)
	
	intervals={}
	for i in range(len(interval_raw)):
		chr_raw,start,stop,strand,type=interval_raw[i].strip().split()
		if chr_prefix == True:
			chr = chr_raw.replace("chr", '', 1)
		else:
			chr = chr_raw
		interval_region=[int(start),int(stop)]
		if chr in intervals:
			intervals[chr].append(interval_region)
		else:
			intervals[chr]=[interval_region]
	
	return intervals	

	
def CNV_list(CNVlist,chr_prefix):
	CNVs={}
	for i in range(len(CNVlist)):
		chr_raw,start,stop,strand,type=CNVlist[i].strip().split()
		if chr_prefix == True:
			chr = chr_raw.replace("chr", '', 1)
		else:
			chr = chr_raw
		CNV_region=[int(start),int(stop)]
		if chr in CNVs:
			CNVs[chr].append(CNV_region)
		else:
			CNVs[chr]=[CNV_region]
	return CNVs


	#excludes CNVs from the intervals for random generation 
def random_intervals(interval_list, CNVs):
	unsorted_rand_intervals = deepcopy(interval_list)	#make copy of intervals that does not reference original
	
	for chr in CNVs:
		if  chr not in unsorted_rand_intervals.keys():
			continue
		for CNV_start, CNV_stop in CNVs[chr]:
			for interval in unsorted_rand_intervals[chr]:
				interval_start = interval[0]
				interval_stop = interval[1]
				if (CNV_start <= interval_start) and (CNV_stop >= interval_stop):
					unsorted_rand_intervals[chr].remove([interval_start, interval_stop])
				elif (interval_start <= CNV_start <= interval_stop) or (interval_start <= CNV_stop <= interval_stop): #if ranges overlap
					if CNV_start <= interval_start:
						interval_start = CNV_stop + 1
					else:
						interval[1] = CNV_start - 1
						if CNV_stop <= interval_stop:
							unsorted_rand_intervals[chr].append([CNV_stop+1, interval_stop])
						
	rand_intervals = {}
	for chr in unsorted_rand_intervals:
		rand_intervals[chr] = sorted(unsorted_rand_intervals[chr])
	
	return rand_intervals

	
def intervals_chr_length(interval_list):	
	chr_interval_len={}
	chroms = str(list(range(1,23)))
	for chr in interval_list:	
		if chr in chroms:
			chr_interval_len[chr] = interval_list[chr][-1][1]
		else:
			continue
	
	chr_lst=[]
	for chr in chr_interval_len:							#account for different lengths of chromosomes so that random takes this into account
		chr_rep=[str(chr)] * int(round(chr_interval_len[str(chr)]/10000,0)) #10000 needs to be smaller than the smallest total chromosome exome intervals 
		for rep in chr_rep:
			chr_lst.append(rep)
	
	return chr_interval_len, chr_lst
	
def intervals_rand_sites(no_of_samples,window_size,chr_interval_len,chr_lst,interval_list,chr_prefix):
	seed(None)
	count=1
	windows = []
	
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(1,chr_interval_len[chr]-window_size))
		past_start = False
		window_count = 0
		
		for start,stop in interval_list[chr]:
			if (past_start == False) and (point < start):
				break
				
			elif start<=point<=stop:	#start found
				window_count = (stop - point) + 1
				past_start = True

			elif past_start: 
				window_count += (stop - start) + 1
				
			if window_count >= window_size:
				window_end = stop - (window_count - window_size)
				if chr_prefix == True:
					windows.append("chr" + str(chr) +" "+ str(point) +" "+ str(window_end) +" "+ str(count))
				else:
					windows.append(str(chr) +" "+ str(point) +" "+ str(window_end) +" "+ str(count))
				count+=1
				break
				
	return windows

