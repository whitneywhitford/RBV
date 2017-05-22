#Adapated from http://genometoolbox.blogspot.co.nz/2014/06/generate-random-genomic-positions.html accessed 10/03/2017


import sys
from copy import deepcopy
from operator import itemgetter



#define functions

#WGS no intervals functions
def gaps(gapfile):
	# Import chromosome gaps
	gap_raw=open(gapfile).readlines()
	gaps={}
	for i in range(len(gap_raw)):
		chr,start,stop,type=gap_raw[i].strip().split()
		gap_region=[int(start),int(stop)]
		if chr in gaps:
			gaps[chr].append(gap_region)
		else:
			gaps[chr]=[gap_region]
	
	return gaps

def gaps_total(gaplist,CNVs):
	for i in range(len(CNVs)):
		chr,start,stop,type=CNVs[i].strip().split()
		CNV_region=[int(start),int(stop)]
		if chr in gaplist:
			gaplist[chr].append(CNV_region)
		else:
			gaplist[chr]=[CNV_region]
	
	return gaplist

#chr_length modified to fit Readbal
def gaps_chr_length(ref_fai):	
	reffile = open(ref_fai).readlines()
	chr_len={}
	chroms = str(list(range(1,23)))
	for i in range(len(reffile)):
		raw_chr,length,bite_index,bases_pl,bites_pl=reffile[i].strip().split()
		if raw_chr in chroms:
			chr = str(raw_chr)
			chr_len[chr]=int(length)
		else:
			continue
	
	chr_lst=[]
	for chr in chr_len:							#account for different lengths of chromosomes so that random takes this into account
		chr_rep=[str(chr)] * int(round(chr_len[str(chr)]/10000000,0))
		for rep in chr_rep:
			chr_lst.append(rep)
	
	return chr_len, chr_lst


def gaps_rand_sites(no_of_samples,size,chr_len,chr_lst,gaps):
	from random import choice,uniform,seed
	seed(None)
	count=1
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(0,chr_len[chr]-size))

		
		# Exclude inaccessible regions
		include="T"
		for start,stop in gaps[chr]:
			if (point <= start) and ((point+size)>= stop):
				include="F"
				break #if encounters any gap, window won't be included
			if (start<=point<=stop or start<=(point+size)<=stop):
				include="F"
				break #if encounters any gap, window won't be included
		
		# Return points in accessible regions
		if include=="T":
			yield count,int(chr),point,(point+size)
			count+=1


def gaps_make_random(no_of_samples,size,chr_len,chr_lst,gaps):
	from operator import itemgetter
	
	a=gaps_rand_sites(no_of_samples,size,chr_len,chr_lst,gaps)
	b=sort_coords(a)
	
	return b

		
		
#intervals functions
def intervals(intervalfile):
	interval_raw=open(intervalfile).readlines()
	intervals={}
	for i in range(len(interval_raw)):
		chr,start,stop,type=interval_raw[i].strip().split()
		interval_region=[int(start),int(stop)]
		if chr in intervals:
			intervals[chr].append(interval_region)
		else:
			intervals[chr]=[interval_region]
	
	return intervals	

	
def CNV_list(CNVfile):
	CNV_raw=open(CNVfile).readlines()
	CNVs={}
	for i in range(len(CNV_raw)):
		chr,start,stop,type=CNV_raw[i].strip().split()
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
	
def intervals_rand_sites(no_of_samples,window_size,chr_interval_len,chr_lst,interval_list):
	from random import choice,uniform,seed
	seed(None)
	count=1
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(1,chr_interval_len[chr]-window_size))
		past_start = False
		window_count = 0
		
		for start,stop in interval_list[chr]:
			if start<=point<=stop:	#start found
				window_count = (stop - point) + 1
				past_start = True

			elif past_start: 
				window_count += (stop - start) + 1
				
			if window_count >= window_size:
				window_end = stop - (window_count - window_size)
				yield count,int(chr),point,window_end
				count+=1
				break
					
def intervals_make_random(no_of_samples,size,chr_interval_len,chr_lst,interval_list):
	from operator import itemgetter
	
	a=intervals_rand_sites(no_of_samples,size,chr_interval_len,chr_lst,interval_list)
	b=sort_coords(a)
	
	return b

	
#shared functions
def make_bed(gen_out):
	windows = []
	for count,chr,start,end in gen_out:
		windows.append(str(chr) +" "+ str(start) +" "+ str(end) +" "+ str(count))
	return windows

def sort_coords(coords,cols=itemgetter(1,2)):
	unsorted=[]
	coord=iter(coords)
	item=next(coord,None)
	
	while item:
		unsorted.append(item)
		item=next(coord,None)
	
	sorted_regions=sorted(unsorted, key=cols)
	
	for i in range(len(sorted_regions)):
		yield sorted_regions[i]

		
#main functions
def gaps_windows(ref_fai, gaps, no_of_samples, size):
	chr_len,chr_lst=gaps_chr_length(ref_fai)
	
	a=gaps_make_random(no_of_samples,size,chr_len,chr_lst,gaps)
	random_windows = make_bed(a)
	return random_windows
	
def intervals_window(interval_list, no_of_samples, size):	
	chr_interval_len,chr_lst=intervals_chr_length(interval_list)
	
	a=intervals_make_random(no_of_samples,size,chr_interval_len,chr_lst,interval_list)
	random_windows = make_bed(a)
	return random_windows


#if __name__=='__main__':
#	main(ref_fai, gapfile, no_of_samples, size)