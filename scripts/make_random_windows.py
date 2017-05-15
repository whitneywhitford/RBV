#modified from http://genometoolbox.blogspot.co.nz/2014/06/generate-random-genomic-positions.html


#TO DO: include the CNVs to be investigated to be included in gaps?



import sys
from operator import itemgetter



#define functions

#WGS no intervals functions
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
			if (start<=point<=stop and start<=(point+size)<=stop):
				include="F"
		
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
	gap_raw=open(intervalfile).readlines()
	gaps={}
	for i in range(len(gap_raw)):
		chr,start,stop,type=gap_raw[i].strip().split()
		gap_region=[int(start),int(stop)]
		if chr in gaps:
			gaps[chr].append(gap_region)
		else:
			gaps[chr]=[gap_region]
	
	return gaps	

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
				window_start = point
				window_count = (stop - point) + 1
				past_start = True
			elif past_start: 
				window_count += (stop - start) + 1
				
			if window_count >= window_size:
				window_end = stop - (window_count - window_size)
				yield count,int(chr),window_start,window_end
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