#modified from http://genometoolbox.blogspot.co.nz/2014/06/generate-random-genomic-positions.html


#TO DO: include the CNVs to be investigated to be included in gaps?



import sys
from operator import itemgetter



#define functions

#WGS functions

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

def WGS_intervals(intervalfile):
	# Import chromosome gaps
	interval_raw=open(intervalfile).readlines()
	intervals={}
	for i in range(len(interval_raw)):
		chr,start,stop,type=gap_raw[i].strip().split()
		interval_region=[int(start),int(stop)]
		if chr in gaps:
			intervals[chr].append(interval_region)
		else:
			intervals[chr]=[interval_region]
	
	return intervals


#chr_length modified to fit Readbal
def WGS_chr_length(ref_fai):	
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


def WGS_rand_sites(no_of_samples,size,chr_len,chr_lst,gaps):
	from random import choice,uniform,seed
	seed(None)
	count=1
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(0+size/2,chr_len[chr]-size/2))
		
		rand_start = point-size/2
		rand_stop = point+size/2
		
		# Exclude inaccessible regions
		include="T"
		for start,stop in gaps[chr]:
			if (start<=rand_start<=stop or start<=rand_stop<=stop) and ((rand_start not in range(start,stop+1)) or (rand_stop not in range(start,stop+1))):
				include="F"
		
		# Return points in accessible regions
		if include=="T":
			yield count,int(chr),rand_start,rand_stop
			count+=1


def WGS_make_random(no_of_samples,size,chr_len,chr_lst,gaps):
	from operator import itemgetter
	
	a=WGS_rand_sites(no_of_samples,size,chr_len,chr_lst,gaps)
	b=sort_coords(a)
	
	return b

		
		
#WES functions
def WES_intervals(intervalfile):
	interval_raw=open(intervalfile).readlines()
	intervals={}
	for i in range(len(interval_raw)):
		chr,start,stop,type=interval_raw[i].strip().split()
		interval_region=[int(start),int(stop)]
		if chr in intervals:
			intervals[chr].extend(range(int(start),(int(stop))+1))
		else:
			intervals[chr]=range(int(start),(int(stop)+1))
			
	return intervals
	
	
def WES_chr_length(intervals):	
	chr_exome_len={}
	for chr in intervals:
		chr_exome_len[chr]=int(len(intervals[chr]))
	
	chr_lst=[]
	for chr in chr_exome_len:							#account for different lengths of chromosomes so that random takes this into account
		chr_rep=[str(chr)] * int(round(chr_exome_len[str(chr)]/100,0)) #100 needs to be smaller than the smallest total chromosome exome intervals 
		for rep in chr_rep:
			chr_lst.append(rep)
	
	return chr_exome_len, chr_lst
	
def WES_rand_sites(no_of_samples,size,chr_exome_len,chr_lst,intervals):
	from random import choice,uniform,seed
	seed(None)
	count=1
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(0+size/2,chr_exome_len[chr]-size/2))
		
		chr_intervals = intervals[chr]
		
		print chr, point, chr_intervals[point]
		
		start = chr_intervals[point-size/2]
		stop = chr_intervals[point+size/2]
		
		yield count,int(chr),start,stop
		count+=1


def WES_make_random(no_of_samples,size,chr_exome_len,chr_lst,intervals):
	from operator import itemgetter
	
	a=WES_rand_sites(no_of_samples,size,chr_exome_len,chr_lst,intervals)
	b=sort_coords(a)
	
	return b
	

def make_bed(gen_out):
	windows = []
	for count,chr,start,end in gen_out:
		windows.append(str(chr) +" "+ str(start) +" "+ str(end) +" "+ str(count))
	return windows


#shared functions
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
		

def WGS(ref_fai, gaps, no_of_samples, size, iteration_number, intervals):
	chr_len,chr_lst=WGS_chr_length(ref_fai)
	
	a=WGS_make_random(no_of_samples,size,chr_len,chr_lst,gaps)
	random_windows = make_bed(a)
	return random_windows
	
def WES(intervals, no_of_samples, size):
	chr_exome_len,chr_lst=WES_chr_length(intervals)
	
	a=WES_make_random(no_of_samples,size,chr_exome_len,chr_lst,intervals)
	random_windows = make_bed(a)
	return random_windows


#if __name__=='__main__':
#	main(ref_fai, gapfile, no_of_samples, size)