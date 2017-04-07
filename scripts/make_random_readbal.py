#modified from http://genometoolbox.blogspot.co.nz/2014/06/generate-random-genomic-positions.html


#TO DO: include the CNVs to be investigated to be included in gaps?



import sys
from operator import itemgetter

# Define parameters
ref_fai = sys.argv[1]
gapfile = sys.argv[2]
no_of_samples=int(sys.argv[3])
size=int(sys.argv[4])

#define functions	
#chr_length modified to fit Readbal
def chr_length():	
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
	
def chr_gaps():
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


def rand_sites(no_of_samples,size,chr_len,chr_lst,gaps):
	from random import choice,uniform,seed
	seed(None)
	count=1
	# Choose chromosome and point on chromosome
	while count<=no_of_samples:
		chr=choice(chr_lst)
		point=int(uniform(0+size/2,chr_len[chr]-size/2))
		
		# Exclude inaccessible regions
		include="T"
		for start,stop in gaps[chr]:
			if start<=point<=stop:
				include="F"
		
		# Return points in accessible regions
		if include=="T":
			yield count,int(chr),point-size/2,point+size/2
			count+=1


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


def make_random(no_of_samples,size,chr_len,chr_lst,gaps):
	from operator import itemgetter
	
	a=rand_sites(no_of_samples,size,chr_len,chr_lst,gaps)
	b=sort_coords(a)
	
	return b


def make_bed(gen_out, out_file_dir):
	out=open(out_file_dir, "w")
	for count,chr,start,end in gen_out:
		print >> out, str(chr) +"\t"+ str(start) +"\t"+ str(end) +"\t"+ str(count)


def main(ref_fai, gapfile, no_of_samples, size):
	chr_len,chr_lst=chr_length()
	gaps=chr_gaps()
	
	a=make_random(no_of_samples,size,chr_len,chr_lst,gaps)
	make_bed(a, "permutation_"+".bed")
	print "Permutation "+" complete."


if __name__=='__main__':
	main(ref_fai, gapfile, no_of_samples, size)