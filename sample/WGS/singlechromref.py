with open('/gpfs1m/projects/nesi00322/Readbal/EAGLE_sims/human_g1k_v37_chrnumberonly.fasta') as b:
    ref_in = b.read().splitlines()

a = 0

for line in ref_in:
	if line.startswith(">"):
		chrom = line.strip(">")
		out = open(chrom + ".fa", 'w')
		out.write(str(line) + "\n")
	else:
		out.write(str(line) + "\n")

out.close()