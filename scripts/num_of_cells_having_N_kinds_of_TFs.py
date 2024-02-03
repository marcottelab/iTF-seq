#!/usr/bin/env python3

import sys

num_gene_dict = {}
with open(sys.argv[1]) as INPUT:
	for line in INPUT:
		cb, num_gene = line.strip().split()
		if num_gene not in num_gene_dict:
			num_gene_dict[num_gene] = 0
		num_gene_dict[num_gene] += 1

'''
for num_gene in sorted(num_gene_dict, key = int):
	print(num_gene, num_gene_dict[num_gene], sep = "\t")
'''

for i in range(25+1):
	if str(i) in num_gene_dict:
		print(str(i), num_gene_dict[str(i)], sep="\t")
	else:
		print(str(i), "0", sep="\t")
