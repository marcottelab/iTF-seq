#!/usr/bin/env python3
# Usage: [THIS SCRIPT] [d1.cb_gene_UMIcount.mat.QCed_CBs.txt] [threshold (int)]

import sys

num_of_genes_in_cb = {}
with open(sys.argv[1]) as INPUT:
	header_words = INPUT.readline().strip().split("\t")
	for line in INPUT:
		words = line.strip().split("\t")
		cb = words[0]
		if cb not in num_of_genes_in_cb:
			num_of_genes_in_cb[cb] = 0
		for i in range(1, len(header_words)):
			if int(words[i]) >= int(sys.argv[2]):
				num_of_genes_in_cb[cb] += 1

for cb in num_of_genes_in_cb:
	print(cb, num_of_genes_in_cb[cb], sep = "\t")
