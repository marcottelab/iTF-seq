#!/usr/bin/env python3
# Usage: [THIS SCRIPT] [d1.cb_gene_UMIcount.mat.QCed_CBs.txt] [threshold (int)]

import sys

genes_in_cb = {}
with open(sys.argv[1]) as INPUT:
	header_words = INPUT.readline().strip().split("\t")
	for line in INPUT:
		words = line.strip().split("\t")
		cb = words[0]
		if cb not in genes_in_cb:
			genes_in_cb[cb] = []
		for i in range(1, len(header_words)):
			if int(words[i]) >= int(sys.argv[2]):
				genes_in_cb[cb].append(header_words[i])

gene2num_of_cells = {}
for cb in genes_in_cb:
	if len(genes_in_cb[cb]) == 1:
		gene = genes_in_cb[cb][0]
		if gene not in gene2num_of_cells:
			gene2num_of_cells[gene] = 0
		gene2num_of_cells[gene] += 1

for gene in gene2num_of_cells:
	print(gene, gene2num_of_cells[gene], sep="\t")
