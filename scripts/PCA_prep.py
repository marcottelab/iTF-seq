#!/usr/bin/env python3

import glob
import os

os.system("mkdir -p PCA")
file_list = glob.glob("*.tsv")

target_genes = []
with open("../common_genes.txt") as REF:
	for line in REF:
		target_genes.append(line.strip())

output1 = open("PCA/iTF-expressing_cells_and_ctrl.common_genes.sorted.tsv", "w")
output2 = open("PCA/CBs_per_iTF.tsv", "w")
cb_gene_exp = {}
for file_ in file_list:
	with open(file_) as INPUT:
		cb_list = INPUT.readline().strip().replace('"', '').split("\t")
		iTF = file_.split("/")[1].split(".")[0]
		print(iTF, "\t".join(cb_list), sep="\t", file=output2)
		for cb in cb_list:
			cb_gene_exp[cb] = {}
		for line in INPUT:
			words = line.strip().replace('"', '').split("\t")
			gene = words[0]
			if gene in target_genes:
				for i in range(1, len(words)):
					cb = cb_list[i-1]
					cb_gene_exp[cb][gene] = words[i]

print("CB", "\t".join(sorted(target_genes)), sep="\t", file=output1)
for cb in sorted(cb_gene_exp):
	output_str = f"{cb}\t"
	for gene in sorted(cb_gene_exp[cb]):
		output_str += (cb_gene_exp[cb][gene] + "\t")
	print(output_str.rstrip(), file=output1)
output1.close()
output2.close()
