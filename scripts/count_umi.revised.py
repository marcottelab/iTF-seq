#!/usr/bin/env python3

from collections import defaultdict
days = ["1", "3", "5"]

for day in days:
	print(day)
	count1 = defaultdict(dict)
	count2 = defaultdict(dict)
	output1 = open(f"d{day}.cb_gene_UMIcount.including_no_align.revised.txt", "w")
	output2 = open(f"d{day}.cb_gene_UMIcount.matched_align_only.revised.txt", "w")
	with open(f"d{day}.cb_gene_UMI.txt") as INPUT:
		for line in INPUT:
			words = line.strip().split("\t")
			cb, gene, umi, gene_CR = words[0], words[1], words[2], words[3]
			if gene_CR == "NA":
				if gene not in count1[cb]:
					count1[cb][gene] = set()
				count1[cb][gene].add(umi)
			else:
				#if gene != gene_CR:
				#	print(line.strip())
				if gene not in count1[cb]:
					count1[cb][gene] = set()
				if gene not in count2[cb]:
					count2[cb][gene] = set()
				count1[cb][gene].add(umi)
				count2[cb][gene].add(umi)
	for cb in sorted(count1):
		for gene in sorted(count1[cb]):
			print(cb, gene, len(count1[cb][gene]), sep="\t", file=output1)
	for cb in sorted(count2):
		for gene in sorted(count2[cb]):
			print(cb, gene, len(count2[cb][gene]), sep="\t", file=output2)
	output1.close()
	output2.close()
