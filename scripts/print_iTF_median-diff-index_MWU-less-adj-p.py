#!/usr/bin/env python3
#Description: Print the [iTF, median of differentiation index, one-sided MWU test adjusted p-value] table.
#Usage: [THIS SCRIPT]

for day in ["1", "3", "5"]:
	ordered = []
	with open(f"day{day}.PCA.z-values.for_figures.txt") as FOR_ORDER:
		for line in FOR_ORDER:
			ordered.append(line.strip().split("\t")[0])
	median_index = {}
	with open(f"day{day}.PCA.z-values.txt") as INDEX:
		INDEX.readline()
		for line in INDEX:
			words = line.strip().split("\t")
			tf, value = words[0], float(words[1])
			median_index[tf] = value
	adjp = {}
	with open(f"mwu_test_results.less.all_TF.txt") as PVALUE:
		for line in PVALUE:
			words = line.strip().split("\t")
			if words[0] == day:
				tf, value = words[1], words[2]
				adjp[tf] = value
	for itf in ordered:
		if itf != "Control":
			print(f"{day}\t{itf}\t{median_index[itf]}\t{adjp[itf]}")
		else:
			print(f"{day}\t{itf}\t{median_index[itf]}\tNA")
