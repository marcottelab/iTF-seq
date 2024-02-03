#!/usr/bin/env python3

days = ["1", "3", "5"]

for day in days:
	cb_to_zval = {}
	with open(f"day{day}.PCA.z-values.CB.txt") as INPUT1:
		INPUT1.readline()
		for line in INPUT1:
			words = line.strip().split("\t")
			cb, zval = words[1], words[2]
			cb_to_zval[cb] = zval

	output = open(f"d{day}.metadata.g1k_mt10_umi10k.umi3.PCA_z-values.txt", "w")

	with open(f"d{day}.metadata.g1k_mt10_umi10k.umi3.txt") as INPUT2:
		print(INPUT2.readline().strip(), "PCA_zvalues", sep="\t", file=output)
		for line in INPUT2:
			cb = line.strip().split("\t")[0]
			if cb in cb_to_zval:
				print(line.strip(), cb_to_zval[cb], sep="\t", file=output)
			else:
				print(line.strip(), "NA", sep="\t", file=output)
	output.close()
