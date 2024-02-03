#!/usr/bin/env python3
# Usage: [THIS SCRIPT] [Cell barcodes after QC by Seurat] [cb_gene_UMIcount.mat.txt]

import sys

cb_list = []
with open(sys.argv[1]) as CB_LIST:
	CB_LIST.readline()
	for line in CB_LIST:
		cb = line.strip().split("\t")[0]
		cb_list.append(cb)

output = open(sys.argv[2].replace("txt", "QCed_CBs.txt"), "w")
with open(sys.argv[2]) as INPUT:
	header_words = INPUT.readline().strip().replace("_UMIcount", "").split("\t")
	new_header = ""
	for i in range(1, len(header_words)):
		new_header += (header_words[i] + "\t")
	print("cell_barcode", new_header.strip(), sep = "\t", file = output)

	for line in INPUT:
		words = line.strip().split()
		count = 0
		out_str = ""
		for i in range(1, len(header_words)):
			count += int(words[i])
			out_str += (words[i] + "\t")
		if count == 0:
			print("ZERO LINE", line.strip())
		else:
			if words[0] in cb_list:
				print(words[0], out_str.strip(), sep = "\t", file = output)
output.close()
