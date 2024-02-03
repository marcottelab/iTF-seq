#!/usr/bin/env python3

import sys
import numpy as np

tra_dict = {}

with open(sys.argv[1]) as BLAST_OUT:
	for line in BLAST_OUT:
		words = line.strip().split("\t")
		read_info, TF = words[0], words[1].split("_")[0]
		if TF == "GSC":
			TF = "Gsc"
		elif TF == "Neud4":
			TF = "Dpf1"
		elif TF == "Rbm9":
			TF = "Rbfox2"

		gene = "NA"
		for tag in read_info.split("_"):
			if tag.startswith("GN:"):
				gene = tag.split(":")[-1]
			elif tag.startswith("CB:"):
				cb = tag.split(":")[-1]
			elif tag.startswith("UB:"):
				umi = tag.split(":")[-1]
		transcript = cb + "_" + umi
		if transcript not in tra_dict:
				tra_dict[transcript] = [0, 0, 0]
		if gene == TF:
			tra_dict[transcript][0] += 1
		else:
			if gene != "NA":
				tra_dict[transcript][1] += 1
			else:
				tra_dict[transcript][2] += 1

for transcript in tra_dict:
	scores = tra_dict[transcript]
	qc_score = scores[1] / (scores[0] + scores[1] + scores[2])
	if qc_score >= 0.5:
		cb, umi = transcript.split("_")
		print(cb, umi, sep="\t")
