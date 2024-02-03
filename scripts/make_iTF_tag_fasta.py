#!/usr/bin/env python3
# This is an example script to change a text file into FASTA format.

with open("iTF_tags.txt") as INPUT:
	INPUT.readline()
	for line in INPUT:
		_, gene, fw, rv = line.strip().split("\t")
		print(f">{gene}_forward\n{fw}")
		print(f">{gene}_reverse\n{rv}")
