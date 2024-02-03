#!/usr/bin/env python3

import sys

tran2tf = {}

with open(sys.argv[1]) as INPUT:
	for line in INPUT:
		words = line.strip().split("\t")
		qname, tf = words[0], words[1].split("_")[0]
		for tag in qname.split("_"):
			if tag.startswith("CB:"):
				cb = tag.split(":")[-1]
			elif tag.startswith("UB:"):
				umi = tag.split(":")[-1]
		transcript = cb + "_" + umi
		if transcript not in tran2tf:
			tran2tf[transcript] = set()
		tran2tf[transcript].add(tf)

print("total transcripts:", len(tran2tf))
for transcript in tran2tf:
	if len(tran2tf[transcript]) > 1:
		print(transcript, tran2tf[transcript])
