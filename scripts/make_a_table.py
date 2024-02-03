#!/usr/bin/env python3

import sys

thr = int(sys.argv[1])

tf2num_cells = {}
for day in ["1", "3", "5"]:
	with open(f"d{day}.num_of_singlets_for_each_iTF.UMI_th{thr}.txt") as INPUT:
		for line in INPUT:
			[tf, num_cells] = line.strip().split("\t")
			if tf not in tf2num_cells:
				tf2num_cells[tf] = {}
			tf2num_cells[tf][day] = int(num_cells)

print("iTF", "day1", "day3", "day5", sep="\t")

for tf in sorted(tf2num_cells):
	out_str = (tf + "\t")
	for day in ["1", "3", "5"]:
		if day not in tf2num_cells[tf]:
			out_str += "0\t"
		else:
			out_str += f"{tf2num_cells[tf][day]}\t"
	print(out_str.rstrip())
