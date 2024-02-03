#!/usr/bin/env python3

target_TFs = ["Control", "Gata3", "Fosl1", "Tbx5", "Cited1", "Id2", "Snai1"]

print("CB\tday_iTF\tiTF_of_interest\tday_iTF_of_interest\tExample_Control_only")
with open("int.SCT.metadata.tsv") as INPUT:
	INPUT.readline()
	for line in INPUT:
		words = line.strip().split("\t")
		cb, day, iTF = words[0], words[1], words[4]
		if iTF == "NO_iTF_TAGS":
			iTF = "Control"
		day_iTF = f"{day}_{iTF}"
		if iTF == "Control":
			print(cb, day_iTF, iTF, day_iTF, day_iTF, sep="\t")
		elif iTF in target_TFs:
			print(cb, day_iTF, iTF, day_iTF, "others", sep="\t")
		else:
			print(cb, day_iTF, "others", "others", "others", sep="\t")
