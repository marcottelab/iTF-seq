#!/usr/bin/env python3

days = ["1", "3", "5"]

for day in days:	
	output = open(f"d{day}.cb_gene_UMI.txt", "w")
	black_list = set()
	with open(f"day{day}.blast_out.transcripts_with_multi_iTFs.txt") as INPUT0:
		INPUT0.readline()
		for line in INPUT0:
			black_list.add(line.strip().split(" ")[0])
	
	with open(f"day{day}.blacklist.diff_align.txt") as INPUT1:
		INPUT1.readline()
		for line in INPUT1:
			words = line.strip().split("\t")
			transcript = words[0] + "_" + words[1]
			black_list.add(transcript)

	with open(f"day{day}.blast_out.Dlx4.txt") as INPUT2:
		for line in INPUT2:
			words = line.strip().split("\t")
			read_info = words[0]
			gene_by_tag = words[1].split("_")[0]
			if gene_by_tag == "GSC":
				gene_by_tag = "Gsc"
			elif gene_by_tag == "Neud4":
				gene_by_tag = "Dpf1"
			elif gene_by_tag == "Rbm9":
				gene_by_tag = "Rbfox2"

			geneID_by_CR = "NA"
			gene_by_CR = "NA"
			for tag in read_info.split("_"):
				if tag.startswith("CB:Z:"):
					cb = tag.split(":")[-1]
				elif tag.startswith("UB:Z:"):
					umi = tag.split(":")[-1]
				elif tag.startswith("GX:Z:"):
					geneID_by_CR = tag.split(":")[-1]
				elif tag.startswith("GN:Z:"):
					gene_by_CR= tag.split(":")[-1]
			transcript = cb + "_" + umi
			if transcript not in black_list:
				print(cb, gene_by_tag, umi, gene_by_CR, geneID_by_CR, sep="\t", file=output)
	output.close()
