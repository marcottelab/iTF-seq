#!/usr/bin/env python3
# Usage: [THIS SCRIPT]

for day in ["1", "3", "5"]:
	cb_by_CR = set()
	with open(f"cells_by_CR.day{day}.tsv") as INPUT0:
		INPUT0.readline()
		for line in INPUT0:
			cb = line.strip().split("\t")[0].replace('"', '')
			cb_by_CR.add(cb)
		
	target_cb = set()
	zero_by_thr = set()
	with open(f"QCed/d{day}.CB_num_of_gene.UMI_th3.txt") as INPUT1:
		for line in INPUT1:
			cb, num_gene = line.strip().split("\t")
			if num_gene == "1":
				target_cb.add(cb)
			elif num_gene == "0":
				zero_by_thr.add(cb)
			
	without_any_tags = set()
	with open(f"d{day}.g1k_mt10_umi10k.without_iTF_tags.txt") as INPUT3:
		for line in INPUT3:
			cb = line.strip()
			without_any_tags.add(cb)

	output = open(f"d{day}.metadata.g1k_mt10_umi10k.umi3.temp.txt", "w")
	written_cb = set()
	with open(f"d{day}.cb_gene_UMIcount.mat.QCed_CBs.txt") as INPUT2:
		header_words = INPUT2.readline().strip().split("\t")
		temp = len(header_words) - 1
		zero_list = []
		for i in range(temp):
			zero_list.append("0")
		temp_words = [i + "_iTF" for i in header_words[1:]]
		print("CB", "detected_iTF", "\t".join(temp_words), "small_ctrl", "big_ctrl", sep="\t", file=output)

		for line in INPUT2:
			words = line.strip().split("\t")
			cb = words[0]
			if cb not in target_cb and cb not in (without_any_tags | zero_by_thr):
				print(cb, "NOT_USED", "\t".join(zero_list), 0, 0, sep="\t", file=output)
				written_cb.add(cb)
			elif cb in target_cb and cb in (without_any_tags | zero_by_thr):
				print(cb, "ERROR: Controls but target as well?!", sep="\t")
			elif cb in target_cb:
				count = 0
				for i in range(1, len(header_words)):
					if int(words[i]) >= 3: # UMI treshold
						temp_list = zero_list[:]
						temp_list[i-1] = "1"
						print(cb, header_words[i], "\t".join(temp_list), 0, 0, sep="\t", file=output)
						count += 1
				if count != 1:
					print(cb, "ERROR: having more than 1 iTF", sep="\t")
				written_cb.add(cb)

	for cb in (cb_by_CR - written_cb):
		if cb in without_any_tags:
			print(cb, "NO_iTF_TAGS", "\t".join(zero_list), 1, 1, sep="\t", file=output)	
		elif cb in zero_by_thr:
			print(cb, "ZERO_BY_THR", "\t".join(zero_list), 0, 1, sep="\t", file=output)
		else:
			print(cb, "NOT_USED", "\t".join(zero_list), 0, 0, sep="\t", file=output)

	output.close()
