#!/usr/bin/env python3
# Description: From the QC-passed cell list and number of genes per cell (nUMI >= 1, dX.CB_num_of_gene.UMI_th1.txt), print control cells.
# Usage: [THIS SCRIPT]

for day in ["1", "3", "5"]:
	targets = set()
	with open(f"d{day}.g1k_mt10_umi10k.tsv") as WHOLE_CELLS:
		WHOLE_CELLS.readline()
		for line in WHOLE_CELLS:
			cb = line.strip().split("\t")[0]
			targets.add(cb)
	with open(f"QCed/d{day}.CB_num_of_gene.UMI_th1.txt") as CELLS_WITH_ITF_TAGS:
		for line in CELLS_WITH_ITF_TAGS:
			cb, num_gene = line.strip().split("\t")
			if int(num_gene) > 0:
				targets.remove(cb)

	with open(f"d{day}.g1k_mt10_umi10k.without_iTF_tags.txt", "w") as OUTPUT:
		for cb in sorted(list(targets)):
			print(cb, file=OUTPUT)
