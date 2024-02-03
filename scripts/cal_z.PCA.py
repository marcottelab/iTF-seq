#!/usr/bin/env python3

import numpy as np

gene_to_cb = {}
with open("PCA/CBs_per_iTF.tsv") as CB_LIST:
	for line in CB_LIST:
		words = line.strip().split("\t")
		gene = words[0]
		gene_to_cb[gene] = words[1:]

cb_to_values = {}
with open("PCA/PCA_result.tsv") as VALUES:
	VALUES.readline()
	for line in VALUES:
		words = line.strip().split("\t")
		cb = words[0]
		values = []
		for i in range(1, len(words)):
			values.append(float(words[i]))
		cb_to_values[cb] = values

ctrl_center = []
for i in range(50):
	values = []
	for cb in gene_to_cb["NO_iTF_TAGS"]:
		values.append(cb_to_values[cb][i])
	ctrl_center.append(np.mean(values))
a = np.array(ctrl_center)

ctrl_cells = {}
ctrl_dists = []
for cb in gene_to_cb["NO_iTF_TAGS"]:
	b = np.array(cb_to_values[cb])
	dist = np.linalg.norm(a-b)
	ctrl_dists.append(dist)

ctrl_avg = np.mean(ctrl_dists)
ctrl_std = np.std(ctrl_dists)

print("TF\tmedian\tvalues")
for gene in sorted(gene_to_cb):
	z_list = []
	for cb in gene_to_cb[gene]:
		b = np.array(cb_to_values[cb])
		dist = np.linalg.norm(a-b)
		z = (dist - ctrl_avg) / ctrl_std 
		z_list.append(z)
	z_median = np.median(z_list)
	value_str = ""
	for z in z_list:
		value_str += f"{str(z)}\t"
	if gene  == "NO_iTF_TAGS":
		name = "Control"
	else:
		name = gene
	print(name, z_median, value_str.rstrip(), sep="\t")
