#!/usr/bin/env python3

import scipy.stats as stats
import statsmodels.stats.multitest as multitest

for day in ["1", "3", "5"]:
	total_dict = {}

	with open(f"day{day}.PCA.z-values.for_figures.txt") as INPUT:
		for line in INPUT:
			words = line.strip().split("\t")
			total_dict[words[0]] = []
			for word in words[1:]:
				total_dict[words[0]].append(float(word))

	p_value_list = []
	info_list = []

	for tf in total_dict:
		if tf == "Control":
			continue
		_, p_value = stats.mannwhitneyu(total_dict["Control"], total_dict[tf], alternative="less")
		p_value_list.append(p_value)
		info_list.append(tf)

	_, p_value_corrected, _, _ = multitest.multipletests(p_value_list, alpha=0.05, method="fdr_bh")

	for i in range(len(info_list)):
		print(day, info_list[i], str(p_value_corrected[i]), sep="\t")
