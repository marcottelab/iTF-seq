#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

for day in ["1", "3", "5"]:
	total_dict = {}

	with open(f"day{day}.PCA.z-values.for_figures.txt") as INPUT:
		for line in INPUT:
			words = line.strip().split("\t")
			total_dict[words[0]] = []
			for word in words[1:]:
				total_dict[words[0]].append(float(word))

	box_plot_data = []
	colors = []
	for tf in total_dict:
		box_plot_data.append(total_dict[tf])
		if tf == "Control":
			colors.append("gray")
		elif tf == "Gcm1" and day != "5":
			colors.append("gray")
		elif tf == "Otx2":
			colors.append("gray")
		else:
			colors.append("white")

	figure(figsize=(15, 7), dpi=300)
	plt.ylim(-2, 10)
	plt.xticks(rotation=90)
	box = plt.boxplot(box_plot_data, patch_artist=True, labels=total_dict.keys(),
				flierprops={'marker': 'o', 'markersize': 1}, widths=0.8)
	for patch, color in zip(box["boxes"], colors):
		patch.set_facecolor(color)	
	plt.savefig(f"day{day}.z-value.boxplots.pdf")		
