#!/usr/bin/env python3

import sys

with open(sys.argv[1]) as INPUT:
	for line in INPUT:
		if "Dlx4_forward" in line:
			words = line.strip().split("\t")
			if int(words[3]) == 30:
				print(line.strip())
		else:
			print(line.strip())
