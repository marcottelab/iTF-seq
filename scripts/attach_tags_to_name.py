#!/usr/bin/env python3
import sys

with open(sys.argv[1]) as INPUT:
	for line in INPUT:
		if line.startswith(">"):
			print(line.strip().replace("\t", "_"))
		else:
			print(line.strip())
