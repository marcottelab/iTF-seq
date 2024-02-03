#!/usr/bin/env python3

import os

for i in [1, 3, 5]:
	for j in [1, 2, 3, 4, 5]:
		cmd = f"./how_many_TFs_in_a_cell.py d{i}.cb_gene_UMIcount.mat.QCed_CBs.txt {j} > QCed/d{i}.CB_num_of_gene.UMI_th{j}.txt"
		cmd2 = f"./num_of_cells_having_N_kinds_of_TFs.py QCed/d{i}.CB_num_of_gene.UMI_th{j}.txt > QCed/d{i}.num_of_cells_having_N_kinds_of_TFs.UMI_th{j}.txt"
		print(cmd)
		os.system(cmd)
		print(cmd2)
		os.system(cmd2)
		
		cmd = f"./print_num_of_singlets_for_each_iTF.py d{i}.cb_gene_UMIcount.mat.QCed_CBs.txt {j} | sort -k2 -nr > d{i}.num_of_singlets_for_each_iTF.UMI_th{j}.txt"
		print(cmd)
		os.system(cmd)
