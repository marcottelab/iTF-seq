#!/usr/bin/env python3

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

data = pd.read_table("PCA/iTF-expressing_cells_and_ctrl.common_genes.sorted.tsv", header=0)

pca = PCA(n_components=50)
principal_components = pca.fit_transform(data.iloc[:, 1:])
principal_df = pd.DataFrame(data=principal_components)
output_df = pd.concat([data["CB"], principal_df], axis=1)
output_df.to_csv("PCA/PCA_result.tsv", sep="\t", index=False)
