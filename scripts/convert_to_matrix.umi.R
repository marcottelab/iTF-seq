library(tidyverse)

dayX = read_tsv("d1.cb_gene_UMIcount.including_no_align.revised.txt", col_names=FALSE)
names(dayX) <- c("A0_cb", "gene", "UMIcount")
data <- pivot_wider(dayX, names_from=gene, values_from=UMIcount, names_glue="{gene}_{.value}")
data <- data[ , order(names(data))]
write_tsv(data, "d1.cb_gene_UMIcount.mat.txt", na="0")

dayX = read_tsv("d3.cb_gene_UMIcount.including_no_align.revised.txt", col_names=FALSE)
names(dayX) <- c("A0_cb", "gene", "UMIcount")
data <- pivot_wider(dayX, names_from=gene, values_from=UMIcount, names_glue="{gene}_{.value}")
data <- data[ , order(names(data))]
write_tsv(data, "d3.cb_gene_UMIcount.mat.txt", na="0")

dayX = read_tsv("d5.cb_gene_UMIcount.including_no_align.revised.txt", col_names=FALSE)
names(dayX) <- c("A0_cb", "gene", "UMIcount")
data <- pivot_wider(dayX, names_from=gene, values_from=UMIcount, names_glue="{gene}_{.value}")
data <- data[ , order(names(data))]
write_tsv(data, "d5.cb_gene_UMIcount.mat.txt", na="0")
