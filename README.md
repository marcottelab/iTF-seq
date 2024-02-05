# Systematic mapping of TF-mediated cell fate changes by a pooled induction coupled with scRNA-seq and multi-omics approaches
Transcriptional regulation controls cellular functions through interactions between transcription factors (TFs) and their chromosomal targets. However, understanding the fate conversion potential of multiple TFs in an inducible manner remains limited. To address this, ***iTF-seq was introduced as a method for identifying individual TFs that can alter cell fate towards specific lineages at a single-cell level. iTF-seq enables time-course monitoring of transcriptome changes, and with biotinylated individual TFs, it provides a multi-omics approach to understanding the mechanisms behind TF-mediated cell fate changes.*** Our iTF-seq study in mouse embryonic stem cells identified multiple TFs that trigger rapid transcriptome changes indicative of differentiation within a day of induction. Moreover, cells expressing these potent TFs often show a slower cell cycle and increased cell death. Further analysis using bioChIP-seq revealed that Gcm1 and Otx2 act as pioneer factors and activators by increasing gene accessibility and activating the expression of lineage specification genes during cell fate conversion. ***iTF-seq has utility in both mapping cell fate conversion and understanding cell fate conversion mechanisms.***

# iTF-seq pipeline
This is the overall explanation of the iTF-seq pipeline, which includes shell commands, Python scripts, and R scripts. We uploaded those scripts to the ["scripts" directory](./scripts) in this repository. We will focus on newly generated scripts rather than how to run already well-known programs.
>[!Note]
>We used the expression "dayX" in the following code blocks to describe a code line that was repeatedly used for day 1, day 3, and day 5 data.

## 0. Validation of iTF tag sequences
>[!Important]
>This is a prerequisite step.

You need to check whether your tag sequences can be found from reference sequences. We made a FASTA file with our tag sequences and ran NCBI BLAST+[^1] agaisnt references.
```bash
# Example codes
# Make a FASTA file with your tag sequences.
./make_iTF_tag_fasta.py > iTF_tags.fasta 

makeblastdb -in genome.fa -parse_seqids -dbtype nucl
blastn -db genome.fa -query iTF_tags.fa -out blast_result.txt -perc_identity 80 -outfmt 6 -num_threads 4
# You may adjust perc_ientity, outfmt, and num_threads based on your needs.

gffread genes.gtf -g genome.fa -w transcripts.fa
makeblastdb -in transcripts.fa -parse_seqids -dbtype nucl
blastn -db transcripts.fa -query iTF_tags.fa -out blast_result.transcripts.txt -outfmt 6 -num_threads 4
```

## 1. Mapping and cell quality control
### Mapping of scRNA-seq reads (alignment)
We ran 10x Genomics Cell Ranger Count[^2] and used the *outs/filtered_feature_bc_matrix* directory for the following steps. If you use other aligners, create your own Seurat[^3] Objects with similar feature-barcode matrices.
### Cell Quality Control
We used the Seurat R package to check nFeature_RNA(The number of unique genes detected in each cell), nCount_RNA(the total number of molecules detected within a cell), and percent.mt(The percentage of reads that map to the mitochondrial genome).  
Our final criteria: nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000

## 2. Detection of cells overexpressing single iTF and control cells
### BAM to FASTA conversion
Reads are filtered by samtools[^4] based on [barcoded BAM tags by Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam). The *outs/possorted_genome_bam.bam* from Cell Ranger was used. Reads without either cell barcodes or unique molecular identifiers (UMIs) were removed.
```bash
samtools collate -u --threads 4 outs/possorted_genome_bam.bam dayX.collated.bam
samtools view --tag CB dayX.collated.bam -b | samtools view --tag UB -b -o dayX.reads_with_CB_UB.bam
samtools fasta dayX.reads_with_CB_UB.bam -T CB,UB,GX,GN,xf --threads 4
./attach_tags_to_name.py dayX.reads_with_CB_UB.fasta > dayX.reads_with_CB_UB.v2.fasta
```
### Run BLASTN against iTF tags
BLAST was performed to detect reads aligned with the tags.
```bash
makeblastdb -in iTF_tags.fasta -dbtype nucl -parse_seqids
blastn -db iTF_tags.fasta -query dayX.reads_with_CB_UB.v2.fasta -num_threads 4 -perc_identity 85 -evalue 1e-7 -max_target_seqs 5 -outfmt 6 -out dayX.blast_out.txt
```
### Remove problematic transcripts from BLASTN results
Collected reads with the tags were organized to the transcript level. The alignment results based on our tags were also compared to the alignment by Cell Ranger. When more than half of a transcript's reads aligned to two different genes by the two methods, that transcript was excluded. We did not exclude the cases where a gene was annotated as "NA" by Cell Ranger but aligned with our tags and used the alignment results for further analysis.
```bash
# Double-check for the Dlx4 forward tag
grep Dlx4_forward day?.blast_out.txt | cut -f4 | sort | uniq 
./remove_Dlx4_shorter_than_30bp.py dayX.blast_out.txt > dayX.blast_out.Dlx4.txt

# Print transcripts that have multiple iTF tags
./check_transcripts.py dayX.blast_out.Dlx4.txt > dayX.blast_out.transcripts_with_multi_iTFs.txt

# Compare the alignment by Cell Ranger and by our tags.
# If a transcript has equal to or more than half of its reads aligned to two different genes by the two methods, print it.
./print_blacklist_transcripts_by_qc_score.py dayX.blast_out.Dlx4.txt > dayX.blacklist.diff_align.txt
./organize_blast_output.py
# outputs: dX.cb_gene_UMI.txt
./count_umi.revised.py
# outputs: dX.cb_gene_UMIcount.including_no_align.revised.txt,
           dX.cb_gene_UMIcount.matched_align_only.revised.txt
```
### Selection of UMI thresholds and preparation of metadata file for Seurat
To reduce false positive events in detecting iTFs, we exploited UMIs. We applied different thresholds for the minimum number of UMI (1-5) and counted cells with single iTF detection. The control cells were identified by excluding cells expressing any iTFs with a minimum UMI threshold of 1.
```bash
Rscript convert_to_matrix.umi.R
# outputs: dX.cb_gene_UMIcount.mat.txt
Rscript print_QCed_cell_list.R
# outputs: dX.g1k_mt10_umi10k.tsv # These are all QC-passed cells regardless of iTF tags.
./trim_matrix.py dX.g1k_mt10_umi10k.tsv dX.cb_gene_UMIcount.mat.txt
# outputs: dX.cb_gene_UMIcount.mat.QCed_CBs.txt

# Testing UMI threshold 1 - 5. (equal to or more than)
mkdir QCed
./cmds.py # This script runs the following three scripts.
  a. how_many_TFs_in_a_cell.py
     outputs: QCed/dX.CB_num_of_gene.UMI_thX.txt
  b. num_of_cells_having_N_kinds_of_TFs.py
     outputs: QCed/dX.num_of_cells_having_N_kinds_of_TFs.UMI_thX.txt
  c. print_num_of_singlets_for_each_iTF.py
     outputs: dX.num_of_singlets_for_each_iTF.UMI_thX.txt

# Print cell barcodes of control cells
./find_QCpassed_no_iTF_cells.py 
# outputs: dX.g1k_mt10_umi10k.without_iTF_tags.txt
```
With *QCed/dX.num_of_cells_having_N_kinds_of_TFs.UMI_thX.txt*, you can determine the optimal UMI threshold. In our case, we draw *number of cells* x *number of genes in cells* x *UMI threshold* chart. In our case, as the number of cells with a single iTF peaked in the minimum UMI 3 threshold, this threshold was adopted to designate iTF cells. 
```bash
./make_a_table.py 3 > num_of_cells_with_one_iTF.UMI_th3.txt
./print_metadata_for_seurat.py # This file requires cells_by_CR.dayX.tsv, simple metadata files from Seurat Object before the cell quality control steps.
# outputs: dX.metadata.g1k_mt10_umi10k.umi3.txt
```
Explanation about the 2nd column ("detected_iTF") of the metadata files.
* Name of TFs: iTF cells with that TF
* "NO_iTF_TAGS": cells without any iTF tags (We chose these as controls).
* "ZERO_BY_THR": cells without iTF because of the UMI threshold.
* "NOT USED": The others, like poor-quality cells, cells having multiple iTFs. Literally not used for analyses.

## 3. Calculation of the differentiation index
> [!important]
> ### Make a text file with your genes of interest
> For example, for "all genes", we used genes that are commonly detected for all time points, excluding our 80 TFs of interest. You can use other gene sets related to your study. The format is a text file having one gene name per line.

```bash
# This R script prints "scale.data" after applying SCTransform.
# Also, it requires outs/filtered_feature_bc_matrix, dX.metadata.g1k_mt10_umi10k.umi3.txt(metadata made from previous steps)
Rscript print_expression_values.SCT.dayX.R 
# Outputs: dayX/TF.tsv files and dayX/NO_iTF_TAGS.tsv
```
### Calculation
Repeat the following code block for each time point.
```bash
# Principal component analysis
cd dayX
./PCA_prep.py
# outputs: PCA/iTF-expressing_cells_and_ctrl.common_genes.sorted.tsv
           PCA/CBs_per_iTF.tsv
./PCA.py
# output: PCA/PCA_result.tsv
./cal_z.PCA.py > ../dayX.PCA.z-values.txt
# If you need each cell's score, use cal_z.PCA.with_CB.py
./cal_z.PCA.with_CB.py > ../dayX.PCA.z-values.CB.txt
```
```bash
# cd .. (Not in each day's folder)
./add_z-values_to_metadata.py
# outputs: dX.metadata.g1k_mt10_umi10k.umi3.PCA_z-values.txt
```
### Drawing boxplots
```bash
# Before running the following Python scripts,
# Sort dayX.PCA.z-values.txt based on the median of z-values and save as dayX.PCA.z-values.for_figures.txt
# Then...
./boxplots.py
# outputs: dayX.z-value.boxplot.pdf
```
### Predict potent iTFs by comparing their differentiation indexes to controls.
```bash
./mwu_test_with_cell_diff_index.py > mwu_test_results.less.all_TF.txt
./print_iTF_median-diff-index_MWU-less-adj-p.py > iTF_median-diff-index_MWU-less-adjp.tsv
```
You may use [Morpheus](https://software.broadinstitute.org/morpheus)[^5] to conduct clusterings or make a similarity matrix.

## 4. Visualization (UMAP plots)
```bash
Rscript umap_plots.dayX.R
# This code requires
    a. filtered_feature_bc_matrix from Cell Ranger
    b. metadata files from previous steps
# outputs: seurat_data.dayX.rds,
           overview UMAP plot
           UMAP plot for expression of iTF (endo + ecto expression),
           UMAP plot for the location of iTF cells (by tags),
```
### Integrated version of multiple time points
This is an example code. You would need to edit this script for your study.
```bash
Rscript integrate.SCT.R
# This code also requires
    a. filtered_feature_bc_matrix from Seurat
    b. metadata files from 3_Detecting_cells_having_single_iTF_and_controls
# outputs: int.SCT.after_integration.rds,
           int.SCT.after_clustering.rds,
           int.SCT.metadata.tsv,
           integrated_plot.pdf
```
If you need UMAP focusing on specific TFs, you may refer to the following code block.
```bash
./print_meta_for_some_TFs.py > metadata_for_some_TFs.tsv
Rscript someTFs.visulazation.R
```

## References
[^1]: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421. https://doi.org/10.1186/1471-2105-10-421.
[^2]: Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, et al. 2017. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8: 14049. https://doi.org/10.1038/ncomms14049.
[^3]: Hao Y, Hao S, Andersen-Nissen E, Mauck WM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. 2021. Integrated analysis of multimodal single-cell data. Cell 184: 3573-3587.e29. https://www.sciencedirect.com/science/article/pii/S0092867421005833.
[^4]: Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, et al. 2021. Twelve years of SAMtools and BCFtools. Gigascience 10: giab008. https://doi.org/10.1093/gigascience/giab008.
[^5]: Morpheus, https://software.broadinstitute.org/morpheus
