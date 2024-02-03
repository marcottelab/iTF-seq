# Systematic mapping of TF-mediated cell fate changes by a pooled induction coupled with scRNA-seq and multi-omics approaches
Transcriptional regulation controls cellular functions through interactions between transcription factors (TFs) and their chromosomal targets. However, understanding the fate conversion potential of multiple TFs in an inducible manner remains limited. To address this, iTF-seq was introduced as a method for identifying individual TFs that can alter cell fate towards specific lineages at a single-cell level. iTF-seq enables time-course monitoring of transcriptome changes, and with biotinylated individual TFs, it provides a multi-omics approach to understanding the mechanisms behind TF-mediated cell fate changes. Our iTF-seq study in mouse embryonic stem cells identified multiple TFs that trigger rapid transcriptome changes indicative of differentiation within a day of induction. Moreover, cells expressing these potent TFs often show a slower cell cycle and increased cell death. Further analysis using bioChIP-seq revealed that Gcm1 and Otx2 act as pioneer factors and activators by increasing gene accessibility and activating the expression of lineage specification genes during cell fate conversion. iTF-seq has utility in both mapping cell fate conversion and understanding cell fate conversion mechanisms.

# iTF-seq pipeline
This is the overall explanation of the iTF-seq pipeline, which includes Bash, Python, and R scripts. We will focus on newly generated scripts rather than how to run already well-known programs.

## Validation of iTF tag sequences
>[!Important]
>This is a prerequisite step.

You need to check whether your tag sequences can be found from reference sequences. We used NCBI BLAST+[^1].
```bash
# Example codes
makeblastdb -in genome.fa -parse_seqids -dbtype nucl
blastn -db genome.fa -query tags.fa -out blast_result.txt -perc_identity 80 -outfmt 6 -num_threads 4
gffread genes.gtf -g genome.fa -w transcripts.fa
makeblastdb -in transcripts.fa -parse_seqids -dbtype nucl
blastn -db transcripts.fa -query tags.fa -out blast_result.transcripts.txt -outfmt 6 -num_threads 4
```

## Mapping of scRNA-seq reads and cell quality control
### Mapping (alignment)
We used Cell Ranger Count[^2]. Among the outputs, we used the *outs/filtered_feature_bc_matrix* directory for the next step.
### Cell Quality Control
We used the Seurat R package[^3] to check nFeature_RNA, nCount_RNA, and percent.mt.  
Our final criteria: nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000

## Detection of cells overexpressing single iTF and control cells
### BAM to FASTA conversion
Reads are filtered by samtools[^4] based on [barcoded BAM tags by Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam). The *outs/possorted_genome_bam.bam* from Cell Ranger was used.
```bash
samtools collate -u --threads 4 outs/possorted_genome_bam.bam dayX.collated.bam
samtools view --tag CB dayX.collated.bam -b | samtools view --tag UB -b -o dayX.reads_with_CB_UB.bam
samtools fasta dayX.reads_with_CB_UB.bam -T CB,UB,GX,GN,xf --threads 4
./attach_tags_to_name.py dayX.reads_with_CB_UB.fasta > dayX.reads_with_CB_UB.v2.fasta
```
### Run BLASTN against iTF tags
```bash
./make_iTF_tag_fasta.py > iTF_tags.fasta
makeblastdb -in iTF_tags.fasta -dbtype nucl -parse_seqids
blastn -db iTF_tags.fasta -query dayX.reads_with_CB_UB.v2.fasta -num_threads 4 -perc_identity 85 -evalue 1e-7 -max_target_seqs 5 -outfmt 6 -out dayX.blast_out.txt
```
### Remove problematic transcripts from BLASTN results
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






## References
[^1]: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421. https://doi.org/10.1186/1471-2105-10-421.
[^2]: Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, et al. 2017. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8: 14049. https://doi.org/10.1038/ncomms14049.
[^3]: Hao Y, Hao S, Andersen-Nissen E, Mauck WM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. 2021. Integrated analysis of multimodal single-cell data. Cell 184: 3573-3587.e29. https://www.sciencedirect.com/science/article/pii/S0092867421005833.
[^4]: Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, et al. 2021. Twelve years of SAMtools and BCFtools. Gigascience 10: giab008. https://doi.org/10.1093/gigascience/giab008.
