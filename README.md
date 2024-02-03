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

## Mapping and cell quality control
### Mapping (alignment)
We used Cell Ranger Count[^2]. Among the outputs, we used *outs/filtered_feature_bc_matrix* directory for the next step.
### Cell Quality Control
We used the Seurat R package[^3] to check nFeature_RNA, nCount_RNA, and percent.mt.  
Our final criteria: nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000



## References
[^1]: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421. https://doi.org/10.1186/1471-2105-10-421.
[^2]: Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, et al. 2017. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8: 14049. https://doi.org/10.1038/ncomms14049.
[^3]: Hao Y, Hao S, Andersen-Nissen E, Mauck WM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. 2021. Integrated analysis of multimodal single-cell data. Cell 184: 3573-3587.e29. https://www.sciencedirect.com/science/article/pii/S0092867421005833.
