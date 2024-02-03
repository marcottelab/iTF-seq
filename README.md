# Systematic mapping of TF-mediated cell fate changes by a pooled induction coupled with scRNA-seq and multi-omics approaches
Transcriptional regulation controls cellular functions through interactions between transcription factors (TFs) and their chromosomal targets. However, understanding the fate conversion potential of multiple TFs in an inducible manner remains limited. To address this, iTF-seq was introduced as a method for identifying individual TFs that can alter cell fate towards specific lineages at a single-cell level. iTF-seq enables time-course monitoring of transcriptome changes, and with biotinylated individual TFs, it provides a multi-omics approach to understanding the mechanisms behind TF-mediated cell fate changes. Our iTF-seq study in mouse embryonic stem cells identified multiple TFs that trigger rapid transcriptome changes indicative of differentiation within a day of induction. Moreover, cells expressing these potent TFs often show a slower cell cycle and increased cell death. Further analysis using bioChIP-seq revealed that Gcm1 and Otx2 act as pioneer factors and activators by increasing gene accessibility and activating the expression of lineage specification genes during cell fate conversion. iTF-seq has utility in both mapping cell fate conversion and understanding cell fate conversion mechanisms.

# iTF-seq pipeline
This is the overall explanation of the iTF-seq pipeline, which includes Bash, Python, and R scripts. I will focus on newly generated scripts rather than how to run already well-known programs.

## Validation of iTF tag sequences
[!Important]
>This is a prerequisite step.
You need to check whether your tag sequences can be found from reference sequences.
```bash
# Example codes
makeblastdb -in genome.fa -parse_seqids -dbtype nucl
blastn -db genome.fa -query tags.fa -out blast_result.txt -perc_identity 80 -outfmt 6 -num_threads 4

gffread genes.gtf -g genome.fa -w transcripts.fa
makeblastdb -in transcripts.fa -parse_seqids -dbtype nucl
blastn -db transcripts.fa -query tags.fa -out blast_result.transcripts.txt -outfmt 6 -num_threads 4
```
