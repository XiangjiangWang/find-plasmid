# Setup environment
## Bowtie2
`conda create -n bt2 bowtie2`
## samtools
[Build from source](http://www.htslib.org/download/)
## pysam
pysam is acting weird so have to create a separate env
`conda create -n ngs pysam=0.20 tqdm`


# Align to human genome (roughly) and filter out aligned reads
```
bowtie2 -p 10 --very-fast -x GRCh38_noalt_as/GRCh38_noalt_as -1 U2OS-LacO-TetO_R1_paired.fastq.gz -2 U2OS-LacO-TetO_R2_paired.fastq.gz -S mu.sam
```
```
373498309 reads; of these:
  373498309 (100.00%) were paired; of these:
    6300223 (1.69%) aligned concordantly 0 times
    259419471 (69.46%) aligned concordantly exactly 1 time
    107778615 (28.86%) aligned concordantly >1 times
    ----
    6300223 pairs aligned concordantly 0 times; of these:
      911113 (14.46%) aligned discordantly 1 time
    ----
    5389110 pairs aligned 0 times concordantly or discordantly; of these:
      10778220 mates make up the pairs; of these:
        5939657 (55.11%) aligned 0 times
        1605143 (14.89%) aligned exactly 1 time
        3233420 (30.00%) aligned >1 times
99.20% overall alignment rate
```

```
samtools view -@ 10 -b -o mu.bam mu.sam
samtools sort -@ 10 -o mu.sorted.bam mu.bam
samtools index -@ 10 mu.sorted.bam
python filter.py
samtools collate -@ 10 -o filtered_collate.bam filtered.bam
samtools fastq -@ 10 -1 filtered_R1.fastq.gz -2 filtered_R2.fastq.gz -0 filtered_0.fastq -s filtered_single.fastq filtered_collate.bam
```
```
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 7797180 reads
filtered_0.fastq and filtered_single.fastq are empty
```


# Align to human genome again (failed, too many candidates)
```
bowtie2 -p 10 --local --very-sensitive-local -x GRCh38_noalt_as/GRCh38_noalt_as -1 filtered_R1.fastq.gz -2 filtered_R2.fastq.gz -S final.sam
```
```
3898590 reads; of these:
  3898590 (100.00%) were paired; of these:
    983501 (25.23%) aligned concordantly 0 times
    78991 (2.03%) aligned concordantly exactly 1 time
    2836098 (72.75%) aligned concordantly >1 times
    ----
    983501 pairs aligned concordantly 0 times; of these:
      25586 (2.60%) aligned discordantly 1 time
    ----
    957915 pairs aligned 0 times concordantly or discordantly; of these:
      1915830 mates make up the pairs; of these:
        1270057 (66.29%) aligned 0 times
        206205 (10.76%) aligned exactly 1 time
        439568 (22.94%) aligned >1 times
83.71% overall alignment rate
```

```
samtools view -@ 10 -b -o final.bam final.sam
samtools sort -@ 10 -o final.sorted.bam final.bam
samtools index -@ 10 final.sorted.bam
python analysis.py
samtools index -@ 10 candidates.bam
```


# Align to plasmid (laco-i-scei-teto_curated.fasta)
```
bowtie2-build plasmid/laco-i-scei-teto_curated.fasta plasmid/plasmid
```
```
bowtie2 -p 10 --local --very-sensitive-local -x plasmid/plasmid -1 filtered_R1.fastq.gz -2 filtered_R2.fastq.gz -S final.sam
```
```
3898590 reads; of these:
  3898590 (100.00%) were paired; of these:
    3804477 (97.59%) aligned concordantly 0 times
    18951 (0.49%) aligned concordantly exactly 1 time
    75162 (1.93%) aligned concordantly >1 times
    ----
    3804477 pairs aligned concordantly 0 times; of these:
      1904 (0.05%) aligned discordantly 1 time
    ----
    3802573 pairs aligned 0 times concordantly or discordantly; of these:
      7605146 mates make up the pairs; of these:
        7598652 (99.91%) aligned 0 times
        2773 (0.04%) aligned exactly 1 time
        3721 (0.05%) aligned >1 times
2.55% overall alignment rate
```

```
samtools view -@ 10 -b -o final.bam final.sam
samtools sort -@ 10 -o final.sorted.bam final.bam
samtools index -@ 10 final.sorted.bam
python analysis.py
samtools index -@ 10 candidates.bam
```

## Results
Found 8 candidate sequences (eyeballing soft-clip.txt and candidates.bam on IGV)

### Local BLAST against plasmid
```
makeblastdb -in plasmid/laco-i-scei-teto_curated.fasta -dbtype nucl -out plasmid/plasmid
```
```
blastn -db plasmid/plasmid -query BLAST/blast_in.fasta -out BLAST/blast_out.txt -task blastn-short
```
4 sequences are from plasmid

### BLAST the other 4 sequences against human ([NCBI online BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome))
 - Not found in Human RefSeq representative genome using blastn (somewhat similar), try other databases
 - Database: Whole-genome shotgun contigs (wgs), Organism: human (taxid: 9606)
 - Program: Somewhat similar sequences (blastn)
 - Results: 44A9DJ97016-Alignment.txt
```
>1 (plasmid)
ATTACCCTGTTATCCCTACCCTCGAGGGATCGACTCTAGAGGCGCCGA
>2 (found contig)
CCCCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGTCAGCTGACGCGT
>3 (plasmid)
CCGGGTCGAGGTAGGCGTGTACGGTGGGAGGCCTATATAAGCAGAGCT
>4 (plasmid)
GAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAA
>5 (not found)
GTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAA
>6 (not found)
TCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACTGGGAGTAGTC
>7 (found contig)
TCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTAT
>8 (plasmid)
TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAA


# examples for candidate reads
>NB551658:237:HCW3VBGXN:4:21604:16146:2857 right soft (6)
TCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACTGGGAGTAGTC
>NB551658:237:HCW3VBGXN:1:12101:7542:15042 left soft (2)
GCGTCAGCTGACTAGAGGATCCCCGGGTACCGAGCTCGAATTCGGGG
>NB551658:237:HCW3VBGXN:2:22307:13446:1253 left soft (3, plasmid)
AGCTCTGCTTATATAGGCCTCCCACCGTACACGCCTACCTCGACCCGG
>NB551658:237:HCW3VBGXN:3:13404:1730:16508 left soft (7)
ATACGCCTATTTTTATAGGTTAATGTCATGATAATAATGGTTTCTTAGA


# these two are plasmid break junction
>NB551658:237:HCW3VBGXN:3:23401:10026:19181 left soft (4)
CACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTC TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCT
>NB551658:237:HCW3VBGXN:4:22610:19235:8689 right soft (8)
TCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTC TTGAGATCCTTTTTTTCTGCGCGTAATCT
```
