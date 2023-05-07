# Setup environment
## Bowtie2 and blast
```
conda create -n bt2 bowtie2 blast --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority
```
## samtools
[Build from source](http://www.htslib.org/download/)
## pysam
pysam is acting weird so have to create a separate env
```
conda create -n ngs pysam=0.20 tqdm --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority
```


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


# de novo assembly of filtered reads
```
abyss-pe name=plasmid k=36 B=2G in='../filtered_R1.fastq.gz ../filtered_R2.fastq.gz' j=10
```
cannot find a plasmid contig, since only ~2.5% of reads are plasmid. 
```
abyss-pe name=plasmid k=18 B=2G in='final.sorted.bam' j=10
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
```

## Find soft clipped sequences
```
python find-soft-clip.py
```
Found 8 sequences (eyeballing soft-clip.txt)

### Local BLAST against plasmid
```
makeblastdb -in plasmid/laco-i-scei-teto_curated.fasta -dbtype nucl -out plasmid/plasmid
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
AGCTCTGCTTATATAGGCCTCCCACCGTACACGCCTACCTCGACCCGG
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

>NB551658:237:HCW3VBGXN:3:22403:13997:6467 right soft (1, plasmid, left is repeat)
ATCACTGATAGGGAGTGGTAAACTCGAC ATTACCCTGTTATCCCTACCCTCGAGGGATCGACTCTAGAGGCGCCGA
NB551658:237:HCW3VBGXN:3:22403:13997:6467 99  LacO-I-SceI-TetO_curated  1071  1 76M = 1252  257 CCGGTGTCTTCTATGGAGGTCAAAACAGCGTGGATGGCGTCTCCAGGCGATCTGACGGTTCACTAAACGAGCTCTG  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  AS:i:152  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:76 YS:i:56 YT:Z:CP
NB551658:237:HCW3VBGXN:3:22403:13997:6467 147 LacO-I-SceI-TetO_curated  1252  1 28M48S  = 1071  -257  ATCACTGATAGGGAGTGGTAAACTCGACATTACCCTGTTATCCCTACCCTCGAGGGATCGACTCTAGAGGCGCCGA  EEE<EEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:56 XS:i:152  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:28 YS:i:152  YT:Z:CP


>NB551658:237:HCW3VBGXN:2:21206:11540:9879 left soft (2, found contig)
GCTGACTAGAGGATCCCCGGGTACCGAGCTCGAATTCGGGG CCGCGGAGGCTGGATCGGTCCCGGTGTCTTCTATG
NB551658:237:HCW3VBGXN:2:21206:11540:9879 99  LacO-I-SceI-TetO_curated  1051  1 41S35M  = 1342  407 GCTGACTAGAGGATCCCCGGGTACCGAGCTCGAATTCGGGGCCGCGGAGGCTGGATCGGTCCCGGTGTCTTCTATG  AAAAAEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEE  AS:i:70 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:35 YS:i:150  YT:Z:CP
NB551658:237:HCW3VBGXN:2:21206:11540:9879 147 LacO-I-SceI-TetO_curated  1342  1 75M = 1051  -407  GATAGGGAGTGGTAAACTCGACTTTCACTTTTCTCTATCACTGATAGGGAGTGGTAAACTCGACTTTCACTTTTC EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA AS:i:150  XS:i:150  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:70 YT:Z:CP


>NB551658:237:HCW3VBGXN:4:11612:14985:9924 left soft (3, plasmid, right is repeat)
TAGGCCTCCCACCGTACACGCCTACCTCGACCCGG TCGACGACTTTCACTTTTCTCTATCACTGATAGGGAGTGGT
NB551658:237:HCW3VBGXN:4:11612:14985:9924 99  LacO-I-SceI-TetO_curated  5013  11  35S41M  = 5284  382 TAGGCCTCCCACCGTACACGCCTACCTCGACCCGGTCGACGACTTTCACTTTTCTCTATCACTGATAGGGAGTGGT  AAAAAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  AS:i:82 XS:i:152  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:41 YS:i:152  YT:Z:CP
NB551658:237:HCW3VBGXN:4:11612:14985:9924 147 LacO-I-SceI-TetO_curated  5284  11  76M = 5013  -382  ATCCCTACCCTCGAGGGATCGACTCTAGAGGCGCCGAATTCCACAAATTGTTATCCGCTCACAATTCCACATGTGG  EEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:152  XS:i:116  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:76 YS:i:82 YT:Z:CP


>NB551658:237:HCW3VBGXN:3:11503:3566:12828 right soft (4, plasmid, start/end of fasta)
CTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTC TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTG
NB551658:237:HCW3VBGXN:3:11503:3566:12828 83  LacO-I-SceI-TetO_curated  17215 41  39M37S  = 17097 -194  CTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTG  EEEEEAE/EEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:78 XS:i:74 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:39 YS:i:152  YT:Z:CP
NB551658:237:HCW3VBGXN:3:11503:3566:12828 163 LacO-I-SceI-TetO_curated  17097 41  76M = 17215 194 TTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGA  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  AS:i:152  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:76 YS:i:78 YT:Z:CP


>NB551658:237:HCW3VBGXN:1:12304:22582:2673 right soft (5, not found)
ACTTTTCTCTATCACTGATAGGGAGTGGTAAACTCGA GTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCA
NB551658:237:HCW3VBGXN:1:12304:22582:2673 99  LacO-I-SceI-TetO_curated  1233  1 75M = 1368  211 TCGACTTTCACTTTTCTCTATCACTGATAGGGAGTGGTAAACTCGACTTTCACTTTTCTCTATCACTGATAGGGA AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEE6AEE/EEEE<EE AS:i:150  XS:i:150  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:74 YT:Z:CP
NB551658:237:HCW3VBGXN:1:12304:22582:2673 147 LacO-I-SceI-TetO_curated  1368  1 37M39S  = 1233  -211  ACTTTTCTCTATCACTGATAGGGAGTGGTAAACTCGAGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCA  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:74 XS:i:74 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:37 YS:i:150  YT:Z:CP


>NB551658:237:HCW3VBGXN:4:23402:26687:8949 right soft (6, not found)
GCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAG TCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACTG
NB551658:237:HCW3VBGXN:4:23402:26687:8949 69  LacO-I-SceI-TetO_curated  790 0 * = 790 0 AAAGGTTGGCTATAAAGAGGTCATCAGTATATGAAACAGCCCC AAAAAEAE/E/EEEEEAEE/EEEAEEE/EEEAAEE/EEAE6/6 YT:Z:UP
NB551658:237:HCW3VBGXN:4:23402:26687:8949 137 LacO-I-SceI-TetO_curated  790 22  36M40S  = 790 0 GCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGTCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACTG  AAAAAAEEAAAEAAE/EEEAEEAAEEEEEEEEAEEEEEEEAEEEEEEEAEEEEEEA<EEAAAEAE/AAAE/EAE/A  AS:i:72 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YT:Z:UP


>NB551658:237:HCW3VBGXN:3:13508:15800:9326 left soft (7, found contig)
TATAGGTTAATGTCATGATAATAATGGTTTCTTAGA CGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCC
NB551658:237:HCW3VBGXN:3:13508:15800:9326 99  LacO-I-SceI-TetO_curated  16089 41  36S40M  = 16323 346 TATAGGTTAATGTCATGATAATAATGGTTTCTTAGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCC  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEAEEA/EEEE  AS:i:80 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YS:i:152  YT:Z:CP
NB551658:237:HCW3VBGXN:3:13508:15800:9326 147 LacO-I-SceI-TetO_curated  16323 41  76M = 16089 -346  CTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCG  AEEEEEE/EEEEEEEEAEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:152  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:76 YS:i:80 YT:Z:CP


>NB551658:237:HCW3VBGXN:1:11310:9141:1212 right soft (8, plasmid, start/end of fasta)
AGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTC TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAA
NB551658:237:HCW3VBGXN:1:11310:9141:1212  83  LacO-I-SceI-TetO_curated  17218 36  36M40S  = 17090 -204  AGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAA  AEEE6EEEAEEEEEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  AS:i:72 XS:i:80 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YS:i:145  YT:Z:CP
NB551658:237:HCW3VBGXN:1:11310:9141:1212  163 LacO-I-SceI-TetO_curated  17090 36  76M = 17218 204 GACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTTATTTTTAATTTAAAAGGATCTAGGTGAAGATCC  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA/AEEEEEEAEEEEEEEEEE/EEE  AS:i:145  XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:42C33  YS:i:72 YT:Z:CP
```

## Find unmapped mates of soft clipped sequence 2, 5, 6 and 7
```
python find-unmapped-mate.py
```

### BLAST each group of sequences against GRCh38 primary
tsv output:
```
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/ACTCGA GTGTGT.fasta" -out "unmapped-mates/ACTCGA GTGTGT.tsv" -outfmt 6 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/ATGCAG TCGAGG.fasta" -out "unmapped-mates/ATGCAG TCGAGG.tsv" -outfmt 6 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/CTTAGA CGTCAG.fasta" -out "unmapped-mates/CTTAGA CGTCAG.tsv" -outfmt 6 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/TCGGGG CCGCGG.fasta" -out "unmapped-mates/TCGGGG CCGCGG.tsv" -outfmt 6 -task blastn -num_threads 10
```
text output:
```
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/ACTCGA GTGTGT.fasta" -out "unmapped-mates/ACTCGA GTGTGT.txt" -outfmt 0 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/ATGCAG TCGAGG.fasta" -out "unmapped-mates/ATGCAG TCGAGG.txt" -outfmt 0 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/CTTAGA CGTCAG.fasta" -out "unmapped-mates/CTTAGA CGTCAG.txt" -outfmt 0 -task blastn -num_threads 10
blastn -db GRCh38_primary/GRCh38_primary -query "unmapped-mates/TCGGGG CCGCGG.fasta" -out "unmapped-mates/TCGGGG CCGCGG.txt" -outfmt 0 -task blastn -num_threads 10
```
Output columns:
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```
Found two high confident junctions:
```
Soft clip junction 6 "ATGCAG TCGAGG" at chr11 5225524:5225455, eval=1.12E-11

NB551658:237:HCW3VBGXN:3:12504:7475:1216  73  LacO-I-SceI-TetO_curated  789 22  37M39S  = 789 0 CGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGTCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACT  AAAAAEEEEEEEEEEE/EEEEEEEAEEEEEEEEEEEEAEEEEEEAEEAEEEEEEEEEEAE6EEE/EAEEEEEEAE<  AS:i:74 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:37 YT:Z:UP
NB551658:237:HCW3VBGXN:3:12504:7475:1216  133 LacO-I-SceI-TetO_curated  789 0 * = 789 0 ACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCTAATAAAGGAAATTTATTTTCATTGCAATAGTGTGTTGGA //AAAEEE6AEEEEEEEEEEEEEEEEEEEEEAEE/EEEAEE/EEEEEAEAEEEEEEAEEEEE/EEEEEEEEEEEE YT:Z:UP

@NB551658:237:HCW3VBGXN:3:12504:7475:1216
CGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGTCGAGGGATCTCCATAAGAGAAGAGGGACAGCTATGACT
@NB551658:237:HCW3VBGXN:3:12504:7475:1216
ACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCTAATAAAGGAAATTTATTTTCATTGCAATAGTGTGTTGGA

>NC_000011.10 Homo sapiens chromosome 11, GRCh38.p13 Primary Assembly
Length=135086622

 Score = 73.4 bits (80),  Expect = 1e-11
 Identities = 58/70 (83%), Gaps = 0/70 (0%)
 Strand=Plus/Minus

Query  3        ATCATGAAGCCCCTTGAGCATCTGACTTCTGGCTAATAAAGGAAATTTATTTTCATTGCA  62
                || ||||||  |||||||||||||  ||||| ||||||||  | ||||||||||||||||
Sbjct  5225524  ATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCA  5225465

Query  63       ATAGTGTGTT  72
                ||  ||| ||
Sbjct  5225464  ATGATGTATT  5225455
```
```
Soft clip junction 2 "TCGGGG CCGCGG" at chr11 5226513:5226589, eval=2.53E-07; 5226517:5226590, eval=8.82E-07

NB551658:237:HCW3VBGXN:4:13502:7832:5881  101 LacO-I-SceI-TetO_curated  1051  0 * = 1051  0 TCCATATAACATGAATTTTACAATAGCGAAAAAGAAAGAACAATCAAGGGTCCCCAAACTCACCCTGAAGTTCTCA  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EAEEEEEEEEEEEE  YT:Z:UP
NB551658:237:HCW3VBGXN:4:13502:7832:5881  153 LacO-I-SceI-TetO_curated  1051  28  25S51M  = 1051  0 CCGGGTACCGAGCTCGAATTCGGGGCCGCGGAGGCTGGATCGGTCCCGGTGTCTTCTATGGAGGTCAAAACAGCGT  EEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEAAAAA  AS:i:102  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:51 YT:Z:UP

@NB551658:237:HCW3VBGXN:4:13502:7832:5881
TCCATATAACATGAATTTTACAATAGCGAAAAAGAAAGAACAATCAAGGGTCCCCAAACTCACCCTGAAGTTCTCA
@NB551658:237:HCW3VBGXN:4:13502:7832:5881
ACGCTGTTTTGACCTCCATAGAAGACACCGGGACCGATCCAGCCTCCGCGGCCCCGAATTCGAGCTCGGTACCCGG

>NC_000011.10 Homo sapiens chromosome 11, GRCh38.p13 Primary Assembly
Length=135086622

 Score = 59.0 bits (64),  Expect = 3e-07
 Identities = 61/77 (79%), Gaps = 5/77 (6%)
 Strand=Plus/Plus

Query  5        TATAACATGAATTTTACAATAGCGAAAAAG---AAAGAACA-ATCAAGGGTCCCCAA-AC  59
                ||| ||||||| || || ||||  || |||   |||||| | |||||| |||||  | ||
Sbjct  5226513  TATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGAC  5226572

Query  60       TCACCCTGAAGTTCTCA  76
                |||||||||||||||||
Sbjct  5226573  TCACCCTGAAGTTCTCA  5226589


NB551658:237:HCW3VBGXN:4:13607:12002:2387 101 LacO-I-SceI-TetO_curated  1051  0 * = 1051  0 TAACATGAATTTTACAATAGCGAAAAAGAAAGAACAATCAAGGGTCCCCAAACTCACCCTGAAGTTCTCAGCTCTA  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEE<EEEEEEEEAEEEEEEEEE  YT:Z:UP
NB551658:237:HCW3VBGXN:4:13607:12002:2387 153 LacO-I-SceI-TetO_curated  1051  24  34S42M  = 1051  0 AGAGGATCCCCGGGTACCGAGCTCGAATTCGGGGCCGCGGAGGCTGGATCGGTCCCGGTGTCTTCTATGGAGGTCA  EEEEEEEEEAEEEEEEEEEEEAEEEEEEEEEEEEEEAEEEEEAEEEEEAAAEEEEEEAEEEEEEEEEEEEEAA6AA  AS:i:84 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:42 YT:Z:UP

>NC_000011.10 Homo sapiens chromosome 11, GRCh38.p13 Primary Assembly
Length=135086622

 Score = 58.1 bits (63),  Expect = 9e-07
 Identities = 59/74 (80%), Gaps = 5/74 (7%)
 Strand=Plus/Plus

Query  3        ACATGAATTTTACAATAGCGAAAAAG---AAAGAACA-ATCAAGGGTCCCCAA-ACTCAC  57
                ||||||| || || ||||  || |||   |||||| | |||||| |||||  | ||||||
Sbjct  5226517  ACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCAC  5226576

Query  58       CCTGAAGTTCTCAG  71
                ||||||||||||||
Sbjct  5226577  CCTGAAGTTCTCAG  5226590
```
Structure in genome:
```

```