import pysam
from tqdm import tqdm

junctions = [
    'TCGGGG CCGCGG', 
    'ACTCGA GTGTGT', 
    'ATGCAG TCGAGG', 
    'CTTAGA CGTCAG'
]

for j in junctions:
    seq = j.replace(' ', '')
    qname = []
    samfile1 = pysam.AlignmentFile('soft-clip.sam', 'r')
    for read in tqdm(samfile1.fetch(until_eof=True)):
        if seq in read.query_sequence and read.mate_is_unmapped:
            qname.append(read.query_name)

    samfile1 = pysam.AlignmentFile('soft-clip.sam', 'r')
    with open(f'unmapped-mates/{j}.fasta', 'w') as f:
        for read in tqdm(samfile1.fetch(until_eof=True)):
            if read.query_name in qname and read.is_unmapped and len(read.query_sequence)>50:
                f.write(f'>{read.query_name}\n{read.query_sequence}\n')
