import pysam
from tqdm import tqdm

def revcomp(seq):
    '''
    give the reverse complement of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        reverse complement of input
    '''
    if isinstance(seq, str):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', 
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4'}
        bases = [complement[base] for base in seq]
        bases.reverse()
        return ''.join(bases)
    elif isinstance(seq, list):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', 
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4', 
        0:1, 1:0, 2:3, 3:2, 4:4}
        bases = [complement[base] for base in seq]
        bases.reverse()
        return bases
    else:
        return -1


TH = 0.3 # relative length of soft clipped sequence

samfile1 = pysam.AlignmentFile('final.sorted.bam', 'rb')
samfile2 = pysam.AlignmentFile('candidates.bam', 'wb', template=samfile1)

fa = '' # fasta
sc = [] # soft clipped sequences
for read in tqdm(samfile1.fetch(until_eof=True)):
    if read.is_unmapped:
        continue
    qlen = read.query_length
    alen = read.query_alignment_length
    astart = read.query_alignment_start
    astop = read.query_alignment_end

    if astart > TH*qlen:
        # left soft
        samfile2.write(read)
        s = read.query_sequence[:astart]
        sc.append(revcomp(s))
        fa += f'>{read.query_name} left soft\n{s} {read.query_sequence[astart:]}\n'
    elif astop < (1-TH)*qlen:
        # right soft
        samfile2.write(read)
        s = read.query_sequence[astop:]
        sc.append(s)
        fa += f'>{read.query_name} right soft\n{read.query_sequence[:astop]} {s}\n'

sc = '\n'.join(sorted(sc))

with open('candidates.fasta', 'w') as f:
    f.write(fa)
with open('soft-clip.txt', 'w') as f:
    f.write(sc)