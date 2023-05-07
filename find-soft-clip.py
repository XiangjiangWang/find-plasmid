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
samfile2 = pysam.AlignmentFile('soft-clip.sam', 'w', template=samfile1)

qname = {} # names of read pairs to output
for read in tqdm(samfile1.fetch(until_eof=True)):
    if read.is_unmapped:
        continue
    qlen = read.query_length
    alen = read.query_alignment_length
    astart = read.query_alignment_start
    astop = read.query_alignment_end

    # cautious when both reads have soft clipped sequences
    if astart > TH*qlen:
        # left soft
        s = read.query_sequence[:astart]
        fa = f'>{read.query_name} left soft\n{s} {read.query_sequence[astart:]}\n' # add a space between human and plasmid
        qname[read.query_name] = [None, None, revcomp(s), fa] # r1, r2, soft clipped sequence, fasta record
    elif astop < (1-TH)*qlen:
        # right soft
        s = read.query_sequence[astop:]
        fa = f'>{read.query_name} right soft\n{read.query_sequence[:astop]} {s}\n' # add a space between human and plasmid
        qname[read.query_name] = [None, None, s, fa] # r1, r2, soft clipped sequence, fasta record

# sort read pairs with soft clipped sequences
qname = {k: v for k, v in sorted(qname.items(), key=lambda x: x[1][2])}

# fetch read pairs in qname
samfile1.reset()
for read in tqdm(samfile1.fetch(until_eof=True)):
    try:
        v = qname[read.query_name]
        if read.is_read1:
            v[0] = read
        elif read.is_read2:
            v[1] = read
    except KeyError:
        continue
samfile1.close()

# write read pairs
all_sc = ''
all_fa = ''
for n, (r1, r2, sc, fa) in qname.items():
    if r1 is not None:
        samfile2.write(r1)
    else:
        print(f'r1 does not exist: {n}')
    if r2 is not None:
        samfile2.write(r2)
    else:
        print(f'r2 does not exist: {n}')
    all_sc += sc+'\n'
    all_fa += fa
samfile2.close()

with open('soft-clip.fasta', 'w') as f:
    f.write(all_fa)
with open('soft-clip.txt', 'w') as f:
    f.write(all_sc)