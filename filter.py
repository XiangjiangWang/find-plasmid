import pysam
from tqdm import tqdm

samfile1 = pysam.AlignmentFile('mu.sorted.bam', 'rb')
samfile2 = pysam.AlignmentFile('filtered.bam', 'wb', template=samfile1)

for read in tqdm(samfile1.fetch(until_eof=True)):
    if read.is_proper_pair:
        # concordant
        continue
    if read.is_mapped and read.mate_is_mapped:
        # discordant
        continue
    elif read.is_mapped and read.mate_is_unmapped:
        samfile2.write(read)
    else:
        # read.is_unmapped
        samfile2.write(read)
samfile1.close()
samfile2.close()
