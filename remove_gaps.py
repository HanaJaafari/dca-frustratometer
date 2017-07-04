#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Removes sequence with too many gaps from alignments')
parser.add_argument('-g','--gapPercent', type=float, default=0.2,
                    help='The maximum proportion of gaps on any sequence in the alignment')
parser.add_argument('--remove',choices=['gap_columns','all_insertions','None'],default='all_insertions',
                    help='Remove all insertions, only empty columns or keep all columns')
parser.add_argument('alignment',help='Input alignment')
parser.add_argument('output', help='Output File')


args = parser.parse_args()
print(args)

gapPercent=args.gapPercent
pfam_Alignment=args.alignment
out_file=args.output
Remove_gap_only_columns=True if args.remove=='gap_columns' else False
Remove_all_insertions=True if args.remove=='all_insertions' else False

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

#Parse alignment
alignment=AlignIO.read(pfam_Alignment,'fasta')
print 'Starting sequences: ', len(alignment)

#Calculate length of query
Query_len=len(str(alignment[0].seq).replace('.',''))

#Calculate max number of allowed gaps
max_gaps=gapPercent*Query_len

#Select the sequences
selected_sequences=[]
for r in alignment:
    if r.seq.count('-')<max_gaps:
        selected_sequences+=[r]

print 'Selected sequences: ', len(selected_sequences)
new_alignment = MultipleSeqAlignment(selected_sequences)


if Remove_gap_only_columns:
    l=len(selected_sequences)
    ee=0
    final_alignment=[]
    for i in range(new_alignment.get_alignment_length()):
        s=new_alignment[:,i]
        if s.count('.')==l:
            ee+=1
        else:
            if type(final_alignment)==list:
                final_alignment=new_alignment[:,i:i+1]
            else:
                final_alignment+=new_alignment[:,i:i+1]
    print 'Removed gap only columns: ', ee
    print AlignIO.write(final_alignment,out_file,'fasta')
    print final_alignment
elif Remove_all_insertions:
    import re
    handle=open(out_file,'w+')
    for a in new_alignment:
        handle.write('>%s\n'%a.description)
        handle.write('%s\n'%re.sub('[a-z,.]','',str(a.seq)))
    handle.close()
else:
    final_alignment=new_alignment
    print AlignIO.write(final_alignment,out_file,'fasta')
    print final_alignment


