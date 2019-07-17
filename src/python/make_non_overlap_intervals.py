'''
To handle GATK Depthh of Coverage quirk that merges adjacent intervals, we leave a gap of 1bp in between windows.
So that we can feed the intervals to xHMM
'''

import sys
import csv
import Bio.SeqIO


REF_FASTA = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7/3D7_V3.fasta"
WIN_SIZE = 300


if len(sys.argv) < 2:
    raise ValueError("Must specify arguments  <bed output file to write to>")


bed_output = sys.argv[1]


ref = Bio.SeqIO.to_dict(Bio.SeqIO.parse(REF_FASTA, "fasta"))

with open(bed_output, 'w') as fh_out:
    writer = csv.DictWriter(fh_out, fieldnames=["chr", "start", "end"], delimiter="\t")
    for rec in Bio.SeqIO.parse(REF_FASTA, "fasta"):
        for start_0based in range(0, len(rec.seq), WIN_SIZE):
            end_1based = start_0based + WIN_SIZE - 2
            if end_1based > len(rec.seq):
                end_1based =  len(rec.seq)
            outrow = {"chr": rec.id, "start": start_0based, "end": end_1based}
            writer.writerow(outrow)
            
            
            
    
    
        