'''
Converts 8 column blast format for BLAT output into putative binding site bed file.
Only exact and full matches allowed for now.
'''

import sys
import csv
import Bio.SeqIO


REF_FASTA = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7/3D7_V3.fasta"



if len(sys.argv) < 3:
    raise ValueError("Must specify arguments <blat output> <primer fasta> <bed output file to write to>")



blat_output_tsv = sys.argv[1]
primer_fasta = sys.argv[2]
bed_output = sys.argv[3]

rec_iter = Bio.SeqIO.parse(primer_fasta, "fasta")
primers = Bio.SeqIO.to_dict(rec_iter)

ref = Bio.SeqIO.to_dict(Bio.SeqIO.parse(REF_FASTA, "fasta"))

primer_start_chr_pos_0based = []
with open(blat_output_tsv, 'r') as fh_in:
    reader = csv.DictReader(fh_in, fieldnames=["query_id", "subject_id", "perc_identity", "align_len", "mismatches", "gap_opens", 
                                        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score" ], delimiter="\t")
    
    
    
    for inrow in reader:
        query_id = inrow["query_id"]
        mismatches = int(inrow["mismatches"])
        gap_opens = int(inrow["gap_opens"])
        perc_identity = float(inrow["perc_identity"])
        align_len = int(inrow["align_len"])
        
        if query_id not in primers:
            raise ValueError("Can't find primer " + str(query_id))
        
        query_primer = primers[query_id]
        if (len(query_primer.seq) != align_len or mismatches > 0 or gap_opens > 0 or perc_identity < 100):
            continue
        
        
        primer_start_chr_pos_0based.append((inrow["subject_id"], int(inrow["s_start"]) - 1))


primer_start_chr_pos_0based.sort(key=lambda tup: (tup[0], tup[1]))
print (primer_start_chr_pos_0based)

last_primer_start_chr  = None
last_primer_start_pos_0based = None
last_target_start_pos_0based = None
last_target_end_pos_1based = None

with open(bed_output, 'w') as fh_out:
    writer = csv.DictWriter(fh_out, fieldnames=["chr", "start", "end"], delimiter="\t")
    for primer_start_chr, primer_start_pos_0based in primer_start_chr_pos_0based:

        if last_primer_start_chr:

            # Bug: https://groups.google.com/a/broadinstitute.org/forum/#!topic/xhmm-users/h-ppA_QXvG4
            # Need to make > 1bp gap between intervals
            if last_primer_start_chr == primer_start_chr:
                last_target_end_pos_1based = primer_start_pos_0based - 2
            else:
                ref_chr = ref[last_primer_start_chr].seq
                last_target_end_pos_1based = len(ref_chr)
        
            # Bug: https://groups.google.com/a/broadinstitute.org/forum/#!topic/xhmm-users/h-ppA_QXvG4
            if last_target_end_pos_1based - last_target_start_pos_0based  <= 1:  # Need to make > 1bp gap between intervals
                continue
            
            
            
            
            outrow = {"chr": last_primer_start_chr, "start": last_target_start_pos_0based, "end": last_target_end_pos_1based}
            writer.writerow(outrow)
        
        if last_primer_start_chr != primer_start_chr:
            last_target_end_pos_1based = None
            
        last_primer_start_chr = primer_start_chr
        last_target_start_pos_0based = primer_start_pos_0based
    
        
            
            
    
    
        