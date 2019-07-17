'''
SAMPLE  CNV     INTERVAL        KB      CHR     MID_BP  TARGETS NUM_TARG        Q_EXACT Q_SOME  Q_NON_DIPLOID   Q_START Q_STOP  MEAN_RD MEAN_ORIG_RD
PF1155-CW       DUP     Pf3D7_05_v3:209401-217798       8.40    Pf3D7_05_v3     213599  443..464        22      3       99      99      4       3       3.84    180.30
PF1155-CW       DUP     Pf3D7_05_v3:499801-506998       7.20    Pf3D7_05_v3     503399  1306..1329      24      8       99      99      9       10      3.23    218.70
PF1155-CW       DUP     Pf3D7_05_v3:1063501-1065298     1.80    Pf3D7_05_v3     1064399 2925..2930      6       14      55      55      10      3       3.33    197.67
PF1128-CW       DUP     Pf3D7_05_v3:113101-115198       2.10    Pf3D7_05_v3     114149  184..190        7       13      71      71      2       2       3.20    252.76
PF1128-CW       DEL     Pf3D7_05_v3:611101-618598       7.50    Pf3D7_05_v3     614849  1635..1650      16      3       99      99      4       2       -2.76   180.50
'''
import csv

XHMM_CNVC = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/output/xhmm/xhmm_test.xcnv"
XHMM_CNVC_WINDOWIZED = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/output/xhmm/xhmm_test.xcnv.win300.tsv"
with  open(XHMM_CNVC, 'r') as fh_in, open(XHMM_CNVC_WINDOWIZED, 'w') as fh_out:
    reader = csv.DictReader(fh_in, delimiter = "\t")
    writer = csv.DictWriter(fh_out, delimiter="\t", fieldnames=["sample", "chrom", "window", "copy_number"])
    writer.writeheader()

    for inrow in reader:
        sample = inrow["SAMPLE"]
        interval = inrow["INTERVAL"]

        copy_number = 1
        if inrow["CNV"] == "DUP":
            copy_number = 2
        else:
            copy_number =  0

        chromo, pos_range = interval.split(":")
        range_start, range_end = pos_range.split("-")
        range_start = int(range_start)
        range_end = int(range_end)
        win_start = 300 * (range_start // 300 )
        for i in range(win_start, int(range_end), 300):
            outrow = {"sample": sample,
                      "chrom": chromo,
                      "window": i,
                      "copy_number": copy_number}
            writer.writerow(outrow)
