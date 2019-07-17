# Snakemake file for generating CNV calls in whole genome shotgun, whole genome amplified samples used in lancet paper for ppq mfl resistance

import pandas as pd

# PPQ_FOFN = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/resource/ppq_resistance.fofn"
PPQ_OUTPUT_DIR='/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/input/ppq_resistance'
DEPTH_COV_OUTPUT_DIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/depth_of_cov'

REFERENCE_FASTA='/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7/3D7_V3.fasta'

INTERVAL_LIST='/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7_v3.window300bp.bed'

PPQ_SAMPLES_LIST = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/resource/ppq_samples.txt"
PYSAMDIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/ppq_resistance/pysamstats_gc'
HMMDIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/ppq_resistance/gaussian_hmm_cnv'
# Additional piperaquine samples that show copy number variation
ppq_resist_samples = pd.read_csv(PPQ_SAMPLES_LIST, header=None, names = ["sample"])
ppq_resist_samples["library"] = "unknown"
ppq_resist_samples["lane"] = "unknown"

samples = pd.concat([ppq_resist_samples])


rule all:
    input:
        expand(DEPTH_COV_OUTPUT_DIR + '/{sample}.{library}.{lane}.win300.sample_interval_summary', zip, sample=samples["sample"], library=samples["library"], lane=samples["lane"]),
        expand(PYSAMDIR + '/{sample}.{chrom}.pysamstats_gc.txt', sample=samples["sample"], chrom=["Pf3D7_10_v3", "Pf3D7_10_v3", "Pf3D7_07_v3", "Pf3D7_03_v3", "Pf3D7_13_v3", "Pf3D7_11_v3", "Pf3D7_14_v3", "Pf3D7_09_v3", "Pf3D7_01_v3", "Pf3D7_12_v3", "Pf3D7_08_v3", "Pf3D7_06_v3", "Pf3D7_04_v3", "Pf3D7_02_v3", "M76611", "PFC10_API_IRAB"]),
        expand(HMMDIR + '/{sample}.{chrom}.cnv.tsv', sample=samples["sample"], chrom=["Pf3D7_10_v3", "Pf3D7_10_v3", "Pf3D7_07_v3", "Pf3D7_03_v3", "Pf3D7_13_v3", "Pf3D7_11_v3", "Pf3D7_14_v3", "Pf3D7_09_v3", "Pf3D7_01_v3", "Pf3D7_12_v3", "Pf3D7_08_v3", "Pf3D7_06_v3", "Pf3D7_04_v3", "Pf3D7_02_v3", "M76611", "PFC10_API_IRAB"])



rule depth_of_cov:
    input:
        bam = PPQ_OUTPUT_DIR + '/{sample}.bam'
    output:
        DEPTH_COV_OUTPUT_DIR + '/{sample}.{library}.{lane}.win300.sample_interval_summary'
    shell:
        """
        set +u
        set -e
        set -x



        echo SAMPLE={wildcards.sample}  LIBRARY={wildcards.library} LANE={wildcards.lane} BAM={input.bam}

        samtools index {input.bam}

        java -jar /software/vertres/bin-external/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
                 -T DepthOfCoverage \
                 -R {REFERENCE_FASTA} \
                 -o {DEPTH_COV_OUTPUT_DIR}/{wildcards.sample}.{wildcards.library}.{wildcards.lane}.win300 \
                 -I {input.bam} \
                 -L {INTERVAL_LIST} \
                  --countType COUNT_FRAGMENTS         \
                 --includeRefNSites \
                 --minBaseQuality 0 \
                 --minMappingQuality 20 \
                 --start 1 --stop 1000 --nBins 999 \
                 --omitDepthOutputAtEachBase \
                 --omitLocusTable \
                 --omitPerSampleStats
        set -u
        """

rule pysam:
    input:
        bam = PPQ_OUTPUT_DIR + '/{sample}.bam'
    output:
        PYSAMDIR + '/{sample}.{chrom}.pysamstats_gc.txt'
    shell:
        """
        /nfs/team112/pythonbrew/pythons/Python-2.7.2.farm3/bin/pysamstats -t coverage_gc -c {wildcards.chrom} \
          --output {PYSAMDIR}/{wildcards.sample}.{wildcards.chrom}.pysamstats_gc.txt  \
          -f {REFERENCE_FASTA} \
          {input.bam}
        """


rule hmm:
    input:
        bam = PPQ_OUTPUT_DIR + '/{sample}.bam'
    output:
        HMMDIR + '/{sample}.{chrom}.cnv.tsv'
    shell:
        """
        python $HOME/gitrepo/pf_swga_cnv/src/python/gaussian_hmm_cnv.py -pysam_gc {PYSAMDIR}/{wildcards.sample}.{wildcards.chrom}.pysamstats_gc.txt -out_tsv {HMMDIR}/{wildcards.sample}.{wildcards.chrom}.cnv.tsv
        """
