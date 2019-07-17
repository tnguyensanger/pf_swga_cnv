# Snakemake file for generating CNV calls in whole genome shotgun, whole genome amplified samples

import pandas as pd

WGS_FOFN = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/resource/pf_swga_cnv_wgs_indel_realign_bqsr.output.fofn"
BQSR_OUTPUT_DIR='/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/pf_swga_cnv_wgs_indel_realign_bqsr'
DEPTH_COV_OUTPUT_DIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/depth_of_cov'

# WGA samples that are paired with the sWGA samples that are under investigation
wga_pair_samples = pd.read_table(WGS_FOFN).set_index("sample", drop=False)


samples = pd.concat([wga_pair_samples])

rule all:
    output:
        expand(DEPTH_COV_OUTPUT_DIR + '/{sample}.{library}.{lane}', zip, sample=samples["sample"], library=samples["library"], lane=samples["lane"])



rule depth_of_cov:
    input:
        bam = BQSR_OUTPUT_DIR + '/{sample}.bam'
    output:
        DEPTH_COV_OUTPUT_DIR + '/{sample}.{library}.{lane}.sample_interval_statistics',
        DEPTH_COV_OUTPUT_DIR + '/{sample}.{library}.{lane}.sample_interval_summary'
    shell:
        """
        set -e
        set -x

        # vrpipe-fileinfo --setup pf_swga_plex_indel_realign_bqsr --metadata sample,library,lane --step "gatk_print_reads_gatk3_v2|recalibrated_bam_files" > /lustre/scratch118/malaria/team112/pipelines/setups/pf_swga/resource/pf_swga_plex_indel_realign_bqsr.output.fofn
        # vrpipe-fileinfo --setup pf_swga_cnv_3D7_indel_realign_bqsr --metadata sample,library,lane --step "gatk_print_reads_gatk3_v2|recalibrated_bam_files" > /lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_3D7/resource/pf_swga_cnv_3D7_indel_realign_bqsr.output.fofn
        # vrpipe-fileinfo --setup pf_swga_cnv_wgs_indel_realign_bqsr --metadata sample,library,lane --step "gatk_print_reads_gatk3_v2|recalibrated_bam_files" > /lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/resource/pf_swga_cnv_wgs_indel_realign_bqsr.output.fofn


        REFERENCE_FASTA=/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7/3D7_V3.fasta

        BAM={wildcards.bam}
        SAMPLE={wildcards.sample}
        LIBRARY={wildcards.library}
        LANE={wildcards.lane}

        echo SAMPLE=$SAMPLE  LIBRARY=$LIBRARY LANE=$LANE BAM=$BAM

        java -jar /software/vertres/bin-external/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
                 -T DepthOfCoverage \
                 -R $REFERENCE_FASTA \
                 -o {DEPTH_COV_OUTPUT_DIR}/$SAMPLE.$LIBRARY.$LANE \
                 -I $BAM \
                 --countType COUNT_FRAGMENTS \
                 --includeRefNSites \
                 --minBaseQuality 0 \
                 --minMappingQuality 20 \
                 --start 1 --stop 7000 --nBins 6999
        """
