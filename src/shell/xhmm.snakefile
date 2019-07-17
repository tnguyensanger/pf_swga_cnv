# Snakemake file for generating CNV calls in whole genome shotgun, whole genome amplified samples used in lancet paper for ppq mfl resistance

import pandas as pd

DEPTH_COV_OUTPUT_DIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/depth_of_cov'

XHMM_OUTPUT_DIR = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv_wgs/output/xhmm'

REFERENCE_FASTA = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7/3D7_V3.fasta'

INTERVAL_LIST = '/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/resource/3D7_v3.window300bp.bed'



rule all:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.xcnv.win300.tsv',
        XHMM_OUTPUT_DIR + '/xhmm_test.centered.RD.txt'


rule find_sample_interval_summary_list:
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.sample_interval_summary.list.txt'
    shell:
        """

        # Get all the WGS windowed
        find {DEPTH_COV_OUTPUT_DIR} -name "*.win300.sample_interval_summary"  > {XHMM_OUTPUT_DIR}/xhmm_test.sample_interval_summary.list.txt

        """

rule merge_sample_interval_summary_list:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.sample_interval_summary.list.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.RD.txt'
    shell:
        """

        # Merge the gatk sample interval summaries
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --mergeGATKdepths -o {XHMM_OUTPUT_DIR}/xhmm_test.RD.txt \
        --GATKdepthsList {XHMM_OUTPUT_DIR}/xhmm_test.sample_interval_summary.list.txt

        """


rule interval_gc:
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.locus_GC.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.extreme_gc_targets.txt',
    shell:
        """

        # Optionally, run GATK to calculate the per-target GC content and create a list of the targets with extreme GC content:
        java -jar /software/vertres/bin-external/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
                 -T GCContentByInterval \
                 -R {REFERENCE_FASTA} \
                 -o {XHMM_OUTPUT_DIR}/xhmm_test.locus_GC.txt \
                 -L {INTERVAL_LIST}

        cat {XHMM_OUTPUT_DIR}/xhmm_test.locus_GC.txt | awk '{{if ($2 < 0.1 || $2 > 0.9) print $1}}' > {XHMM_OUTPUT_DIR}/xhmm_test.extreme_gc_targets.txt

        """


rule filter_mean_rd:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.locus_GC.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.extreme_gc_targets.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.RD.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.filtered_centered.RD.txt'
    shell:
        """

        # Filters samples and targets and then mean-centers the targets:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --matrix -r {XHMM_OUTPUT_DIR}/xhmm_test.RD.txt --centerData --centerType target \
        -o {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt \
        --outputExcludedTargets {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeTargets {XHMM_OUTPUT_DIR}/xhmm_test.extreme_gc_targets.txt \
        --minTargetSize 10 --maxTargetSize 10000 \
        --minMeanTargetRD 10 --maxMeanTargetRD 1000 \
        --minMeanSampleRD 25 --maxMeanSampleRD 1000 \
        --maxSdSampleRD 700

        """


rule mean_rd:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.locus_GC.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.centered.RD.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.centered.RD.txt.filtered_targets.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.centered.RD.txt.filtered_samples.txt'
    shell:
        """

        # Does NOT Filters samples and targets but does then mean-centers the targets:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --matrix -r {XHMM_OUTPUT_DIR}/xhmm_test.RD.txt --centerData --centerType target \
        -o {XHMM_OUTPUT_DIR}/xhmm_test.centered.RD.txt \
        --outputExcludedTargets {XHMM_OUTPUT_DIR}/xhmm_test.centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples {XHMM_OUTPUT_DIR}/xhmm_test.centered.RD.txt.filtered_samples.txt

        """



rule pca_on_filtered_mean_rd:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.filtered_centered.RD.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.RD_PCA.PC_LOADINGS.txt'
    shell:
        """

        # Runs PCA on mean-centered data:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --PCA -r {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt --PCAfiles {XHMM_OUTPUT_DIR}/xhmm_test.RD_PCA

        """


rule norm_by_pc:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.filtered_centered.RD.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.RD_PCA.PC_LOADINGS.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.txt'
    shell:
        """

        # Normalizes mean-centered data using PCA information:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --normalize -r {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt --PCAfiles {XHMM_OUTPUT_DIR}/xhmm_test.RD_PCA \
        --normalizeOutput {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.txt \
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

        """

rule zscore:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.txt'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt'
    shell:
        """

        # Filters and z-score centers (by sample) the PCA-normalized data:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --matrix -r {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.txt --centerData --centerType sample --zScoreData \
        -o {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt \
        --outputExcludedTargets {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --outputExcludedSamples {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        --maxSdTargetRD 30

        """



rule filter_orig_rd:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.RD.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt',
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.same_filtered.RD.txt'
    shell:
        """

        # Filters original read-depth data to be the same as filtered, normalized data:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --matrix -r {XHMM_OUTPUT_DIR}/xhmm_test.RD.txt \
        --excludeTargets {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt.filtered_targets.txt \
        --excludeTargets {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --excludeSamples {XHMM_OUTPUT_DIR}/xhmm_test.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeSamples {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        -o {XHMM_OUTPUT_DIR}/xhmm_test.same_filtered.RD.txt

        """




rule call_cnv:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt',
        XHMM_OUTPUT_DIR + '/xhmm_test.same_filtered.RD.txt',
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.xcnv'
    shell:
        """

        # Discovers CNVs in normalized data:
        /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --discover -p {XHMM_OUTPUT_DIR}/xhmm_test.params.txt \
        -r {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt -R {XHMM_OUTPUT_DIR}/xhmm_test.same_filtered.RD.txt \
        -c {XHMM_OUTPUT_DIR}/xhmm_test.xcnv -a {XHMM_OUTPUT_DIR}/xhmm_test.aux_xcnv -s {XHMM_OUTPUT_DIR}/xhmm_test


        """



rule output_cnv_csv:
    input:
        XHMM_OUTPUT_DIR + '/xhmm_test.xcnv'
    output:
        XHMM_OUTPUT_DIR + '/xhmm_test.xcnv.win300.tsv'
    shell:
        """
        python $HOME/gitrepo/pf_swga_cnv/src/python/xhmm_interval_to_win.py {XHMM_OUTPUT_DIR}/xhmm_test.xcnv {XHMM_OUTPUT_DIR}/xhmm_test.xcnv.win300.tsv

        """



# # Genotypes discovered CNVs in all samples:
# /nfs/team112/tn6/programs/xhmm/statgen-xhmm-cc14e528d909/xhmm --genotype -p {XHMM_OUTPUT_DIR}/xhmm_test.params.txt \
# -r {XHMM_OUTPUT_DIR}/xhmm_test.PCA_normalized.filtered.sample_zscores.RD.txt -R {XHMM_OUTPUT_DIR}/xhmm_test.same_filtered.RD.txt \
# -g {XHMM_OUTPUT_DIR}/xhmm_test.xcnv -F $REFERENCE_FASTA \
# -v {XHMM_OUTPUT_DIR}/xhmm_test.vcf




# /nfs/team112/tn6/programs/plinkseq/plinkseq-0.10/pseq  . loc-intersect --group refseq --locdb ./RefSeq.locdb --file {INTERVAL_LIST} --out {XHMM_OUTPUT_DIR}/xhmm_test.annotated_targets.refseq
