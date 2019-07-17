
import os
import os
import glob
import re
import pandas
import argparse

import numpy as np

from hmmlearn import hmm

PF_FIELD_GDNA_COV_DIR = "/lustre/scratch118/malaria/team112/personal/tn6/pf_swga_cnv/output"


def get_norm_gc(pysam_gc_stats_tsv):
    """
    Takes pysam gc coverage statistics as input.
    For each window of 300bp, bins the GC content.
    Normalies that coverage in that window by the average content in all windows with the same GC bin.
    """

    # Per-position GC content, read coverage.
    # We don't care about properly paired reads here, since there could be
    # inversions, deletions, insertions that alter orientation of read pairs.
    pysam_gc = pandas.read_csv(pysam_gc_stats_tsv, sep="\t")
    # Indicate the window start position
    pysam_gc["window"] =  300 * (pysam_gc.pos // 300)

    # Aggregate reads by window
    pysam_gc_win = pysam_gc.groupby("window").agg({'gc': {"gc_win_med": "median", 'gc_win_mean': "mean"},
                                              'reads_all': {"reads_all_win_med": "median",
                                                            'reads_all_win_mean': "mean"}})
    pysam_gc_win.columns = pysam_gc_win.columns.droplevel()
    # Label each window with the ave gc_bin content
    pysam_gc_win["gc_bin"] = 10 * (pysam_gc_win["gc_win_mean"] //10)

    # Find the mean reads in each GC Bin
    pysam_gc_per_pos_per_win = pysam_gc.reset_index().set_index("window").join(pysam_gc_win[["gc_bin"]])

    pysam_by_gc_bin = pysam_gc_per_pos_per_win.groupby("gc_bin").agg({'reads_all': "mean"})
    pysam_by_gc_bin = pysam_by_gc_bin.rename(columns={'reads_all': 'reads_all_gc_bin_mean'})

    pysam_per_win_per_gc_bin = pysam_gc_win.reset_index().set_index("gc_bin").join(pysam_by_gc_bin, rsuffix="_median")
    pysam_per_win_per_gc_bin = pysam_per_win_per_gc_bin.reset_index()
    pysam_per_win_per_gc_bin["window"] = pysam_per_win_per_gc_bin["window"].astype(float)
    pysam_per_win_per_gc_bin = pysam_per_win_per_gc_bin.sort_values("window")

    pysam_per_win_per_gc_bin["norm_reads_all"] = pysam_per_win_per_gc_bin["reads_all_win_mean"] / pysam_per_win_per_gc_bin["reads_all_gc_bin_mean"]



    return pysam_per_win_per_gc_bin



def fit_hmm(depth_normed,  # normalised coverage array
            transition_probability,  # probability of state transition
            variance,  # variance per copy
            variance_fixed,  # variance for the zero copy number state
            min_swga_copy_number=0,  # minimum copy number to consider in the model
            max_swga_copy_number=5,  # maximum copy number to consider in the model
            n_iter=0,  # number of iterations to perform when fitting the model
            params='st',  # parameters that can be changed through fitting.  Default start prob, trans matrix.
            init_params='',  # parameters that are initialised from the data
           ):
    """
    Fits a Gaussian Hidden Markov Model.
    """

    n_states = max_swga_copy_number - min_swga_copy_number

    # construct the transition matrix
    transmat = np.zeros((n_states, n_states))
    transmat[:] = transition_probability
    transmat[np.diag_indices(n_states)] = 1-((n_states-1)*transition_probability)


    # construct means and covariance
    means = np.array([[n] for n in range(min_swga_copy_number, max_swga_copy_number)])
    covars = np.array([[variance*n + variance_fixed] for n in range(min_swga_copy_number, max_swga_copy_number)])

    # setup HMM
    model = hmm.GaussianHMM(n_states,
                        covariance_type='diag',
                        n_iter=n_iter,
                        transmat_prior=transmat,
                        params=params,
                        init_params=init_params)
    model.means_ = means
    model.covars_ = covars

    # fit HMM
    obs = np.column_stack([depth_normed])
#     obs = depth_normed
    model.fit(obs)

    # predict hidden states
    h = model.predict(obs)

    return h



if __name__ == "__main__":
    TRANSITION_PROB = 0.0001
    VARIANCE_FIXED = 0.017120675716  # Based on variance in 3D7 which has no copy number variation


    parser = argparse.ArgumentParser(description='Calculates normalized coverage by GC content, and predicts the copy number using a Gaussian HMM')
    parser.add_argument('-pysam_gc', help="Sample pysam coverage by gc")
    parser.add_argument('-out_tsv', help="Output copy number tsv per window")
    args = parser.parse_args()

    cov_stats = get_norm_gc(pysam_gc_stats_tsv=args.pysam_gc)

    cov_stats["copy_number"] = fit_hmm(depth_normed =cov_stats["norm_reads_all"].values,
                          transition_probability=TRANSITION_PROB,
                          variance=cov_stats["norm_reads_all"].var(),
                          variance_fixed=VARIANCE_FIXED)

    cov_stats.to_csv(path_or_buf=args.out_tsv, sep='\t', na_rep='', header=True, index=False)
