"""
This script takes aligned WGS and amplicon data and produces a plot.
"""
import sys
import os
import pysam
import vcf
import string
import itertools
import argparse
import numpy as np
from operator import itemgetter
from itertools import groupby
from collections import Counter, defaultdict
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from sonLib.bioio import popenCatch, fastaRead, system
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from lib.general_lib import format_ratio

bases = {"A", "T", "G", "C", "a", "t", "g", "c"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--wgs_bam", required=True)
    parser.add_argument("--consensus_vcf", default="/hive/users/ifiddes/amplicon_analysis_with_new_suns/all_samples_reference_free/index/Notch2NL_SUN_UniqueIndels_ConsensusRef.vcf.gz")
    parser.add_argument("--index_path", default="/hive/users/ifiddes/amplicon_analysis_with_new_suns/all_samples_reference_free/index/notch2_aligned_consensus.fasta")
    parser.add_argument("--ranges", default="/hive/users/ifiddes/amplicon_analysis_with_new_suns/all_samples_reference_free/index/amplicons_ispcr_consensus.bed")
    parser.add_argument("--outdir", default="plots")
    return parser.parse_args()


def analyze_amplicon(bam_path, vcf_handle, ranges):
    results = {}
    for r in ranges:
        range_results = defaultdict(list)
        name = r[3]
        start = int(r[1])
        end = int(r[2])
        vcf_recs = list(vcf_handle.fetch(r[0], start, end))
        for vcf_rec in vcf_recs:
            if vcf_rec.is_indel:
                continue
            pos_str = "{0}:{1}-{1}".format(vcf_rec.CHROM, vcf_rec.POS)
            mpileup_rec = popenCatch("samtools mpileup -q 20 -Q 20 -r {} {}".format(pos_str, bam_path))
            try:
                mpileup_rec = mpileup_rec.split()[:-1]
                pile_up_result = Counter(x.upper() for x in mpileup_rec[4] if x in bases)
                depth = len(mpileup_rec[4])
            except IndexError:
                pile_up_result = Counter()
                depth = 0
            sample_dict = {s.sample: s.gt_bases for s in vcf_rec.samples}
            for s in vcf_rec.samples:
                if len([x for x in sample_dict.itervalues() if x == s.gt_bases]) != 1:
                    continue
                c = format_ratio(pile_up_result[s.gt_bases], depth)
                range_results[s.sample].append(c)
        avg_results = {}
        for para, vals in range_results.iteritems():
            vals = np.array(vals)
            avg_results[para] = [min(vals), np.mean(vals), max(vals)]
        results[r[3]] = avg_results
    return results


def analyze_wgs(bam_path, vcf_handle):
    wgs_results = defaultdict(list)
    vcf_recs = list(vcf_handle.fetch('Notch2NL_consensus', 0, 1000000))
    for vcf_rec in vcf_recs:
        if vcf_rec.is_indel:
            continue
        pos_str = "{0}:{1}-{1}".format(vcf_rec.CHROM, vcf_rec.POS)
        mpileup_rec = popenCatch("samtools mpileup -q 20 -Q 20 -r {} {}".format(pos_str, bam_path))
        mpileup_rec = mpileup_rec.split()[:-1]
        pile_up_result = Counter(x.upper() for x in mpileup_rec[4] if x in bases)
        sample_dict = {s.sample: s.gt_bases for s in vcf_rec.samples}
        for s in vcf_rec.samples:
            if len([x for x in sample_dict.itervalues() if x == s.gt_bases]) != 1:
                continue
            c = 1.0 * pile_up_result[s.gt_bases] / len(mpileup_rec[4])
            wgs_results[s.sample].append([vcf_rec.POS, c])
    return wgs_results


def make_plot(wgs_results, out_path, sample):
    paralogs = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D']
    fig, plots = plt.subplots(5, sharey=True, sharex=True)
    plt.yticks((0, 0.1, 0.2, 0.3, 0.4))
    plt.ylim((0, 0.4))
    plt.xticks((0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))
    plt.xlim((0, 100143))
    plt.xlabel("Alignment position")
    for i, (p, para) in enumerate(zip(plots, paralogs)):
        p.set_title(para)
        wgs = wgs_results[para]
        xvals, yvals = zip(*wgs)
        p.vlines(xvals, np.zeros(len(xvals)), yvals, color=sns.color_palette()[0], alpha=0.7, linewidth=0.8)
        # mark the zeros
        zero_wgs = [[x, y + 0.02] for x, y in wgs if y == 0]
        if len(zero_wgs) > 0:
            z_xvals, z_yvals = zip(*zero_wgs)
            p.vlines(z_xvals, np.zeros(len(z_xvals)), z_yvals, color=sns.color_palette()[2], alpha=0.7, linewidth=0.8)
    plt.tight_layout(pad=2.5, h_pad=0.3)
    fig.suptitle("{} combined amplicon results".format(sample))
    zero_line = matplotlib.lines.Line2D([], [], color=sns.color_palette()[2])
    reg_line = matplotlib.lines.Line2D([], [], color=sns.color_palette()[0])
    fig.legend(handles=(reg_line, zero_line), labels=["WGS SUN Fraction", "WGS Missing SUN"], loc="upper right")
    fig.text(0.01, 0.5, 'SUN fraction of reads', va='center', rotation='vertical')
    plt.savefig(out_path, format="pdf")
    plt.close()


def main():
    args = parse_args()
    ranges = [x.split() for x in open(args.ranges) if not x.startswith("#")]
    vcf_handle = vcf.Reader(open(args.consensus_vcf))
    wgs_results = analyze_wgs(args.wgs_bam, vcf_handle)
    out_path = os.path.join(args.outdir, "{}.wgs_amplicon_combined.pdf".format(args.sample))
    make_plot(wgs_results, out_path, args.sample)


if __name__ == "__main__":
    main()