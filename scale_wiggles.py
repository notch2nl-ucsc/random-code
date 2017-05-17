import pysam
from tools.procOps import *
from tools.fileOps import *
from tools.mathOps import *
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-bam', required=True)
    parser.add_argument('--tgt-bam', required=True)
    parser.add_argument('--scaling-region', default='chr1:147136783-147902914')
    parser.add_argument('--chrom-sizes', default='/cluster/data/hg38/chrom.sizes')
    parser.add_argument('--out-diff', required=True)
    parser.add_argument('--out-ref')
    parser.add_argument('--out-tgt', required=True)
    return parser.parse_args()


def find_read_depth(bam, region):
    s = pysam.Samfile(bam)
    reads = [x for x in s.fetch(region=region)]
    return len(reads)


def create_chr1_wiggle(bam, wig, scale_factor=None):
    cmd = ['bamCoverage', '-b', bam, '-o', wig,
           '--numberOfProcessors', '20', '--region', 'chr1', '--smoothLength', '250']
    if scale_factor is not None:
        cmd.extend(['--scaleFactor', scale_factor])
    run_proc(cmd)


def subtract_wiggles(ref_wiggle, tgt_scaled_wiggle):
    tmp = get_tmp_file(suffix='bw')
    cmd = ['wiggletools', 'diff', ref_wiggle, tgt_scaled_wiggle]
    run_proc(cmd, stdout=tmp)
    return tmp


def convert_wig_bigwig(subtracted_wiggle, chrom_sizes, out_bw):
    cmd = ['wigToBigWig', subtracted_wiggle, chrom_sizes, out_bw]
    run_proc(cmd)


if __name__ == '__main__':
    args = parse_args()
    ref_depth = find_read_depth(args.ref_bam, args.scaling_region)
    tgt_depth = find_read_depth(args.tgt_bam, args.scaling_region)
    scaling_factor = format_ratio(ref_depth, tgt_depth)
    print 'Scaling factor: {:.5}'.format(scaling_factor)
    if args.out_ref is None:
        out_ref = get_tmp_file()
    else:
        out_ref = args.out_ref
    if not os.path.exists(out_ref):
        create_chr1_wiggle(args.ref_bam, out_ref)
    create_chr1_wiggle(args.tgt_bam, args.out_tgt, scaling_factor)
    subtracted_wiggle = subtract_wiggles(out_ref, args.out_tgt)
    convert_wig_bigwig(subtracted_wiggle, args.chrom_sizes, args.out_diff)
    os.remove(subtracted_wiggle)
    if args.out_ref is None:
        os.remove(out_ref)
