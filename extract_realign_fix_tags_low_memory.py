"""
1. Find all barcodes that overlap the Notch loci
2. Walk the full BAM, including unmapped reads, for tags that match these barcodes
3. Align these reads to the consensus
4. Modify the resulting alignment to retain the barcodes
"""
import os
import sys
import pysam
import argparse
from collections import defaultdict
from tools.procOps import run_proc
from tools.fileOps import get_tmp_file
from tools.bio import reverse_complement


class ReadHolder(object):
    """holds a read"""
    def __init__(self, read):
        self.qual = read.qual
        self.bx = read.get_tag('BX')
        self.name = read.qname
        self.is_read2 = read.is_read2
        if self.is_read2:
            self.seq = reverse_complement(read.seq)
        else:
            self.seq = read.seq


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regions', nargs='+', default=['chr1:119982766-120091937',
                                                         'chr1:120702227-120803079',
                                                         'chr1:146147022-146250287',
                                                         'chr1:148559151-148669602',
                                                         'chr1:149371037-149472386'])
    parser.add_argument('--consensus_ref', required=True)
    parser.add_argument('inbam')
    parser.add_argument('outbam')
    parser.add_argument('out_fwd')
    parser.add_argument('--paired', action='store_true')
    return parser.parse_args()


def find_notch_barcodes(inbam, regions):
    """find all unique BX tags within our regions of interest. Also save the records"""
    barcodes = set()
    rec_dict = defaultdict(list)
    bam_handle = pysam.Samfile(inbam)
    for region in regions:
        aln_iter = bam_handle.fetch(region=region, multiple_iterators=True)
        for rec in aln_iter:
            if rec.is_duplicate or rec.is_secondary or rec.is_supplementary:
                continue
            try:
                tag = rec.get_tag('BX')
                barcodes.add(tag)
                rh = ReadHolder(rec)
                rec_dict[rec.qname].append(rh)
            except KeyError:
                continue
    bam_handle.close()
    return barcodes, rec_dict


def write_recs_to_fastq(rec_dict, tmp_fastq, is_paired):
    """write out all of these reads to fastq, maintaining a tag map for future use"""
    tag_map = {}
    with open(tmp_fastq, 'w') as outf:
        for qname, recs in rec_dict.iteritems():
            if is_paired:
                if len(recs) != 2:
                    continue
                l, r = recs
                if l.is_read2:
                    l, r = r, l
                assert l.bx == r.bx and l.name == r.name
                left_str = '@{}/1 {}\n{}\n+\n{}\n'.format(l.name, l.bx, l.seq, l.qual)
                right_str = '@{}/2 {}\n{}\n+\n{}\n'.format(r.name, r.bx, r.seq, r.qual)
                outf.write(left_str)
                outf.write(right_str)
                tag_map[l.name] = l.bx
            else:
                if len(recs) != 1:
                    continue
                rec = recs[0]
                left_str = '@{}/1 {}\n{}\n+\n{}\n'.format(rec.name, rec.bx, rec.seq, rec.qual)
                outf.write(left_str)
                tag_map[rec.name] = rec.bx
    return tag_map


def walk_full_bam(inbam, barcodes, out_fwd):
    """walk the full bam, including unmapped reads, and extract all reads that match our barcodes"""
    bam_handle = pysam.Samfile(inbam)
    with open(out_fwd, 'w') as out_fwd_fh:
        for rec in bam_handle.fetch(until_eof=True):
            if rec.is_duplicate or rec.is_secondary or rec.is_supplementary:
                continue
            try:
                tag = rec.get_tag('BX')
            except KeyError:
                continue
            if tag in barcodes:
                bx = rec.get_tag('BX')
                r = '/1' if rec.is_read1 else '/2'
                seq = rec.seq if rec.is_read1 else reverse_complement(rec.seq)
                out_fwd_fh.write('@{}{} {}\n{}\n+\n{}\n'.format(rec.qname, r, bx, seq, rec.qual))


def align_fastq(tmp_fastq, bwa_index, tmp_bam, is_paired):
    """run the BWA alignment"""
    bwa_cmd = [['bwa', 'mem', bwa_index, tmp_fastq, '-R', '@RG \tID:NOTCH2NL\tSM:NOTCH2NL', '-t', '8'],
               ['samtools', 'view', '-b', '-'],
               ['sambamba', 'sort', '-o', tmp_bam, '/dev/stdin']]
    if is_paired:
        bwa_cmd[0].append('-p')
    run_proc(bwa_cmd)


def add_tag_to_bam(tmp_bam, out_bam, tag_map):
    sam_handle = pysam.Samfile(tmp_bam)
    out_handle = pysam.Samfile(out_bam, 'wb', template=sam_handle)
    for rec in sam_handle:
        rec.set_tag('BX', tag_map[rec.qname])
        out_handle.write(rec)


def main():
    args = parse_args()
    barcodes, notch_recs = find_notch_barcodes(args.inbam, args.regions)
    tmp_fastq = get_tmp_file()
    tmp_bam = get_tmp_file(suffix='bam')
    tag_map = write_recs_to_fastq(notch_recs, tmp_fastq, args.paired)
    align_fastq(tmp_fastq, args.consensus_ref, tmp_bam, args.paired)
    add_tag_to_bam(tmp_bam, args.outbam, tag_map)
    run_proc(['samtools', 'index', args.outbam])
    os.remove(tmp_bam)
    os.remove(tmp_fastq)
    walk_full_bam(args.inbam, barcodes, args.out_fwd)


if __name__ == '__main__':
    main()
