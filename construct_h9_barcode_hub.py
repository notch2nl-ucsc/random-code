"""
Takes the H9 wildtpye assembly and associated transcript assembly, as well as a paired set of deletion mutant BAMs/fastqs

Constructs a assembly hub out of the wildtype, creating transcript tracks.
For each BAM/fastq pair, produces a single bigwig track.


"""

import argparse
import numpy as np
import seaborn as sns
from tools.fileOps import *
from tools.procOps import *
from tools.transcripts import *
from tools.bio import *
from pyfasta import Fasta
import pysam
from collections import *
from itertools import *


hub_str = '''hub h9_hub
shortLabel {0}
longLabel {0}
genomesFile genomes.txt
email NoEmail

'''


genomes_str = '''genome Notch2NL_scaffolds
twoBitPath scaffolds/scaffolds.2bit
trackDb scaffolds/trackDb.txt
organism Notch2NL_scaffolds
description Notch2NL_scaffolds
scientificName Notch2NL_scaffolds
defaultPos A1:1-10000000

'''


annot_str = '''track notch2nl
shortLabel NOTCH2NL
longLabel NOTCH2NL
bigDataUrl annotations.bb
type bigBed 12
visibility pack
priority 1

'''

mappability_str = '''track mappability
type bigWig
shortLabel 150mer mappability
longLabel 150mer mappability
visibility full
bigDataUrl scaffolded_map_150.bw
viewLimits 0:1
maxHeightPixels 30:30:30
priority 2

'''


wig_str = '''track {name}
type bigWig
shortLabel {name}
longLabel {name}
visibility full
bigDataUrl {path}
viewLimits 0:20
maxHeightPixels 40:40:40
color {color}
priority {priority}

'''


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scaffold-fa', default='/hive/users/ifiddes/notch2nl_berkeley_data/E2del19N_E2del68_combined_longranger/E2del68_E2del19N_combined/new-assembly/scaffolded.fa')
    parser.add_argument('--transcript-fa', default='/hive/users/ifiddes/notch2nl_berkeley_data/E2del19N_E2del68_combined_longranger/E2del68_E2del19N_combined/new-assembly/transcripts.fa')
    parser.add_argument('--mappability-bw', default='/hive/users/ifiddes/notch2nl_berkeley_data/E2del19N_E2del68_combined_longranger/E2del68_E2del19N_combined/new-assembly/scaffolded_map_150.bw')
    parser.add_argument('--names', required=True, nargs='+', help='Names of base-level directories to dig through. '
        'Expects each directory have $name/$name.consensus_mapped.sorted.bam and $name/$name.fastq')
    parser.add_argument('--out-hub', default='hub')
    parser.add_argument('--hub-name', default='H9 wildtype hub')
    return parser.parse_args()


def get_colors():
    colors = sns.color_palette()
    for c in colors:
        c = map(int, (np.array(c) * 256).round())
        yield ','.join(map(str, c))


def add_transcripts(scaffold_fa, scaffold_sizes, transcript_fa, out_tx):
    transcripts = Fasta(transcript_fa)
    scaffold = Fasta(scaffold_fa)
    beds = []
    txs = []
    assert len(scaffold) == 10, scaffold.keys()
    for n, seq in scaffold.iteritems():
        if n not in ['N1', 'N2']:
            tx_name = 'NOTCH2NL-' + '-'.join(n)
        else:
            tx_name = 'NOTCH2-{}'.format(n[-1])
        tx = transcripts[tx_name]
        with TemporaryFilePath() as tmp_tx, TemporaryFilePath() as tmp_fa, TemporaryFilePath() as tmp_psl:
            write_fasta(tmp_tx, tx_name, str(tx))
            write_fasta(tmp_fa, n, str(seq))
            cmd = ['blat', '-noHead', tmp_fa, tmp_tx, tmp_psl]
            run_proc(cmd)
            cmd = ['pslToBed', tmp_psl, '/dev/stdout']
            r = call_proc_lines(cmd)
            tx = sorted([Transcript(x.split()) for x in r], key=len)[-1]
            assert len(tx.exon_intervals) == 5
            txs.append(tx)
    assert len(txs) == 10
    with TemporaryFilePath() as tmp_bed:
        with open(tmp_bed, 'w') as outf:
            for tx in txs:
                print_row(outf, tx.get_bed())
        cmd = ['bedSort', tmp_bed, tmp_bed]
        run_proc(cmd)
        cmd = ['bedToBigBed', tmp_bed, scaffold_sizes, out_tx]
        run_proc(cmd)
    return txs


def construct_wiggles(scaffolds_bam, reads_fq, scaffold_fa, out_wig):
    with TemporaryFilePath() as tmp_fq, TemporaryFilePath() as tmp_sam:
        cmd = ['samtools', 'fastq', scaffolds_bam]
        run_proc(cmd, stdout=tmp_fq)
        cmd = ['bwa', 'mem', '-t', '30', scaffold_fa, reads_fq]
        run_proc(cmd, stdout=tmp_sam)

        # load barcode map
        bcode_map = {}
        for x in open(reads_fq):
            if x.startswith('@'):
                name, bcode = x.split()
                name = name[1:-2]
                bcode_map[name] = bcode

        # divide reads by barcode
        sh = pysam.Samfile(tmp_sam)
        s = defaultdict(list)
        for x in sh:
            if not x.is_unmapped:
                bcode = bcode_map[x.qname]
                s[bcode].append(x)

        # find haplotype with highest average mapq
        by_avg = {}
        for bcode, reads in s.iteritems():
            if len(reads) < 10:
                continue
            by_tgt = defaultdict(list)
            for r in reads:
                by_tgt[r.reference_name].append([r, r.mapq])
            avgs = {}
            for rname, rname_vals in by_tgt.iteritems():
                rname_reads, mapqs = zip(*rname_vals)
                avgs[rname] = np.mean(mapqs)
            ordered = sorted(avgs.iteritems(), key=lambda x:x[1])
            best_name, best_score = ordered[-1]
            if (len(ordered) > 1 and best_score > ordered[-2][1]) or len(ordered) == 1:
                mappings = [x for x in reads if x.reference_name == best_name]
                by_avg[bcode] = mappings

        # place these separated alignments into a BAM
        with TemporaryFilePath() as tmp_bam, TemporaryFilePath() as tmp_sorted_bam:
            with pysam.Samfile(tmp_bam, 'wb', template=sh) as outf:
                for read_set in by_avg.itervalues():
                    for read in read_set:
                        outf.write(read)

            # convert this to a coverage plot
            cmd = ['sambamba', 'sort', '-o', tmp_sorted_bam, tmp_bam]
            run_proc(cmd)
            cmd = ['samtools', 'index', tmp_sorted_bam]
            run_proc(cmd)
            # get total size
            f = Fasta(scaffold_fa)
            tot = sum(len(x) for x in f.itervalues())
            cmd = ['bamCoverage', '-b', tmp_sorted_bam, '-o', out_wig, '--binSize', '5', '--normalizeTo1x', tot]
            #cmd = ['bamCoverage', '-b', tmp_sorted_bam, '-o', out_wig, '--binSize', '5', '--normalizeUsingRPKM']
            run_proc(cmd)


def construct_multi_region(scaffold_sizes, txs, e2_bed, e2_e5_bed, total_bed):
    # first, a BED for the sizes
    with open(total_bed, 'w') as outf:
        for name, end in iter_lines(scaffold_sizes):
            print_row(outf, [name, 0, end])
    # now, exon2
    with open(e2_bed, 'w') as outf:
        for tx in txs:
            assert len(tx.exon_intervals) == 5, (txs, tx)
            e2 = tx.exon_intervals[1]
            print_row(outf, [e2.chromosome, e2.start - 1500, e2.stop + 1500])
    # now, exon2-5
    with open(e2_e5_bed, 'w') as outf:
        for tx in txs:
            e2 = tx.exon_intervals[1]
            e5 = tx.exon_intervals[4]
            print_row(outf, [e2.chromosome, e2.start - 5000, e5.stop + 5000])



if __name__ == '__main__':
    args = parse_args()
    # create hub
    scaffolds_dir = os.path.join(args.out_hub, 'scaffolds')
    tx_fa = os.path.join(scaffolds_dir, 'annotations.bb')
    scaffold_2bit = os.path.join(scaffolds_dir, 'scaffolds.2bit')
    scaffold_sizes = os.path.join(scaffolds_dir, 'scaffolds.chrom.sizes')
    mappability_bw = os.path.join(scaffolds_dir, 'scaffolded_map_150.bw')
    e2_bed = os.path.join(args.out_hub, 'exon2_regions.bed')
    e2_e5_bed = os.path.join(args.out_hub, 'exon2_exon5_regions.bed')
    total_bed = os.path.join(args.out_hub, 'total_region.bed')
    ensure_dir(scaffolds_dir)

    cmd = ['faToTwoBit', args.scaffold_fa, scaffold_2bit]
    run_proc(cmd)
    cmd = ['twoBitInfo', scaffold_2bit, scaffold_sizes]
    run_proc(cmd)

    with open(os.path.join(args.out_hub, 'hub.txt'), 'w') as outf:
        outf.write(hub_str.format(args.hub_name))

    with open(os.path.join(args.out_hub, 'genomes.txt'), 'w') as outf:
        outf.write(genomes_str)

    with open(os.path.join(scaffolds_dir, 'trackDb.txt'), 'w') as outf:
        outf.write(annot_str)
        outf.write(mappability_str)
        os.link(args.mappability_bw, mappability_bw)
        txs = add_transcripts(args.scaffold_fa, scaffold_sizes, args.transcript_fa, tx_fa)
        color_iter = cycle(get_colors())
        for i, (name, color) in enumerate(zip(*[args.names, color_iter])):
            bam = os.path.join(name, name + '.consensus_mapped.sorted.bam')
            fq = os.path.join(name, name + '.fastq')
            out_wig = os.path.join(scaffolds_dir, name + '.bw')
            construct_wiggles(bam, fq, args.scaffold_fa, out_wig)
            outf.write(wig_str.format(name=name, path=os.path.basename(out_wig), color=color, priority=10 + i))

    construct_multi_region(scaffold_sizes, txs, e2_bed, e2_e5_bed, total_bed)
