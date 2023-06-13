#! /usr/bin/env python

"""Estimate reasonable bin sizes from BAM read counts or depths."""
import logging
import os
import sys
import pysam
from io import StringIO
from itertools import islice
import numpy as np
import pandas as pd
from pathlib import PurePath
import argparse
from functools import wraps

def midsize_file(fnames):
    """Select the median-size file from several given filenames.
    If an even number of files is given, selects the file just below the median.
    """
    assert fnames, 'No files provided to calculate the median size.'
    return sorted(fnames, key=lambda f: os.stat(f).st_size)[(len(fnames) - 1) // 2]


def do_autobin(bam_fname, access=None, bp_per_bin=37500, min_size=20, max_size=50000, fasta=None):
    """Quickly calculate reasonable bin sizes from BAM read counts.
    Parameters
    ----------
    bam_fname : string
        BAM filename.
    access : GenomicArray
        Sequencing-accessible regions of the reference genome
    bp_per_bin : int
        Desired number of sequencing read nucleotide bases mapped to each bin.
    Returns
    -------
    An integer:
        avg. bin size
    """

    # Closes over bp_per_bin
    def depth2binsize(depth, min_size, max_size):
        if depth:
            bin_size = int(round(bp_per_bin / depth))
            if bin_size < min_size:
                logging.info("Limiting est. bin size %d to given min. %d",
                             bin_size, min_size)
                bin_size = min_size
            elif bin_size > max_size:
                logging.info("Limiting est. bin size %d to given max. %d",
                             bin_size, max_size)
                bin_size = max_size
            return bin_size

    ensure_bam_index(bam_fname)
    rc_table = idxstats(bam_fname, drop_unmapped=True, fasta=fasta)
    read_len = get_read_length(bam_fname, fasta=fasta)

    if access is not None and len(access):
        rc_table = update_chrom_length(rc_table, access)
    tgt_depth = average_depth(rc_table, read_len)

    # Clip bin sizes to specified ranges
    tgt_bin_size = depth2binsize(tgt_depth, min_size, max_size)

    # print(bam_fname, read_len, tgt_depth)
    return tgt_bin_size


def ensure_bam_index(bam_fname):
    """Ensure a BAM file is indexed, to enable fast traversal & lookup.
    For MySample.bam, samtools will look for an index in these files, in order:
    - MySample.bam.bai
    - MySample.bai
    """
    if PurePath(bam_fname).suffix == ".cram":
      if os.path.isfile(bam_fname + '.crai'):
        # MySample.cram.crai
        bai_fname = bam_fname + '.crai'
      else:
        # MySample.crai
        bai_fname = bam_fname[:-1] + 'i'
      if not is_newer_than(bai_fname, bam_fname):
         logging.info("Indexing CRAM file %s", bam_fname)
         pysam.index(bam_fname)
         bai_fname = bam_fname + '.crai'
      assert os.path.isfile(bai_fname), "Failed to generate cram index " + bai_fname
    else:
        if os.path.isfile(bam_fname + '.bai'):
            # MySample.bam.bai
            bai_fname = bam_fname + '.bai'
        else:
            # MySample.bai
            bai_fname = bam_fname[:-1] + 'i'
        if not is_newer_than(bai_fname, bam_fname):
            logging.info("Indexing BAM file %s", bam_fname)
            pysam.index(bam_fname)
            bai_fname = bam_fname + '.bai'
        assert os.path.isfile(bai_fname), "Failed to generate bam index " + bai_fname
    return bai_fname


def is_newer_than(target_fname, orig_fname):
    """Compare file modification times."""
    if not os.path.isfile(target_fname):
        return False
    return (os.stat(target_fname).st_mtime >= os.stat(orig_fname).st_mtime)


def idxstats(bam_fname, drop_unmapped=False, fasta=None):
    """Get chromosome names, lengths, and number of mapped/unmapped reads.
    Use the BAM index (.bai) to get the number of reads and size of each
    chromosome. Contigs with no mapped reads are skipped.
    """
    handle = StringIO(pysam.idxstats(bam_fname, split_lines=False, reference_filename=fasta))
    table = pd.read_csv(handle, sep='\t', header=None,
                        names=['chromosome', 'length', 'mapped', 'unmapped'])
    if drop_unmapped:
        table = table[table.mapped != 0].drop('unmapped', axis=1)
    return table


def get_read_length(bam, span=1000, fasta=None):
    """Get (median) read length from first few reads in a BAM file.
    Illumina reads all have the same length; other sequencers might not.
    Parameters
    ----------
    bam : str or pysam.Samfile
        Filename or pysam-opened BAM file.
    n : int
        Number of reads used to calculate median read length.
    """
    was_open = False
    if isinstance(bam, str):
        bam = pysam.Samfile(bam, 'rb', reference_filename=fasta)
    else:
        was_open = True
    lengths = [read.query_length for read in islice(bam, span) if read.query_length > 0]
    if was_open:
        bam.seek(0)
    else:
        bam.close()
    return np.median(lengths)


def update_chrom_length(rc_table, regions):
    if regions is not None and len(regions):
        chrom_sizes = region_size_by_chrom(regions)
        rc_table = rc_table.merge(chrom_sizes, on='chromosome', how='inner')
        rc_table['length'] = rc_table['length_y']
        rc_table = rc_table.drop(['length_x', 'length_y'], axis=1)
    return rc_table


def region_size_by_chrom(regions):
    chromgroups = regions.groupby('chromosome', sort=False)
    # sizes = chromgroups.apply(total_region_size) # XXX
    sizes = [total_region_size(g) for _key, g in chromgroups]
    return pd.DataFrame({'chromosome': regions.chromosome.drop_duplicates(),
                         'length': sizes})


def total_region_size(regions):
    """Aggregate area of all genomic ranges in `regions`."""
    return (regions.end - regions.start).sum()

    
def average_depth(rc_table, read_length):
    """Estimate the average read depth across the genome.
    Returns
    -------
    float
        Median of the per-chromosome mean read depths, weighted by chromosome
        size.
    """
    mean_depths = read_length * rc_table.mapped / rc_table.length
    return weighted_median(mean_depths, rc_table.length)


def on_weighted_array(f):
    """Ensure `a` and `w` are equal-length numpy arrays with no NaN values.
    For weighted descriptives -- `a` is the array of values, `w` is weights.
    1. Drop any cells in `a` that are NaN from both `a` and `w`
    2. Replace any remaining NaN cells in `w` with 0.
    """
    @wraps(f)
    def wrapper(a, w, **kwargs):
        if len(a) != len(w):
            raise ValueError("Unequal array lengths: a=%d, w=%d" % (len(a), len(w)))
        if not len(a):
            return np.nan
        a = np.asfarray(a)
        w = np.asfarray(w)
        # Drop a's NaN indices from both arrays
        a_nan = np.isnan(a)
        if a_nan.any():
            a = a[~a_nan]
            if not len(a):
                return np.nan
            w = w[~a_nan]
        # Fill w's NaN indices
        w_nan = np.isnan(w)
        if w_nan.any():
            w[w_nan] = 0.0
        return f(a, w, **kwargs)
    return wrapper


@on_weighted_array
def weighted_median(a, weights):
    """Weighted median of a 1-D numeric array."""
    order = a.argsort()
    a = a[order]
    weights = weights[order]
    midpoint = 0.5 * weights.sum()
    if (weights > midpoint).any():
        # Any point with the majority of total weight must be the median
        return a[weights.argmax()]
    cumulative_weight = weights.cumsum()
    midpoint_idx = cumulative_weight.searchsorted(midpoint)
    if (midpoint_idx > 0 and
        cumulative_weight[midpoint_idx-1] - midpoint < sys.float_info.epsilon):
        # Midpoint of 2 array values
        return a[midpoint_idx-1 : midpoint_idx+1].mean()
    return a[midpoint_idx]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parameters for automatically determining bin size")
    parser.add_argument('bams', nargs='+', help="Sample bams")
    parser.add_argument('-g', '--access', metavar="FILENAME", 
            help="Sequencing-accessible genomic regions")
    parser.add_argument('-f', '--fasta', metavar="FILENAME", help="Reference genome file")
    parser.add_argument('-l', '--log', metavar="FILENAME", help="Logger file")
    parser.add_argument('-o', '--output', metavar="FILENAME", help="Output file name")
    parser.add_argument('-b', '--bp-per-bin', type=float, default=37500., 
            help="""Desired average number of sequencing read bases mapped to each
                    bin. [Default: %(default)s]""")
    parser.add_argument('--max-size', metavar="BASES", type=int, default=1000000,
            help="Maximum size of target bins. [Default: %(default)s]")
    parser.add_argument('--min-size', metavar="BASES", type=int, default=20,
            help="Minimum size of target bins. [Default: %(default)s]")
    args = parser.parse_args()
    logging.basicConfig(filename=args.log, filemode="w", 
                        format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                        datefmt="%d-%M-%Y %H:%M:%S", level=logging.DEBUG)
    bam_fname = midsize_file(args.bams)
    access_fname = args.access
    access = pd.read_csv(access_fname, sep='\t', header=None, names=['chromosome', 'start', 'end'])
    bin_size = do_autobin(bam_fname=bam_fname, access=access, bp_per_bin=args.bp_per_bin, 
                max_size=args.max_size, min_size=args.min_size)
    bin_size = (bin_size // 100) * 100
    with open(args.output, 'w') as f:
        print(bin_size, file=f)