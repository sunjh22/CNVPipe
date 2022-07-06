# /usr/bin/env python

# Calculate read coverage in specified bins using pysam module

import sys
import math
import pysam
import numpy as np
import pandas as pd
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from io import StringIO

NULL_LOG2_COVERAGE = -20.0

def region_depth_count(bamfile, chrom, start, end, min_mapq):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    def filter_read(read):
        """True if the given read should be counted towards coverage."""
        return not (read.is_duplicate
                    or read.is_secondary
                    or read.is_unmapped
                    or read.is_qcfail
                    or read.mapq < min_mapq)
    
    if isinstance(bamfile, str):
        bamfile = pysam.Samfile(bamfile, 'rb')
    
    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            bases += sum([1 for p in read.positions if start <= p < end])
    depth = bases / (end - start) if end > start else 0
    row = (chrom, start, end,
           math.log(depth, 2) if depth else NULL_LOG2_COVERAGE,
           depth)
    return count, row


def interval_coverages_pileup(bed_fname, bam_fname, min_mapq, procs=1, fasta=None):
    """Calculate log2 coverages in the BAM file at each interval."""
    if procs == 1:
        table = bedcov(bed_fname, bam_fname, min_mapq, fasta)
    else:
        chunks = []
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bed_chunk, bam_fname, min_mapq, fasta)
                         for bed_chunk in to_chunks(bed_fname))
            for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
                chunks.append(table)
                rm(bed_chunk_fname)
        table = pd.concat(chunks, ignore_index=True)

    # User-supplied bins might be zero-width or reversed -- skip those
    spans = table.end - table.start
    ok_idx = (spans > 0)
    table = table.assign(depth=0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, 'depth'] = (table.loc[ok_idx, 'baseCount']
                                  / spans[ok_idx])
    ok_idx = (table['depth'] > 0)
    table.loc[ok_idx, 'log2'] = np.log2(table.loc[ok_idx, 'depth'])
    return table


def bedcov(bed_fname, bam_fname, min_mapq, fasta=None):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude low-MAPQ reads
    cmd = [bed_fname, bam_fname]
    if min_mapq and min_mapq > 0:
        cmd.extend(['-Q', bytes(min_mapq)])
    if fasta:
        cmd.extend(['--reference', fasta])
    try:
        raw = pysam.bedcov(*cmd, split_lines=False)
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
    if not raw:
        raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
    #columns = detect_bedcov_columns(raw)
    columns = ['chrom', 'start', 'end', 'gcContent', 'baseCount']
    table = pd.read_csv(StringIO(raw), sep='\t', names=columns, usecols=columns)
    return table


def gc_correction(table):
    # Use lowess to correct read depth based on GC content
    gcContent = table['gcContent']
    depth = table['log2']
    corrected = lowess(gcContent, depth, f=0.05, return_sorted = False)
    




'''
def detect_bedcov_columns(text):
    """Determine which 'bedcov' output columns to keep.

    Format is the input BED plus a final appended column with the count of
    basepairs mapped within each row's region. The input BED might have 3
    columns (regions without names), 4 (named regions), or more (arbitrary
    columns after 'gene').
    """
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount < 3:
        raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
    if tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    if tabcount == 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
    # Input BED has arbitrary columns after 'gene' -- ignore them
    fillers = ["_%d" % i for i in range(1, tabcount - 3)]
    return ['chromosome', 'start', 'end', 'gene'] + fillers + ['basecount']
'''

bed_fname = sys.argv[1]     # say bin.boundary.test.bed
bam_fname = sys.argv[2]     # say test0.bam
outputFile = sys.argv[3]    # say test0.coverage.txt
#fasta = '/home/sunjh/data3/refs/hg38/analysisSet/hg38.analysisSet.fa'
'''
with open(bed_fname, 'r') as f:
    for x in f:
        x = x.strip().split('\t')
        chrom = x[0]
        start = int(x[1])
        end = int(x[2])
        count, row = region_depth_count(bam_fname, chrom, start, end, min_mapq=1)
        print(count, row)
'''

table = interval_coverages_pileup(bed_fname, bam_fname, min_mapq=1)
table.to_csv(outputFile, sep='\t', float_format="%.4f", index=False)