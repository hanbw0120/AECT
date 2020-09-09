#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import time
import subprocess
import sys

import py2bit
import pysam
import multiprocessing
import numpy as np
import argparse

from scipy.stats import binom

from deeptools.utilities import tbitToBamChrName, getGC_content
from deeptools import writeBedGraph, parserCommon, mapReduce
from deeptools import utilities
from deeptools.bamHandler import openBam

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parentParser = parserCommon.getParentArgParse(binSize=True, blackList=False)
    requiredArgs = getRequiredArgs()
    parser = argparse.ArgumentParser(
        parents=[requiredArgs, parentParser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool corrects the GC-bias using the'
        ' method proposed by [Benjamini & Speed (2012). '
        'Nucleic Acids Research, 40(10)]. It will remove reads'
        ' from regions with too high coverage compared to the'
        ' expected values (typically GC-rich regions) and will'
        ' add reads to regions where too few reads are seen '
        '(typically AT-rich regions). '
        'The tool ``computeGCBias`` needs to be run first to generate the '
        'frequency table needed here.',
        usage='An example usage is:\n correctGCBias '
        '-b file.bam --effectiveGenomeSize 2150570000 -g mm9.2bit '
        '--GCbiasFrequenciesFile freq.txt -o gc_corrected.bam '
        '[options]',
        conflict_handler='resolve',
        add_help=False)
    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfile', '-b',
                          metavar='BAM file',
                          help='Sorted BAM file to correct.',
                          required=True)
    required.add_argument('--effectiveGenomeSize',
                          help='The effective genome size is the portion '
                          'of the genome that is mappable. Large fractions of '
                          'the genome are stretches of NNNN that should be '
                          'discarded. Also, if repetitive regions were not '
                          'included in the mapping of reads, the effective '
                          'genome size needs to be adjusted accordingly. '
                          'A table of values is available here: '
                          'http://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html .',
                          default=None,
                          type=int,
                          required=True)

    required.add_argument('--genome', '-g',
                          help='Genome in two bit format. Most genomes can be '
                          'found here: http://hgdownload.cse.ucsc.edu/gbdb/  '
                          'Search for the .2bit ending. Otherwise, fasta '
                          'files can be converted to 2bit using faToTwoBit '
                          'available here: '
                          'http://hgdownload.cse.ucsc.edu/admin/exe/',
                          metavar='two bit file',
                          required=True)

    required.add_argument('--GCbiasFrequenciesFile', '-freq',
                          help='Indicate the output file from '
                          'computeGCBias containing '
                          'the observed and expected read frequencies per GC-'
                          'content.',
                          type=argparse.FileType('r'),
                          metavar='FILE',
                          required=True)

    output = parser.add_argument_group('Output options')
    output.add_argument('--correctedFile', '-o',
                        help='Name of the corrected BAM file. The ending should '
                        'be ".bam"',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        required=True)

    # define the optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    return parser


def getReadGCcontent(tbit, read, fragmentLength, chrNameBit):
    """
    The fragments for forward and reverse reads are defined as follows::

           |- read.pos       |- read.aend
        ---+=================>-----------------------+---------    Forward strand

           |-fragStart                               |-fragEnd

        ---+-----------------------<=================+---------    Reverse strand
                                   |-read.pos        |-read.aend

           |-----------------------------------------|
                            read.tlen

    """
    fragStart = None
    fragEnd = None

    if read.is_paired and read.is_proper_pair and abs(read.tlen) < 2 * fragmentLength:
        if read.is_reverse and read.tlen < 0:
            fragEnd = read.reference_end
            fragStart = read.reference_end + read.template_length
        elif read.template_length >= read.query_alignment_length:
            fragStart = read.pos
            fragEnd = read.pos + read.template_length

    if not fragStart:
        if read.is_reverse:
            fragEnd = read.reference_end
            fragStart = read.reference_end - fragmentLength
        else:
            fragStart = read.pos
            fragEnd = fragStart + fragmentLength
    fragStart = max(0, fragStart)
    try:
        gc = getGC_content(tbit, chrNameBit, fragStart, fragEnd)
    except Exception:
        return None
    if gc is None:
        return None

    # match the gc to the given fragmentLength
    gc = int(np.round(gc * fragmentLength))
    return gc


def numCopiesOfRead(value):
    """
    Based int he R_gc value, decides
    whether to keep, duplicate, triplicate or delete the read.
    It returns an integer, that tells the number of copies of the read
    that should be keep.
    >>> np.random.seed(1)
    >>> numCopiesOfRead(0.8)
    1
    >>> numCopiesOfRead(2.5)
    2
    >>> numCopiesOfRead(None)
    1
    """
    copies = 1
    if value:
        copies = int(value) + (1 if np.random.rand() < value % 1 else 0)
    return copies


def writeCorrectedSam_wrapper(args):
    return writeCorrectedSam_worker(*args)


def writeCorrectedSam_worker(chrNameBam, chrNameBit, start, end,
                             step=None,
                             tag_but_not_change_number=True,  #HAN# modified to True, because we only count the "YC" tag for bedtools counting
                             verbose=True):
    r"""
    Writes a BAM file, deleting and adding some reads in order to compensate
    for the GC bias. **This is a stochastic method.**
    >>> np.random.seed(1)
    >>> test = Tester()
    >>> args = test.testWriteCorrectedSam()
    >>> tempFile = writeCorrectedSam_worker(*args, \
    ... tag_but_not_change_number=True, verbose=False)
    >>> try:
    ...     import StringIO
    ... except ImportError:
    ...     from io import StringIO
    >>> ostdout = sys.stdout
    >>> import tempfile
    >>> sys.stdout = tempfile.TemporaryFile()
    >>> idx = pysam.index(tempFile)
    >>> sys.stdout = ostdout
    >>> bam = pysam.Samfile(tempFile)
    >>> [dict(r.tags)['YN'] for r in bam.fetch(args[0], 200, 250)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1]
    >>> res = os.remove(tempFile)
    >>> res = os.remove(tempFile+".bai")
    >>> tempFile = \
    ... writeCorrectedSam_worker(*test.testWriteCorrectedSam_paired(),\
    ... tag_but_not_change_number=True, verbose=False)
    >>> sys.stdout = tempfile.TemporaryFile()
    >>> idx = pysam.index(tempFile)
    >>> sys.stdout = ostdout
    >>> bam = pysam.Samfile(tempFile)
    >>> [dict(r.tags)['YN'] for r in bam.fetch('chr2L', 0, 50)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> res = os.remove(tempFile)
    >>> res = os.remove(tempFile+".bai")
    """
    global R_gc
    fragmentLength = len(R_gc) - 1

    if verbose:
        print("Sam for %s %s %s " % (chrNameBit, start, end))
    i = 0

    tbit = py2bit.open(global_vars['2bit'])

    bam = openBam(global_vars['bam'])
    tempFileName = utilities.getTempFileName(suffix='.bam')

    outfile = pysam.Samfile(tempFileName, 'wb', template=bam)
    startTime = time.time()
    matePairs = {}
    read_repetitions = 0
    removed_duplicated_reads = 0

    # cache data
    # r.flag & 4 == 0 is to filter unmapped reads that
    # have a genomic position
    
    reads = [r for r in bam.fetch(chrNameBam, start, end)
             if r.pos > start and r.flag & 4 == 0 ]
    # and r.flag & 256 == 0 and r.flag & 512 == 0 and r.flag & 1024 == 0

    r_index = -1
    for read in reads:
#        if read.pos <= start or read.is_unmapped:
#            continue
        r_index += 1
        copies = None
        gc = None

        # check if a mate has already been procesed
        # to apply the same correction
        try:
            copies = matePairs[read.qname]['copies']
            gc = matePairs[read.qname]['gc']
            del(matePairs[read.qname])
        except:
            # this exception happens when a mate is
            # not present. This could
            # happen because of removal of the mate
            # by some filtering
            gc = getReadGCcontent(tbit, read, fragmentLength,
                                  chrNameBit)
            if gc:
                copies = numCopiesOfRead(float(1) / R_gc[gc])
            else:
                copies = 1
        # is this read in the same orientation and position as the previous?
        if gc and r_index > 0 and read.pos == reads[r_index - 1].pos \
                and read.is_reverse == reads[r_index - 1].is_reverse \
                and read.pnext == reads[r_index - 1].pnext:
            read_repetitions += 1
            if read_repetitions >= global_vars['max_dup_gc'][gc]:
                copies = 0  # in other words do not take into account this read
                removed_duplicated_reads += 1
        else:
            read_repetitions = 0

        readName = read.qname
        # Each tag is a tuple of (tag name, value, type)
        # Note that get_tags() returns ord(type) rather than type and this must
        # be fixed!
        # It turns out that the "with_value_type" option only started working in
        # pysam-0.8.4, so we can't reliably add tags on earlier versions without
        # potentially creating BAM files that break HTSJDK/IGV/etc.

        if gc:
            count_correct = float(round(float(1) / R_gc[gc], 4))
            
        else:
            count_correct = 0.0

        if read.is_paired and read.is_proper_pair \
                and not read.mate_is_unmapped \
                and not read.is_reverse:
            matePairs[readName] = {'copies': copies,
                                   'gc': gc}
        
        
        
        """
        outfile.write(read)
        """
        if tag_but_not_change_number:
            read.qname = readName + ':gc:' + str(count_correct) #HAN#
            outfile.write(read)
            continue

        if verbose:
            if i % 500000 == 0 and i > 0:
                endTime = time.time()
                print("{},  processing {} ({:.1f} per sec) reads "
                      "@ {}:{}-{}".format(multiprocessing.current_process().name,
                                          i, i / (endTime - startTime),
                                          chrNameBit, start, end))
        i += 1

    outfile.close()
    if verbose:
        endTime = time.time()
        print("{},  processing {} ({:.1f} per sec) reads "
              "@ {}:{}-{}".format(multiprocessing.current_process().name,
                                  i, i / (endTime - startTime),
                                  chrNameBit, start, end))
        percentage = float(removed_duplicated_reads) * 100 / len(reads) \
            if len(reads) > 0 else 0
        print("duplicated reads removed %d of %d (%.2f) " %
              (removed_duplicated_reads, len(reads), percentage))

    return tempFileName


def getFragmentFromRead(read, defaultFragmentLength, extendPairedEnds=True):
    """
    The read has to be pysam object.

    The following values are defined (for forward reads)::


             |--          -- read.tlen --              --|
             |-- read.alen --|
        -----|===============>------------<==============|----
             |               |            |
          read.pos      read.aend      read.pnext


          and for reverse reads


             |--             -- read.tlen --           --|
                                         |-- read.alen --|
        -----|===============>-----------<===============|----
             |                           |               |
          read.pnext                   read.pos      read.aend

    this is a sketch of a pair-end reads

    The function returns the fragment start and end, either
    using the paired end information (if available) or
    extending the read in the appropriate direction if this
    is single-end.

    Parameters
    ----------
    read : pysam read object


    Returns
    -------
    tuple
        (fragment start, fragment end)

    """
    # convert reads to fragments

    # this option indicates that the paired ends correspond
    # to the fragment ends
    # condition read.tlen < maxPairedFragmentLength is added to avoid read pairs
    # that span thousands of base pairs

    if extendPairedEnds is True and read.is_paired and 0 < abs(read.tlen) < 1000:
        if read.is_reverse:
            fragmentStart = read.pnext
            fragmentEnd = read.aend
        else:
            fragmentStart = read.pos
            # the end of the fragment is defined as
            # the start of the forward read plus the insert length
            fragmentEnd = read.pos + read.tlen
    else:
        if defaultFragmentLength <= read.aend - read.pos:
            fragmentStart = read.pos
            fragmentEnd = read.aend
        else:
            if read.is_reverse:
                fragmentStart = read.aend - defaultFragmentLength
                fragmentEnd = read.aend
            else:
                fragmentStart = read.pos
                fragmentEnd = read.pos + defaultFragmentLength

    return fragmentStart, fragmentEnd


def run_shell_command(command):
    """
    Runs the given shell command. Report
    any errors found.
    """
    try:
        subprocess.check_call(command, shell=True)

    except subprocess.CalledProcessError as error:
        sys.stderr.write('Error{}\n'.format(error))
        exit(1)
    except Exception as error:
        sys.stderr.write('Error: {}\n'.format(error))
        exit(1)


def main(args=None):
    args = process_args(args)
    global F_gc, N_gc, R_gc

    data = np.loadtxt(args.GCbiasFrequenciesFile.name)

    F_gc = data[:, 0]
    N_gc = data[:, 1]
    R_gc = data[:, 2]

    global global_vars
    global_vars = {}
    global_vars['2bit'] = args.genome
    global_vars['bam'] = args.bamfile

    # compute the probability to find more than one read (a redundant read)
    # at a certain position based on the gc of the read fragment
    # the binomial function is used for that
    max_dup_gc = [binom.isf(1e-7, F_gc[x], 1.0 / N_gc[x])
                  if F_gc[x] > 0 and N_gc[x] > 0 else 1
                  for x in range(len(F_gc))]

    global_vars['max_dup_gc'] = max_dup_gc

    tbit = py2bit.open(global_vars['2bit'])
    bam, mapped, unmapped, stats = openBam(args.bamfile, returnStats=True, nThreads=args.numberOfProcessors)

    global_vars['genome_size'] = sum(tbit.chroms().values())
    global_vars['total_reads'] = mapped
    global_vars['reads_per_bp'] = \
        float(global_vars['total_reads']) / args.effectiveGenomeSize

    # apply correction
    print("applying correction")
    # divide the genome in fragments containing about 4e5 reads.
    # This amount of reads takes about 20 seconds
    # to process per core (48 cores, 256 Gb memory)
    chunkSize = int(4e5 / global_vars['reads_per_bp'])

    # chromSizes: list of tuples
    chromSizes = [(bam.references[i], bam.lengths[i])
                  for i in range(len(bam.references))]

    regionStart = 0
    if args.region:
        chromSizes, regionStart, regionEnd, chunkSize = \
            mapReduce.getUserRegion(chromSizes, args.region,
                                    max_chunk_size=chunkSize)

    print("genome partition size for multiprocessing: {}".format(chunkSize))
    print("using region {}".format(args.region))
    mp_args = []
    bedGraphStep = args.binSize
    chrNameBitToBam = tbitToBamChrName(list(tbit.chroms().keys()), bam.references)
    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    print(chrNameBitToBam, chrNameBamToBit)
    c = 1
    for chrom, size in chromSizes:
        start = 0 if regionStart == 0 else regionStart
        for i in range(start, size, chunkSize):
            try:
                chrNameBamToBit[chrom]
            except KeyError:
                print("no sequence information for ")
                "chromosome {} in 2bit file".format(chrom)
                print("Reads in this chromosome will be skipped")
                continue
            length = min(size, i + chunkSize)
            mp_args.append((chrom, chrNameBamToBit[chrom], i, length,
                            bedGraphStep))
            c += 1

    pool = multiprocessing.Pool(args.numberOfProcessors)

    if args.correctedFile.name.endswith('bam'):
        if len(mp_args) > 1 and args.numberOfProcessors > 1:
            print(("using {} processors for {} "
                   "number of tasks".format(args.numberOfProcessors,
                                            len(mp_args))))

            res = pool.map_async(
                writeCorrectedSam_wrapper, mp_args).get(9999999)
        else:
            res = list(map(writeCorrectedSam_wrapper, mp_args))

        if len(res) == 1:
            command = "cp {} {}".format(res[0], args.correctedFile.name)
            run_shell_command(command)
        else:
            print("concatenating (sorted) intermediate BAMs")
            header = pysam.Samfile(res[0])
            of = pysam.Samfile(args.correctedFile.name, "wb", template=header)
            header.close()
            for f in res:
                f = pysam.Samfile(f)
                for e in f.fetch(until_eof=True):
                    of.write(e)
                f.close()
            of.close()

        print("indexing BAM")
        pysam.index(args.correctedFile.name)

        for tempFileName in res:
            os.remove(tempFileName)

    if args.correctedFile.name.endswith('bg') or \
            args.correctedFile.name.endswith('bw'):
        
        print('Not supported! Please use BAM file as output.')


if __name__ == "__main__":
    main()
