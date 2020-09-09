import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        usage="%prog -b <input bam> -d <reference bed> -g <genome2bit>",
        description='GC-bias adjustment for cell-free DNA based TSS coverage profile',
        epilog = "Written by Han Bowei (hanbw0120@foxmail.com), 2020\n"
    )

    parser.add_argument('--bam', '-b', type=str, help='Input bam file', required=True)
    parser.add_argument('--bed', '-d', type=str, help='Reference bed file', required=True)
    parser.add_argument('--effectiveGenomeSize', '-s', type=int, default=2827437033, help='effectiveGenomeSize for deepTools (default(hg19): 2827437033)')
    parser.add_argument('--genome2bit', '-g', type=str, required=True,
                        help='Genome 2bit file, download from http://hgdownload.cse.ucsc.edu/gbdb/')
    parser.add_argument('--fragment_size', '-f', type=int, default=167, help='Fragment size of cfDNA')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads (default: 1)')
    parser.add_argument('--mode', '-m', type=str, default="SE", help='Single-end(SE) or Paired-end(PE) (default: SE)')

    args = parser.parse_args()


    os.system("computeGCBias -b %s --effectiveGenomeSize %s -g %s -l %s --GCbiasFrequenciesFile %s.freq.txt -p %s"
              % (args.bam, args.effectiveGenomeSize, args.genome2bit, args.fragment_size, args.bam, args.threads))
    os.system("python ./correctGCBias_for_bedtools.py -b %s --effectiveGenomeSize %s -g %s --GCbiasFrequenciesFile %s.freq.txt -o %s.gc.bam -p %s"
              % (args.bam, args.effectiveGenomeSize, args.genome2bit, args.bam, args.bam, args.threads))
    os.system("python ./computeBEDcov.py -b %s.gc.bam -d %s -m %s" % (args.bam, args.bed, mode))