#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pybedtools


def parse_arguments():
    parser = argparse.ArgumentParser(add_help=False)

    # define the arguments
    parser.add_argument('--bamfile', '-b',
                          metavar='BAM file',
                          help='Sorted BAM file to correct.',
                          required=True)

    parser.add_argument('--bedfile', '-d',
                          metavar='BED file',
                          help='BED file for counting.',
                          required=True)

    parser.add_argument('--mode', '-m',
                          help='single end (SE) or paired end (PE) reads?',
                          default='SE')
    
    parser.add_argument('--out_prefix', '-o',
                        help='Prefix of the corrected file.',
                        default="*default",
                        type=str)

    parser.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
    
    args = parser.parse_args()

    return args


def main(args=None):
    args = parse_arguments()
    
    # output file name
    if args.out_prefix == "*default":
        bam_file_name = os.path.basename(args.bamfile)
        out_file_raw = bam_file_name.replace(".gctmp.bam", "") + ".raw.bed"
        out_file_gc = bam_file_name.replace(".gctmp.bam", "") + ".gc.bed"
    else:
        out_file_raw = args.out_prefix + ".raw.bed"
        out_file_gc = args.out_prefix + ".gc.bed"
    
    # run bedtools
    bam_in = pybedtools.BedTool(args.bamfile)
    bed_in = pybedtools.BedTool(args.bedfile)
    
    bam_intersect = bam_in.intersect(bed_in, bed=True, wb =True) 
    
#    # save all bed lines
#    bed_list = []
#    for bed_line in bed_in:
#        bed_list.append(str(bed_line))
#    
    # single-end mode
    if args.mode == 'SE':
        bed_depth_gc = {}
        bed_depth_raw = {}
        for bam_line in bam_intersect:
            bam_str = str(bam_line)
            # read_id = split('\t')[0].split(':gc:')[0]
            try:
                read_gc = float(bam_str.split('\t')[3].split(':gc:')[1])
                read_bed = bam_str.split(',\t')[-1]
            except IndexError:
                print(bam_str)
                print('No gc score has been add!')
                continue
                
            if read_bed in bed_depth_gc:
                bed_depth_gc[read_bed] += read_gc
                bed_depth_raw[read_bed] += 1
            else:
                bed_depth_gc[read_bed] = read_gc 
                bed_depth_raw[read_bed] = 1
    
    # paired-end mode
    elif args.mode == 'PE':
        bed_depth_gc = {}
        bed_depth_raw = {}
        bed_read_id = {}

        for bam_line in bam_intersect:
            bam_str = str(bam_line)
            read_id = bam_str.split('\t')[3].split(':gc:')[0]
            read_gc = float(bam_str.split('\t')[3].split(':gc:')[1])
            read_bed = bam_str.split(',\t')[-1]
            
            if read_bed in bed_depth_gc:
                if read_id in bed_read_id[read_bed]:
                    pass
                
                else:
                    bed_read_id[read_bed].append(read_id)
                    bed_depth_gc[read_bed] += read_gc
                    bed_depth_raw[read_bed] += 1
                    
            else:
                bed_read_id[read_bed] = [read_id]
                bed_depth_gc[read_bed] = read_gc 
                bed_depth_raw[read_bed] = 1
     
    else:
        raise Exception("please choose mode 'SE' or 'PE'!")
        
    with open(out_file_gc, "w") as gc_out:
        with open(out_file_raw, "w") as raw_out:
            for bed_line in bed_in:
                bed_line_str = str(bed_line)
                if bed_line_str in bed_depth_gc:
                    gc_out.write(bed_line_str.replace("\n", "\t") + str(round(bed_depth_gc[bed_line_str], 3)) + "\n")
                    raw_out.write(bed_line_str.replace("\n", "\t") + str(bed_depth_raw[bed_line_str]) + "\n")
                else:
                    gc_out.write(bed_line_str.replace("\n", "\t") + '0.000\n')
                    raw_out.write(bed_line_str.replace("\n", "\t") + '0\n')
            

if __name__ == "__main__":
    main()
