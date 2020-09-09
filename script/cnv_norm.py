import sys

tss_exp_input = sys.argv[1]
tss_cnv_input = sys.argv[2]
tss_exp_output = tss_exp_input.replace("gc.bed", "gc_cnv.bed")

def calcLocalMean(position, chrom, tss_cnv_input):
    LOG2 = open(tss_cnv_input, "r")
    header = LOG2.readline()
    norm = None
    for line in LOG2.readlines():
        info = line.split("\t")
        if chrom == info[1] and position > int(info[2]) and position < int(info[3]):
            norm = float(info[6].rstrip())
            break
    LOG2.close()

    if norm:
        return norm
    else:
        return 1


with open(tss_exp_input, "r") as tss_exp_in:
    with open(tss_exp_output, "w") as tss_exp_out:
        for tss_exp_line in tss_exp_in:
            tss_exp_line_list = tss_exp_line.rstrip().split("\t")
            norrate = calcLocalMean( (int(tss_exp_line_list[1]) + 1000) , tss_exp_line_list[0], tss_cnv_input)
            norexp = float(tss_exp_line_list[8]) / norrate
            norexp_str = str(format(norexp, '.3f'))

            tss_exp_out.write(tss_exp_line.rstrip() + "\t" + norexp_str + "\n")


