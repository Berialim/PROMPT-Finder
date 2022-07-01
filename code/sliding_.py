# -*- codeing = utf-8 -*-
# @Time : 2021/6/15 9:57 AM
# @Author : Berial
# @File: sliding_window.py
# @Software: PyCharm
##work for antisense chromatin RNA-seq
from statsmodels.distributions.empirical_distribution import ECDF
import argparse
import concurrent.futures
import random
import time
import os

parser = argparse.ArgumentParser(description='process sliding window for chromatin RNA-seq')
parser.add_argument("-f", "--forward", type=str, metavar="forward", required=True, help="chromatin forward strand bedgrph file")
parser.add_argument("-r", "--reverse", type=str, metavar="reverse", required=True, help="chromatin reverse strand bedgrph file")
parser.add_argument("-b", "--binsize", type=int, metavar="binsize", default=200, help="binsize for sliding window default = 200bp")
parser.add_argument("-s", "--sliding", type=int, metavar="sliding", default=10, help="sliding length for sliding window. default = 10bp")
parser.add_argument("-o", "--output", type=str, metavar="output", default="result.bed", help="the output file name. default = result.bed")
parser.add_argument("-c", "--cutoff", type=float, metavar="cutoff", default=0.05, help="the pvalue cutoff of probability distribution. default = 0.05")
parser.add_argument("-p", "--threads", type=int, metavar="threads", default=23, help="the maximum threads allow to use. default = 23")
parser.add_argument("-sp", "--species", type=str, metavar="species", default="hg38", help="the specie of data")

args = parser.parse_args()
_forward_file = args.forward
_reverse_file = args.reverse
_bin = args.binsize
_sliding = args.sliding
_output_file = args.output.strip(".bed")
_threads = args.threads

start = time.perf_counter()
##read genome size
with open("/media/hp/disk4/bioinfo_student/reference/{}/sizes.genome".format(args.species), "r") as f1:
    _size = {}
    for line in f1:
        chr, len = line.strip().split("\t")
        if "_" not in chr:
            if int(len) > 2000000:
                _size[chr] = int(len)

def add_key(dict, k):
    if k not in dict.keys():
        dict[k] = ""


# generate random files
def rando_bin(_chr, _strand):
    print("start {} {}".format(_chr, _strand))
    # read known site
    with open("/media/hp/disk4/bioinfo_student/reference/{}/RefSeq_all.txt".format(args.species), "r") as f2:
        _rna = {}
        for line in f2:
            if "{}\t".format(_chr) in line:
                sep = line.strip().split("\t")
                chr, strand, start, end = sep[2:6]
                end = int(end)
                start = int(start)
                if _strand == "+":
                    end = end + 10000
                    ASstart = start - 20000
                    ASend = start
                elif _strand == "-":
                    start = start - 10000
                    ASstart = end
                    ASend = end +20000
                if strand == _strand:
                    for i in range(start, end):
                        add_key(_rna, i)
                else:
                    for i in range(ASstart, ASend):
                        add_key(_rna, i)
    _rna = _rna.keys()
    ## read bg file
    _bg = {}
    if _strand == "+":
        file = _forward_file
    else:
        file = _reverse_file
    with open(file, "r") as f1:
        for line in f1:
            if line.startswith(_chr + "\t"):
                chr, start, end, count = line.strip().split("\t")
                start = int(start)
                end = int(end)
                count = int(count)
                if chr == _chr and count != 0:
                    for i in range(start, end):
                        _bg[i] = count
                if chr != _chr:
                    break
    # random results
    _random = []
    size = _size[_chr]
    chr = _chr
    for _ in range(0, 10000):
        total = 0
        _bin_start = random.randint(0, size - _bin)
        _bin_end = _bin_start + _bin
        while _bin_start in _rna or _bin_end in _rna:
            _bin_start = random.randint(0, size - _bin)
            _bin_end = _bin_start + _bin
        for i in range(_bin_start, _bin_end):
            total += _bg.get(i, 0)
        _random.append(total)
    ecdf = ECDF(_random)
    for i in range(0, 1000000000000000):
        if ecdf(i) >= (1 - args.cutoff):
            p_cut = i
            break
    with open(_output_file, "a") as fo:
        chr = _chr
        start = 0
        end = _size[chr]
        _start = -1
        _end = -1
        _total = []
        for start1 in range(start, (end - _bin), _sliding):
            if _start == -1:
                _start = start1
            if not _total:
                for i in range(start1, (start1 + _bin)):
                    _total.append(_bg.get(i, 0))
            else:
                del _total[:_sliding]
                for i in range((start1 + _bin - _sliding), (start1 + _bin)):
                    _total.append(_bg.get(i, 0))
            if 0 in _total:
                if _end != -1:
                    fo.write("\t".join([chr, str(_start), str(_end), ".", ".", _strand]) + "\n")
                    _start = -1
                    _end = -1
                    continue
                else:
                    _start = -1
                    continue
            total = sum(_total)
            if total >= p_cut and start1 + _bin < end:
                _end = start1 + _bin
            elif total >= p_cut and start1 + _bin >= end:
                _end = end
                fo.write("\t".join([chr, str(_start), str(_end), ".", ".", _strand]) + "\n")
            elif total < p_cut and _end != -1:
                fo.write("\t".join([chr, str(_start), str(_end), ".", ".", _strand]) + "\n")
                _start = -1
                _end = -1
            elif total < p_cut and _end == -1:
                _start = -1
    print("finished {} {}".format(_chr, _strand))


def main(i):
    rando_bin(i, "+")
    rando_bin(i, "-")


if __name__ == "__main__":
    chr = []
    for i in range(1, 23):
        chr.append("chr" + str(i))
    chr.append("chrX")
    with concurrent.futures.ProcessPoolExecutor(max_workers=_threads) as executor:
        executor.map(main, chr)
    os.system("sort -k1,1 -k2,2n {}>{}".format(_output_file, args.output))
    os.system("rm {}".format(_output_file))
    end = time.perf_counter()
    print("all done")
    time = end - start
    hour = int(time / 3600)
    min = int(time / 60) - hour * 60
    seconds = int(time - hour * 3600 - min * 60)
    print("time: {}hour(s) {}minute(s) {}second(s)".format(hour, min, seconds))