# -*- codeing = utf-8 -*-
# @Time : 2021/6/21 4:16 PM
# @Author : Berial
# @File: merge_.py
# @Software: PyCharm
# -*- codeing = utf-8 -*-
# @Time : 2021/6/15 2:17 PM
# @Author : Berial
# @File: merge_final.py.py
# @Software: PyCharm
import argparse
import os
import re
import numpy as np

print("usage in my: "
      "python merge_.py  -f ~/reference/hg38/RefSeq.bed -i result.bed -e ~/reference/expression_list.txt -c ~/reference/sizes.genome")
parser = argparse.ArgumentParser(
    description='merge the sliding window result and remove windows located on expressed transcripts ')
parser.add_argument("-e", "--expressed", type=str, metavar="", required=True,
                    help="the list of expressed gene in the experiment cell")
parser.add_argument("-f", "--bedfile", type=str, metavar="", required=True, help="the bedfile of human")
parser.add_argument("-i", "--input", type=str, metavar="", required=True, help="the result of sliding window bed file")
parser.add_argument("-c", "--chromatinsize", type=str, metavar="", required=True,
                    help="the file contain the chromatin size (size.genome)")
parser.add_argument("-b", "--before", type=int, metavar="", default=50000,
                    help="the annotate region before TSS default = 50000bp")
parser.add_argument("-a", "--after", type=int, metavar="", default=2000,
                    help="the annotate region after TSS default = 2000bp")
parser.add_argument("-g", "--gap", type=int, metavar="", default=400,
                    help="the gap when we merge the peaks whose distance is less than. default = 400bp")
parser.add_argument("-o", "--output", type=str, metavar="", default="merge.bed",
                    help="the output file name. default = merge.bed")

args = parser.parse_args()


def add_dict(dict, k, v):
    if k in dict.keys():
        dict[k] += [v]
    else:
        dict[k] = [v]



# select expressed transcripts
# focus on the NM : mRNA
_expressed_bed = ""
_list = []

# read the expressed transcript id list analyzed from chromatin RNA-seq
with open(args.expressed, "r") as f1:
    for line in f1:
        _name = re.match(r".*([NX]._.*?\.\d+)", line)
        if _name:
            _name = _name.group(1)
        _list.append(_name)
# read the chromatin size from size.genomes file as a dictionary chromatin: size
with open(args.chromatinsize) as f3:
    _chr = {}  ##all chr length
    for line in f3:
        chr, length = line.strip().split("\t")
        if "_" in chr:
            continue
        _chr[chr] = int(length)
# read all known transcripts bed file and focus on the NM: mRNA transcripts
# merge all
_TSS = {}
_TSS_with_strand = {}  # all expressed transcripts TSS = chr,strand = [TSS...]
with open(args.bedfile, "r") as f2:
    with open("temp_high.bed", "w") as fo:
        for line in f2:
            # find the expressed transcripts
            _name = re.match(r".*([NX]._.*?\.\d+)", line)
            chr, start, end, name, score, strand = line.strip().split("\t")[:6]
            if _name:
                _name = _name.group(1)
            # get the _TSS dict : "chr_strand_TSS" : name_length(longest)
            if "_" in chr:
                continue
            start = int(start)
            end = int(end)
            if strand == "+":
                TSS = start
                length = int(end) - int(start)
            elif strand == "-":
                TSS = end
                length = int(end) - int(start)
            if _name in _list:
                fo.write(line)
                _key1 = "{},{}".format(chr, strand)
                add_dict(_TSS_with_strand, _key1, TSS)
            _key = "{},{},{}".format(chr, strand, TSS)
            _value = "{},{}".format(name, length)
            if _key not in _TSS.keys():
                _TSS[_key] = _value
            else:
                if "NM" not in name:
                    continue
                name1, length1 = _TSS[_key].split(",")
                if "NM" not in name1:
                    _TSS[_key] = _value
                length1 = int(length1)
                if length1 < length:
                    _TSS[_key] = _value
for k, v in _TSS_with_strand.items():
    _TSS_with_strand[k] = set(v)
del _list
with open("try1.txt", "w") as fo:
    for k,v in _TSS.items():
        fo.write("{},{}\n".format(k, v))

# get the upstream and down stream region of transcripts TSS(each TSS select longest transcripts as its
# transcriptID)(b-TSS-a) _TSS_s  name : TSS
with open("temp_TSS.bed", "w") as fo:
    _TSS_s = {}
    for k, v in _TSS.items():
        chr, strand, TSS = k.split(",")
        name, length = v.split(",")
        if "NM" not in name:
            continue
        TSS = int(TSS)
        _TSS_s[name] = TSS
        _key = "{},{}".format(chr, strand)
        if strand == "+":
            start = max(0, TSS - args.before)
            end = min(_chr[chr], TSS + args.after)
            strand1 = "-"
        elif strand == "-":
            end = min(_chr[chr], TSS + args.before)
            start = max(0, TSS - args.after)
            strand1 = "+"
        fo.write("\t".join([chr, str(start), str(end), name, ".", strand1]) + "\n")

# use bedtools intersect to get the annotations located in the b-TSS-a region
os.system("bedtools intersect -s -a temp_TSS.bed -b {}|sort -k1,1 -k2,2n  > temp_intersected.bed".format(args.input))
# merge annotations and remove
with open("temp_intersected.bed") as f1:
    _peak_list = {}  # the dictionary "chr,name,strand" : []
    #  peak relative located to TSS
    for line in f1:
        chr, start, end, name, score, strand = line.strip().split("\t")
        start = int(start)
        end = int(end)
        _key = chr + "," + name + "," + strand
        TSS = _TSS_s.get(name, -1)
        if TSS == -1:
            continue
        start = start - TSS
        end = end - TSS
        if _key not in _peak_list.keys():
            _peak_list[_key] = [start, end]
        else:
            _peak_list[_key] += [start, end]


def rev(x):
    y = []
    for i in x:
        y.insert(0, -i)
    return y


def rm_noice(x):
    rem = []
    for i in range(1, len(x), 2):
        if x[i] <= 400:
            rem += [i - 1, i]
    return [x[i] for i in range(0, len(x)) if i not in rem]


def merge1(x):  # merge the peaks gap less than the args.gap
    rem = []
    for i in range(2, len(x), 2):
        if x[i] - x[i - 1] <= args.gap:
            rem += [i - 1, i]
    return [x[i] for i in range(0, len(x)) if i not in rem]


def rm_far(list):  # remove merged peaks if the distance of start to TSS is longer than 2000
    # and remove peaks that end before TSS
    need = []
    for i in range(0, len(list), 2):
        if list[i] <= 2000 and list[i + 1] > 0:
            need += list[i:i + 2]
    return need


def merge_list(list, strand):
    if strand == "+":
        return rm_far(merge1(rm_noice(list)))
    else:
        a = rm_noice(rev(list))
        b = rm_far(merge1(a))
        return rev(b)


def distance(list):  # select the longest peak
    temp = []
    for i in range(0, len(list), 2):
        temp.append(list[i + 1] - list[i])
    site = temp.index(max(temp)) * 2
    return [list[site], list[site + 1]]


_peak = {}  # the dictionary "chr,name,strand":[start,end]
for k, v in _peak_list.items():
    chr, name, strand = k.split(",")
    merged = merge_list(v, strand)
    if not merged:
        continue
    merged = distance(merged)
    relative_site = np.array(merged) + _TSS_s.get(name)
    _peak[k] = relative_site.tolist()

print("_peak length={}".format(len(_peak)))


def AST_same(site, _peak):
    _same = {}
    for k,v in _peak.items():
        chr, name, strand = k.split(",")
        start, end = v
        if site == "ASS":
            _key = "{},{},{}".format(chr, strand, [start, end][strand == "-"])
            _value = [name,  [start, end][strand == "+"]]
        elif site == "AES":
            _key = "{},{},{}".format(chr, strand, [start, end][strand == "+"])
            _value = [name, [start, end][strand == "-"]]
        add_dict(_same, _key, _value)
    return _same

def closest(_same_ASS):
    peak = {}
    for k,v in _same_ASS.items():
        chr, strand, ASS = k.split(",")
        distance1 = 10**5
        ASS = int(ASS)
        for i in v:
            name, AES = i
            AES = int(AES)
            distance = abs(_TSS_s.get(name) - ASS)
            if distance < distance1:
                ASS1 = ASS
                AES1 = AES
                distance1 = distance
                name1 = name
        start = [ASS1, AES1][strand == "-"]
        end = [ASS1, AES1][strand == "+"]
        _key = "{},{},{}".format(chr, name1, strand)
        peak[_key] = [start, end]
    return peak

_same_ASS = AST_same("ASS", _peak)
_peak1 = closest(_same_ASS)
del _same_ASS
print("_peak1 length{}".format(len(_peak1)))
with open("try.txt", "w") as fo:
    for k,v in _peak1.items():
        fo.write("{},{}\n".format(k, v))

_RNU = {} # chr, strand : [TSS..] a list of all RNU TSS in this chromatin, strand
bed_name = args.bedfile.split("/")[-1]
with open(args.chromatinsize.replace("sizes.genome", "RNU.bed")) as f:
    for line in f:
        chr, start, end, name, score, strand = line.strip().split("\t")
        _key = "{},{}".format(chr, strand)
        TSS = [int(start), int(end)][strand == "-"]
        add_dict(_RNU, _key, TSS)

with open("temp_merge.bed", "w") as fo:
    for k, v in _peak1.items():
        chr, name, strand = k.split(",")
        start = v[0]
        end = v[1]
        sit = list(range(start, end))
        TSS = int([start, end][strand == "-"])
        if end - start > 1000:
            _key = "{},{}".format(chr, strand)
            i = ""
            if _key in _TSS_with_strand.keys() and _key in _RNU.keys():
                i = list(set(_RNU[_key]) & set(sit))
                if i:
                    continue
                i = list(_TSS_with_strand[_key] & set(sit)) + list(set(_RNU[_key]) & set(sit))
            elif _key in _TSS_with_strand.keys():
                i = list(_TSS_with_strand[_key] & set(sit))
            elif _key in _RNU.keys():
                i = list(set(_RNU[_key]) & set(sit))
                if i:
                    continue
            if not i:
                fo.write("\t".join([chr, str(start), str(end), "RT_" + name, ".", strand]) + "\n")
            else:
                distance1 = []
                for ii in i:
                    distance = abs(TSS - ii)
                    distance1.append(distance)
                i = i[distance1.index(min(distance1))]
                distance2 = min(distance1)
                if distance2 < 1000:
                    continue
                if strand == "+":
                    fo.write("\t".join([chr, str(start), str(i), "RT_" + name, ".", strand]) + "\n")
                elif strand == "-":
                    fo.write("\t".join([chr, str(i), str(end), "RT_" + name, ".", strand]) + "\n")
os.system("bedtools intersect -s -v -a temp_merge.bed -b temp_high.bed > {}".format(args.output))