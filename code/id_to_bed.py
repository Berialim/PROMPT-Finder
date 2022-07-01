# -*- codeing = utf-8 -*-
# @Time : 2021/7/14 9:59 AM
# @Author : Berial
# @File: id_to_bed.py
# @Software: PyCharm
import sys
if len(sys.argv) <2:
    print("usage\n"
          "python id_to_bed.py input.txt output.bed NM/NR..... gene/transcript")
    sys.exit(0)
if len(sys.argv) < 4:
    print("usage\n"
          "python id_to_bed.py input.txt output.bed NM/NR..... gene/transcript")
    need = "NM"
    type = "transcript"
else:
    need = sys.argv[3]
    type = sys.argv[4]
list = []
with open(sys.argv[1]) as f1:
    for line in f1:
        list.append(line.strip("\n"))

with open("/media/hp/disk4/bioinfo_student/reference/hg38/RefSeq_all.txt") as f2:
    with open(sys.argv[2], "w") as fo:
        for line in f2:
            sep = line.strip("\n").split("\t")
            name, chr, strand, start, end = sep[1:6]
            name2 = sep[12]
            if need not in name or "_" in chr:
                continue
            id = [name, name2][type == "gene"]
            if id in list:
                fo.write("\t".join([chr, start, end, name, ".", strand]) + "\n")