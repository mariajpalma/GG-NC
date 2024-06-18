
from collections import defaultdict
pairs = defaultdict(float)

for i in range(1,22):
	with open("hgdp.tgp.gwaspy.merged.chr"+str(i)+".merged.qced.maf0.01", "r") as chr:
		for line in chr:
			line=line.split("\t")
			if float(line[4]) > 2:
				pairs[line[0],line[1]]+=float(line[4])

with open ("IBD_cM_sum_NoDivision_segmentSize2", "w") as nf:
	for key, value in pairs.items():
		nf.write(str(key))
		nf.write("\t")
		nf.write(str(value))
		nf.write("\n")

