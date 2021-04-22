import sys

stdin_fileno = sys.stdin

f=open(sys.argv[1],'r')
snplist=dict()
for l in f.readlines():
	if l.startswith("#"):
		continue
	l=l.strip().split("\t")
	ref=l[3]
	alt=l[4].split(",")
	for a in alt:
		snpid="_".join([l[0],l[1],ref,a])
		if not snpid in snplist:
			snplist[snpid]=1
f.close()

def check_snp_list(l):
	snpid="_".join([l[3],l[7],l[8],l[5]])
	if snpid in snplist:
		print(l[0])

ref=sys.argv[2]
alt=sys.argv[3]

for line in stdin_fileno:
	line=line.strip().split("\t")
	if line[8]==ref and line[5]==alt:
		check_snp_list(line)