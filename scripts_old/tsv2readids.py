import sys
 
stdin_fileno = sys.stdin

f=open(sys.argv[1],'r')
snplist=dict()
for l in f.readlines():
	if l.startswith("#"):
		continue
	l=l.strip().split("\t")
	snpid="_".join([l[0],l[1],l[3],l[4]])
	if not snpid in snplist:
		snplist[snpid]=1
f.close()

def check_snp_list(l):
	snpid="_".join([l[3],l[7],l[8],l[5]])
	if snpid in snplist:
		print(l[0])
 

for line in stdin_fileno:
	line=line.strip().split("\t")
	if line[8]=="T" and line[5]=="C":
		check_snp_list(line)
	if line[8]=="A" and line[5]=="G":
		check_snp_list(line)

